//! `syn2bani db` — sketch database management.
//!
//! `build`/`add` write `.s2ba` sketches; `search` compares two sketch sets
//! through the shared two-stage pipeline ([`crate::cli::compare`]): screen,
//! then refine survivors with the validated chain-restricted MLE estimator.
//!
//! `add` reads the enzyme panel recorded in an existing sketch of the
//! database and digests new genomes with exactly that panel — the legacy
//! version hardcoded BcgI, silently producing sketches inconsistent with a
//! panel-built DB.

use anyhow::{Context, Result};
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::ani::resolve_enzymes;
use crate::cli::compare::{self, Verdict};
use crate::cli::search::sketch_dir_paths;
use crate::core::chain_ani::{self, ChainAniConfig};
use crate::core::screen::ScreenConfig;
use crate::enzyme::EnzymeRegistry;
use crate::io::read_sketch;

/// Build a sketch database from a set of genomes.
pub fn run_db_build(
    genomes: &[PathBuf],
    output: &Path,
    enzyme: &str,
    threads: usize,
    parallel: bool,
    multi_enzyme: bool,
    enzyme_list: Option<&str>,
) -> Result<()> {
    crate::cli::sketch::run_sketch(
        genomes,
        output,
        enzyme,
        threads,
        parallel,
        multi_enzyme,
        enzyme_list,
    )
}

/// The enzyme panel recorded in a database, taken from its first readable
/// sketch. `None` when the database is empty (caller picks the default).
fn db_panel(database: &Path) -> Result<Option<Vec<String>>> {
    let Ok(paths) = sketch_dir_paths(database) else {
        return Ok(None);
    };
    let Some(first) = paths.first() else {
        return Ok(None);
    };
    let sk = read_sketch(first).with_context(|| format!("reading {}", first.display()))?;
    if sk.enzymes.is_empty() {
        anyhow::bail!(
            "{} is a v1 sketch with no enzyme table; cannot infer the database panel",
            first.display()
        );
    }
    Ok(Some(sk.enzymes.iter().map(|e| e.name.clone()).collect()))
}

/// Add genomes to an existing sketch database, digesting them with the
/// database's own recorded panel (or the default panel for an empty DB).
pub fn run_db_add(
    genomes: &[PathBuf],
    database: &Path,
    threads: usize,
    parallel: bool,
) -> Result<()> {
    fs::create_dir_all(database)?;
    let panel = db_panel(database)?;
    let spec = panel.unwrap_or_else(|| "BcgI,AlfI,AloI,FalI".split(',').map(String::from).collect());
    crate::cli::sketch::run_sketch(
        genomes,
        database,
        &spec.join(","),
        threads,
        parallel,
        false,
        None,
    )
}

/// Remove genome sketches from a database by ID.
pub fn run_db_remove(genome_ids: &[String], database: &Path) -> Result<()> {
    for id in genome_ids {
        let path = database.join(format!("{}.s2ba", id));
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("Failed to remove: {}", path.display()))?;
        }
    }
    Ok(())
}

/// List all entries in a sketch database.
pub fn run_db_list(database: &Path) -> Result<()> {
    for entry in fs::read_dir(database)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "s2ba").unwrap_or(false) {
            let sketch = read_sketch(&path)
                .with_context(|| format!("Failed to read: {}", path.display()))?;
            println!(
                "{}\t{}\t{}\t{:.4}",
                sketch.genome_id,
                sketch.chromosomes.len(),
                sketch.metadata.tag_count,
                sketch.metadata.gc_content
            );
        }
    }
    Ok(())
}

/// Search query sketches against a sketch database.
///
/// Both sides are directories of `.s2ba` files. Every sketch is loaded once;
/// pairs go screen → refine with the same estimator and columns as `ani`
/// (plus a trailing `flag`), filtered by `min_ani` on the gated estimate.
#[allow(clippy::too_many_arguments)]
pub fn run_db_search(
    queries: &Path,
    database: &Path,
    output: Option<&Path>,
    threads: usize,
    parallel: bool,
    min_ani: f64,
    screen: ScreenConfig,
    refine_min_approx: f64,
    verbose: bool,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;
    // Sketches are self-describing; the CLI panel only seeds geometry for
    // enzymes an input might still add via FASTA. Both sides are sketches
    // here, so the default panel is just a base for the union.
    let registry = EnzymeRegistry::new();
    let enzymes = resolve_enzymes(&registry, "BcgI,AlfI,AloI,FalI")?;
    let geometry = chain_ani::geometry_from(&enzymes);
    let cfg = ChainAniConfig::default();

    let q_paths = sketch_dir_paths(queries)?;
    let db_paths = sketch_dir_paths(database)?;
    let q_set = compare::load_genomes(&q_paths, &enzymes, &pool, screen.window)?;
    let db_set = compare::load_genomes(&db_paths, &enzymes, &pool, screen.window)?;
    let geometry = {
        let mut g = geometry;
        for d in q_set.iter().chain(db_set.iter()) {
            for t in &d.digest.tags {
                g.entry(t.enzyme.clone()).or_insert((32, 6));
            }
        }
        g
    };

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(fs::File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    writeln!(out, "{}\tflag", compare::ani_header(false, verbose))?;

    let pairs: Vec<(usize, usize)> = (0..q_set.len())
        .flat_map(|qi| (0..db_set.len()).map(move |di| (qi, di)))
        .collect();
    let hits: Vec<Option<String>> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(qi, di)| {
                let q = &q_set[qi];
                let d = &db_set[di];
                let (verdict, approx) = compare::screen_pair(q, d, &screen);
                match verdict {
                    Verdict::Reject => None,
                    Verdict::Refine => {
                        if refine_min_approx > 0.0 && approx < refine_min_approx {
                            return None;
                        }
                        let res = compare::refine_pair(q, d, &geometry, &cfg);
                        if !(res.ani_gated.is_finite() && res.ani_gated >= min_ani) {
                            return None;
                        }
                        let mut row = compare::ani_row(
                            &q.digest.genome_id,
                            &d.digest.genome_id,
                            &res,
                            None,
                            verbose,
                        );
                        row.push('\t');
                        row.push_str(compare::flag_str(&res));
                        Some(row)
                    }
                }
            })
            .collect()
    });
    let mut n_hits = 0usize;
    for row in hits.iter().flatten() {
        writeln!(out, "{row}")?;
        n_hits += 1;
    }
    out.flush()?;
    eprintln!(
        "db search: {} hits ({} queries x {} DB sketches, min-ani {})",
        n_hits,
        q_set.len(),
        db_set.len(),
        min_ani
    );
    Ok(())
}

/// Merge multiple sketch databases into one.
pub fn run_db_merge(databases: &[PathBuf], output: &Path) -> Result<()> {
    fs::create_dir_all(output)?;
    for db in databases {
        for entry in fs::read_dir(db)? {
            let entry = entry?;
            let src = entry.path();
            if src.extension().map(|e| e == "s2ba").unwrap_or(false) {
                let dst = output.join(src.file_name().unwrap());
                fs::copy(&src, &dst).with_context(|| {
                    format!("Failed to copy: {} -> {}", src.display(), dst.display())
                })?;
            }
        }
    }
    Ok(())
}
