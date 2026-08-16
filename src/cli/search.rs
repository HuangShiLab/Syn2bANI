//! `syn2bani search` — query genomes against a sketch database.
//!
//! Two stages (see [`crate::cli::compare`]): the recall-first containment
//! screen runs against every DB sketch; survivors are refined with the
//! validated chain-restricted MLE estimator and reported with `ani`-style
//! columns plus a trailing `flag`, best hit first per query.
//!
//! This replaces the legacy v7 path, which digested queries BcgI-only,
//! lumped DB tags enzyme-agnostically, and returned GBRT-debiased ANIs that
//! underestimated genuine hits by 8–11 points (a self-hit scored 0.90).
//!
//! The database is a directory of `.s2ba` sketches (from `syn2bani sketch`)
//! or, via `--rl`, a file listing sketch paths. Sketches carry their own
//! enzyme panel (format v2); queries are digested with `--enzymes`, whose
//! default matches the `sketch` default. A panel mismatch does not crash —
//! only enzymes present in both sides contribute — but a DB sketched with
//! fewer enzymes loses sensitivity, so a warning is printed.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::ani::resolve_enzymes;
use crate::cli::compare::{self, LoadedGenome, Verdict};
use crate::core::chain_ani::{self, ChainAniConfig};
use crate::core::screen::ScreenConfig;
use crate::enzyme::EnzymeRegistry;

/// Collect the `.s2ba` paths of a sketch directory, sorted for determinism.
pub(crate) fn sketch_dir_paths(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut paths: Vec<PathBuf> = std::fs::read_dir(dir)
        .with_context(|| format!("reading sketch database {}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|ext| ext == "s2ba"))
        .collect();
    paths.sort();
    if paths.is_empty() {
        anyhow::bail!("{} contains no .s2ba sketches", dir.display());
    }
    Ok(paths)
}

/// Warn (once) when the DB panel differs from the query panel: only enzymes
/// on both sides contribute, so a mismatched DB silently loses sensitivity.
fn warn_panel_mismatch(db: &[LoadedGenome], query_panel: &[String]) {
    let mut q = query_panel.to_vec();
    q.sort();
    let mut n_mismatch = 0usize;
    for g in db {
        let mut p: Vec<String> = g.screen_keys.iter().map(|(n, _)| n.clone()).collect();
        p.sort();
        if p != q {
            n_mismatch += 1;
        }
    }
    if n_mismatch > 0 {
        eprintln!(
            "warning: {n_mismatch}/{} DB sketches were built with a different enzyme panel \
             than --enzymes ({q:?}); only shared enzymes contribute to these comparisons",
            db.len()
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run_search(
    query: &[PathBuf],
    ql: Option<&Path>,
    database: Option<&Path>,
    rl: Option<&Path>,
    output: Option<&Path>,
    threads: usize,
    parallel: bool,
    min_ani: f64,
    enzymes_spec: &str,
    screen: ScreenConfig,
    refine_min_approx: f64,
    verbose: bool,
) -> Result<()> {
    let query_paths: Vec<PathBuf> = match (ql, query.is_empty()) {
        (Some(q), true) => compare::read_path_list(q)?,
        (None, false) => query.to_vec(),
        (Some(_), false) => anyhow::bail!("pass queries either positionally or via --ql, not both"),
        (None, true) => anyhow::bail!("no queries given; pass paths or use --ql"),
    };
    let db_paths: Vec<PathBuf> = match (database, rl) {
        (Some(dir), None) => sketch_dir_paths(dir)?,
        (None, Some(list)) => compare::read_path_list(list)?,
        (Some(_), Some(_)) => {
            anyhow::bail!("pass the database either as a directory or via --rl, not both")
        }
        (None, None) => anyhow::bail!("no database given; pass a sketch directory or use --rl"),
    };

    let pool = crate::cli::build_pool(parallel, threads)?;
    let registry = EnzymeRegistry::new();
    let enzymes = resolve_enzymes(&registry, enzymes_spec)?;
    if enzymes.is_empty() {
        anyhow::bail!("no enzymes selected");
    }
    let geometry = chain_ani::geometry_from(&enzymes);
    let cfg = ChainAniConfig::default();

    // Load the whole DB once; the legacy path rebuilt every DB tag set per
    // query, which is why it needed ~9 minutes for 100 queries × 5000 sketches.
    let db = compare::load_genomes(&db_paths, &enzymes, &pool, screen.window)?;
    let queries = compare::load_genomes(&query_paths, &enzymes, &pool, screen.window)?;
    warn_panel_mismatch(
        &db,
        &enzymes.iter().map(|e| e.name.clone()).collect::<Vec<_>>(),
    );
    let mut all: Vec<&LoadedGenome> = queries.iter().collect();
    all.extend(db.iter());
    let geometry = {
        let mut g = geometry;
        for d in &all {
            for t in &d.digest.tags {
                g.entry(t.enzyme.clone())
                    .or_insert_with(|| crate::core::chain_ani::geometry_for_name(&t.enzyme));
            }
        }
        g
    };

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    writeln!(out, "{}\tflag", compare::ani_header(false, verbose))?;

    // Screen + refine all query × DB pairs in parallel; keep per-query order.
    let pairs: Vec<(usize, usize)> = (0..queries.len())
        .flat_map(|qi| (0..db.len()).map(move |di| (qi, di)))
        .collect();
    let hits: Vec<Option<String>> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(qi, di)| {
                let q = &queries[qi];
                let d = &db[di];
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

    // Group by query (order preserved) and sort each query's hits by the
    // gated estimate, best first — a search report, not a distance dump.
    let mut per_query: Vec<Vec<&String>> = vec![Vec::new(); queries.len()];
    for ((qi, _), row) in pairs.iter().zip(hits.iter()) {
        if let Some(row) = row {
            per_query[*qi].push(row);
        }
    }
    let ani_of = |row: &str| -> f64 {
        // ani_gated is the second-to-last column before the appended flag.
        row.rsplit('\t')
            .nth(2)
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(f64::NAN)
    };
    let mut n_hits = 0usize;
    for rows in &mut per_query {
        rows.sort_by(|a, b| {
            ani_of(b)
                .partial_cmp(&ani_of(a))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for row in rows.iter() {
            writeln!(out, "{row}")?;
            n_hits += 1;
        }
    }
    out.flush()?;
    eprintln!(
        "search: {} hits ({} queries x {} DB sketches, min-ani {})",
        n_hits,
        queries.len(),
        db.len(),
        min_ani
    );
    Ok(())
}
