//! `syn2bani triangle` — all-vs-all pairwise ANI over a genome set.
//!
//! Two stages (see [`crate::cli::compare`]): the recall-first containment
//! screen rejects pairs that are certainly below the detection floor;
//! survivors are refined with the validated chain-restricted MLE estimator.
//!
//! Output:
//! - `--edge-list`: one `ani`-style row (plus a trailing `flag`) per pair
//!   that passed the screen. Pairs certainly below the floor are omitted —
//!   the legacy path wrote `0.0000` for them, which is indistinguishable
//!   from a measured zero.
//! - default matrix: full symmetric matrix with `NaN` for screened-out
//!   pairs and 100 on the diagonal.
//!
//! This replaces the legacy v7 `TagMatcher`/`AniCalculator` path, on which
//! every reported ≥80% hit on dereplicated GTDB samples was spurious.

use anyhow::Result;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::ani::resolve_enzymes;
use crate::cli::compare::{self, Verdict};
use crate::core::chain_ani::{self, ChainAniConfig};
use crate::core::screen::ScreenConfig;
use crate::enzyme::EnzymeRegistry;

#[allow(clippy::too_many_arguments)]
pub fn run_triangle(
    genomes: &[PathBuf],
    ql: Option<&Path>,
    output: Option<&Path>,
    edge_list: bool,
    threads: usize,
    parallel: bool,
    enzymes_spec: &str,
    screen: ScreenConfig,
    refine_min_approx: f64,
    verbose: bool,
) -> Result<()> {
    let paths: Vec<PathBuf> = match (ql, genomes.is_empty()) {
        (Some(q), true) => compare::read_path_list(q)?,
        (None, false) => genomes.to_vec(),
        (Some(_), false) => {
            anyhow::bail!("pass genomes either positionally or via --ql, not both")
        }
        (None, true) => anyhow::bail!("no genomes given; pass paths or use --ql"),
    };
    if paths.len() < 2 {
        anyhow::bail!("triangle needs at least two genomes");
    }

    let pool = crate::cli::build_pool(parallel, threads)?;
    let registry = EnzymeRegistry::new();
    let enzymes = resolve_enzymes(&registry, enzymes_spec)?;
    if enzymes.is_empty() {
        anyhow::bail!("no enzymes selected");
    }
    let geometry = chain_ani::geometry_from(&enzymes);
    let cfg = ChainAniConfig::default();

    // Digest (or read sketches) once; the legacy path also kept every full
    // sequence in RAM for the whole run, which `Digest` does not.
    let loaded = compare::load_genomes(&paths, &enzymes, &pool, screen.window)?;
    let geometry = compare::geometry_union(&loaded, &geometry);
    let n = loaded.len();

    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();

    // Screen all pairs; refine survivors. `None` = certainly below the floor.
    let results: Vec<Option<String>> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(i, j)| {
                let q = &loaded[i];
                let r = &loaded[j];
                let (verdict, approx) = compare::screen_pair(q, r, &screen);
                match verdict {
                    Verdict::Reject => None,
                    Verdict::Refine => {
                        if refine_min_approx > 0.0 && approx < refine_min_approx {
                            return None;
                        }
                        let res = compare::refine_pair(q, r, &geometry, &cfg);
                        let mut row = compare::ani_row(
                            &q.digest.genome_id,
                            &r.digest.genome_id,
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
    let n_refined = results.iter().flatten().count();

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    if edge_list {
        // Omit pairs the estimator itself could not measure (non-finite
        // ani_gated = below its detection floor); a distance edge list with
        // millions of NaN rows is noise. The matrix mode below keeps them as
        // explicit NaN cells instead.
        writeln!(out, "{}\tflag", compare::ani_header(false, verbose))?;
        for row in results.iter().flatten() {
            let ani = row
                .rsplit('\t')
                .nth(2)
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(f64::NAN);
            if ani.is_finite() {
                writeln!(out, "{row}")?;
            }
        }
    } else {
        // Full symmetric matrix. Screened-out pairs are NaN, never 0.0000.
        let ani_of = |row: &str| -> f64 {
            row.rsplit('\t')
                .nth(2)
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(f64::NAN)
        };
        let mut matrix = vec![vec![f64::NAN; n]; n];
        for ((i, j), row) in pairs.iter().zip(results.iter()) {
            if let Some(row) = row {
                let a = ani_of(row);
                matrix[*i][*j] = a;
                matrix[*j][*i] = a;
            }
        }

        write!(out, "\t")?;
        for g in &loaded {
            write!(out, "{}\t", g.digest.genome_id)?;
        }
        writeln!(out)?;

        for i in 0..n {
            write!(out, "{}\t", loaded[i].digest.genome_id)?;
            for j in 0..n {
                if i == j {
                    write!(out, "100.0")?;
                } else if matrix[i][j].is_finite() {
                    write!(out, "{:.4}", matrix[i][j])?;
                } else {
                    write!(out, "NaN")?;
                }
                if j < n - 1 {
                    write!(out, "\t")?;
                }
            }
            writeln!(out)?;
        }
    }
    out.flush()?;
    eprintln!(
        "triangle: {}/{} pairs passed the screen and were refined",
        n_refined,
        pairs.len()
    );
    Ok(())
}
