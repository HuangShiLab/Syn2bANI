//! `syn2bani dist` — pairwise ANI over query × reference sets.
//!
//! Two stages (see [`crate::cli::compare`]): a recall-first per-enzyme
//! tag-window containment screen rejects pairs that are certainly below the
//! detection floor; survivors are refined with the validated chain-restricted
//! MLE estimator and reported with the same columns as `ani`, plus a trailing
//! `flag` column. Screened-out pairs are omitted (they would be
//! BELOW_DETECTION noise, not signal).
//!
//! This replaces the legacy v7 `TagMatcher`/`AniCalculator` path, whose
//! exact-tag `min_af=0.1` screen rejected 94% of true ≥80% ANI pairs and
//! whose GBRT debias was miscalibrated near identity.

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
pub fn run_dist(
    positional: &[PathBuf],
    ql: Option<&Path>,
    rl: Option<&Path>,
    enzymes_spec: &str,
    threads: usize,
    parallel: bool,
    verbose: bool,
    min_ani: f64,
    screen: ScreenConfig,
    refine_min_approx: f64,
    output: Option<&Path>,
) -> Result<()> {
    // Same input convention as `ani`: either both list files, or the
    // positional "queries... reference" form.
    let (query, reference): (Vec<PathBuf>, Vec<PathBuf>) = match (ql, rl) {
        (Some(q), Some(r)) => (compare::read_path_list(q)?, compare::read_path_list(r)?),
        (None, None) => {
            if positional.len() < 2 {
                anyhow::bail!(
                    "need at least one query and one reference; pass \
                     `<queries...> <reference>` or use --ql/--rl"
                );
            }
            let (qs, rs) = positional.split_at(positional.len() - 1);
            (qs.to_vec(), rs.to_vec())
        }
        _ => anyhow::bail!("--ql and --rl must be given together"),
    };
    if query.is_empty() || reference.is_empty() {
        anyhow::bail!("empty query or reference list");
    }

    let pool = crate::cli::build_pool(parallel, threads)?;
    let registry = EnzymeRegistry::new();
    let enzymes = resolve_enzymes(&registry, enzymes_spec)?;
    if enzymes.is_empty() {
        anyhow::bail!("no enzymes selected");
    }
    let geometry = chain_ani::geometry_from(&enzymes);
    let cfg = ChainAniConfig::default();

    // Every genome is loaded/digested exactly once — the legacy path
    // re-parsed and re-digested every reference per query (O(Q×R) I/O).
    let all_paths: Vec<PathBuf> = query.iter().chain(reference.iter()).cloned().collect();
    let loaded = compare::load_genomes(&all_paths, &enzymes, &pool, screen.window)?;
    let geometry = compare::geometry_union(&loaded, &geometry);
    let (q_sets, r_sets) = loaded.split_at(query.len());

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    // `ani` columns plus a trailing `flag` (appended, so shared columns keep
    // their positions and stay byte-comparable with `ani` output).
    writeln!(out, "{}\tflag", compare::ani_header(false, verbose))?;

    let pairs: Vec<(usize, usize)> = (0..q_sets.len())
        .flat_map(|qi| (0..r_sets.len()).map(move |ri| (qi, ri)))
        .collect();

    let collected: Vec<Option<String>> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(qi, ri)| {
                let q = &q_sets[qi];
                let r = &r_sets[ri];
                let (verdict, approx) = compare::screen_pair(q, r, &screen);
                match verdict {
                    Verdict::Reject => None,
                    Verdict::Refine => {
                        if refine_min_approx > 0.0 && approx < refine_min_approx {
                            return None;
                        }
                        let res = compare::refine_pair(q, r, &geometry, &cfg);
                        if min_ani > 0.0
                            && !(res.ani_gated.is_finite() && res.ani_gated >= min_ani)
                        {
                            return None;
                        }
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

    let mut n_reported = 0usize;
    for row in collected.into_iter().flatten() {
        writeln!(out, "{row}")?;
        n_reported += 1;
    }
    out.flush()?;
    eprintln!(
        "dist: {}/{} pairs passed the screen and were refined",
        n_reported,
        pairs.len()
    );
    Ok(())
}
