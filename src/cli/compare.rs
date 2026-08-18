//! Shared two-stage pairwise machinery for the database-scale subcommands
//! (`dist`, `search`, `triangle`, `db search`).
//!
//! Stage 1 (screen): per-enzyme tag-window containment ([`crate::core::screen`])
//! rejects pairs that are certainly below the ANI detection floor.
//! Stage 2 (refine): screen survivors go through the validated
//! chain-restricted MLE estimator ([`crate::core::chain_ani::compute`]) — the
//! exact code path `ani` uses — and are reported with `ani`-style columns.
//!
//! The legacy v7 `TagMatcher`/`AniCalculator` (GBRT debias) path is
//! deliberately not used here: it was measured to produce 100% spurious
//! `triangle` hits and to underestimate genuine `search` hits by 8–11 ANI
//! points (see DB_SCALE_BENCHMARK.md in the paper repo).

use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::ani::{digest_all, load_sketch, Digest};
use crate::core::chain_ani::{self, ChainAniConfig, ChainAniResult, Geometry};
use crate::core::screen::{self, ScreenConfig};
use crate::enzyme::EnzymeConfig;

/// A genome ready for pairwise work: the estimator digest plus its
/// precomputed screen keys.
pub(crate) struct LoadedGenome {
    pub(crate) digest: Digest,
    /// Per-enzyme sorted-unique screen keys, from [`screen::keys_per_enzyme`].
    pub(crate) screen_keys: Vec<(String, Vec<u64>)>,
}

/// Load (digest or read sketch) every genome once, in parallel, and build
/// screen keys. `.s2ba` inputs skip digestion; their recorded enzyme panel is
/// honoured via the tags themselves.
pub(crate) fn load_genomes(
    paths: &[PathBuf],
    enzymes: &[EnzymeConfig],
    pool: &rayon::ThreadPool,
    screen_window: usize,
) -> Result<Vec<LoadedGenome>> {
    pool.install(|| {
        paths
            .par_iter()
            .map(|p| {
                let digest = if p.extension().is_some_and(|e| e == "s2ba") {
                    load_sketch(p)?
                } else {
                    digest_all(p, enzymes)?
                };
                let screen_keys = screen::keys_per_enzyme(&digest.tags, screen_window);
                Ok(LoadedGenome { digest, screen_keys })
            })
            .collect::<Result<Vec<_>>>()
    })
}

/// Geometry taken from the loaded genomes' own enzyme tables, so sketched
/// inputs do not depend on `--enzymes` matching the sketch panel.
pub(crate) fn geometry_union(genomes: &[LoadedGenome], base: &Geometry) -> Geometry {
    let mut g = base.clone();
    for d in genomes {
        for t in &d.digest.tags {
            g.entry(t.enzyme.clone())
                .or_insert_with(|| chain_ani::geometry_for_name(&t.enzyme));
        }
    }
    g
}

/// The outcome of screening one pair.
pub(crate) enum Verdict {
    /// Certainly below the detection floor; skip.
    Reject,
    /// Cannot rule out a detectable ANI; run the MLE estimator.
    Refine,
}

/// Screen one pair. `approx` is the crude containment-implied ANI, returned
/// for the optional second-tier gate only — never reported as an estimate.
pub(crate) fn screen_pair(
    q: &LoadedGenome,
    r: &LoadedGenome,
    cfg: &ScreenConfig,
) -> (Verdict, f64) {
    let (shared, min_total) = screen::shared_keys(&q.screen_keys, &r.screen_keys);
    let approx = screen::approx_ani(shared, min_total, cfg.window);
    let verdict = if cfg.passes(shared, min_total) {
        Verdict::Refine
    } else {
        Verdict::Reject
    };
    (verdict, approx)
}

/// Run the validated estimator on one pair. This is the identical call `ani`
/// makes, so refined numbers are interchangeable with `ani` output.
pub(crate) fn refine_pair(
    q: &LoadedGenome,
    r: &LoadedGenome,
    geometry: &Geometry,
    cfg: &ChainAniConfig,
) -> ChainAniResult {
    chain_ani::compute(
        &q.digest.tags,
        &r.digest.tags,
        geometry,
        q.digest.total_length,
        r.digest.total_length,
        &q.digest.contig_lens,
        &r.digest.contig_lens,
        cfg,
    )
}

/// The `ani` TSV header. Column order is frozen; new columns must be appended.
pub(crate) fn ani_header(calibrate: bool, verbose: bool) -> String {
    let mut h = String::from("query\treference\tani\tani_uniform\taf_query\taf_reference\tstd_err");
    if calibrate {
        h.push_str("\tani_cal");
    }
    h.push_str("\tsynteny_blocks\tsynteny_score\tbreakpoint_count");
    if verbose {
        h.push_str(
            "\thet_shape\tretention\tani_from_loss\tani_from_hist\tenzyme_spread\tenzyme_chi2\tper_enzyme\tn_anchors\tn_chains\tn_tags\tmax_block_anchors\tmean_block_anchors\tflag",
        );
    }
    // Appended last so every pre-existing column keeps its position.
    h.push_str("\tani_gated\tgate\tani_upper95");
    h
}

/// One `ani` TSV row (without trailing newline). Extracted verbatim from
/// `ani` so `dist`/`search`/`triangle` report identical numbers and columns.
pub(crate) fn ani_row(
    qid: &str,
    rid: &str,
    res: &ChainAniResult,
    cal: Option<f64>,
    verbose: bool,
) -> String {
    let mut line = format!(
        "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.5}",
        qid,
        rid,
        res.ani_het * 100.0,
        res.ani * 100.0,
        res.af_query,
        res.af_reference,
        res.std_err * 100.0
    );
    if let Some(c) = cal {
        line.push_str(&format!("\t{:.4}", c));
    }
    line.push_str(&format!(
        "\t{}\t{:.4}\t{}",
        res.synteny_blocks, res.synteny_score, res.breakpoint_count
    ));
    if verbose {
        line.push_str(&format!(
            "\t{:.3}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            res.het_shape,
            res.retention,
            res.ani_from_loss * 100.0,
            res.ani_from_hist * 100.0,
            res.agreement.spread * 100.0,
            res.agreement.reduced_chi2,
            if res.agreement.fits.is_empty() {
                "-".to_string()
            } else {
                res.agreement
                    .fits
                    .iter()
                    .map(|f| format!("{}:{:.2}", f.enzyme, f.ani * 100.0))
                    .collect::<Vec<_>>()
                    .join(",")
            },
            res.n_anchors,
            res.n_chains,
            res.n_tags_in_chains,
            res.max_block_anchors,
            res.mean_block_anchors,
            flag_str(res),
        ));
    }
    let gate = if !res.ani_gated.is_finite() {
        "none"
    } else if res.gate_fallback {
        "uniform_fallback"
    } else if res.het_shape.is_finite() {
        "gamma"
    } else {
        "uniform"
    };
    line.push_str(&format!(
        "\t{:.4}\t{}\t{:.4}",
        res.ani_gated * 100.0,
        gate,
        res.ani_upper95 * 100.0
    ));
    line
}

/// The reliability flag, same wording as `ani --verbose`.
pub(crate) fn flag_str(res: &ChainAniResult) -> &'static str {
    if res.below_detection {
        "BELOW_DETECTION"
    } else if res.unreliable {
        "INCONSISTENT"
    } else {
        "ok"
    }
}

/// Read a file of paths, one per line, ignoring blanks.
pub(crate) fn read_path_list(path: &Path) -> Result<Vec<PathBuf>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading path list {}", path.display()))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(PathBuf::from)
        .collect())
}

#[cfg(test)]
mod output_tests {
    use crate::core::chain_ani::ChainAniResult;
    use crate::core::mle::EnzymeAgreement;

    fn dummy() -> ChainAniResult {
        ChainAniResult {
            ani: 0.92,
            ani_from_loss: 0.90,
            ani_from_hist: 0.93,
            std_err: 0.01,
            inconsistent: false,
            ani_gated: 0.91,
            gate_fallback: false,
            unreliable: false,
            af_query: 0.5,
            af_reference: 0.5,
            n_chains: 10,
            n_anchors: 100,
            n_tags_in_chains: 200,
            synteny_blocks: 1,
            synteny_score: 1.0,
            breakpoint_count: 0,
            max_block_anchors: 10,
            mean_block_anchors: 10.0,
            ani_het: 0.91,
            het_shape: 2.0,
            retention: 0.6,
            below_detection: false,
            ani_upper95: 0.95,
            agreement: EnzymeAgreement::default(),
            strata: Vec::new(),
            chains: Vec::new(),
        }
    }

    #[test]
    fn header_and_row_column_counts_match() {
        let res = dummy();
        for verbose in [false, true] {
            let h = super::ani_header(false, verbose);
            let row = super::ani_row("q", "r", &res, None, verbose);
            assert_eq!(
                h.split('\t').count(),
                row.split('\t').count(),
                "verbose={verbose}: header {h} vs row {row}"
            );
        }
        // dist/search/triangle append `flag` after the row; the new
        // ani_upper95 column must sit inside the shared row, before it.
        let row = super::ani_row("q", "r", &res, None, false);
        assert!(row.ends_with("95.0000"), "row {row}");
    }
}

#[cfg(test)]
mod tests {
    /// Screen calibration dump (dev tool): prints per-pair screen statistics
    /// as TSV to compare against `ani` truth. Not a real test — run manually:
    ///
    /// SYN2BANI_CAL_QL=q.txt SYN2BANI_CAL_RL=r.txt \
    ///   cargo test --release screen_calibration_dump -- --ignored --nocapture
    #[test]
    #[ignore]
    fn screen_calibration_dump() {
        let ql = std::env::var("SYN2BANI_CAL_QL").expect("SYN2BANI_CAL_QL");
        let rl = std::env::var("SYN2BANI_CAL_RL").expect("SYN2BANI_CAL_RL");
        let spec = std::env::var("SYN2BANI_CAL_ENZYMES")
            .unwrap_or_else(|_| "BcgI,AlfI,AloI,FalI".to_string());
        let qlist = super::read_path_list(std::path::Path::new(&ql)).unwrap();
        let rlist = super::read_path_list(std::path::Path::new(&rl)).unwrap();
        let registry = crate::enzyme::EnzymeRegistry::new();
        let enzymes = crate::cli::ani::resolve_enzymes(&registry, &spec).unwrap();
        let window: usize = std::env::var("SYN2BANI_CAL_WINDOW")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(crate::core::screen::SCREEN_WINDOW);
        let pool = crate::cli::build_pool(true, 0).unwrap();
        let mut paths = qlist.clone();
        paths.extend(rlist.clone());
        let loaded = super::load_genomes(&paths, &enzymes, &pool, window).unwrap();
        let (qs, rs) = loaded.split_at(qlist.len());
        println!("query\treference\tn_keys_q\tn_keys_r\tshared\tmin_total\tcontainment\tapprox_ani");
        for q in qs {
            let nq: usize = q.screen_keys.iter().map(|(_, v)| v.len()).sum();
            for r in rs {
                let nr: usize = r.screen_keys.iter().map(|(_, v)| v.len()).sum();
                let (shared, min_total) =
                    crate::core::screen::shared_keys(&q.screen_keys, &r.screen_keys);
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.4}",
                    q.digest.genome_id,
                    r.digest.genome_id,
                    nq,
                    nr,
                    shared,
                    min_total,
                    shared as f64 / min_total.max(1) as f64,
                    crate::core::screen::approx_ani(shared, min_total, window) * 100.0,
                );
            }
        }
    }
}
