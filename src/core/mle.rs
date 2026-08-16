//! Chain-restricted, per-enzyme-stratified ANI estimation by maximum likelihood.
//!
//! # Why not a regression
//!
//! `raw_ani` (the mean tag identity) and `mash_ani` (the containment-derived
//! estimate) are not two independent features to be combined by a learned
//! model. They are two moments of a single likelihood, and treating them as
//! features is what makes a calibration model degenerate once multi-enzyme
//! exact matching pins `raw_ani` near 1.0.
//!
//! For a Type IIB enzyme the tag contains its own recognition site, so a query
//! tag only exists if that site survived in the query. Given a query tag inside
//! a homologous (chained) region, its outcome against the reference is:
//!
//! ```text
//! found with m body mismatches   P_m(a) = C(b, m) (1-a)^m a^(k - m)     m <= tol
//! not found                      P_miss(a) = 1 - sum_{m<=tol} P_m(a)
//! ```
//!
//! where `a` is ANI, `k` is the tag length, `b = k - site_len` is the mutable
//! body, and `tol` is the accepted mismatch budget. Mismatches can only be
//! observed in the body: a mutation in the recognition site removes the tag
//! from that genome entirely, which is why the site contributes `a^site_len`
//! (folded into `a^(k-m)`) but never shows up as an observable mismatch.
//!
//! Maximising the sum of these log-probabilities over a single scalar `a`, with
//! one stratum per enzyme so each contributes its own `k` and `b`, uses the
//! mismatch histogram and the loss rate together. The relative weighting falls
//! out of the Fisher information rather than being fitted, and `tol = 0`
//! degrades gracefully to a pure loss-rate estimator instead of going flat.
//!
//! Counts must come from inside chained regions. Genome-wide counts fold the
//! accessory-genome fraction into the loss rate, which biases ANI by
//! `ln(shared_fraction) / k` — the very bias a constant offset appears to fix
//! and then fails to fix whenever shared content varies between pairs.

/// Observed outcomes for one enzyme inside the chained regions.
#[derive(Debug, Clone)]
pub struct EnzymeStratum {
    /// Enzyme name, for diagnostics only.
    pub enzyme: String,
    /// Effective tag length in bases (the full tag must match, site included).
    pub tag_len: usize,
    /// Mutable body length = tag_len - exact site - d2 - d3.
    pub body_len: usize,
    /// Number of 2-of-4 IUPAC degenerate site positions (Y R S W K M). Such a
    /// position survives a mutation w.p. `a + (1-a)/3`, and the site-preserving
    /// case shows up as a mismatch in the histogram.
    pub d2: usize,
    /// Number of 3-of-4 IUPAC degenerate site positions (V H D B). Such a
    /// position survives a mutation w.p. `a + 2(1-a)/3`.
    pub d3: usize,
    /// `hist[m]` = number of query tags matched with exactly `m` mismatches.
    /// Length is `tol + 1`.
    pub hist: Vec<u64>,
    /// Query tags inside a chained region with no acceptable match.
    pub n_miss: u64,
}

impl EnzymeStratum {
    pub fn total(&self) -> u64 {
        self.hist.iter().sum::<u64>() + self.n_miss
    }

    pub fn n_found(&self) -> u64 {
        self.hist.iter().sum()
    }
}

/// Outcome of the MLE fit.
#[derive(Debug, Clone)]
pub struct MleResult {
    /// Maximum-likelihood ANI.
    pub ani: f64,
    /// ANI implied by the loss rate alone (containment-only information).
    pub ani_from_loss: f64,
    /// ANI implied by the mismatch histogram alone (identity-only information).
    pub ani_from_hist: f64,
    /// Total query tags entering the fit.
    pub n_tags: u64,
    /// Approximate standard error from the observed Fisher information.
    pub std_err: f64,
    /// True when the two partial estimators disagree by more than ~5 standard
    /// errors. That is a data problem (repeats, contamination, a mis-built
    /// chain), not a divergence estimate — treat `ani` as untrustworthy.
    pub inconsistent: bool,
}

const A_LO: f64 = 0.50;
const A_HI: f64 = 0.999_999;

/// Gamma shape search bounds. Real substitution-rate heterogeneity fits shapes
/// of roughly 0.1-2; the upper bound only has to be far enough above that to
/// look uniform.
const ALPHA_LO: f64 = 0.1;
const ALPHA_HI: f64 = 200.0;
/// Stand-in for the alpha -> infinity (uniform-rate) limit.
const ALPHA_UNIFORM: f64 = 1.0e6;
/// Chi-square 95th percentile, 1 degree of freedom.
const LRT_CRIT: f64 = 3.841;

/// Natural log of `C(n, m)` via lgamma, safe for the sizes involved here.
fn ln_binom(n: usize, m: usize) -> f64 {
    if m > n {
        return f64::NEG_INFINITY;
    }
    ln_gamma(n as f64 + 1.0) - ln_gamma(m as f64 + 1.0) - ln_gamma((n - m) as f64 + 1.0)
}

/// Lanczos approximation of ln(gamma(x)) for x > 0.
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 8] = [
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection formula; not needed for our inputs but keeps this total.
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = 0.999_999_999_999_809_93;
    let t = x + 7.5;
    for (i, g) in G.iter().enumerate() {
        a += g / (x + i as f64 + 1.0);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// log P(found with exactly `m` body mismatches) for one stratum.
fn ln_p_found(a: f64, s: &EnzymeStratum, m: usize) -> f64 {
    if s.d2 == 0 && s.d3 == 0 {
        // No degenerate site positions: the original binomial form, kept as the
        // dedicated branch so fully-specific enzymes take the bit-identical
        // numerical path they always have.
        if m > s.body_len || m > s.tag_len {
            return f64::NEG_INFINITY;
        }
        let ln_a = a.ln();
        let ln_1ma = (1.0 - a).ln();
        ln_binom(s.body_len, m) + (m as f64) * ln_1ma + ((s.tag_len - m) as f64) * ln_a
    } else {
        ln_p_found_deg(a, s, m)
    }
}

/// log-add: `ln(exp(x) + exp(y))`.
fn ln_add(x: f64, y: f64) -> f64 {
    if x == f64::NEG_INFINITY {
        y
    } else if y == f64::NEG_INFINITY {
        x
    } else if x >= y {
        x + (y - x).exp().ln_1p()
    } else {
        y + (x - y).exp().ln_1p()
    }
}

/// log P(found with exactly `m` mismatches) for a stratum whose recognition
/// site contains IUPAC degenerate positions — exact convolution over the
/// position classes (cheap: d2 + d3 <= 3 for every enzyme in the panel).
///
/// Identity framework, mutant base uniform over the 3 alternatives. Per
/// position:
///
/// - exact site (e positions): survives w.p. `a`, never a mismatch;
/// - body (b positions): survives, mismatch w.p. `1-a`;
/// - d2 class (2-of-4): no mutation w.p. `a`, site-preserving mutation w.p.
///   `q2 = (1-a)/3` (observed as a mismatch), site-killing w.p. `2(1-a)/3`
///   (tag gone, lands in the miss count);
/// - d3 class (3-of-4): same with `q3 = 2(1-a)/3` preserving and `(1-a)/3`
///   killing.
///
/// ```text
/// P_m(a) = a^e * sum over j2,j3 with j2+j3<=m of
///   C(d2,j2) q2^j2 a^(d2-j2) * C(d3,j3) q3^j3 a^(d3-j3)
///   * C(b, m-j2-j3) (1-a)^(m-j2-j3) a^(b-(m-j2-j3))
/// ```
fn ln_p_found_deg(a: f64, s: &EnzymeStratum, m: usize) -> f64 {
    if m > s.body_len + s.d2 + s.d3 || m > s.tag_len {
        return f64::NEG_INFINITY;
    }
    let e = s.tag_len.saturating_sub(s.body_len + s.d2 + s.d3);
    let b = s.body_len;
    let ln_a = a.ln();
    let ln_1ma = (1.0 - a).ln();
    let ln_q2 = ln_1ma - 3.0_f64.ln();
    let ln_q3 = ln_1ma + (2.0_f64 / 3.0).ln();
    let mut acc = f64::NEG_INFINITY;
    for j2 in 0..=s.d2.min(m) {
        for j3 in 0..=s.d3.min(m - j2) {
            let mb = m - j2 - j3;
            if mb > b {
                continue;
            }
            let term = ln_binom(s.d2, j2)
                + j2 as f64 * ln_q2
                + (s.d2 - j2) as f64 * ln_a
                + ln_binom(s.d3, j3)
                + j3 as f64 * ln_q3
                + (s.d3 - j3) as f64 * ln_a
                + ln_binom(b, mb)
                + mb as f64 * ln_1ma
                + (b - mb) as f64 * ln_a
                + e as f64 * ln_a;
            acc = ln_add(acc, term);
        }
    }
    acc
}

/// Negative log-likelihood of the whole observation set at ANI `a`.
fn nll(a: f64, strata: &[EnzymeStratum]) -> f64 {
    let a = a.clamp(1e-9, 1.0 - 1e-12);
    let mut total = 0.0;
    for s in strata {
        let mut p_found_sum = 0.0;
        for m in 0..s.hist.len() {
            let lp = ln_p_found(a, s, m);
            if lp.is_finite() {
                p_found_sum += lp.exp();
                if s.hist[m] > 0 {
                    total -= s.hist[m] as f64 * lp;
                }
            }
        }
        if s.n_miss > 0 {
            let p_miss = (1.0 - p_found_sum).max(1e-300);
            total -= s.n_miss as f64 * p_miss.ln();
        }
    }
    total
}

/// Tag-count-weighted mean of `P(found)` across strata at ANI `a`.
///
/// Doubles as a usability signal: when this is small, tag loss carries nearly
/// all the information and the observed mismatch histogram is only the
/// truncated tail of the real distribution.
pub fn expected_retention(a: f64, strata: &[EnzymeStratum]) -> f64 {
    let total: f64 = strata.iter().map(|s| s.total() as f64).sum();
    if total <= 0.0 {
        return f64::NAN;
    }
    strata
        .iter()
        .map(|s| {
            let p: f64 = (0..s.hist.len())
                .map(|m| ln_p_found(a, s, m))
                .filter(|lp| lp.is_finite())
                .map(f64::exp)
                .sum();
            s.total() as f64 * p
        })
        .sum::<f64>()
        / total
}

/// Golden-section minimisation of a unimodal function on `[lo, hi]`.
fn minimize<F: Fn(f64) -> f64>(f: F, mut lo: f64, mut hi: f64) -> f64 {
    const INV_PHI: f64 = 0.618_033_988_749_894_9;
    let mut c = hi - INV_PHI * (hi - lo);
    let mut d = lo + INV_PHI * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..200 {
        if (hi - lo).abs() < 1e-10 {
            break;
        }
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - INV_PHI * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + INV_PHI * (hi - lo);
            fd = f(d);
        }
    }
    0.5 * (lo + hi)
}

/// ANI from the loss rate alone: solve `sum_e Q_e * P_found_e(a) = sum_e M_e`.
///
/// This is the containment-only estimator, stratified so each enzyme uses its
/// own tag length. It is what `mash_ani` should be once restricted to chains.
pub fn ani_from_loss_rate(strata: &[EnzymeStratum]) -> f64 {
    let observed: f64 = strata.iter().map(|s| s.n_found() as f64).sum();
    let total: f64 = strata.iter().map(|s| s.total() as f64).sum();
    if total <= 0.0 || observed <= 0.0 {
        return f64::NAN;
    }
    minimize(
        |a| {
            let expected: f64 = strata
                .iter()
                .map(|s| {
                    let p: f64 = (0..s.hist.len())
                        .map(|m| ln_p_found(a, s, m))
                        .filter(|lp| lp.is_finite())
                        .map(f64::exp)
                        .sum();
                    s.total() as f64 * p
                })
                .sum();
            (expected - observed).abs()
        },
        A_LO,
        A_HI,
    )
}

/// ANI from the mismatch histogram alone, conditioned on being found.
///
/// This is the information `raw_ani` tries to use, but corrected for
/// truncation: the histogram is renormalised by `P(found)` so the fact that
/// tags with more than `tol` mismatches are invisible no longer compresses the
/// estimate toward 1.
pub fn ani_from_histogram(strata: &[EnzymeStratum]) -> f64 {
    if strata.iter().all(|s| s.n_found() == 0) {
        return f64::NAN;
    }
    minimize(
        |a| {
            let mut total = 0.0;
            for s in strata {
                if s.n_found() == 0 {
                    continue;
                }
                let lps: Vec<f64> = (0..s.hist.len()).map(|m| ln_p_found(a, s, m)).collect();
                let p_found: f64 = lps.iter().filter(|lp| lp.is_finite()).map(|lp| lp.exp()).sum();
                if p_found <= 0.0 {
                    return f64::INFINITY;
                }
                let ln_p_found_total = p_found.ln();
                for (m, &count) in s.hist.iter().enumerate() {
                    if count > 0 && lps[m].is_finite() {
                        total -= count as f64 * (lps[m] - ln_p_found_total);
                    }
                }
            }
            total
        },
        A_LO,
        A_HI,
    )
}

/// Fit ANI by maximum likelihood over all strata.
pub fn estimate(strata: &[EnzymeStratum]) -> MleResult {
    let strata: Vec<EnzymeStratum> = strata.iter().filter(|s| s.total() > 0).cloned().collect();
    let n_tags: u64 = strata.iter().map(|s| s.total()).sum();
    if strata.is_empty() || n_tags == 0 {
        return MleResult {
            ani: f64::NAN,
            ani_from_loss: f64::NAN,
            ani_from_hist: f64::NAN,
            n_tags: 0,
            std_err: f64::NAN,
            inconsistent: true,
        };
    }

    let ani = minimize(|a| nll(a, &strata), A_LO, A_HI);

    // Observed Fisher information by central differences on the NLL.
    let h = 1e-5_f64.min((1.0 - ani).max(1e-9) / 4.0);
    let f0 = nll(ani, &strata);
    let fp = nll(ani + h, &strata);
    let fm = nll(ani - h, &strata);
    let curvature = (fp - 2.0 * f0 + fm) / (h * h);
    let std_err = if curvature > 0.0 {
        (1.0 / curvature).sqrt()
    } else {
        f64::NAN
    };

    let ani_from_loss = ani_from_loss_rate(&strata);
    let ani_from_hist = ani_from_histogram(&strata);

    // Consistency check: the loss rate and the histogram carry independent
    // information about the same quantity, so a large gap means the model
    // assumptions are violated rather than that the genomes are divergent.
    //
    // Only meaningful while enough of the mismatch distribution falls inside
    // the tolerance window. Once retention is low the surviving histogram is
    // just the truncated tail — it stops being an independent estimator and
    // comparing against it produces false alarms, while the joint fit is
    // already carried almost entirely by the loss rate.
    let p_found = expected_retention(ani, &strata);
    let hist_is_informative = p_found >= 0.2;
    let inconsistent = match (
        ani_from_loss.is_finite(),
        ani_from_hist.is_finite(),
        hist_is_informative,
    ) {
        (true, true, true) => {
            let tol = if std_err.is_finite() {
                (5.0 * std_err).max(0.01)
            } else {
                0.01
            };
            (ani_from_loss - ani_from_hist).abs() > tol
        }
        _ => false,
    };

    MleResult {
        ani,
        ani_from_loss,
        ani_from_hist,
        n_tags,
        std_err,
        inconsistent,
    }
}

/// One enzyme's independent fit.
#[derive(Debug, Clone)]
pub struct EnzymeFit {
    pub enzyme: String,
    pub ani: f64,
    pub std_err: f64,
    pub n_tags: u64,
}

/// Agreement between enzymes, as a check that does not share a denominator.
#[derive(Debug, Clone, Default)]
pub struct EnzymeAgreement {
    pub fits: Vec<EnzymeFit>,
    /// Inverse-variance weighted mean of the per-enzyme estimates.
    pub weighted_mean: f64,
    /// Largest minus smallest per-enzyme ANI, in the same units as `ani`.
    pub spread: f64,
    /// Cochran's Q divided by its degrees of freedom. 1.0 means the enzymes
    /// differ no more than their own standard errors allow; large values mean
    /// they are measuring genuinely different things.
    pub reduced_chi2: f64,
    /// Degrees of freedom, i.e. usable enzymes minus one.
    pub dof: usize,
}

/// Fit each enzyme on its own and measure how well they agree.
///
/// # Why this is worth having
///
/// [`MleResult::inconsistent`] compares `ani_from_loss` against
/// `ani_from_hist`, and those are computed over **the same** chain-restricted
/// tag set. That check sees the two signals disagree; it cannot see both being
/// wrong in the same direction, which is exactly what a tag sample biased toward
/// conserved regions produces. On GTDB that blind spot showed up as pairs
/// flagged `ok` scoring *worse* than pairs flagged `INCONSISTENT`.
///
/// Enzymes give an independent handle because their tag sets are disjoint and
/// their recognition sites differ in composition — site GC runs from 33% (FalI)
/// to 80% (BslFI) — so they sample different sequence contexts. If divergence
/// were uniform they would all measure the same number and disagree only by
/// sampling noise. Under mosaic divergence they need not, and the overdispersion
/// is the signal.
///
/// The homogeneous fit is used per enzyme: with one stratum and a handful of
/// histogram bins the shape parameter is not identifiable, and the quantity of
/// interest is disagreement *between* enzymes rather than each one's absolute
/// accuracy.
pub fn enzyme_agreement(strata: &[EnzymeStratum]) -> EnzymeAgreement {
    let mut fits = Vec::new();
    for s in strata {
        if s.total() == 0 || s.n_found() == 0 {
            continue;
        }
        let one = [s.clone()];
        let r = estimate(&one);
        if r.ani.is_finite() && r.std_err.is_finite() && r.std_err > 0.0 {
            fits.push(EnzymeFit {
                enzyme: s.enzyme.clone(),
                ani: r.ani,
                std_err: r.std_err,
                n_tags: s.total(),
            });
        }
    }

    if fits.len() < 2 {
        return EnzymeAgreement {
            weighted_mean: fits.first().map(|f| f.ani).unwrap_or(f64::NAN),
            spread: 0.0,
            reduced_chi2: f64::NAN,
            dof: 0,
            fits,
        };
    }

    let wsum: f64 = fits.iter().map(|f| 1.0 / (f.std_err * f.std_err)).sum();
    let mean: f64 = fits
        .iter()
        .map(|f| f.ani / (f.std_err * f.std_err))
        .sum::<f64>()
        / wsum;
    // Cochran's Q: the classic heterogeneity statistic from meta-analysis.
    let q: f64 = fits
        .iter()
        .map(|f| {
            let d = f.ani - mean;
            d * d / (f.std_err * f.std_err)
        })
        .sum();
    let dof = fits.len() - 1;
    let lo = fits.iter().map(|f| f.ani).fold(f64::INFINITY, f64::min);
    let hi = fits.iter().map(|f| f.ani).fold(f64::NEG_INFINITY, f64::max);

    EnzymeAgreement {
        weighted_mean: mean,
        spread: hi - lo,
        reduced_chi2: q / dof as f64,
        dof,
        fits,
    }
}

/// Build a stratum from an enzyme's geometry and observed counts.
///
/// `site_len` is the total constrained site length (exact + degenerate), i.e.
/// the pre-decomposition value; this constructor is the d2 = d3 = 0 special
/// case of [`stratum_deg`].
pub fn stratum(
    enzyme: &str,
    tag_len: usize,
    site_len: usize,
    hist: Vec<u64>,
    n_miss: u64,
) -> EnzymeStratum {
    stratum_deg(enzyme, tag_len, site_len, 0, 0, hist, n_miss)
}

/// Build a stratum with the full IUPAC-aware geometry.
///
/// `exact_site` counts only single-base anchor positions; `d2`/`d3` count the
/// 2-of-4 and 3-of-4 degenerate anchor positions. The body is whatever is
/// left: spacer, flanks, and any unconstrained (N) anchor positions.
pub fn stratum_deg(
    enzyme: &str,
    tag_len: usize,
    exact_site: usize,
    d2: usize,
    d3: usize,
    hist: Vec<u64>,
    n_miss: u64,
) -> EnzymeStratum {
    // The 2-bit packing used for tags caps comparisons at 32 bases, so a longer
    // tag contributes only its first 32 bases to the identity test. The site
    // sits at the start of the tag for every enzyme in the panel, so the cap
    // leaves exact_site/d2/d3 untouched.
    let tag_len = tag_len.min(32);
    let body_len = tag_len.saturating_sub(exact_site + d2 + d3);
    EnzymeStratum {
        enzyme: enzyme.to_string(),
        tag_len,
        body_len,
        d2,
        d3,
        hist,
        n_miss,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate the exact expected counts for a known ANI, then check that the
    /// MLE recovers it. This is the property that matters: no calibration.
    fn synth(ani: f64, tag_len: usize, site_len: usize, tol: usize, n: u64) -> EnzymeStratum {
        let body = tag_len - site_len;
        let mut hist = vec![0u64; tol + 1];
        let mut found = 0.0;
        for m in 0..=tol {
            let ln_c = ln_binom(body, m);
            let p = (ln_c + m as f64 * (1.0 - ani).ln() + (tag_len - m) as f64 * ani.ln()).exp();
            hist[m] = (p * n as f64).round() as u64;
            found += p;
        }
        let n_miss = ((1.0 - found) * n as f64).round() as u64;
        stratum("synth", tag_len, site_len, hist, n_miss)
    }

    #[test]
    fn recovers_known_ani_single_enzyme() {
        for &truth in &[0.85, 0.90, 0.95, 0.98, 0.99, 0.995] {
            let s = synth(truth, 32, 6, 2, 1_000_000);
            let out = estimate(&[s]);
            assert!(
                (out.ani - truth).abs() < 5e-4,
                "truth {truth}, got {} (loss {}, hist {})",
                out.ani,
                out.ani_from_loss,
                out.ani_from_hist
            );
        }
    }

    #[test]
    fn recovers_known_ani_multi_enzyme() {
        // Enzymes with different tag lengths must not bias the joint fit; a
        // single-k estimator does bias here, which is the point of stratifying.
        for &truth in &[0.90, 0.95, 0.99] {
            let strata = vec![
                synth(truth, 32, 6, 2, 300_000),
                synth(truth, 27, 7, 2, 300_000),
                synth(truth, 25, 6, 2, 300_000),
            ];
            let out = estimate(&strata);
            assert!(
                (out.ani - truth).abs() < 5e-4,
                "truth {truth}, got {}",
                out.ani
            );
            assert!(!out.inconsistent, "unexpected inconsistency at {truth}");
        }
    }

    #[test]
    fn tolerance_zero_degrades_to_loss_rate() {
        // With tol = 0 the histogram carries no shape information, so the fit
        // must fall back on the loss rate rather than going flat.
        let truth = 0.97;
        let s = synth(truth, 32, 6, 0, 1_000_000);
        let out = estimate(&[s]);
        assert!(
            (out.ani - truth).abs() < 1e-3,
            "truth {truth}, got {}",
            out.ani
        );
    }

    #[test]
    fn single_k_estimator_is_biased_on_mixed_enzymes() {
        // Sanity-check the motivation: pooling enzymes of different tag length
        // and inverting with one k really does shift the answer.
        let truth = 0.98;
        let strata = vec![
            synth(truth, 32, 6, 2, 300_000),
            synth(truth, 27, 7, 2, 300_000),
            synth(truth, 25, 6, 2, 300_000),
        ];
        let found: f64 = strata.iter().map(|s| s.n_found() as f64).sum();
        let total: f64 = strata.iter().map(|s| s.total() as f64).sum();
        let naive = (found / total).powf(1.0 / 32.0);
        let fitted = estimate(&strata).ani;
        assert!(
            (fitted - truth).abs() < (naive - truth).abs(),
            "stratified fit {fitted} should beat naive single-k {naive} (truth {truth})"
        );
    }

    #[test]
    fn flags_inconsistent_observations() {
        // Every surviving tag is a perfect match (histogram says "identical")
        // while a fifth of tags are missing (loss rate says "diverged"). Under
        // the model those cannot both hold: real divergence that deletes tags
        // also leaves 1- and 2-mismatch tags behind. Retention is high here, so
        // the histogram is informative and the check must fire.
        let s = stratum("bad", 32, 6, vec![500_000, 0, 0], 200_000);
        let out = estimate(&[s]);
        assert!(out.inconsistent, "should flag: {out:?}");
    }

    #[test]
    fn does_not_flag_when_histogram_is_uninformative() {
        // At low retention the surviving histogram is only the truncated tail,
        // so disagreement with the loss rate is expected rather than a fault.
        let truth = 0.85;
        let strata = vec![synth(truth, 32, 6, 2, 1_000_000)];
        let out = estimate(&strata);
        assert!(expected_retention(out.ani, &strata) < 0.2);
        assert!(!out.inconsistent, "false alarm at low retention: {out:?}");
        assert!((out.ani - truth).abs() < 5e-4, "got {}", out.ani);
    }

    #[test]
    fn empty_input_is_not_a_panic() {
        let out = estimate(&[]);
        assert!(out.ani.is_nan());
        assert_eq!(out.n_tags, 0);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Gamma rate heterogeneity
// ═══════════════════════════════════════════════════════════════════════════

/// Fit of the rate-heterogeneous model.
#[derive(Debug, Clone)]
pub struct HetResult {
    /// Mean nucleotide identity, `(1 + d/alpha)^-alpha`. Comparable to what an
    /// alignment-based method reports: the expected fraction of identical bases.
    pub ani: f64,
    /// Mean divergence `d` (expected substitutions per site).
    pub divergence: f64,
    /// Gamma shape `alpha`. Small = strongly heterogeneous (conserved core plus
    /// divergent regions); large = approaching a uniform rate.
    pub shape: f64,
    /// ANI the homogeneous model would report for the same counts, for contrast.
    pub ani_homogeneous: f64,
    pub n_tags: u64,
    /// True when the likelihood-ratio test accepted the extra parameter. When
    /// false, `ani` falls back to the homogeneous fit.
    pub heterogeneity_supported: bool,
    /// Chi-square statistic of the test, 1 degree of freedom.
    pub lrt: f64,
}

/// Log P(tag found with exactly `m` body mismatches) under gamma-distributed
/// regional rates.
///
/// Substitutions in a tag are Poisson with rate `r*d` per site, and because rate
/// variation is regional (kb scale) while a tag is ~30 bp, every site in a tag
/// shares one `r`. Mixing `r ~ Gamma(alpha, alpha)` turns the per-tag mismatch
/// count from Poisson into a negative binomial:
///
/// ```text
/// ln P(m) = m ln(b d) - ln m! + a ln a - lnG(a) + lnG(a+m) - (a+m) ln(a + d k)
/// ```
///
/// At `m = 0` this reduces to `(1 + d k / alpha)^-alpha`, the gamma-mixed
/// survival probability. As `alpha -> inf` the whole thing tends to the
/// homogeneous Poisson limit.
///
/// # Known residual bias
///
/// This still reads high under strong rate heterogeneity. Measured against exact
/// ground truth on mosaic simulations (`prototype/simulate_mosaic.py`), the error
/// is +0.05 at 98% ANI, +0.4 to +1.3 at 95%, and +1.8 to +5.5 at 90%, growing as
/// the fitted shape falls. It is what the +3% seen on GTDB is made of.
///
/// The likely cause is ascertainment: a query tag only exists if its recognition
/// site survived, so fast-evolving regions are depleted from the tag set before
/// any matching happens, and the rates among *observed* tags are not the
/// genome-wide `Gamma(alpha, alpha)`. Two attempts at a conjugate tilt
/// (`Gamma(alpha, alpha + d*s)`) both made things worse — one overshot to a
/// -1.8 bias, the other diverged entirely — so the derivation is not right yet
/// and the correction is deliberately not applied. Do not re-attempt it without
/// checking against `simulate_mosaic.py`, which reproduces the failure in about
/// a minute.
fn ln_p_found_het(d: f64, alpha: f64, s: &EnzymeStratum, m: usize) -> f64 {
    if m > s.body_len + s.d2 + s.d3 {
        return f64::NEG_INFINITY;
    }
    let k = s.tag_len as f64;
    // Degenerate positions extend the mismatch channel: given regional rate r,
    // body substitutions come at rate r*d per site, while a 2-of-4 (3-of-4)
    // position produces an observed, site-preserving mismatch at rate
    // r*d/3 (2*r*d/3). The mismatch count is therefore the gamma-mixed Poisson
    // (NB) with effective length b_eff, and the site-killing channels
    // (rate r*d*(2*d2/3 + d3)) contribute a closed-form survival factor
    // alpha * ln(alpha / (alpha + d*(2*d2/3 + d3))). At d2 = d3 = 0 this is
    // exactly the fully-specific form below.
    let b_eff = if s.d2 == 0 && s.d3 == 0 {
        s.body_len as f64
    } else {
        s.body_len as f64 + s.d2 as f64 / 3.0 + 2.0 * s.d3 as f64 / 3.0
    };
    let mf = m as f64;
    let bd = b_eff * d;
    if bd <= 0.0 && m > 0 {
        return f64::NEG_INFINITY;
    }
    let ln_bd = if m == 0 { 0.0 } else { bd.ln() };
    let mut lp = mf * ln_bd - ln_gamma(mf + 1.0) + alpha * alpha.ln() - ln_gamma(alpha)
        + ln_gamma(alpha + mf)
        - (alpha + mf) * (alpha + d * k).ln();
    if s.d2 + s.d3 > 0 {
        let kill = 2.0 * s.d2 as f64 / 3.0 + s.d3 as f64;
        lp += alpha * (alpha / (alpha + d * kill)).ln();
    }
    lp
}

/// Negative log-likelihood of the observations under `(d, alpha)`.
fn nll_het(d: f64, alpha: f64, strata: &[EnzymeStratum]) -> f64 {
    if !(d > 0.0) || !(alpha > 0.0) || !d.is_finite() || !alpha.is_finite() {
        return f64::INFINITY;
    }
    let mut total = 0.0;
    for s in strata {
        let mut p_found = 0.0;
        for m in 0..s.hist.len() {
            let lp = ln_p_found_het(d, alpha, s, m);
            if lp.is_finite() {
                p_found += lp.exp();
                if s.hist[m] > 0 {
                    total -= s.hist[m] as f64 * lp;
                }
            }
        }
        if s.n_miss > 0 {
            let p_miss = (1.0 - p_found).max(1e-300);
            total -= s.n_miss as f64 * p_miss.ln();
        }
    }
    if total.is_finite() {
        total
    } else {
        f64::INFINITY
    }
}

/// Mean identity under gamma-mixed rates: `E[e^{-r d}] = (1 + d/alpha)^-alpha`.
///
/// This is the quantity an alignment-based method measures — the expected
/// fraction of identical bases — so it is what should be compared against
/// skani or FastANI.
pub fn het_ani(d: f64, alpha: f64) -> f64 {
    (1.0 + d / alpha).powf(-alpha)
}

/// Fit mean divergence and rate heterogeneity jointly.
///
/// # Why this exists
///
/// The homogeneous fit assumes one divergence for the whole homologous region.
/// Real genome pairs are a mosaic: conserved core under purifying selection
/// alongside far more divergent segments. Tags survive preferentially in the
/// conserved parts, so the surviving mismatch histogram looks tighter than the
/// loss rate implies, and a single-rate fit lands between the two signals and
/// reads high. Measured against skani on 13 Enterobacteriaceae pairs, the gap
/// between the two partial estimators predicts that bias with r = 0.96 — it is
/// a measurement of heterogeneity, not noise.
///
/// With `tol = 2` the observable categories are m = 0, 1, 2 and miss: three
/// degrees of freedom for two parameters, so both are identified.
pub fn estimate_heterogeneous(strata: &[EnzymeStratum]) -> HetResult {
    let strata: Vec<EnzymeStratum> = strata.iter().filter(|s| s.total() > 0).cloned().collect();
    let n_tags: u64 = strata.iter().map(|s| s.total()).sum();
    let homogeneous = estimate(&strata);
    if strata.is_empty() || n_tags == 0 {
        return HetResult {
            ani: f64::NAN,
            divergence: f64::NAN,
            shape: f64::NAN,
            ani_homogeneous: homogeneous.ani,
            n_tags: 0,
            heterogeneity_supported: false,
            lrt: f64::NAN,
        };
    }

    // Profile over log(alpha) on a coarse grid, minimising over d at each point,
    // then refine around the best grid cell. The surface is well behaved but not
    // separable, so a plain 1-D method on either axis alone is not enough.
    let d_lo = 1e-6;
    let d_hi = 1.0;
    let mut best = (f64::INFINITY, 0.01, 1.0);
    let mut scan = |lo_ln: f64, hi_ln: f64, steps: usize, best: &mut (f64, f64, f64)| {
        for i in 0..=steps {
            let ln_a = lo_ln + (hi_ln - lo_ln) * i as f64 / steps as f64;
            let alpha = ln_a.exp();
            let d = minimize(|d| nll_het(d, alpha, &strata), d_lo, d_hi);
            let v = nll_het(d, alpha, &strata);
            if v < best.0 {
                *best = (v, d, alpha);
            }
        }
    };
    // Bounded to shapes that mean something biologically. Gamma shapes fitted
    // to real substitution-rate variation sit around 0.1-2; anything below that
    // is the optimiser running to a boundary on data that cannot identify it.
    scan(ALPHA_LO.ln(), ALPHA_HI.ln(), 48, &mut best);
    let span = (best.2.ln() - 0.35, best.2.ln() + 0.35);
    scan(span.0, span.1, 24, &mut best);

    let (nll_best, d, alpha) = best;

    // Likelihood-ratio test against this model's own alpha -> infinity limit,
    // which is a proper nested null with one degree of freedom.
    //
    // Near-identical genomes are almost all zero-mismatch tags with a tiny loss
    // rate, so alpha is unidentified there and an unconstrained fit will happily
    // run to the boundary and over-correct. Spending the second parameter only
    // when the data pay for it keeps those pairs on the homogeneous estimate.
    let d_null = minimize(|d| nll_het(d, ALPHA_UNIFORM, &strata), 1e-6, 1.0);
    let nll_null = nll_het(d_null, ALPHA_UNIFORM, &strata);
    let lrt = (2.0 * (nll_null - nll_best)).max(0.0);
    let supported = lrt > LRT_CRIT && alpha < ALPHA_HI * 0.99;

    let ani = if supported {
        het_ani(d, alpha)
    } else {
        homogeneous.ani
    };

    HetResult {
        ani,
        divergence: if supported { d } else { -homogeneous.ani.ln() },
        shape: if supported { alpha } else { f64::INFINITY },
        ani_homogeneous: homogeneous.ani,
        n_tags,
        heterogeneity_supported: supported,
        lrt,
    }
}

#[cfg(test)]
mod agreement_tests {
    use super::*;

    fn synth(ani: f64, tag_len: usize, site_len: usize, tol: usize, n: u64, name: &str)
        -> EnzymeStratum {
        let body = tag_len - site_len;
        let mut hist = vec![0u64; tol + 1];
        let mut found = 0.0;
        for m in 0..=tol {
            let p = (ln_binom(body, m) + m as f64 * (1.0 - ani).ln()
                + (tag_len - m) as f64 * ani.ln())
            .exp();
            hist[m] = (p * n as f64).round() as u64;
            found += p;
        }
        let n_miss = ((1.0 - found) * n as f64).round() as u64;
        let mut st = stratum(name, tag_len, site_len, hist, n_miss);
        st.enzyme = name.to_string();
        st
    }

    #[test]
    fn enzymes_measuring_the_same_ani_are_not_overdispersed() {
        // Every enzyme sees the same divergence, so disagreement should be
        // explainable by sampling noise alone: reduced chi-square near 1.
        let strata = vec![
            synth(0.97, 32, 6, 2, 200_000, "A"),
            synth(0.97, 27, 7, 2, 200_000, "B"),
            synth(0.97, 28, 6, 2, 200_000, "C"),
        ];
        let ag = enzyme_agreement(&strata);
        assert_eq!(ag.dof, 2);
        assert!(ag.spread < 2e-3, "spread {} too large", ag.spread);
        assert!(
            ag.reduced_chi2 < 25.0,
            "reduced chi2 {} should be small when all enzymes agree",
            ag.reduced_chi2
        );
    }

    #[test]
    fn enzymes_measuring_different_ani_are_flagged() {
        // Enzymes sampling regions of different conservation. The spread must
        // show up, and the overdispersion must dwarf the agreeing case.
        let agree = enzyme_agreement(&vec![
            synth(0.97, 32, 6, 2, 200_000, "A"),
            synth(0.97, 27, 7, 2, 200_000, "B"),
            synth(0.97, 28, 6, 2, 200_000, "C"),
        ]);
        let disagree = enzyme_agreement(&vec![
            synth(0.99, 32, 6, 2, 200_000, "A"),
            synth(0.97, 27, 7, 2, 200_000, "B"),
            synth(0.95, 28, 6, 2, 200_000, "C"),
        ]);
        assert!(
            disagree.spread > 0.03,
            "expected a visible spread, got {}",
            disagree.spread
        );
        assert!(
            disagree.reduced_chi2 > agree.reduced_chi2 * 100.0,
            "overdispersion should dominate: {} vs {}",
            disagree.reduced_chi2,
            agree.reduced_chi2
        );
    }

    #[test]
    fn a_single_enzyme_has_no_agreement_to_report() {
        let ag = enzyme_agreement(&vec![synth(0.97, 32, 6, 2, 100_000, "only")]);
        assert_eq!(ag.dof, 0);
        assert!(ag.reduced_chi2.is_nan());
        assert_eq!(ag.spread, 0.0);
    }
}

#[cfg(test)]
mod het_tests {
    use super::*;

    /// Counts generated by the heterogeneous model itself, to check recovery.
    fn synth_het(d: f64, alpha: f64, tag_len: usize, site_len: usize, tol: usize, n: u64)
        -> EnzymeStratum {
        let s = stratum("h", tag_len, site_len, vec![0; tol + 1], 0);
        let mut hist = vec![0u64; tol + 1];
        let mut found = 0.0;
        for m in 0..=tol {
            let p = ln_p_found_het(d, alpha, &s, m).exp();
            hist[m] = (p * n as f64).round() as u64;
            found += p;
        }
        let n_miss = ((1.0 - found) * n as f64).round() as u64;
        stratum("h", tag_len, site_len, hist, n_miss)
    }

    #[test]
    fn recovers_mean_identity_under_heterogeneity() {
        for &(d, alpha) in &[(0.02, 0.5), (0.05, 1.0), (0.10, 2.0), (0.02, 100.0)] {
            let truth = het_ani(d, alpha);
            let strata = vec![
                synth_het(d, alpha, 32, 6, 2, 400_000),
                synth_het(d, alpha, 27, 7, 2, 400_000),
            ];
            let out = estimate_heterogeneous(&strata);
            if !out.heterogeneity_supported {
                // alpha=100 is effectively uniform; the gate correctly declines.
                assert!(alpha > 10.0, "declined a genuinely heterogeneous fit: {out:?}");
                continue;
            }
            assert!(
                (out.ani - truth).abs() < 2e-3,
                "d={d} alpha={alpha}: truth {truth:.5}, got {:.5} (shape {:.3})",
                out.ani,
                out.shape
            );
        }
    }

    #[test]
    fn homogeneous_fit_reads_high_under_heterogeneity() {
        // The motivating failure: with strongly heterogeneous rates, assuming a
        // single rate overestimates identity. This is the real-genome bias.
        let (d, alpha) = (0.08, 0.6);
        let truth = het_ani(d, alpha);
        let strata = vec![synth_het(d, alpha, 32, 6, 2, 800_000)];
        let out = estimate_heterogeneous(&strata);
        assert!(out.heterogeneity_supported, "gate should accept: {out:?}");
        assert!(
            out.ani_homogeneous > truth + 0.005,
            "expected homogeneous fit to read high: truth {truth:.4}, homog {:.4}",
            out.ani_homogeneous
        );
        assert!(
            (out.ani - truth).abs() < (out.ani_homogeneous - truth).abs(),
            "heterogeneous fit should be closer: het {:.4} homog {:.4} truth {truth:.4}",
            out.ani,
            out.ani_homogeneous
        );
    }

    #[test]
    fn large_shape_matches_the_homogeneous_limit() {
        // alpha large = nearly uniform rates, so the two models must agree.
        let (d, alpha) = (0.01, 1500.0);
        let strata = vec![synth_het(d, alpha, 32, 6, 2, 500_000)];
        let out = estimate_heterogeneous(&strata);
        assert!(
            (out.ani - out.ani_homogeneous).abs() < 3e-3,
            "het {:.5} vs homog {:.5} (supported={})",
            out.ani,
            out.ani_homogeneous,
            out.heterogeneity_supported
        );
    }
}

#[cfg(test)]
mod deg_tests {
    use super::*;

    /// Deterministic LCG (same multiplier as the chain_ani tests) so the
    /// simulation needs no rand dependency.
    struct Lcg(u64);

    impl Lcg {
        /// Uniform f64 in [0, 1).
        fn f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Forward-simulate tag outcomes under the identity framework the design
    /// note specifies — independent of the likelihood code under test. Per tag,
    /// per position: a mutation happens w.p. `1-a`, the mutant base uniform
    /// over the 3 alternatives. An exact-site mutation kills the tag; a body
    /// mutation is a mismatch; a 2-of-4 (3-of-4) degenerate position preserves
    /// the site — showing up as a mismatch — for 1 (2) of the 3 alternatives
    /// and kills the tag otherwise. A surviving tag with more than `tol`
    /// mismatches is not found and joins the miss count.
    fn simulate_counts(
        ani: f64,
        exact: usize,
        d2: usize,
        d3: usize,
        body: usize,
        tol: usize,
        n: u64,
        seed: u64,
    ) -> (Vec<u64>, u64) {
        let mut rng = Lcg(seed);
        let mut hist = vec![0u64; tol + 1];
        let mut n_miss = 0u64;
        for _ in 0..n {
            let mut m = 0usize;
            let mut alive = true;
            for _ in 0..exact {
                if rng.f64() >= ani {
                    alive = false;
                }
            }
            for _ in 0..d2 {
                if alive && rng.f64() >= ani {
                    if rng.f64() < 1.0 / 3.0 {
                        m += 1;
                    } else {
                        alive = false;
                    }
                }
            }
            for _ in 0..d3 {
                if alive && rng.f64() >= ani {
                    if rng.f64() < 2.0 / 3.0 {
                        m += 1;
                    } else {
                        alive = false;
                    }
                }
            }
            for _ in 0..body {
                if alive && rng.f64() >= ani {
                    m += 1;
                }
            }
            if alive && m <= tol {
                hist[m] += 1;
            } else {
                n_miss += 1;
            }
        }
        (hist, n_miss)
    }

    /// HaeIV: GAY-RTC -> tag 27, exact 4, d2 = 2, body 21.
    /// Hin4I: GAY-VTC -> tag 27, exact 4, d2 = 1, d3 = 1, body 21.
    #[test]
    fn recovers_ani_with_degenerate_site_homogeneous() {
        for &(name, exact, d2, d3) in &[("HaeIV", 4usize, 2usize, 0usize), ("Hin4I", 4, 1, 1)] {
            for &truth in &[0.85, 0.90, 0.95, 0.98] {
                let (hist, n_miss) = simulate_counts(truth, exact, d2, d3, 21, 2, 400_000, 42);
                let good = stratum_deg(name, 27, exact, d2, d3, hist.clone(), n_miss);
                let out = estimate(&[good]);
                assert!(
                    (out.ani - truth).abs() < 2e-3,
                    "{name} truth {truth}, got {} (loss {}, hist {})",
                    out.ani,
                    out.ani_from_loss,
                    out.ani_from_hist
                );

                // The pre-fix geometry treated degenerate positions as exact
                // (site_len = 6). On the same counts it must be further from
                // the truth than the corrected fit.
                let old = stratum(name, 27, exact + d2 + d3, hist, n_miss);
                let out_old = estimate(&[old]);
                let new_err = (out.ani - truth).abs();
                let old_err = (out_old.ani - truth).abs();
                assert!(
                    new_err < old_err,
                    "{name} truth {truth}: corrected err {new_err:.5} should beat old-geometry err {old_err:.5} (old est {})",
                    out_old.ani
                );
            }
        }
    }

    /// Counts generated by the heterogeneous degenerate model itself, to check
    /// recovery of the gamma-mixed NB + survival-factor form.
    fn synth_het_deg(
        d: f64,
        alpha: f64,
        exact: usize,
        d2: usize,
        d3: usize,
        tol: usize,
        n: u64,
    ) -> EnzymeStratum {
        let tag_len = exact + d2 + d3 + 21;
        let s = stratum_deg("h", tag_len, exact, d2, d3, vec![0; tol + 1], 0);
        let mut hist = vec![0u64; tol + 1];
        let mut found = 0.0;
        for m in 0..=tol {
            let p = ln_p_found_het(d, alpha, &s, m).exp();
            hist[m] = (p * n as f64).round() as u64;
            found += p;
        }
        let n_miss = ((1.0 - found) * n as f64).round() as u64;
        stratum_deg("h", tag_len, exact, d2, d3, hist, n_miss)
    }

    #[test]
    fn recovers_mean_identity_with_degenerate_site_het() {
        for &(d, alpha) in &[(0.02, 0.5), (0.05, 1.0), (0.10, 2.0)] {
            let truth = het_ani(d, alpha);
            // HaeIV-like and Hin4I-like strata at the same divergence.
            let strata = vec![
                synth_het_deg(d, alpha, 4, 2, 0, 2, 400_000),
                synth_het_deg(d, alpha, 4, 1, 1, 2, 400_000),
            ];
            let out = estimate_heterogeneous(&strata);
            assert!(out.heterogeneity_supported, "gate should accept: {out:?}");
            assert!(
                (out.ani - truth).abs() < 2e-3,
                "d={d} alpha={alpha}: truth {truth:.5}, got {:.5} (shape {:.3})",
                out.ani,
                out.shape
            );
        }
    }

    #[test]
    fn degenerate_forms_reduce_to_fully_specific_at_zero() {
        // d2 = d3 = 0 must be the same numerical path as before: the het
        // formula's extra term is exactly 0 and b_eff = b exactly.
        let a = stratum("x", 27, 6, vec![100, 20, 3], 50);
        let b = stratum_deg("x", 27, 6, 0, 0, vec![100, 20, 3], 50);
        assert_eq!(a.tag_len, b.tag_len);
        assert_eq!(a.body_len, b.body_len);
        for &d in &[0.01, 0.1] {
            for &alpha in &[0.5, 2.0, 50.0] {
                for m in 0..3 {
                    assert_eq!(
                        ln_p_found_het(d, alpha, &a, m),
                        ln_p_found_het(d, alpha, &b, m),
                        "het form must be bit-identical at d2=d3=0"
                    );
                }
            }
        }
        let out_a = estimate(&[a]);
        let out_b = estimate(&[b]);
        assert_eq!(out_a.ani.to_bits(), out_b.ani.to_bits());
    }
}
