//! Database-scale pre-screen: cheap per-enzyme tag containment used by `dist`,
//! `search` and `triangle` to reject pairs that are certainly below the ANI
//! detection floor before paying for [`crate::core::chain_ani`].
//!
//! The screen is recall-first: a false rejection is a silently lost relative,
//! while a false pass only costs one `chain_ani::compute` call, which itself
//! reports the pair as `BELOW_DETECTION`. It therefore does NOT require exact
//! full-tag matches — at 80% ANI a 32 bp tag survives exact matching with
//! probability ~0.8^32 ≈ 1e-3, so a full-tag screen sits on a knife edge
//! exactly at the boundary where it must not lose pairs.
//!
//! Instead each tag contributes ONE key: the strand-canonical centred
//! [`SCREEN_WINDOW`]-bp window of its packed sequence. Window survival is
//! ~0.8^11 ≈ 0.086 at 80% ANI, two orders of magnitude above the full-tag
//! rate, so a simple shared-key count has a wide, statistically comfortable
//! margin between related (≥80%) and unrelated pairs. The window is centred
//! on the tag body, away from the (constant) recognition-site anchors, and
//! keys are grouped per enzyme so tags from different enzymes never collide.

use crate::core::tag_extractor::{revcomp_packed, GenomeTag};

/// Default width of the centred tag window used for screen keys.
///
/// Calibration (500 GTDB pairs with validated ANI 80-100 vs 250k random
/// GTDB pairs; results/db_scale/DB_REWRITE_VALIDATION.md in the paper repo):
/// short windows are useless as a screen — generic bacterial background
/// homology gives ~0.5 key containment for almost any pair at 11 bp and
/// still ~0.15 at 13 bp. At 18 bp the bulk of random pairs shares 0-2 keys
/// while every calibrated true pair kept >= 6 (80-85% band: >= 29).
pub const SCREEN_WINDOW: usize = 18;

/// Default minimum shared-key count for a pair to pass the screen (AND-ed
/// with [`DEFAULT_MIN_CONTAINMENT`]).
///
/// Calibration: at the default window every one of 500 validated >=80% ANI
/// GTDB pairs had >= 6 shared keys (80-85% band: >= 29), while 83% of random
/// GTDB pairs had fewer than 3 or under the containment floor. Small or
/// heavily fragmented genomes produce few tags, so the absolute count alone
/// cannot carry the gate — hence the containment conjunction. The floor is
/// deliberately loose: at n=2000 GTDB scale the measured false-reject rate
/// on estimator-reportable pairs is ~0.1%, and the residual misses are
/// kilobase-scale shared-island pairs (AF ~0.001) that skani's AF>=15 rule
/// also never reports.
pub const DEFAULT_MIN_SHARED: usize = 3;

/// Default containment floor: shared keys as a fraction of the smaller key
/// set must reach this value AND the absolute count must reach
/// [`DEFAULT_MIN_SHARED`]. Protects small-genome pairs whose absolute counts
/// are structurally low (a 0.5 Mb genome has a few hundred panel tags).
pub const DEFAULT_MIN_CONTAINMENT: f64 = 0.001;

/// Strand-canonical screen key for one tag.
///
/// The window is extracted at the same centred offset from both the forward
/// packing and its reverse complement and the smaller value is taken, so a
/// tag read from either strand — inverted segments, arbitrarily oriented
/// draft contigs — yields the same key.
#[inline]
pub fn screen_key(packed: u64, seq_len: u8, window: usize) -> u64 {
    let len = (seq_len as usize).min(32);
    let window = window.min(32);
    if len <= window {
        return canonical_window(packed, revcomp_packed(packed, seq_len), 0, len);
    }
    let lo = (len - window) / 2;
    let rc = revcomp_packed(packed, seq_len);
    canonical_window(packed, rc, lo, window)
}

#[inline]
fn canonical_window(fwd: u64, rev: u64, lo: usize, width: usize) -> u64 {
    let mask = if width >= 32 {
        u64::MAX
    } else {
        (1u64 << (width * 2)) - 1
    };
    let a = (fwd >> (lo * 2)) & mask;
    let b = (rev >> (lo * 2)) & mask;
    a.min(b)
}

/// Per-enzyme screen keys for a digested genome: sorted, deduplicated.
///
/// Deduplication keeps repeat-expanded tag families (transposons, rRNA
/// operons) from inflating the shared count; what the screen measures is the
/// number of DISTINCT conserved loci.
pub fn keys_per_enzyme(tags: &[GenomeTag], window: usize) -> Vec<(String, Vec<u64>)> {
    let mut by_enzyme: Vec<(String, Vec<u64>)> = Vec::new();
    for t in tags {
        let key = screen_key(t.packed_sequence, t.seq_len, window);
        // Tag lists are sorted by (contig, position), not enzyme, and the
        // panel is small, so a linear probe is cheaper than a hash map.
        match by_enzyme.iter_mut().find(|(name, _)| *name == t.enzyme) {
            Some((_, keys)) => keys.push(key),
            None => by_enzyme.push((t.enzyme.clone(), vec![key])),
        }
    }
    for (_, keys) in &mut by_enzyme {
        keys.sort_unstable();
        keys.dedup();
    }
    by_enzyme
}

/// Shared-key count and the combined smaller-set size between two genomes.
///
/// Only enzymes present in BOTH genomes contribute; a panel mismatch (e.g. a
/// BcgI-only legacy sketch against a 4-enzyme query) shrinks the denominator
/// instead of manufacturing zeros.
pub fn shared_keys(q: &[(String, Vec<u64>)], r: &[(String, Vec<u64>)]) -> (usize, usize) {
    let mut shared = 0usize;
    let mut min_total = 0usize;
    for (name, qk) in q {
        let Some((_, rk)) = r.iter().find(|(rn, _)| rn == name) else {
            continue;
        };
        min_total += qk.len().min(rk.len());
        shared += merge_intersection(qk, rk);
    }
    (shared, min_total)
}

#[inline]
fn merge_intersection(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// Screening configuration.
#[derive(Debug, Clone, Copy)]
pub struct ScreenConfig {
    /// Shared-key count floor...
    pub min_shared: usize,
    /// ...AND containment floor (fraction of the smaller key set); both must
    /// be met. Two conditions because neither separates on its own: the
    /// absolute count fails small/fragmented genomes, the fraction alone
    /// admits too much same-family background.
    pub min_containment: f64,
    /// Centred window width (bp) each tag is reduced to.
    pub window: usize,
}

impl Default for ScreenConfig {
    fn default() -> Self {
        Self {
            min_shared: DEFAULT_MIN_SHARED,
            min_containment: DEFAULT_MIN_CONTAINMENT,
            window: SCREEN_WINDOW,
        }
    }
}

impl ScreenConfig {
    /// The recall-first gate. `true` means "cannot rule out ≥80% ANI — refine
    /// with the MLE estimator"; `false` means "certainly below the floor".
    pub fn passes(&self, shared: usize, min_total: usize) -> bool {
        shared >= self.min_shared
            && min_total > 0
            && (shared as f64) >= self.min_containment * min_total as f64
    }
}

/// Crude ANI implied by exact-window containment: the fraction of the smaller
/// genome's tags whose window survives equals ~a^WINDOW under a uniform
/// divergence model, so a ≈ c^(1/WINDOW).
///
/// This is ONLY used as an optional second-tier gate (`--refine-min-approx`)
/// to bound `chain_ani` calls on huge all-vs-all runs; it is never reported
/// as an ANI estimate.
pub fn approx_ani(shared: usize, min_total: usize, window: usize) -> f64 {
    if min_total == 0 {
        return 0.0;
    }
    (shared as f64 / min_total as f64).powf(1.0 / window.max(1) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::tag_extractor::pack_bytes;

    fn packed_from(seq: &str) -> (u64, u8) {
        let mut buf = [0u8; 32];
        buf[..seq.len()].copy_from_slice(seq.as_bytes());
        (pack_bytes(&buf, seq.len() as u8), seq.len() as u8)
    }

    fn revcomp_str(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            })
            .collect()
    }

    #[test]
    fn key_is_strand_invariant() {
        for seq in [
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTTTGCAGCACAAAAACGTACGTACGTACG",
            "AAAACCCCGGGGTTTTAAAACCCC",
        ] {
            let (f, l) = packed_from(&seq);
            let rc = revcomp_str(&seq);
            let (r, rl) = packed_from(&rc);
            assert_eq!(l, rl);
            assert_eq!(
                screen_key(f, l, 11),
                screen_key(r, rl, 11),
                "strand mismatch for {seq}"
            );
        }
    }

    #[test]
    fn key_changes_with_window_base() {
        // Flipping a base inside the centred window must change the key.
        let (a, l) = packed_from("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let (b, _) = packed_from("ACGTACGTACGTTCGTACGTACGTACGTACGT");
        assert_ne!(screen_key(a, l, 11), screen_key(b, l, 11));
    }

    #[test]
    fn key_ignores_bases_outside_window() {
        // Flipping the FIRST base (outside the centred 11 bp window of a 32 bp
        // tag) must NOT change the key — this looseness is what gives the
        // screen its recall margin.
        let (a, l) = packed_from("ACGTACGTACGTACGTACGTACGTACGTACGT");
        let (b, _) = packed_from("TCGTACGTACGTACGTACGTACGTACGTACGT");
        assert_eq!(screen_key(a, l, 11), screen_key(b, l, 11));
    }

    #[test]
    fn identical_genomes_share_everything() {
        let tags: Vec<GenomeTag> = ["ACGTACGTACGTACGTACGTACGTACGTACGT", "TTGCATGCATGCATGCATGCATGC"]
            .iter()
            .map(|s| {
                let (p, l) = packed_from(s);
                GenomeTag {
                    position: 0,
                    contig_id: 0,
                    sequence: [0u8; 32],
                    packed_sequence: p,
                    seq_len: l,
                    direction: '+',
                    enzyme: "BcgI".to_string(),
                }
            })
            .collect();
        let keys = keys_per_enzyme(&tags, 11);
        let (shared, min_total) = shared_keys(&keys, &keys);
        assert_eq!(shared, min_total);
        assert_eq!(shared, 2);
    }

    #[test]
    fn enzymes_never_collide() {
        let mk = |enzyme: &str| {
            let (p, l) = packed_from("ACGTACGTACGTACGTACGTACGTACGTACGT");
            vec![GenomeTag {
                position: 0,
                contig_id: 0,
                sequence: [0u8; 32],
                packed_sequence: p,
                seq_len: l,
                direction: '+',
                enzyme: enzyme.to_string(),
            }]
        };
        let q = keys_per_enzyme(&mk("BcgI"), 11);
        let r = keys_per_enzyme(&mk("AloI"), 11);
        // Identical sequence, different enzyme: no shared keys, and the
        // denominator only covers enzymes present in both genomes.
        assert_eq!(shared_keys(&q, &r), (0, 0));
    }

    #[test]
    fn threshold_logic() {
        let cfg = ScreenConfig::default(); // shared >= 3 AND cont >= 0.001
        assert!(cfg.passes(3, 1_000));
        assert!(cfg.passes(300, 10_000));
        assert!(!cfg.passes(2, 10_000));
        // Both conditions required: a big absolute count on a huge key set
        // still fails when the fraction is under the floor...
        assert!(!cfg.passes(10, 100_000));
        // ...and a high fraction with too few keys fails too.
        assert!(!cfg.passes(2, 100));
        assert!(!cfg.passes(100, 0));
    }

    #[test]
    fn approx_ani_is_monotone() {
        let a80 = approx_ani(86, 1_000, 11);
        let a90 = approx_ani(310, 1_000, 11);
        assert!((a80 - 0.80).abs() < 1e-3);
        assert!(a90 > a80);
        assert_eq!(approx_ani(5, 0, 11), 0.0);
    }
}
