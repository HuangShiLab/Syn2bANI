//! Structural variation calls derived from the chain-restricted collinear
//! chains of [`crate::core::chain_ani`].
//!
//! This replaces the v7 `StructureAnalyzer` path: the chains come from the
//! gap-penalised, divergence-adaptive DP, so a call here is backed by the same
//! anchor evidence the ANI estimate trusts — not by the old sparse_chain
//! output that re-admitted non-collinear anchors.
//!
//! # Call rules
//!
//! All coordinates are contig-local, half-open, and refer to the chain spans
//! and anchor positions exported in [`ChainBlock`].
//!
//! - **Inversion**: any reverse-orientation chain. Its span is the inverted
//!   segment; the flanking chains on the same query contig supply the phasing
//!   evidence (`support_left` / `support_right`).
//! - **Translocation**: two chains adjacent along a query contig that map to
//!   different reference contigs, or that map to the same reference contig in
//!   an order inconsistent with their query order (sign-aware).
//! - **Insertion / Deletion**, two sources:
//!   (a) *within a chain*: consecutive anchors whose query and reference
//!   offsets disagree by more than `indel_min` bp. The chaining gap penalty
//!   absorbs small indels silently, so only jumps past the threshold surface
//!   here. `Insertion` = the query carries the extra bases.
//!   (b) *between chains*: two adjacent same-contig, same-orientation,
//!   collinear chains bracket an accessory segment; the excess of the query
//!   gap over the reference gap (or vice versa) past `indel_min` is called.
//!   The reference interval between the two anchors must be empty of every
//!   other chain's anchors: if a third chain maps inside it, the junction is
//!   a relocation/order artifact, and reporting `q_gap - r_gap` would invent
//!   a deletion the size of the reference jump (observed on K-12 vs Sakai:
//!   a 12 kb relocated block produced a spurious 961 kb "deletion").
//!
//! All rules operate on the non-shadowed chains only
//! ([`shadowed_mask`]): chains fully contained in a denser chain's query
//! span are secondary mappings of multi-copy repeats (rRNA operons and the
//! like), not rearrangement evidence. Without this filter every such repeat
//! is miscalled as an inversion — including for a genome compared to
//! itself — and each one breaks the query-order adjacency of the backbone
//! chains it sits between, silently suppressing between-chain indel calls.

use crate::core::chain_ani::ChainBlock;

/// Type of structural variation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SvType {
    Inversion,
    Insertion,
    Deletion,
    Translocation,
    Duplication,
}
use crate::utils::fxhash::FastHashMap;

/// One structural variation call between a query and a reference genome.
#[derive(Debug, Clone)]
pub struct SvCall {
    pub sv_type: SvType,
    /// Contig indices into the caller's FASTA record order.
    pub q_contig: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub r_contig: usize,
    pub r_start: usize,
    pub r_end: usize,
    /// Event size in bp (span for inversions, length change for indels,
    /// query-side gap for translocation junctions).
    pub size: usize,
    /// Anchors in the chain flanking the event on the left in query order
    /// (for within-chain indels, the anchors up to and including the left
    /// one). This is the tag-phasing evidence across the breakpoint.
    pub support_left: usize,
    /// Anchors in the chain flanking the event on the right.
    pub support_right: usize,
}

/// Signed reference step between two anchors of a chain, in query order.
/// Positive means collinear progression.
fn r_step(orientation: char, r_from: usize, r_to: usize) -> i64 {
    if orientation == '-' {
        r_from as i64 - r_to as i64
    } else {
        r_to as i64 - r_from as i64
    }
}

/// Mark chains whose entire query span is contained in another chain's span
/// on the same query contig while that chain carries strictly more anchors.
///
/// These are secondary mappings of multi-copy repeats (e.g. rRNA operons):
/// the repeat region chains once into the collinear backbone and once more
/// against the other repeat copy, producing a sparse nested chain. Such a
/// chain is not rearrangement evidence — a genuine inversion breaks the
/// backbone, so nothing spans it — yet it pollutes every downstream rule:
/// it is miscalled as an inversion, and sitting between two backbone chains
/// in query order it breaks their adjacency, which silently suppresses
/// between-chain indel calls (observed on H. pylori 26695 vs a 36 kb cagPAI
/// deletion: the two rRNA shadow chains hid a clean 36,154 bp insertion).
/// Ties in anchor count are broken by chain index so that duplicate chains
/// with identical spans shadow exactly one of the pair.
fn shadowed_mask(chains: &[ChainBlock]) -> Vec<bool> {
    let n = chains.len();
    let mut shadowed = vec![false; n];
    for i in 0..n {
        let c = &chains[i];
        for (j, p) in chains.iter().enumerate() {
            if i == j || p.q_contig != c.q_contig {
                continue;
            }
            let contains = p.q_start <= c.q_start && p.q_end >= c.q_end;
            let dominates =
                p.n_anchors > c.n_anchors || (p.n_anchors == c.n_anchors && j < i);
            if contains && dominates {
                shadowed[i] = true;
                break;
            }
        }
    }
    shadowed
}

/// Detect structural variations from the final adaptive-pass chains.
///
/// `indel_min` is the minimum offset disagreement (bp) reported as an indel;
/// smaller jumps are indistinguishable from the gap arithmetic the chaining
/// DP already tolerates.
pub fn detect(chains: &[ChainBlock], indel_min: usize) -> Vec<SvCall> {
    let mut calls = Vec::new();
    let shadowed = shadowed_mask(chains);

    // (a) Within-chain indels: anchor-to-anchor offset jumps. Shadow chains
    // are skipped: their anchor offsets describe a repeat's secondary
    // mapping, not the comparison's backbone.
    for (c, &sh) in chains.iter().zip(&shadowed) {
        if sh {
            continue;
        }
        for i in 1..c.anchors.len() {
            let (q0, r0) = c.anchors[i - 1];
            let (q1, r1) = c.anchors[i];
            let dq = q1 as i64 - q0 as i64;
            let dr = r_step(c.orientation, r0, r1);
            let diff = dq - dr;
            let sv_type = if diff > indel_min as i64 {
                SvType::Insertion
            } else if diff < -(indel_min as i64) {
                SvType::Deletion
            } else {
                continue;
            };
            calls.push(SvCall {
                sv_type,
                q_contig: c.q_contig,
                q_start: q0,
                q_end: q1,
                r_contig: c.r_contig,
                r_start: r0.min(r1),
                r_end: r0.max(r1),
                size: diff.unsigned_abs() as usize,
                support_left: i,
                support_right: c.anchors.len() - i,
            });
        }
    }

    // Per query contig, order non-shadowed chains along the query and
    // inspect adjacencies.
    let mut by_q: FastHashMap<usize, Vec<usize>> = FastHashMap::default();
    for (i, c) in chains.iter().enumerate() {
        if shadowed[i] {
            continue;
        }
        by_q.entry(c.q_contig).or_default().push(i);
    }

    // Anchor reference positions per reference contig, for the relocation
    // guard in branch (b): a between-chain indel is only credible when no
    // third chain anchors inside the reference interval the junction spans.
    let mut r_index: FastHashMap<usize, Vec<(usize, usize)>> = FastHashMap::default();
    for (i, c) in chains.iter().enumerate() {
        if shadowed[i] {
            continue;
        }
        for &(_, r) in &c.anchors {
            r_index.entry(c.r_contig).or_default().push((r, i));
        }
    }
    for v in r_index.values_mut() {
        v.sort_unstable();
    }
    // True when some chain other than `excl` anchors strictly inside (lo, hi).
    let interval_occupied = |r_contig: usize, lo: usize, hi: usize, excl: [usize; 2]| -> bool {
        r_index.get(&r_contig).is_some_and(|v| {
            let start = v.partition_point(|&(r, _)| r <= lo);
            v[start..]
                .iter()
                .take_while(|&&(r, _)| r < hi)
                .any(|&(_, i)| i != excl[0] && i != excl[1])
        })
    };

    for idxs in by_q.values_mut() {
        idxs.sort_by_key(|&i| (chains[i].q_start, chains[i].q_end));

        // Every reverse-orientation chain is an inversion call.
        for (pos, &i) in idxs.iter().enumerate() {
            let c = &chains[i];
            if c.orientation != '-' {
                continue;
            }
            calls.push(SvCall {
                sv_type: SvType::Inversion,
                q_contig: c.q_contig,
                q_start: c.q_start,
                q_end: c.q_end,
                r_contig: c.r_contig,
                r_start: c.r_start,
                r_end: c.r_end,
                size: c.q_end - c.q_start,
                support_left: pos
                    .checked_sub(1)
                    .map(|p| chains[idxs[p]].n_anchors)
                    .unwrap_or(0),
                support_right: idxs
                    .get(pos + 1)
                    .map(|&j| chains[j].n_anchors)
                    .unwrap_or(0),
            });
        }

        // Adjacent-chain junctions.
        for w in idxs.windows(2) {
            let (a, b) = (&chains[w[0]], &chains[w[1]]);
            let (a_q, a_r) = *a.anchors.last().unwrap();
            let (b_q, b_r) = *b.anchors.first().unwrap();

            if a.r_contig != b.r_contig {
                // Query-adjacent chains landing on different reference contigs:
                // a translocation (or contig-break) junction. The query gap
                // between them, if any, is accessory sequence at the junction.
                calls.push(SvCall {
                    sv_type: SvType::Translocation,
                    q_contig: a.q_contig,
                    q_start: a_q,
                    q_end: b_q,
                    r_contig: b.r_contig,
                    r_start: b_r,
                    r_end: b_r,
                    size: b_q.saturating_sub(a_q),
                    support_left: a.n_anchors,
                    support_right: b.n_anchors,
                });
                continue;
            }

            if a.orientation != b.orientation {
                // An inversion boundary; the reverse chain itself carries the
                // Inversion call, so there is nothing to add here.
                continue;
            }

            let orient = a.orientation;
            if r_step(orient, a_r, b_r) < 0 {
                // Same reference contig but the reference order contradicts the
                // query order: an intra-contig rearrangement junction.
                calls.push(SvCall {
                    sv_type: SvType::Translocation,
                    q_contig: a.q_contig,
                    q_start: a_q,
                    q_end: b_q,
                    r_contig: a.r_contig,
                    r_start: a_r.min(b_r),
                    r_end: a_r.max(b_r),
                    size: b_q.saturating_sub(a_q),
                    support_left: a.n_anchors,
                    support_right: b.n_anchors,
                });
                continue;
            }

            // (b) Collinear neighbours: the inter-chain segment is accessory.
            // Compare its length on the two genomes — but only when no third
            // chain anchors inside the reference interval. If one does, the
            // junction is a relocation, not an indel, and the offset
            // difference measures the jump distance, not an indel size.
            let q_gap = b_q as i64 - a_q as i64;
            let r_gap = r_step(orient, a_r, b_r);
            let (lo, hi) = (a_r.min(b_r), a_r.max(b_r));
            if interval_occupied(a.r_contig, lo, hi, [w[0], w[1]]) {
                continue;
            }
            let diff = q_gap - r_gap;
            let sv_type = if diff > indel_min as i64 {
                SvType::Insertion
            } else if diff < -(indel_min as i64) {
                SvType::Deletion
            } else {
                continue;
            };
            calls.push(SvCall {
                sv_type,
                q_contig: a.q_contig,
                q_start: a_q,
                q_end: b_q,
                r_contig: a.r_contig,
                r_start: a_r.min(b_r),
                r_end: a_r.max(b_r),
                size: diff.unsigned_abs() as usize,
                support_left: a.n_anchors,
                support_right: b.n_anchors,
            });
        }
    }

    calls.sort_by_key(|s| (s.q_contig, s.q_start, s.q_end));
    calls
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A chain with evenly spaced anchors. `r_positions` are taken in query
    /// order, so a reverse chain lists them descending.
    fn block(
        q_contig: usize,
        r_contig: usize,
        orientation: char,
        q_start: usize,
        r_positions: &[usize],
        spacing: usize,
    ) -> ChainBlock {
        let anchors: Vec<(usize, usize)> = r_positions
            .iter()
            .enumerate()
            .map(|(i, &r)| (q_start + i * spacing, r))
            .collect();
        let q_end = q_start + (r_positions.len() - 1) * spacing;
        ChainBlock {
            q_contig,
            r_contig,
            orientation,
            q_start,
            q_end,
            r_start: *r_positions.iter().min().unwrap(),
            r_end: *r_positions.iter().max().unwrap(),
            n_anchors: anchors.len(),
            anchors,
        }
    }

    fn forward(n: usize, q0: usize, r0: usize, spacing: usize) -> ChainBlock {
        let rpos: Vec<usize> = (0..n).map(|i| r0 + i * spacing).collect();
        block(0, 0, '+', q0, &rpos, spacing)
    }

    #[test]
    fn inversion_chain_is_called() {
        let fwd = forward(10, 0, 0, 1_000);
        let inv = block(
            0,
            0,
            '-',
            10_000,
            &[19_000, 18_000, 17_000, 16_000, 15_000, 14_000],
            1_000,
        );
        let calls = detect(&[fwd, inv], 1_000);
        let invs: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Inversion)
            .collect();
        assert_eq!(invs.len(), 1, "{calls:?}");
        assert_eq!(invs[0].q_start, 10_000);
        assert_eq!(invs[0].q_end, 15_000);
        assert_eq!(invs[0].support_left, 10);
        assert_eq!(invs[0].support_right, 0);
    }

    #[test]
    fn translocation_across_reference_contigs() {
        let a = forward(10, 0, 0, 1_000);
        // Same query contig, but the next chain maps to reference contig 1.
        let rpos: Vec<usize> = (0..8).map(|i| i * 1_000).collect();
        let b = block(0, 1, '+', 10_000, &rpos, 1_000);
        let calls = detect(&[a, b], 1_000);
        let t: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Translocation)
            .collect();
        assert_eq!(t.len(), 1, "{calls:?}");
        assert_eq!(t[0].r_contig, 1);
        assert_eq!(t[0].q_start, 9_000);
        assert_eq!(t[0].q_end, 10_000);
    }

    #[test]
    fn intra_contig_order_reversal_is_translocation() {
        let a = forward(10, 0, 0, 1_000);
        // Same reference contig, but this chain sits *before* `a` on the
        // reference while coming after it on the query: order contradiction.
        let rpos: Vec<usize> = (0..10).map(|i| 5_000 + i * 100).collect();
        let b = block(0, 0, '+', 20_000, &rpos, 1_000);
        let calls = detect(&[a, b], 1_000);
        assert!(
            calls.iter().any(|c| c.sv_type == SvType::Translocation),
            "{calls:?}"
        );
    }

    #[test]
    fn within_chain_deletion_and_insertion() {
        // Deletion: reference offset jumps 40 kb between adjacent anchors.
        let mut rpos: Vec<usize> = (0..10).map(|i| i * 1_000).collect();
        for r in rpos.iter_mut().skip(5) {
            *r += 40_000;
        }
        let del = block(0, 0, '+', 0, &rpos, 1_000);
        let calls = detect(std::slice::from_ref(&del), 1_000);
        let d: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Deletion)
            .collect();
        assert_eq!(d.len(), 1, "{calls:?}");
        assert_eq!(d[0].size, 40_000);
        assert_eq!(d[0].support_left, 5);
        assert_eq!(d[0].support_right, 5);

        // Insertion: query gains 12 kb the reference does not have. Build the
        // anchors directly since the helper uses fixed spacing.
        let mut c = forward(10, 0, 0, 1_000);
        for a in c.anchors.iter_mut().skip(5) {
            a.0 += 12_000;
        }
        c.q_end += 12_000;
        let calls = detect(&[c], 1_000);
        let ins: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Insertion)
            .collect();
        assert_eq!(ins.len(), 1, "{calls:?}");
        assert_eq!(ins[0].size, 12_000);
    }

    #[test]
    fn small_offset_jumps_stay_silent() {
        // 400 bp disagreement is under a 1 kb threshold: the chaining gap
        // penalty already absorbs it, so reporting it would be noise.
        let mut rpos: Vec<usize> = (0..10).map(|i| i * 1_000).collect();
        for r in rpos.iter_mut().skip(5) {
            *r += 400;
        }
        let c = block(0, 0, '+', 0, &rpos, 1_000);
        assert!(detect(&[c], 1_000).is_empty());
    }

    #[test]
    fn accessory_segment_between_chains() {
        // Two collinear chains bracket a 30 kb query segment that has only
        // 5 kb of reference between the flanking anchors: a 25 kb insertion.
        let a = forward(10, 0, 0, 1_000);
        let rpos: Vec<usize> = (0..8).map(|i| 14_000 + i * 1_000).collect();
        let b = block(0, 0, '+', 39_000, &rpos, 1_000);
        let calls = detect(&[a, b], 1_000);
        let ins: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Insertion)
            .collect();
        assert_eq!(ins.len(), 1, "{calls:?}");
        // q gap = 39_000 - 9_000 = 30_000; r gap = 14_000 - 9_000 = 5_000.
        assert_eq!(ins[0].size, 25_000);
        assert_eq!(ins[0].support_left, 10);
        assert_eq!(ins[0].support_right, 8);
    }

    #[test]
    fn identical_genomes_call_nothing() {
        let c = forward(20, 0, 0, 1_000);
        assert!(detect(&[c], 1_000).is_empty());
    }

    #[test]
    fn shadowed_repeat_chain_is_not_an_inversion() {
        // Backbone chain spans 0-50 kb; a sparse reverse chain nested inside
        // it (a repeat's secondary mapping, rRNA-operon style) must not be
        // called an inversion.
        let backbone = forward(51, 0, 0, 1_000);
        let rpos: Vec<usize> = (0..6).map(|i| 25_000 - i * 1_000).collect();
        let repeat = block(0, 0, '-', 20_000, &rpos, 1_000);
        let calls = detect(&[backbone, repeat], 1_000);
        assert!(
            !calls.iter().any(|c| c.sv_type == SvType::Inversion),
            "nested sparse reverse chain is a shadow, not an inversion: {calls:?}"
        );
    }

    #[test]
    fn shadow_chain_does_not_suppress_between_chain_indel() {
        // H. pylori 26695 vs a 36 kb cagPAI deletion: two backbone chains
        // bracket the deletion, and a shadow repeat chain nested in the left
        // one sits between them in query order. The indel must still be
        // called through the shadow.
        let left = forward(10, 0, 0, 1_000); // q/r 0-9 kb
        // Right chain: query resumes 36 kb past the deletion point, the
        // reference continues collinearly right after the left chain.
        let rpos: Vec<usize> = (0..8).map(|i| 10_000 + i * 1_000).collect();
        let right = block(0, 0, '+', 46_000, &rpos, 1_000);
        // Shadow: sparse reverse chain nested inside the left backbone.
        let rpos_s: Vec<usize> = (0..4).map(|i| 7_000 - i * 1_000).collect();
        let shadow = block(0, 0, '-', 4_000, &rpos_s, 1_000);
        let calls = detect(&[left, right, shadow], 1_000);
        let ins: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Insertion)
            .collect();
        assert_eq!(ins.len(), 1, "{calls:?}");
        // q gap = 46_000 - 9_000 = 37_000; r gap = 10_000 - 9_000 = 1_000.
        assert_eq!(ins[0].size, 36_000);
    }

    #[test]
    fn shadow_chain_does_not_block_translocation_junction() {
        // A shadow chain nested in chain A must not break the A->B
        // cross-contig translocation junction.
        let a = forward(10, 0, 0, 1_000);
        let rpos_b: Vec<usize> = (0..8).map(|i| i * 1_000).collect();
        let b = block(0, 1, '+', 10_000, &rpos_b, 1_000);
        let rpos_s: Vec<usize> = (0..4).map(|i| 7_000 - i * 1_000).collect();
        let shadow = block(0, 0, '-', 4_000, &rpos_s, 1_000);
        let calls = detect(&[a, b, shadow], 1_000);
        let t: Vec<_> = calls
            .iter()
            .filter(|c| c.sv_type == SvType::Translocation)
            .collect();
        assert_eq!(t.len(), 1, "{calls:?}");
    }

    #[test]
    fn relocation_junction_is_not_a_giant_indel() {
        // A small query block relocated far away on the reference: chain B
        // (q 10-18 kb) maps to r 100-108 kb, while chain C covers the
        // reference interval between A and B. The A->B junction must NOT be
        // called a ~91 kb deletion; B->C is the translocation junction.
        let a = forward(10, 0, 0, 1_000);
        let rpos_b: Vec<usize> = (0..9).map(|i| 100_000 + i * 1_000).collect();
        let b = block(0, 0, '+', 10_000, &rpos_b, 1_000);
        let c = forward(10, 20_000, 10_000, 1_000);
        let calls = detect(&[a, b, c], 1_000);
        assert!(
            !calls
                .iter()
                .any(|s| matches!(s.sv_type, SvType::Insertion | SvType::Deletion)),
            "relocation junction must not produce an indel: {calls:?}"
        );
        assert!(
            calls.iter().any(|s| s.sv_type == SvType::Translocation),
            "the order-contradiction junction is still reported: {calls:?}"
        );
    }
}
