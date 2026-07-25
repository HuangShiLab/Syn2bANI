use crate::utils::fxhash::FastHashMap;

use crate::core::tag_extractor::{GenomeTag, TagSet};
use crate::core::synteny_builder::{SyntenyBlock, SyntenyBuilder};
use crate::parallel::simd::diff_count_u64;

/// Configuration controlling tag matching behavior.
#[derive(Debug, Clone)]
pub struct MatchConfig {
    pub allow_near_match: bool,
    pub near_match_tolerance: usize,
}

impl Default for MatchConfig {
    fn default() -> Self {
        Self {
            allow_near_match: true,
            near_match_tolerance: 2,
        }
    }
}

/// A pair of matched query and reference tags.
#[derive(Debug, Clone)]
pub struct MatchedPair {
    pub query_tag: GenomeTag,
    pub ref_tag: GenomeTag,
    pub hamming_distance: usize,
    pub local_ani: f64,
    pub gap_diff: isize,
}

/// Result of matching two tag sets.
#[derive(Debug, Clone)]
pub struct MatchResult {
    pub matched_pairs: Vec<MatchedPair>,
    pub unmatched_query: Vec<GenomeTag>,
    pub unmatched_ref: Vec<GenomeTag>,
    pub synteny_blocks: Vec<SyntenyBlock>,
    pub shared_tag_fraction: f64,
    /// Raw query contig sequences, indexed by `contig_id`.
    /// Carried forward so that chain-interval k-mer ANI can be computed.
    pub q_sequences: Vec<Vec<u8>>,
    /// Raw reference contig sequences, indexed by `contig_id`.
    pub r_sequences: Vec<Vec<u8>>,
}

/// Matches tags between query and reference genomes.
pub struct TagMatcher;

impl TagMatcher {
    /// Match query tags against reference tags and produce a `MatchResult`.
    ///
    /// Builds a hash index of reference tags by packed sequence, then matches query tags
    /// using Hamming distance (via 64-bit XOR + diff-count) for near-match tolerance.
    /// When a reference tag occurs multiple times, the candidate that best extends the
    /// current synteny block (smallest position-consistency penalty) is chosen.
    pub fn match_tag_sets(
        query: &TagSet,
        reference: &TagSet,
        config: &MatchConfig,
    ) -> MatchResult {
        let mut ref_index: FastHashMap<u64, Vec<usize>> = FastHashMap::default();
        for (i, tag) in reference.tags.iter().enumerate() {
            ref_index.entry(tag.packed_sequence).or_default().push(i);
        }

        let mut matched_pairs: Vec<MatchedPair> = Vec::new();
        let mut unmatched_query: Vec<GenomeTag> = Vec::new();
        let mut matched_ref_flags = vec![false; reference.tags.len()];

        for q_tag in &query.tags {
            let last_pair = matched_pairs.last();

            // 1) Exact packed-sequence matches
            let mut best_idx: Option<usize> = None;
            if let Some(ref_indices) = ref_index.get(&q_tag.packed_sequence) {
                best_idx = Self::select_best_candidate(
                    q_tag,
                    ref_indices,
                    reference,
                    &matched_ref_flags,
                    last_pair,
                    0,
                );
            }

            // 2) Near-match fallback
            if best_idx.is_none() && config.allow_near_match {
                let all_indices: Vec<usize> = (0..reference.tags.len()).collect();
                best_idx = Self::select_best_candidate(
                    q_tag,
                    &all_indices,
                    reference,
                    &matched_ref_flags,
                    last_pair,
                    config.near_match_tolerance,
                );
            }

            if let Some(idx) = best_idx {
                let r_tag = &reference.tags[idx];
                let best_dist = hamming_distance(q_tag, r_tag);
                let accept = if !config.allow_near_match {
                    best_dist == 0
                } else {
                    best_dist <= config.near_match_tolerance
                };
                if accept {
                    let tag_len = q_tag.seq_len.max(r_tag.seq_len) as usize;
                    let local_ani = 1.0 - (best_dist as f64 / tag_len.max(1) as f64);
                    matched_pairs.push(MatchedPair {
                        query_tag: q_tag.clone(),
                        ref_tag: r_tag.clone(),
                        hamming_distance: best_dist,
                        local_ani,
                        gap_diff: 0,
                    });
                    matched_ref_flags[idx] = true;
                    continue;
                }
            }

            unmatched_query.push(q_tag.clone());
        }

        // Compute gap_diff for consecutive matched pairs
        for i in 1..matched_pairs.len() {
            let q_gap = matched_pairs[i]
                .query_tag
                .position
                .saturating_sub(matched_pairs[i - 1].query_tag.position);
            let r_gap = matched_pairs[i]
                .ref_tag
                .position
                .saturating_sub(matched_pairs[i - 1].ref_tag.position);
            matched_pairs[i].gap_diff = q_gap as isize - r_gap as isize;
        }

        // Build synteny blocks from matched pairs
        let synteny_blocks = SyntenyBuilder::build_blocks(&matched_pairs);

        let mut unmatched_ref = Vec::new();
        for (i, tag) in reference.tags.iter().enumerate() {
            if !matched_ref_flags[i] {
                unmatched_ref.push(tag.clone());
            }
        }

        let shared_tag_fraction = if query.tags.is_empty() {
            0.0
        } else {
            matched_pairs.len() as f64 / query.tags.len() as f64
        };

        MatchResult {
            matched_pairs,
            unmatched_query,
            unmatched_ref,
            synteny_blocks,
            shared_tag_fraction,
            q_sequences: query.sequences.clone(),
            r_sequences: reference.sequences.clone(),
        }
    }

    /// Choose the best unused reference candidate for a query tag.
    ///
    /// Ties on Hamming distance are broken by a position-consistency penalty that
    /// favours the candidate continuing the current synteny block.  This addresses
    /// repeated tags on the same reference contig: the anchor that yields the longer
    /// shared (collinear) fraction is preferred.
    fn select_best_candidate(
        q_tag: &GenomeTag,
        candidates: &[usize],
        reference: &TagSet,
        matched_ref_flags: &[bool],
        last_pair: Option<&MatchedPair>,
        tolerance: usize,
    ) -> Option<usize> {
        let mut best: Option<(usize, usize, isize)> = None; // (idx, dist, position_penalty)

        for &idx in candidates {
            if matched_ref_flags[idx] {
                continue;
            }
            let r_tag = &reference.tags[idx];
            let dist = hamming_distance(q_tag, r_tag);
            if dist > tolerance {
                continue;
            }

            let penalty = if let Some(lp) = last_pair {
                if lp.query_tag.contig_id == q_tag.contig_id
                    && lp.ref_tag.contig_id == r_tag.contig_id
                {
                    let expected_ref_pos = if q_tag.position >= lp.query_tag.position {
                        lp.ref_tag.position.saturating_add(q_tag.position - lp.query_tag.position)
                    } else {
                        lp.ref_tag.position.saturating_sub(lp.query_tag.position - q_tag.position)
                    };
                    (r_tag.position as isize - expected_ref_pos as isize).abs()
                } else {
                    isize::MAX
                }
            } else {
                0
            };

            let better = match best {
                None => true,
                Some((_, best_dist, best_penalty)) => {
                    if dist != best_dist {
                        dist < best_dist
                    } else if penalty != best_penalty {
                        penalty < best_penalty
                    } else {
                        false
                    }
                }
            };

            if better {
                best = Some((idx, dist, penalty));
            }
        }

        best.map(|(idx, _, _)| idx)
    }
}

/// Compute Hamming distance between two `GenomeTag`s using 64-bit packed sequences.
///
/// Only compares up to `min(seq_len_a, seq_len_b)` bases.
/// Uses XOR + diff-count (one CPU instruction per comparison, cross-platform).
#[inline]
fn hamming_distance(a: &GenomeTag, b: &GenomeTag) -> usize {
    let cmp_len = (a.seq_len as usize).min(b.seq_len as usize);
    let xor = a.packed_sequence ^ b.packed_sequence;
    let mask = if cmp_len >= 32 {
        u64::MAX
    } else {
        (1u64 << (cmp_len * 2)) - 1
    };
    diff_count_u64(xor & mask) as usize
}
