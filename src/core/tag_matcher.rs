use std::collections::HashMap;

use crate::core::tag_extractor::{GenomeTag, TagSet};
use crate::core::synteny_builder::{SyntenyBlock, SyntenyBuilder};

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
}

/// Matches tags between query and reference genomes.
pub struct TagMatcher;

impl TagMatcher {
    /// Match query tags against reference tags and produce a `MatchResult`.
    ///
    /// Builds a hash index of reference tags by sequence, then matches query tags
    /// using Hamming distance for near-match tolerance.
    pub fn match_tag_sets(
        query: &TagSet,
        reference: &TagSet,
        config: &MatchConfig,
    ) -> MatchResult {
        let mut ref_index: HashMap<[u8; 32], Vec<usize>> = HashMap::new();
        for (i, tag) in reference.tags.iter().enumerate() {
            ref_index.entry(tag.sequence).or_default().push(i);
        }

        let mut matched_pairs: Vec<MatchedPair> = Vec::new();
        let mut unmatched_query: Vec<GenomeTag> = Vec::new();
        let mut matched_ref_flags = vec![false; reference.tags.len()];

        for q_tag in &query.tags {
            if let Some(ref_indices) = ref_index.get(&q_tag.sequence) {
                let mut best_idx = None;
                let mut best_dist = usize::MAX;

                for &idx in ref_indices {
                    if matched_ref_flags[idx] {
                        continue;
                    }
                    let r_tag = &reference.tags[idx];
                    let dist = hamming_distance(&q_tag.sequence, q_tag.seq_len, &r_tag.sequence, r_tag.seq_len);
                    if dist < best_dist {
                        best_dist = dist;
                        best_idx = Some(idx);
                    }
                }

                if let Some(idx) = best_idx {
                    let accept = if !config.allow_near_match {
                        best_dist == 0
                    } else {
                        best_dist <= config.near_match_tolerance
                    };
                    if accept {
                        let r_tag = &reference.tags[idx];
                        let tag_len = q_tag.seq_len.max(r_tag.seq_len) as usize;
                        let local_ani =
                            1.0 - (best_dist as f64 / tag_len.max(1) as f64);
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
            } else if config.allow_near_match {
                // Near-match fallback: scan all reference tags for best Hamming match
                let mut best_idx = None;
                let mut best_dist = usize::MAX;
                for (idx, r_tag) in reference.tags.iter().enumerate() {
                    if matched_ref_flags[idx] {
                        continue;
                    }
                    let dist = hamming_distance(&q_tag.sequence, q_tag.seq_len, &r_tag.sequence, r_tag.seq_len);
                    if dist < best_dist {
                        best_dist = dist;
                        best_idx = Some(idx);
                    }
                }
                if let Some(idx) = best_idx {
                    if best_dist <= config.near_match_tolerance {
                        let r_tag = &reference.tags[idx];
                        let tag_len = q_tag.seq_len.max(r_tag.seq_len) as usize;
                        let local_ani =
                            1.0 - (best_dist as f64 / tag_len.max(1) as f64);
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
        }
    }
}

/// Compute Hamming distance between two fixed-length tag sequences.
/// Only compares up to `min(len_a, len_b)` bases.
fn hamming_distance(a: &[u8; 32], len_a: u8, b: &[u8; 32], len_b: u8) -> usize {
    let cmp_len = (len_a as usize).min(len_b as usize);
    a.iter().zip(b.iter()).take(cmp_len).filter(|(x, y)| x != y).count()
}
