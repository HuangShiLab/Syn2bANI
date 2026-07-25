use crate::core::tag_matcher::MatchResult;
use crate::core::gbrt;
use std::collections::HashMap;

/// Weighting strategy for ANI calculation.
#[derive(Debug, Clone, Copy)]
pub enum WeightStrategy {
    /// Uniform weight (1.0) for all matched pairs.
    Uniform,
    /// Higher weight for tags inside long synteny blocks (sqrt of block length).
    Synteny,
    /// Weight by normalized position in the genome.
    Position,
    /// Weight penalized by gap differences between consecutive tags.
    GapAdjusted,
}

/// Configuration for ANI calculation.
#[derive(Debug, Clone)]
pub struct AniConfig {
    pub weight_strategy: WeightStrategy,
    pub min_shared_tags: usize,
    pub min_af: f64,
    pub debias: bool,
    /// Use the embedded GBRT model for debiasing instead of the simple polynomial correction.
    pub use_gbrt_debias: bool,
    /// Use GBRT v3 (trained on GTDB-R207) instead of v2.
    pub use_gbrt_v3: bool,
    /// Use GBRT v3.6 (trained on 622 pairs, 83-100% ANI) instead of earlier versions.
    pub use_gbrt_v3_6: bool,
    /// Use GBRT v4 (clean inference-time features: raw_ani, shared_log, af_q, af_r).
    pub use_gbrt_v4: bool,
    /// Use GBRT v7 (includes mash_ani and chained_kmer_ani features).
    pub use_gbrt_v7: bool,
    /// Report mash_ani as the final ANI instead of the debiased raw ANI.
    pub use_mash_ani: bool,
    /// Empirical calibration offset applied to mash_ani before reporting.
    /// A small negative offset corrects the systematic overestimation observed
    /// for multi-enzyme exact matching on mid-ANI pairs.
    pub mash_calibration_offset: f64,
    /// Use chain-interval k-mer ANI when `use_mash_ani` is true.
    /// This recomputes ANI inside each sparse-chaining synteny block using
    /// short canonical k-mers, recovering SNP signal lost by exact tag matching.
    pub use_chained_kmer: bool,
    /// k-mer length used for chain-interval ANI (15 or 21 recommended).
    pub chained_kmer_size: usize,
}

impl Default for AniConfig {
    fn default() -> Self {
        Self {
            weight_strategy: WeightStrategy::Uniform,
            min_shared_tags: 10,
            min_af: 0.1,
            debias: true,
            use_gbrt_debias: true,
            use_gbrt_v3: false,
            use_gbrt_v3_6: false,
            use_gbrt_v4: false,
            use_gbrt_v7: true,
            use_mash_ani: true,
            mash_calibration_offset: 0.0,
            use_chained_kmer: true,
            chained_kmer_size: 15,
        }
    }
}

/// Result of an ANI calculation.
#[derive(Debug, Clone)]
pub struct AniResult {
    pub ani: f64,
    /// Raw (uncorrected) ANI before debiasing.
    pub raw_ani: f64,
    /// Mash-like ANI estimate from bidirectional tag containment.
    pub mash_ani: f64,
    /// K-mer ANI estimated inside sparse-chaining synteny blocks.
    pub chained_kmer_ani: f64,
    pub af_query: f64,
    pub af_reference: f64,
    pub weighted_ani: f64,
    pub confidence: f64,
    pub local_ani_profile: Vec<f64>,
    /// True if the ANI is below the reliable detection threshold (~83%).
    pub below_detection: bool,
}


/// Compute a per-query-contig Mash-like ANI estimate.
///
/// For each query contig, we find the reference contig that shares the most matched
/// tags and compute a local mash ANI from the per-contig alignment fractions.  The
/// per-contig estimates are weighted by the number of matched tags.  This makes the
/// estimate robust to incomplete genomes / MAGs where contig counts differ, while
/// avoiding the instability of narrow interval-based denominators.
fn compute_contig_mash_ani(match_result: &MatchResult, _avg_tag_len: f64) -> f64 {
    if match_result.matched_pairs.is_empty() {
        return 0.0;
    }

    // Total tags (matched + unmatched) per query/reference contig.
    let mut q_total: HashMap<usize, usize> = HashMap::new();
    let mut r_total: HashMap<usize, usize> = HashMap::new();

    for p in &match_result.matched_pairs {
        *q_total.entry(p.query_tag.contig_id).or_insert(0) += 1;
        *r_total.entry(p.ref_tag.contig_id).or_insert(0) += 1;
    }
    for t in &match_result.unmatched_query {
        *q_total.entry(t.contig_id).or_insert(0) += 1;
    }
    for t in &match_result.unmatched_ref {
        *r_total.entry(t.contig_id).or_insert(0) += 1;
    }

    // Matched-tag counts and average tag lengths per (q_contig, r_contig) pair.
    let mut pair_counts: HashMap<(usize, usize), (usize, f64)> = HashMap::new();
    for p in &match_result.matched_pairs {
        let entry = pair_counts
            .entry((p.query_tag.contig_id, p.ref_tag.contig_id))
            .or_insert((0, 0.0));
        entry.0 += 1;
        entry.1 += (p.query_tag.seq_len as f64 + p.ref_tag.seq_len as f64) / 2.0;
    }

    // For each query contig, pick the reference contig with the most matched tags.
    let mut best_for_query: HashMap<usize, (usize, usize, f64)> = HashMap::new(); // q_cid -> (r_cid, matched, mash)
    for ((q_cid, r_cid), (matched, len_sum)) in pair_counts {
        let q_tot = q_total.get(&q_cid).copied().unwrap_or(1).max(1) as f64;
        let r_tot = r_total.get(&r_cid).copied().unwrap_or(1).max(1) as f64;
        let local_avg_len = len_sum / matched.max(1) as f64;
        if local_avg_len <= 0.0 {
            continue;
        }
        let af_q = matched as f64 / q_tot;
        let af_r = matched as f64 / r_tot;
        if af_q <= 0.0 || af_r <= 0.0 {
            continue;
        }
        let containment_geo = (af_q * af_r).sqrt();
        let local_mash = (1.0 + containment_geo.ln() / local_avg_len).clamp(0.0, 1.0);

        let entry = best_for_query.entry(q_cid).or_insert((r_cid, matched, local_mash));
        if matched > entry.1 {
            *entry = (r_cid, matched, local_mash);
        }
    }

    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;
    for (_, (_r_cid, matched, mash)) in best_for_query {
        weighted_sum += mash * matched as f64;
        total_weight += matched as f64;
    }

    if total_weight > 0.0 {
        weighted_sum / total_weight
    } else {
        0.0
    }
}

/// Compute canonical k-mer hash sets for a DNA sequence slice.
/// Only A/C/G/T k-mers are emitted; windows containing N are skipped.
/// Each k-mer is represented by the minimum of its forward and reverse-complement
/// hash so that strand orientation does not matter.
fn canonical_kmer_set(seq: &[u8], k: usize) -> std::collections::HashSet<u64> {
    if seq.len() < k {
        return std::collections::HashSet::new();
    }
    let mut set = std::collections::HashSet::with_capacity(seq.len().saturating_sub(k).max(1));
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut valid = 0usize;
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    for &base in seq.iter() {
        let fw_bits = match base {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1u64,
            b'G' | b'g' => 2u64,
            b'T' | b't' => 3u64,
            _ => {
                valid = 0;
                forward = 0;
                reverse = 0;
                continue;
            }
        };
        let rc_bits = 3u64 - fw_bits;
        forward = ((forward << 2) | fw_bits) & mask;
        reverse = (reverse >> 2) | (rc_bits << (2 * k - 2));
        valid += 1;
        if valid >= k {
            let canonical = forward.min(reverse);
            set.insert(canonical);
        }
    }
    set
}

/// Reverse-complement a DNA sequence in place.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &base in seq.iter().rev() {
        rc.push(match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        });
    }
    rc
}

/// Compute k-mer ANI inside each sparse-chaining synteny block and return a
/// weighted average.  If raw sequences are unavailable, returns `mash_ani`.
fn compute_chained_kmer_ani(match_result: &MatchResult, mash_ani: f64, k: usize) -> f64 {
    if match_result.synteny_blocks.is_empty()
        || match_result.q_sequences.is_empty()
        || match_result.r_sequences.is_empty()
        || k == 0
    {
        return mash_ani;
    }

    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;

    for block in &match_result.synteny_blocks {
        let q_seq = match match_result.q_sequences.get(block.q_contig_id) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };
        let r_seq = match match_result.r_sequences.get(block.r_contig_id) {
            Some(s) if !s.is_empty() => s,
            _ => continue,
        };

        // Clip intervals to contig bounds.
        let q_lo = block.query_start.min(q_seq.len());
        let q_hi = (block.query_end + 1).min(q_seq.len());
        let r_lo = block.ref_start.min(r_seq.len());
        let r_hi = (block.ref_end + 1).min(r_seq.len());

        if q_hi <= q_lo + k || r_hi <= r_lo + k {
            continue;
        }

        let q_interval = &q_seq[q_lo..q_hi];
        let r_interval = if block.orientation == '-' {
            reverse_complement(&r_seq[r_lo..r_hi])
        } else {
            r_seq[r_lo..r_hi].to_vec()
        };

        let q_kmers = canonical_kmer_set(q_interval, k);
        let r_kmers = canonical_kmer_set(&r_interval, k);
        if q_kmers.is_empty() || r_kmers.is_empty() {
            continue;
        }

        let min_size = q_kmers.len().min(r_kmers.len()) as f64;
        let shared = q_kmers.intersection(&r_kmers).count() as f64;
        let containment = (shared / min_size).min(1.0).max(1e-9);
        let local_ani = (1.0 + containment.ln() / k as f64).clamp(0.0, 1.0);

        // Weight by the number of anchor tags in the chain (more anchors =
        // higher confidence in the interval alignment).
        let weight = block.matched_tags.max(1) as f64;
        weighted_sum += local_ani * weight;
        total_weight += weight;
    }

    if total_weight > 0.0 {
        (weighted_sum / total_weight).clamp(0.0, 1.0)
    } else {
        mash_ani
    }
}

/// Calculates ANI from matched tag pairs.
pub struct AniCalculator;

impl AniCalculator {
    /// Calculate ANI from a `MatchResult` using the given configuration.
    pub fn calculate_ani(match_result: &MatchResult, config: &AniConfig) -> AniResult {
        let total_q = match_result.matched_pairs.len() + match_result.unmatched_query.len();
        let total_r = match_result.matched_pairs.len() + match_result.unmatched_ref.len();

        let af_query = if total_q > 0 {
            match_result.matched_pairs.len() as f64 / total_q as f64
        } else {
            0.0
        };

        let af_reference = if total_r > 0 {
            match_result.matched_pairs.len() as f64 / total_r as f64
        } else {
            0.0
        };

        if match_result.matched_pairs.len() < config.min_shared_tags
            || af_query < config.min_af
            || af_reference < config.min_af
        {
            return AniResult {
                ani: 0.0,
                raw_ani: 0.0,
                mash_ani: 0.0,
                chained_kmer_ani: 0.0,
                af_query,
                af_reference,
                weighted_ani: 0.0,
                confidence: 0.0,
                local_ani_profile: Vec::new(),
                below_detection: true,
            };
        }

        let local_ani_profile: Vec<f64> =
            match_result.matched_pairs.iter().map(|p| p.local_ani).collect();

        let ani = if !local_ani_profile.is_empty() {
            local_ani_profile.iter().sum::<f64>() / local_ani_profile.len() as f64
        } else {
            0.0
        };

        let ani_percent = ani * 100.0;

        // Mash-like ANI from bidirectional tag containment.
        // Uses geometric mean containment and average matched tag length.
        let avg_tag_len = if match_result.matched_pairs.is_empty() {
            0.0
        } else {
            match_result
                .matched_pairs
                .iter()
                .map(|p| (p.query_tag.seq_len as f64 + p.ref_tag.seq_len as f64) / 2.0)
                .sum::<f64>()
                / match_result.matched_pairs.len() as f64
        };
        let _global_mash_ani = if af_query > 0.0 && af_reference > 0.0 && avg_tag_len > 0.0 {
            let containment_geo = (af_query * af_reference).sqrt();
            (1.0 + containment_geo.ln() / avg_tag_len).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let mash_ani = compute_contig_mash_ani(match_result, avg_tag_len);
        let chained_kmer_ani =
            compute_chained_kmer_ani(match_result, mash_ani, config.chained_kmer_size);

        let weights = Self::compute_weights(match_result, config);
        let weighted_ani = if !local_ani_profile.is_empty() {
            local_ani_profile
                .iter()
                .zip(weights.iter())
                .map(|(a, w)| a * w)
                .sum::<f64>()
                / weights.iter().sum::<f64>()
        } else {
            0.0
        };

        let shared_count = match_result.matched_pairs.len();
        let final_ani = if config.use_mash_ani {
            if config.use_chained_kmer {
                chained_kmer_ani
            } else {
                (mash_ani + config.mash_calibration_offset).clamp(0.0, 1.0)
            }
        } else if config.debias {
            if config.use_gbrt_debias {
                gbrt_debias_ani(ani, mash_ani, chained_kmer_ani, af_query, af_reference, shared_count, config)
            } else {
                simple_debias_ani(ani_percent, af_query, af_reference) / 100.0
            }
        } else {
            ani
        };

        let confidence = Self::compute_confidence(match_result, af_query, af_reference);

        AniResult {
            ani: final_ani,
            raw_ani: ani,
            mash_ani,
            chained_kmer_ani,
            af_query,
            af_reference,
            weighted_ani,
            confidence,
            local_ani_profile,
            below_detection: final_ani < 0.83,
        }
    }

    fn compute_weights(match_result: &MatchResult, config: &AniConfig) -> Vec<f64> {
        match config.weight_strategy {
            WeightStrategy::Uniform => vec![1.0; match_result.matched_pairs.len()],
            WeightStrategy::Synteny => {
                let block_map = Self::map_pairs_to_blocks(match_result);
                match_result
                    .matched_pairs
                    .iter()
                    .enumerate()
                    .map(|(i, _)| {
                        if let Some(block_idx) = block_map.get(&i) {
                            if let Some(block) = match_result.synteny_blocks.get(*block_idx) {
                                let len = block.matched_tags.max(1) as f64;
                                len.sqrt()
                            } else {
                                1.0
                            }
                        } else {
                            1.0
                        }
                    })
                    .collect()
            }
            WeightStrategy::Position => match_result
                .matched_pairs
                .iter()
                .map(|p| {
                    let norm_pos =
                        p.query_tag.position as f64 / (p.query_tag.position.max(1) as f64 + 1.0);
                    1.0 + norm_pos.sin()
                })
                .collect(),
            WeightStrategy::GapAdjusted => match_result
                .matched_pairs
                .iter()
                .map(|p| {
                    let gap_penalty = (p.gap_diff.abs() as f64).min(10.0) / 10.0;
                    1.0 - gap_penalty * 0.5
                })
                .collect(),
        }
    }

    fn map_pairs_to_blocks(
        match_result: &MatchResult,
    ) -> crate::utils::fxhash::FastHashMap<usize, usize> {
        let mut map = crate::utils::fxhash::FastHashMap::default();
        for (block_idx, block) in match_result.synteny_blocks.iter().enumerate() {
            for &pair_idx in &block.pair_indices {
                if pair_idx < match_result.matched_pairs.len() {
                    map.insert(pair_idx, block_idx);
                }
            }
        }
        map
    }

    fn compute_confidence(match_result: &MatchResult, af_q: f64, af_r: f64) -> f64 {
        let shared_count = match_result.matched_pairs.len() as f64;
        let af_min = af_q.min(af_r);
        let raw = (1.0 - (-shared_count / 100.0).exp()) * af_min.sqrt();
        raw.min(1.0).max(0.0)
    }
}

/// Simple polynomial ANI debias correction.
fn simple_debias_ani(ani: f64, af_q: f64, af_r: f64) -> f64 {
    let af_min = af_q.min(af_r);
    let correction = 0.02 * (100.0 - ani) * (1.0 - af_min);
    ani + correction
}

/// GBRT-based ANI debias correction.
/// Uses the embedded gradient-boosted regression tree model.
fn gbrt_debias_ani(raw_ani: f64, mash_ani: f64, chained_kmer_ani: f64, af_q: f64, af_r: f64, shared_count: usize, config: &AniConfig) -> f64 {
    // For tests, avoid calling the singleton model (which panics in test cfg).
    #[cfg(test)]
    {
        let _ = (shared_count, config.use_gbrt_v3, config.use_gbrt_v3_6, config.use_gbrt_v4, config.use_gbrt_v7);
        gbrt::simple_debias(raw_ani, af_q, af_r)
    }
    #[cfg(not(test))]
    {
        if config.use_gbrt_v7 {
            gbrt::load_v7_model().predict_runtime_v7(raw_ani, mash_ani, chained_kmer_ani, shared_count, af_q, af_r)
        } else if config.use_gbrt_v4 {
            gbrt::load_v4_model().predict_runtime_v4(raw_ani, shared_count, af_q, af_r)
        } else if config.use_gbrt_v3_6 {
            gbrt::load_v3_6_model().predict_runtime_v3_6(raw_ani, af_q, af_r)
        } else if config.use_gbrt_v3 {
            let total_min = (shared_count as f64 / raw_ani.max(1e-9)).max(1.0) as usize;
            let total_max = total_min;
            let containment = if total_max > 0 {
                shared_count as f64 / total_max as f64
            } else {
                0.0
            };
            gbrt::load_v3_model().predict_runtime(raw_ani, af_q, af_r, shared_count as f64, containment)
        } else {
            let total_min = (shared_count as f64 / raw_ani.max(1e-9)).max(1.0) as usize;
            let total_max = total_min;
            let containment = if total_max > 0 {
                shared_count as f64 / total_max as f64
            } else {
                0.0
            };
            gbrt::model().predict_runtime(raw_ani, af_q, af_r, shared_count as f64, containment)
        }
    }
}
