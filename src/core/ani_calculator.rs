use crate::core::tag_matcher::MatchResult;
use crate::core::gbrt;

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
}

impl Default for AniConfig {
    fn default() -> Self {
        Self {
            weight_strategy: WeightStrategy::Uniform,
            min_shared_tags: 10,
            min_af: 0.1,
            debias: true,
            use_gbrt_debias: true,
        }
    }
}

/// Result of an ANI calculation.
#[derive(Debug, Clone)]
pub struct AniResult {
    pub ani: f64,
    /// Raw (uncorrected) ANI before debiasing.
    pub raw_ani: f64,
    pub af_query: f64,
    pub af_reference: f64,
    pub weighted_ani: f64,
    pub confidence: f64,
    pub local_ani_profile: Vec<f64>,
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
                af_query,
                af_reference,
                weighted_ani: 0.0,
                confidence: 0.0,
                local_ani_profile: Vec::new(),
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

        let final_ani = if config.debias {
            if config.use_gbrt_debias {
                gbrt_debias_ani(ani, af_query, af_reference, total_q, total_r)
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
            af_query,
            af_reference,
            weighted_ani,
            confidence,
            local_ani_profile,
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
    ) -> std::collections::HashMap<usize, usize> {
        let mut map = std::collections::HashMap::new();
        for (block_idx, block) in match_result.synteny_blocks.iter().enumerate() {
            for (pair_idx, pair) in match_result.matched_pairs.iter().enumerate() {
                if pair.query_tag.position >= block.query_start
                    && pair.query_tag.position <= block.query_end
                {
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
fn gbrt_debias_ani(raw_ani: f64, af_q: f64, af_r: f64, total_q: usize, total_r: usize) -> f64 {
    let shared = (raw_ani * (total_q.min(total_r) as f64)).max(1.0) as f64; // approximate shared tags
    let max_tags = total_q.max(total_r).max(1) as f64;
    let containment = shared / max_tags;

    // For tests, avoid calling the singleton model (which panics in test cfg).
    #[cfg(test)]
    {
        gbrt::simple_debias(raw_ani, af_q, af_r)
    }
    #[cfg(not(test))]
    {
        gbrt::model().predict_runtime(raw_ani, af_q, af_r, shared, containment)
    }
}
