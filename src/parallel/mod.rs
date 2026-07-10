use rayon::ThreadPoolBuilder;
use std::path::PathBuf;

pub mod simd;

use crate::core::{
    AniCalculator, AniConfig, MatchConfig, TagExtractor, TagMatcher, TagSet,
};
use crate::enzyme::EnzymeRegistry;
use crate::io::parse_fasta;

/// Result of a single pairwise comparison.
pub struct ComparisonResult {
    pub query_path: PathBuf,
    pub ref_path: PathBuf,
    pub ani: f64,
    pub af_query: f64,
    pub af_reference: f64,
    pub shared_tags: usize,
}

/// Genome-level parallel comparison using a Rayon thread pool.
///
/// Each pair is loaded, digested, matched, and scored independently.
pub fn parallel_compare(
    pairs: &[(PathBuf, PathBuf)],
    config: &AniConfig,
    threads: usize,
) -> Vec<ComparisonResult> {
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build thread pool");

    let registry = EnzymeRegistry::new();
    let default_enz = registry.get("BcgI").unwrap().clone();

    pool.install(|| {
        use rayon::prelude::*;
        pairs
            .par_iter()
            .map(|(q_path, r_path)| {
                let q_records = parse_fasta(q_path).unwrap_or_default();
                let r_records = parse_fasta(r_path).unwrap_or_default();

                if q_records.is_empty() || r_records.is_empty() {
                    return ComparisonResult {
                        query_path: q_path.clone(),
                        ref_path: r_path.clone(),
                        ani: 0.0,
                        af_query: 0.0,
                        af_reference: 0.0,
                        shared_tags: 0,
                    };
                }

                let q_tags: Vec<_> = q_records
                    .iter()
                    .flat_map(|r| TagExtractor::extract_from_sequence(&r.sequence, &default_enz))
                    .collect();
                let q_total_len: usize = q_records.iter().map(|r| r.sequence.len()).sum();
                let q_gc_count: usize = q_records
                    .iter()
                    .map(|r| {
                        r.sequence
                            .iter()
                            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                            .count()
                    })
                    .sum();

                let r_tags: Vec<_> = r_records
                    .iter()
                    .flat_map(|r| TagExtractor::extract_from_sequence(&r.sequence, &default_enz))
                    .collect();
                let r_total_len: usize = r_records.iter().map(|r| r.sequence.len()).sum();
                let r_gc_count: usize = r_records
                    .iter()
                    .map(|r| {
                        r.sequence
                            .iter()
                            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                            .count()
                    })
                    .sum();

                let q_tag_set = TagSet {
                    genome_id: q_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown")
                        .to_string(),
                    chromosome: "all".to_string(),
                    tags: q_tags,
                    total_length: q_total_len,
                    gc_content: q_gc_count as f64 / q_total_len.max(1) as f64,
                };

                let r_tag_set = TagSet {
                    genome_id: r_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown")
                        .to_string(),
                    chromosome: "all".to_string(),
                    tags: r_tags,
                    total_length: r_total_len,
                    gc_content: r_gc_count as f64 / r_total_len.max(1) as f64,
                };

                let match_config = MatchConfig::default();
                let match_result = TagMatcher::match_tag_sets(&q_tag_set, &r_tag_set, &match_config);
                let ani_result = AniCalculator::calculate_ani(&match_result, config);

                ComparisonResult {
                    query_path: q_path.clone(),
                    ref_path: r_path.clone(),
                    ani: ani_result.ani,
                    af_query: ani_result.af_query,
                    af_reference: ani_result.af_reference,
                    shared_tags: match_result.matched_pairs.len(),
                }
            })
            .collect()
    })
}
