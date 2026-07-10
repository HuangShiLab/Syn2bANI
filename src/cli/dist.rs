use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::{self, Write};
use std::fs::File;
use std::path::Path;

use crate::core::{AniCalculator, AniConfig, GenomeTag, MatchConfig, TagExtractor, TagMatcher, TagSet, WeightStrategy};
use crate::enzyme::EnzymeRegistry;
use crate::io::{parse_fasta, TsvFormatter};

/// Handler for the `dist` subcommand.
///
/// Performs a two-pass comparison:
/// 1. Coarse screening via shared tag count / max-containment.
/// 2. Fine ANI calculation using the full tag-matching pipeline.
pub fn run_dist(
    query: &[std::path::PathBuf],
    reference: &[std::path::PathBuf],
    enzyme: &str,
    threads: usize,
    parallel: bool,
    multi_enzyme: bool,
    structural: bool,
    output: Option<&Path>,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let registry = EnzymeRegistry::new();
    let enzymes = if multi_enzyme {
        registry.all().to_vec()
    } else {
        vec![registry
            .get(enzyme)
            .with_context(|| format!("Unknown enzyme: {}", enzyme))?
            .clone()]
    };

    let match_config = MatchConfig::default();
    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 10,
        min_af: 0.1,
        debias: true,
        use_gbrt_debias: true,
    };

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(File::create(path)?)
    } else {
        Box::new(io::stdout())
    };

    TsvFormatter::write_header(&mut writer)?;

    for q_path in query {
        let q_name = q_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");
        let q_records = parse_fasta(q_path)
            .with_context(|| format!("Failed to parse query FASTA: {}", q_path.display()))?;

        let mut all_q_tags: Vec<GenomeTag> = Vec::new();
        let mut q_total_len = 0usize;
        let mut q_gc_count = 0usize;
        for record in &q_records {
            for enz in &enzymes {
                all_q_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz));
            }
            q_total_len += record.sequence.len();
            q_gc_count += record
                .sequence
                .iter()
                .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                .count();
        }

        let q_tag_set = TagSet {
            genome_id: q_name.to_string(),
            chromosome: "all".to_string(),
            tags: all_q_tags,
            total_length: q_total_len,
            gc_content: q_gc_count as f64 / q_total_len.max(1) as f64,
        };

        // Parallelize over references while preserving output order
        let ref_indices: Vec<usize> = (0..reference.len()).collect();
        let results: Vec<_> = pool.install(|| {
            ref_indices
                .into_par_iter()
                .filter_map(|idx| {
                    let r_path = &reference[idx];
                    let r_name = r_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown");
                    let r_records = parse_fasta(r_path).ok()?;

                    let mut all_r_tags: Vec<GenomeTag> = Vec::new();
                    let mut r_total_len = 0usize;
                    let mut r_gc_count = 0usize;
                    for record in &r_records {
                        for enz in &enzymes {
                            all_r_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz));
                        }
                        r_total_len += record.sequence.len();
                        r_gc_count += record
                            .sequence
                            .iter()
                            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                            .count();
                    }

                    let r_tag_set = TagSet {
                        genome_id: r_name.to_string(),
                        chromosome: "all".to_string(),
                        tags: all_r_tags,
                        total_length: r_total_len,
                        gc_content: r_gc_count as f64 / r_total_len.max(1) as f64,
                    };

                    // Pass 1: coarse screening
                    let shared = shared_tag_count(&q_tag_set, &r_tag_set);
                    let max_tags = q_tag_set.tags.len().max(r_tag_set.tags.len()).max(1);
                    let coarse_ani = shared as f64 / max_tags as f64;

                    if coarse_ani < ani_config.min_af {
                        return None;
                    }

                    // Pass 2: fine ANI
                    let match_result = TagMatcher::match_tag_sets(&q_tag_set, &r_tag_set, &match_config);
                    let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);

                    let sv_count = if structural { 0 } else { 0 };

                    Some((idx, r_path.clone(), r_tag_set.genome_id, ani_result, sv_count))
                })
                .collect()
        });

        // Sort by original index to maintain deterministic output order
        let mut sorted = results;
        sorted.sort_by_key(|(idx, _, _, _, _)| *idx);

        for (_idx, r_path, r_genome_id, ani_result, sv_count) in sorted {
            TsvFormatter::write_record(
                &mut writer,
                &q_path.display().to_string(),
                &r_path.display().to_string(),
                &q_tag_set.genome_id,
                &r_genome_id,
                &ani_result,
                sv_count,
            )?;
        }
    }

    Ok(())
}

fn shared_tag_count(q: &TagSet, r: &TagSet) -> usize {
    let q_seqs: HashSet<&[u8]> = q.tags.iter().map(|t| t.sequence.as_slice()).collect();
    r.tags
        .iter()
        .filter(|t| q_seqs.contains(t.sequence.as_slice()))
        .count()
}
