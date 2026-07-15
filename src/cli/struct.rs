use anyhow::{Context, Result};
use std::io::{self, Write};
use std::fs::File;
use std::path::Path;

use crate::core::{
    AniCalculator, AniConfig, MatchConfig, StructureAnalyzer, TagExtractor, TagMatcher, TagSet,
    WeightStrategy,
};
use crate::enzyme::EnzymeRegistry;
use crate::io::{parse_fasta, ExtendedTsvFormatter};

/// Handler for the `struct` subcommand.
///
/// Performs structural variation analysis between query and reference genomes,
/// outputting either PAF or extended TSV.
pub fn run_struct(
    query: &[std::path::PathBuf],
    reference: &[std::path::PathBuf],
    output: Option<&Path>,
    paf: bool,
    rearrangement: bool,
    indel: bool,
) -> Result<()> {
    let registry = EnzymeRegistry::new();
    let default_enz = registry.get("BcgI").unwrap().clone();

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(File::create(path)?)
    } else {
        Box::new(io::stdout())
    };

    for q_path in query {
        let q_records = parse_fasta(q_path)
            .with_context(|| format!("Failed to parse query: {}", q_path.display()))?;

        let mut all_q_tags = Vec::new();
        let mut q_total_len = 0usize;
        let mut q_gc_count = 0usize;
        for record in &q_records {
            all_q_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, &default_enz));
            q_total_len += record.sequence.len();
            q_gc_count += record
                .sequence
                .iter()
                .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                .count();
        }

        let q_tag_set = TagSet {
            genome_id: q_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string(),
            chromosome: "all".to_string(),
            tags: all_q_tags,
            total_length: q_total_len,
            gc_content: q_gc_count as f64 / q_total_len.max(1) as f64,
        };

        for r_path in reference {
            let r_records = parse_fasta(r_path)
                .with_context(|| format!("Failed to parse reference: {}", r_path.display()))?;

            let mut all_r_tags = Vec::new();
            let mut r_total_len = 0usize;
            let mut r_gc_count = 0usize;
            for record in &r_records {
                all_r_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, &default_enz));
                r_total_len += record.sequence.len();
                r_gc_count += record
                    .sequence
                    .iter()
                    .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                    .count();
            }

            let r_tag_set = TagSet {
                genome_id: r_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string(),
                chromosome: "all".to_string(),
                tags: all_r_tags,
                total_length: r_total_len,
                gc_content: r_gc_count as f64 / r_total_len.max(1) as f64,
            };

            let match_config = MatchConfig::default();
            let match_result = TagMatcher::match_tag_sets(&q_tag_set, &r_tag_set, &match_config);

            let mut svs = Vec::new();
            if rearrangement {
                svs.extend(StructureAnalyzer::detect_rearrangements(&match_result.synteny_blocks));
            }
            if indel {
                svs.extend(StructureAnalyzer::detect_indels(&match_result));
            }

            if paf {
                let paf_str = StructureAnalyzer::to_paf(&svs);
                writeln!(writer, "{}", paf_str)?;
            } else {
                ExtendedTsvFormatter::write_header(&mut writer)?;
                let ani_config = AniConfig {
                    weight_strategy: WeightStrategy::Uniform,
                    min_shared_tags: 10,
                    min_af: 0.1,
                    debias: true,
                    use_gbrt_debias: true,
                    use_gbrt_v3: false,
                    use_gbrt_v3_6: true,
                };
                let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);
                let rearrangements = if rearrangement {
                    svs.iter()
                        .filter(|sv| {
                            matches!(
                                sv.sv_type,
                                crate::core::SvType::Inversion | crate::core::SvType::Translocation
                            )
                        })
                        .count()
                } else {
                    0
                };
                let indels = if indel {
                    svs.iter()
                        .filter(|sv| {
                            matches!(
                                sv.sv_type,
                                crate::core::SvType::Insertion | crate::core::SvType::Deletion
                            )
                        })
                        .count()
                } else {
                    0
                };
                ExtendedTsvFormatter::write_record(
                    &mut writer,
                    &q_path.display().to_string(),
                    &r_path.display().to_string(),
                    &q_tag_set.genome_id,
                    &r_tag_set.genome_id,
                    &ani_result,
                    svs.len(),
                    rearrangements,
                    indels,
                    match_result.synteny_blocks.len(),
                )?;
            }
        }
    }

    Ok(())
}
