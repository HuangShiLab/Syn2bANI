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
    multi_enzyme: bool,
    enzymes: Option<&str>,
) -> Result<()> {
    let registry = EnzymeRegistry::new();
    let default_enzyme = "AloI,BslFI";
    let use_multi = multi_enzyme || enzymes.is_some();
    let enzyme_source = enzymes.unwrap_or(default_enzyme);
    let enzyme_list: Vec<_> = if multi_enzyme {
        registry.all().to_vec()
    } else {
        enzyme_source.split(',')
            .map(|name| name.trim())
            .filter(|name| !name.is_empty())
            .map(|name| {
                registry
                    .get(name)
                    .with_context(|| format!("Unknown enzyme: {}", name))
                    .map(|e| e.clone())
            })
            .collect::<Result<Vec<_>, _>>()?
    };

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
        let mut q_seqs: Vec<Vec<u8>> = Vec::with_capacity(q_records.len());
        for (cid, record) in q_records.iter().enumerate() {
            q_seqs.push(record.sequence.clone());
            for enz in &enzyme_list {
                all_q_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz, cid));
            }
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
            sequences: q_seqs,
        };

        for r_path in reference {
            let r_records = parse_fasta(r_path)
                .with_context(|| format!("Failed to parse reference: {}", r_path.display()))?;

            let mut all_r_tags = Vec::new();
            let mut r_total_len = 0usize;
            let mut r_gc_count = 0usize;
            let mut r_seqs: Vec<Vec<u8>> = Vec::with_capacity(r_records.len());
            for (cid, record) in r_records.iter().enumerate() {
                r_seqs.push(record.sequence.clone());
                for enz in &enzyme_list {
                    all_r_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz, cid));
                }
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
                sequences: r_seqs,
            };

            let match_config = if use_multi {
        MatchConfig {
            allow_near_match: false,
            near_match_tolerance: 0,
        }
    } else {
        MatchConfig::default()
    };
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
                    min_shared_tags: 0,
                    min_af: 0.0,
                    debias: true,
                    use_gbrt_debias: true,
                    use_gbrt_v3: false,
                    use_gbrt_v3_6: false,
                    use_gbrt_v4: false,
                    use_gbrt_v7: false,
                    use_mash_ani: true,
                    mash_calibration_offset: 0.0,
                    use_chained_kmer: true,
                    chained_kmer_size: 15,
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
