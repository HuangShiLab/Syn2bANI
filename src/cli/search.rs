use anyhow::{Context, Result};
use std::io::{self, Write};
use std::fs::File;
use std::path::{Path, PathBuf};
use rayon::prelude::*;

use crate::core::{
    AniCalculator, AniConfig, GenomeTag, MatchConfig, TagExtractor, TagMatcher, TagSet,
    WeightStrategy,
};
use crate::enzyme::EnzymeRegistry;
use crate::io::{parse_fasta, read_sketch, TgtSketch, TsvFormatter};

/// Handler for the `search` subcommand.
///
/// Loads a sketch database and searches query genomes against it,
/// filtering results by a minimum ANI threshold.
pub fn run_search(
    query: &[PathBuf],
    database: &Path,
    output: Option<&Path>,
    threads: usize,
    parallel: bool,
    min_ani: f64,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let db_entries: Vec<TgtSketch> = std::fs::read_dir(database)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|ext| ext == "s2ba").unwrap_or(false))
        .filter_map(|e| read_sketch(&e.path()).ok())
        .collect();

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(File::create(path)?)
    } else {
        Box::new(io::stdout())
    };

    TsvFormatter::write_header(&mut writer)?;

    let match_config = MatchConfig::default();
    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 10,
        min_af: 0.0,
        debias: true,
        use_gbrt_debias: true,
    };

    let registry = EnzymeRegistry::new();
    let default_enz = registry.get("BcgI").unwrap().clone();

    for q_path in query {
        let q_records = parse_fasta(q_path)
            .with_context(|| format!("Failed to parse query: {}", q_path.display()))?;

        let mut all_q_tags: Vec<GenomeTag> = Vec::new();
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

        let results: Vec<_> = pool.install(|| {
            db_entries
                .par_iter()
                .filter_map(|db_sketch| {
                    let mut db_tags = Vec::new();
                    for chrom in &db_sketch.chromosomes {
                        for tag in &chrom.tags {
                            let unpacked = crate::io::unpack_sequence(tag.seq);
                            let mut sequence = [0u8; 32];
                            let copy_len = unpacked.len().min(32);
                            sequence[..copy_len].copy_from_slice(&unpacked[..copy_len]);
                            db_tags.push(GenomeTag {
                                position: tag.position as usize,
                                sequence,
                                packed_sequence: tag.seq,
                                seq_len: copy_len as u8,
                                direction: if tag.direction == 0 { '+' } else { '-' },
                                enzyme: "unknown".to_string(),
                            });
                        }
                    }

                    let db_tag_set = TagSet {
                        genome_id: db_sketch.genome_id.clone(),
                        chromosome: "all".to_string(),
                        tags: db_tags,
                        total_length: db_sketch.metadata.total_length as usize,
                        gc_content: db_sketch.metadata.gc_content,
                    };

                    let match_result = TagMatcher::match_tag_sets(&q_tag_set, &db_tag_set, &match_config);
                    let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);

                    if ani_result.ani >= min_ani {
                        Some((db_sketch.genome_id.clone(), ani_result, match_result))
                    } else {
                        None
                    }
                })
                .collect()
        });

        for (db_id, ani_result, _match_result) in results {
            TsvFormatter::write_record(
                &mut writer,
                &q_path.display().to_string(),
                &format!("{}/{}.s2ba", database.display(), db_id),
                &q_tag_set.genome_id,
                &db_id,
                &ani_result,
                0,
            )?;
        }
    }

    Ok(())
}
