use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, Write};

use rayon::prelude::*;

use crate::core::{
    AniCalculator, AniConfig, GenomeTag, MatchConfig, TagExtractor, TagMatcher, TagSet,
    WeightStrategy,
};
use crate::enzyme::EnzymeRegistry;
use crate::io::{
    parse_fasta, read_sketch, write_sketch, ChromSketch, SketchMetadata, SketchTag, TgtSketch,
    TsvFormatter,
};

/// Build a sketch database from a set of genomes.
pub fn run_db_build(
    genomes: &[PathBuf],
    output: &Path,
    enzyme: &str,
    threads: usize,
    parallel: bool,
    multi_enzyme: bool,
) -> Result<()> {
    crate::cli::sketch::run_sketch(genomes, output, enzyme, threads, parallel, multi_enzyme)
}

/// Add genomes to an existing sketch database.
pub fn run_db_add(genomes: &[PathBuf], database: &Path) -> Result<()> {
    fs::create_dir_all(database)?;
    let registry = EnzymeRegistry::new();
    let default_enz = registry.get("BcgI").unwrap().clone();

    for genome in genomes {
        let records = parse_fasta(genome)
            .with_context(|| format!("Failed to parse: {}", genome.display()))?;

        let mut chromosomes = Vec::new();
        let mut total_length = 0u64;
        let mut total_gc = 0.0f64;
        let mut total_tags = 0u64;

        for record in &records {
            let tags = TagExtractor::extract_from_sequence(&record.sequence, &default_enz);
            let chrom_tags: Vec<_> = tags
                .into_iter()
                .map(|tag| SketchTag {
                    position: tag.position as u64,
                    seq: crate::io::pack_sequence(&tag.sequence),
                    direction: if tag.direction == '+' { 0 } else { 1 },
                    enzyme_id: 0,
                })
                .collect();

            let len = record.sequence.len() as u64;
            let gc = crate::utils::gc_content(&record.sequence);
            total_length += len;
            total_gc += gc * len as f64;
            total_tags += chrom_tags.len() as u64;

            chromosomes.push(ChromSketch {
                name: record.id.clone(),
                tags: chrom_tags,
                gc_content: gc,
                length: len,
            });
        }

        let sketch = TgtSketch {
            genome_id: genome
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string(),
            chromosomes,
            metadata: SketchMetadata {
                total_length,
                gc_content: if total_length > 0 {
                    total_gc / total_length as f64
                } else {
                    0.0
                },
                tag_count: total_tags,
            },
            ..Default::default()
        };

        let out_path = database.join(format!("{}.s2ba", sketch.genome_id));
        write_sketch(&sketch, &out_path)
            .with_context(|| format!("Failed to write sketch: {}", out_path.display()))?;
    }

    Ok(())
}

/// Remove genome sketches from a database by ID.
pub fn run_db_remove(genome_ids: &[String], database: &Path) -> Result<()> {
    for id in genome_ids {
        let path = database.join(format!("{}.s2ba", id));
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("Failed to remove: {}", path.display()))?;
        }
    }
    Ok(())
}

/// List all entries in a sketch database.
pub fn run_db_list(database: &Path) -> Result<()> {
    for entry in fs::read_dir(database)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "s2ba").unwrap_or(false) {
            let sketch = read_sketch(&path)
                .with_context(|| format!("Failed to read: {}", path.display()))?;
            println!(
                "{}\t{}\t{}\t{:.4}",
                sketch.genome_id,
                sketch.chromosomes.len(),
                sketch.metadata.tag_count,
                sketch.metadata.gc_content
            );
        }
    }
    Ok(())
}

/// Search query sketches against a sketch database.
///
/// Loads all `.s2ba` files from `queries` and `database`, then computes
/// pairwise ANI for every query–db combination, filtering by `min_ani`.
pub fn run_db_search(
    queries: &Path,
    database: &Path,
    output: Option<&Path>,
    threads: usize,
    parallel: bool,
    min_ani: f64,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let query_entries: Vec<TgtSketch> = fs::read_dir(queries)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|ext| ext == "s2ba").unwrap_or(false))
        .filter_map(|e| read_sketch(&e.path()).ok())
        .collect();

    let db_entries: Vec<TgtSketch> = fs::read_dir(database)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|ext| ext == "s2ba").unwrap_or(false))
        .filter_map(|e| read_sketch(&e.path()).ok())
        .collect();

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(std::fs::File::create(path)?)
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
        use_gbrt_v3: true,
    };

    for q_sketch in &query_entries {
        let q_tags: Vec<GenomeTag> = q_sketch.chromosomes.iter()
            .flat_map(|chrom| {
                chrom.tags.iter().map(|tag| {
                    let unpacked = crate::io::unpack_sequence(tag.seq);
                    let mut sequence = [0u8; 32];
                    let copy_len = unpacked.len().min(32);
                    sequence[..copy_len].copy_from_slice(&unpacked[..copy_len]);
                    GenomeTag {
                        position: tag.position as usize,
                        sequence,
                        packed_sequence: tag.seq,
                        seq_len: copy_len as u8,
                        direction: if tag.direction == 0 { '+' } else { '-' },
                        enzyme: "unknown".to_string(),
                    }
                })
            })
            .collect();

        let q_tag_set = TagSet {
            genome_id: q_sketch.genome_id.clone(),
            chromosome: "all".to_string(),
            tags: q_tags,
            total_length: q_sketch.metadata.total_length as usize,
            gc_content: q_sketch.metadata.gc_content,
        };

        let results: Vec<_> = pool.install(|| {
            db_entries
                .par_iter()
                .filter_map(|db_sketch| {
                    let db_tags: Vec<GenomeTag> = db_sketch.chromosomes.iter()
                        .flat_map(|chrom| {
                            chrom.tags.iter().map(|tag| {
                                let unpacked = crate::io::unpack_sequence(tag.seq);
                                let mut sequence = [0u8; 32];
                                let copy_len = unpacked.len().min(32);
                                sequence[..copy_len].copy_from_slice(&unpacked[..copy_len]);
                                GenomeTag {
                                    position: tag.position as usize,
                                    sequence,
                                    packed_sequence: tag.seq,
                                    seq_len: copy_len as u8,
                                    direction: if tag.direction == 0 { '+' } else { '-' },
                                    enzyme: "unknown".to_string(),
                                }
                            })
                        })
                        .collect();

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
                        Some((db_sketch.genome_id.clone(), ani_result))
                    } else {
                        None
                    }
                })
                .collect()
        });

        for (db_id, ani_result) in results {
            TsvFormatter::write_record(
                &mut writer,
                &format!("{}/{}.s2ba", queries.display(), q_sketch.genome_id),
                &format!("{}/{}.s2ba", database.display(), db_id),
                &q_sketch.genome_id,
                &db_id,
                &ani_result,
                0,
            )?;
        }
    }

    Ok(())
}

/// Merge multiple sketch databases into one.
pub fn run_db_merge(databases: &[PathBuf], output: &Path) -> Result<()> {
    fs::create_dir_all(output)?;
    for db in databases {
        for entry in fs::read_dir(db)? {
            let entry = entry?;
            let src = entry.path();
            if src.extension().map(|e| e == "s2ba").unwrap_or(false) {
                let dst = output.join(src.file_name().unwrap());
                fs::copy(&src, &dst).with_context(|| {
                    format!("Failed to copy: {} -> {}", src.display(), dst.display())
                })?;
            }
        }
    }
    Ok(())
}
