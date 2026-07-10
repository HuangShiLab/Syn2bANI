use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::fs;

use crate::core::TagExtractor;
use crate::enzyme::EnzymeRegistry;
use crate::io::{
    parse_fasta, read_sketch, write_sketch, ChromSketch, SketchMetadata, SketchTag, TgtSketch,
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
