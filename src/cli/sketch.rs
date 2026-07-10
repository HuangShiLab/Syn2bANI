use anyhow::{Context, Result};
use std::path::{Path, PathBuf};

use crate::core::TagExtractor;
use crate::enzyme::EnzymeRegistry;
use crate::io::{
    parse_fasta, ChromSketch, SketchMetadata, SketchTag, TgtSketch, write_sketch,
};

/// Handler for the `sketch` subcommand.
///
/// Builds binary sketch files (`.s2ba`) from one or more input genomes.
pub fn run_sketch(
    genomes: &[PathBuf],
    output: &Path,
    enzyme: &str,
    _threads: usize,
    multi_enzyme: bool,
) -> Result<()> {
    let registry = EnzymeRegistry::new();
    let enzymes = if multi_enzyme {
        registry.all().to_vec()
    } else {
        vec![registry
            .get(enzyme)
            .with_context(|| format!("Unknown enzyme: {}", enzyme))?
            .clone()]
    };

    std::fs::create_dir_all(output)?;

    for genome_path in genomes {
        let records = parse_fasta(genome_path)
            .with_context(|| format!("Failed to parse: {}", genome_path.display()))?;

        let mut chromosomes = Vec::new();
        let mut total_length = 0u64;
        let mut total_gc = 0.0f64;
        let mut total_tags = 0u64;

        for record in &records {
            let mut chrom_tags = Vec::new();
            for (enz_idx, enz) in enzymes.iter().enumerate() {
                let tags = TagExtractor::extract_from_sequence(&record.sequence, enz);
                for tag in tags {
                    chrom_tags.push(SketchTag {
                        position: tag.position as u64,
                        seq: crate::io::pack_sequence(&tag.sequence),
                        direction: if tag.direction == '+' { 0 } else { 1 },
                        enzyme_id: enz_idx as u16,
                    });
                }
            }

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
            genome_id: genome_path
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

        let out_path = output.join(format!("{}.s2ba", sketch.genome_id));
        write_sketch(&sketch, &out_path)
            .with_context(|| format!("Failed to write sketch: {}", out_path.display()))?;
    }

    Ok(())
}
