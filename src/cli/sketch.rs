use anyhow::{Context, Result};
use rayon::prelude::*;
use std::path::{Path, PathBuf};

use crate::core::TagExtractor;
use crate::enzyme::EnzymeRegistry;
use crate::io::{
    parse_fasta, ChromSketch, SketchEnzyme, SketchMetadata, SketchTag, TgtSketch, write_sketch,
};

/// Handler for the `sketch` subcommand.
///
/// Builds binary sketch files (`.s2ba`) from one or more input genomes.
pub fn run_sketch(
    genomes: &[PathBuf],
    output: &Path,
    enzyme: &str,
    threads: usize,
    parallel: bool,
    multi_enzyme: bool,
    enzyme_list: Option<&str>,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let registry = EnzymeRegistry::new();
    // A comma-separated list wins, so a sketch can be built with exactly the
    // panel `ani` will use.
    let enzymes = if let Some(list) = enzyme_list {
        list.split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|name| {
                registry
                    .get(name)
                    .with_context(|| format!("Unknown enzyme: {name}"))
                    .map(|e| e.clone())
            })
            .collect::<Result<Vec<_>>>()?
    } else if multi_enzyme {
        registry.all().to_vec()
    } else {
        vec![registry
            .get(enzyme)
            .with_context(|| format!("Unknown enzyme: {}", enzyme))?
            .clone()]
    };

    // Recorded in the sketch so readers never have to guess the panel.
    let enzyme_table: Vec<SketchEnzyme> = enzymes
        .iter()
        .map(|e| SketchEnzyme {
            name: e.name.clone(),
            tag_length: e.tag_length as u8,
            site_length: (e.left_anchor.len() + e.right_anchor.len()) as u8,
        })
        .collect();

    std::fs::create_dir_all(output)?;

    // Parallel sketch computation, serial I/O for deterministic file output
    let sketches: Vec<_> = pool.install(|| {
        genomes
            .par_iter()
            .filter_map(|genome_path| {
                let records = parse_fasta(genome_path).ok()?;

                let mut chromosomes = Vec::new();
                let mut total_length = 0u64;
                let mut total_gc = 0.0f64;
                let mut total_tags = 0u64;

                for (cid, record) in records.iter().enumerate() {
                    let mut chrom_tags = Vec::new();
                    let mut per_enzyme_tags = Vec::new();
                    for (enz_idx, enz) in enzymes.iter().enumerate() {
                        let tags = TagExtractor::extract_from_sequence(&record.sequence, enz, cid);
                        per_enzyme_tags.push((enz_idx, tags));
                    }
                    for (enz_idx, tags) in per_enzyme_tags {
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
                    enzymes: enzyme_table.clone(),
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

                Some((genome_path.clone(), sketch))
            })
            .collect()
    });

    for (_genome_path, sketch) in sketches {
        let out_path = output.join(format!("{}.s2ba", sketch.genome_id));
        write_sketch(&sketch, &out_path)
            .with_context(|| format!("Failed to write sketch: {}", out_path.display()))?;
    }

    Ok(())
}
