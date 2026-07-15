use anyhow::Result;
use rayon::prelude::*;
use std::io::{self, Write};
use std::fs::File;
use std::path::Path;

use crate::core::{
    AniCalculator, AniConfig, MatchConfig, TagExtractor, TagMatcher, TagSet, WeightStrategy,
};
use crate::enzyme::EnzymeRegistry;
use crate::io::parse_fasta;

/// Handler for the `triangle` subcommand.
///
/// Performs all-to-all pairwise comparisons and outputs either an edge list
/// or a full symmetric matrix.
pub fn run_triangle(
    genomes: &[std::path::PathBuf],
    output: Option<&Path>,
    edge_list: bool,
    threads: usize,
    parallel: bool,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let registry = EnzymeRegistry::new();
    let default_enz = registry.get("BcgI").unwrap().clone();

    let mut tag_sets = Vec::new();
    for path in genomes {
        let records = parse_fasta(path)?;
        let mut all_tags = Vec::new();
        let mut total_len = 0usize;
        let mut gc_count = 0usize;
        for record in &records {
            let tags = TagExtractor::extract_from_sequence(&record.sequence, &default_enz);
            all_tags.extend(tags);
            total_len += record.sequence.len();
            gc_count += record
                .sequence
                .iter()
                .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                .count();
        }
        tag_sets.push((
            path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string(),
            TagSet {
                genome_id: path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string(),
                chromosome: "all".to_string(),
                tags: all_tags,
                total_length: total_len,
                gc_content: gc_count as f64 / total_len.max(1) as f64,
            },
        ));
    }

    let pairs: Vec<_> = (0..tag_sets.len())
        .flat_map(|i| ((i + 1)..tag_sets.len()).map(move |j| (i, j)))
        .collect();

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(File::create(path)?)
    } else {
        Box::new(io::stdout())
    };

    let match_config = MatchConfig::default();
    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 10,
        min_af: 0.0,
        debias: true,
        use_gbrt_debias: true,
        use_gbrt_v3: false,
        use_gbrt_v3_6: true,
    };

    let results: Vec<_> = pool.install(|| {
        if parallel {
            pairs
                .par_iter()
                .map(|&(i, j)| {
                    let match_result =
                        TagMatcher::match_tag_sets(&tag_sets[i].1, &tag_sets[j].1, &match_config);
                    let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);
                    (i, j, ani_result.ani, match_result.matched_pairs.len())
                })
                .collect()
        } else {
            pairs
                .iter()
                .map(|&(i, j)| {
                    let match_result =
                        TagMatcher::match_tag_sets(&tag_sets[i].1, &tag_sets[j].1, &match_config);
                    let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);
                    (i, j, ani_result.ani, match_result.matched_pairs.len())
                })
                .collect()
        }
    });

    if edge_list {
        writeln!(writer, "query\treference\tani\tshared_tags")?;
        for (i, j, ani, shared) in results {
            writeln!(writer, "{}\t{}\t{:.4}\t{}", tag_sets[i].0, tag_sets[j].0, ani, shared)?;
        }
    } else {
        // Full symmetric matrix
        let n = tag_sets.len();
        let mut matrix = vec![vec![0.0f64; n]; n];
        for (i, j, ani, _shared) in results {
            matrix[i][j] = ani;
            matrix[j][i] = ani;
        }

        write!(writer, "\t")?;
        for (name, _) in &tag_sets {
            write!(writer, "{}\t", name)?;
        }
        writeln!(writer)?;

        for i in 0..n {
            write!(writer, "{}\t", tag_sets[i].0)?;
            for j in 0..n {
                if i == j {
                    write!(writer, "100.0")?;
                } else {
                    write!(writer, "{:.4}", matrix[i][j])?;
                }
                if j < n - 1 {
                    write!(writer, "\t")?;
                }
            }
            writeln!(writer)?;
        }
    }

    Ok(())
}
