use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::{self, Write};
use std::fs::File;
use std::path::Path;

use crate::core::{AniCalculator, AniConfig, GenomeTag, MatchConfig, TagExtractor, TagMatcher, TagSet, WeightStrategy};
use crate::enzyme::EnzymeRegistry;
use crate::io::{parse_fasta, TsvFormatter};

/// Write a raw-features training TSV header.
fn write_raw_features_header<W: Write>(writer: &mut W) -> io::Result<()> {
    writeln!(
        writer,
        "query_file\tref_file\tquery_name\tref_name\t\
         raw_ani\tmash_ani\tchained_kmer_ani\taf_q\taf_r\tshared_tags\tcontainment\tdiv_proxy\tref_gc\t\
         corrected_ani"
    )
}

/// Write a raw-features training record.
fn write_raw_features_record<W: Write>(
    writer: &mut W,
    query_file: &str,
    ref_file: &str,
    query_name: &str,
    ref_name: &str,
    raw_ani: f64,
    mash_ani: f64,
    chained_kmer_ani: f64,
    af_q: f64,
    af_r: f64,
    shared_tags: usize,
    containment: f64,
    ref_gc: f64,
    corrected_ani: f64,
) -> io::Result<()> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
        query_file, ref_file, query_name, ref_name,
        raw_ani, mash_ani, chained_kmer_ani, af_q, af_r, shared_tags, containment,
        1.0 - raw_ani, ref_gc, corrected_ani
    )
}

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
    enzymes: Option<&str>,
    structural: bool,
    raw_features: bool,
    use_mash_ani: bool,
    min_af: f64,
    output: Option<&Path>,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;

    let registry = EnzymeRegistry::new();
    let use_multi = multi_enzyme || enzymes.is_some() || enzyme.contains(',');
    let enzyme_list: Vec<_> = if let Some(list) = enzymes {
        list.split(',')
            .map(|name| name.trim())
            .filter(|name| !name.is_empty())
            .map(|name| {
                registry
                    .get(name)
                    .with_context(|| format!("Unknown enzyme: {}", name))
                    .map(|e| e.clone())
            })
            .collect::<Result<Vec<_>, _>>()?
    } else if multi_enzyme {
        registry.all().to_vec()
    } else if enzyme.contains(',') {
        enzyme.split(',')
            .map(|name| name.trim())
            .filter(|name| !name.is_empty())
            .map(|name| {
                registry
                    .get(name)
                    .with_context(|| format!("Unknown enzyme: {}", name))
                    .map(|e| e.clone())
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        vec![registry
            .get(enzyme)
            .with_context(|| format!("Unknown enzyme: {}", enzyme))?
            .clone()]
    };

    let match_config = if use_multi {
        // Multi-enzyme tag sets are much larger and mixes tags from different
        // enzymes. The default near-match fallback scans all reference tags for
        // every unmatched query tag, which is O(n*m) and prohibitively slow.
        // Exact matching is sufficient when many enzymes are combined.
        MatchConfig {
            allow_near_match: false,
            near_match_tolerance: 0,
        }
    } else {
        MatchConfig::default()
    };
    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 10,
        min_af,
        debias: true,
        use_gbrt_debias: true,
        use_gbrt_v3: false,
        use_gbrt_v3_6: false,
        use_gbrt_v4: false,
        use_gbrt_v7: use_mash_ani,  // --mash-ani flag requests GBRT-debiased output; use v7
        use_mash_ani: !use_mash_ani,  // default output is mash-like (chained k-mer)
        mash_calibration_offset: 0.0,
        use_chained_kmer: !use_mash_ani,  // chain-interval k-mer ANI when mash_ani is default
        chained_kmer_size: 15,
    };

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(File::create(path)?)
    } else {
        Box::new(io::stdout())
    };

    if raw_features {
        write_raw_features_header(&mut writer)?;
    } else {
        TsvFormatter::write_header(&mut writer)?;
    }

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
        let mut q_seqs: Vec<Vec<u8>> = Vec::with_capacity(q_records.len());
        for (cid, record) in q_records.iter().enumerate() {
            q_seqs.push(record.sequence.clone());
            let mut record_tags: Vec<GenomeTag> = Vec::new();
            for enz in &enzyme_list {
                record_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz, cid));
            }
            all_q_tags.extend(record_tags);
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
            sequences: q_seqs,
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
                    let mut r_seqs: Vec<Vec<u8>> = Vec::with_capacity(r_records.len());
                    for (cid, record) in r_records.iter().enumerate() {
                        r_seqs.push(record.sequence.clone());
                        let mut record_tags: Vec<GenomeTag> = Vec::new();
                        for enz in &enzyme_list {
                            record_tags.extend(TagExtractor::extract_from_sequence(&record.sequence, enz, cid));
                        }
                        all_r_tags.extend(record_tags);
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
                        sequences: r_seqs,
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

                    Some((idx, r_path.clone(), r_tag_set.genome_id, ani_result, sv_count, r_tag_set.gc_content))
                })
                .collect()
        });

        // Sort by original index to maintain deterministic output order
        let mut sorted = results;
        sorted.sort_by_key(|(idx, _, _, _, _, _)| *idx);

        for (_idx, r_path, r_genome_id, ani_result, sv_count, ref_gc) in sorted {
            if raw_features {
                let shared = ani_result.local_ani_profile.len();
                let max_tags = (shared + ani_result.local_ani_profile.len()).max(1);
                let containment = shared as f64 / max_tags as f64;
                write_raw_features_record(
                    &mut writer,
                    &q_path.display().to_string(),
                    &r_path.display().to_string(),
                    &q_tag_set.genome_id,
                    &r_genome_id,
                    ani_result.raw_ani,
                    ani_result.mash_ani,
                    ani_result.chained_kmer_ani,
                    ani_result.af_query,
                    ani_result.af_reference,
                    shared,
                    containment,
                    ref_gc,
                    ani_result.ani,
                )?;
            } else {
                TsvFormatter::write_record(
                    &mut writer,
                    &q_path.display().to_string(),
                    &r_path.display().to_string(),
                    &q_tag_set.genome_id,
                    &r_genome_id,
                    &ani_result,
                    sv_count,
                )?;

                if ani_result.below_detection {
                    eprintln!(
                        "WARNING: {} vs {} — ANI below detection threshold (~83%). Result ({:.2}%) may be unreliable. Consider using a full-alignment tool (e.g., FastANI) for distant comparisons.",
                        q_tag_set.genome_id, r_genome_id, ani_result.ani * 100.0
                    );
                }
            }
        }
    }

    Ok(())
}

fn shared_tag_count(q: &TagSet, r: &TagSet) -> usize {
    let q_seqs: HashSet<u64> = q.tags.iter().map(|t| t.packed_sequence).collect();
    r.tags
        .iter()
        .filter(|t| q_seqs.contains(&t.packed_sequence))
        .count()
}
