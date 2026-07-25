//! `syn2bani ani` — chain-restricted ANI by maximum likelihood.
//!
//! Deliberately separate from `dist`: no GBRT model, no empirical offset, no
//! `--mash-ani` variants. The estimate comes from fitting
//! [`crate::core::mle`]'s likelihood to chain-restricted tag outcomes, so there
//! is nothing to calibrate.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::core::chain_ani::{self, ChainAniConfig};
use crate::core::{GenomeTag, TagExtractor};
use crate::enzyme::{EnzymeConfig, EnzymeRegistry};

fn resolve_enzymes(registry: &EnzymeRegistry, spec: &str) -> Result<Vec<EnzymeConfig>> {
    spec.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|name| {
            registry
                .get(name)
                .with_context(|| format!("Unknown enzyme: {name}"))
                .map(|e| e.clone())
        })
        .collect()
}

/// Digest one genome with the whole enzyme panel into a single tag list.
fn digest_all(path: &Path, enzymes: &[EnzymeConfig]) -> Result<(Vec<GenomeTag>, usize, String)> {
    let mut tags: Vec<GenomeTag> = Vec::new();
    let mut total_length = 0usize;
    let mut genome_id = String::new();
    for (i, enzyme) in enzymes.iter().enumerate() {
        let set = TagExtractor::extract_from_fasta(path, enzyme)
            .with_context(|| format!("digesting {} with {}", path.display(), enzyme.name))?;
        if i == 0 {
            total_length = set.total_length;
            genome_id = set.genome_id.clone();
        }
        tags.extend(set.tags);
    }
    // Sort by (contig, position) so chaining and window lookups see genome order
    // regardless of which enzyme produced each tag.
    tags.sort_by_key(|t| (t.contig_id, t.position));
    Ok((tags, total_length, genome_id))
}

#[allow(clippy::too_many_arguments)]
pub fn run_ani(
    query: &[PathBuf],
    reference: &[PathBuf],
    enzymes_spec: &str,
    mismatch_tolerance: usize,
    min_chain_anchors: usize,
    max_gap: usize,
    threads: usize,
    parallel: bool,
    verbose: bool,
    output: Option<&Path>,
) -> Result<()> {
    let pool = crate::cli::build_pool(parallel, threads)?;
    let registry = EnzymeRegistry::new();
    let enzymes = resolve_enzymes(&registry, enzymes_spec)?;
    if enzymes.is_empty() {
        anyhow::bail!("no enzymes selected");
    }
    let geometry = chain_ani::geometry_from(&enzymes);
    let cfg = ChainAniConfig {
        mismatch_tolerance,
        min_chain_anchors,
        max_gap,
        ..Default::default()
    };

    // Digest every genome once, in parallel, then compare all query x reference.
    let all_paths: Vec<&PathBuf> = query.iter().chain(reference.iter()).collect();
    let digested: Vec<(Vec<GenomeTag>, usize, String)> = pool.install(|| {
        all_paths
            .par_iter()
            .map(|p| digest_all(p, &enzymes))
            .collect::<Result<Vec<_>>>()
    })?;
    let (q_sets, r_sets) = digested.split_at(query.len());

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    write!(out, "query\treference\tani\taf_query\taf_reference\tstd_err")?;
    if verbose {
        write!(
            out,
            "\tani_from_loss\tani_from_hist\tn_anchors\tn_chains\tn_tags\tflag"
        )?;
    }
    writeln!(out)?;

    for (q_tags, q_len, q_id) in q_sets.iter() {
        let rows: Vec<String> = pool.install(|| {
            r_sets
                .par_iter()
                .map(|(r_tags, r_len, r_id)| {
                    let res =
                        chain_ani::compute(q_tags, r_tags, &geometry, *q_len, *r_len, &cfg);
                    let mut line = format!(
                        "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.5}",
                        q_id,
                        r_id,
                        res.ani * 100.0,
                        res.af_query,
                        res.af_reference,
                        res.std_err * 100.0
                    );
                    if verbose {
                        line.push_str(&format!(
                            "\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                            res.ani_from_loss * 100.0,
                            res.ani_from_hist * 100.0,
                            res.n_anchors,
                            res.n_chains,
                            res.n_tags_in_chains,
                            if res.inconsistent { "INCONSISTENT" } else { "ok" }
                        ));
                    }
                    line
                })
                .collect()
        });
        for line in rows {
            writeln!(out, "{line}")?;
        }
    }
    out.flush()?;
    Ok(())
}
