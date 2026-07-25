//! `syn2bani ani` — chain-restricted ANI by maximum likelihood.
//!
//! Deliberately separate from `dist`: no GBRT model, no empirical offset, no
//! `--mash-ani` variants. The estimate comes from fitting
//! [`crate::core::mle`]'s likelihood to chain-restricted tag outcomes, so there
//! is nothing to calibrate.
//!
//! `ani` is the rate-heterogeneous estimate and is the one to use.
//! `ani_uniform` is the single-rate fit, kept alongside because the difference
//! between them measures how mosaic the two genomes are.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::core::chain_ani::{self, ChainAniConfig};
use crate::core::{GenomeTag, TagExtractor};
use crate::enzyme::{EnzymeConfig, EnzymeRegistry};

/// Longest tag the 2-bit packing can hold. Tags longer than this are truncated,
/// and truncation is not reverse-complement symmetric.
const MAX_PACKABLE_TAG: usize = 32;

fn resolve_enzymes(registry: &EnzymeRegistry, spec: &str) -> Result<Vec<EnzymeConfig>> {
    let requested: Vec<EnzymeConfig> = spec
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|name| {
            registry
                .get(name)
                .with_context(|| format!("Unknown enzyme: {name}"))
                .map(|e| e.clone())
        })
        .collect::<Result<_>>()?;

    // Drop enzymes whose tags cannot be packed losslessly.
    //
    // Tags are matched through a strand-canonical 2-bit packing capped at 32
    // bases. Keeping the first 32 bases of a longer tag is not symmetric under
    // reverse complement: a tag read from the forward strand keeps bases 0..31,
    // while the same locus read from the other strand keeps what corresponds to
    // bases 1..32. The canonical forms then disagree and the tags cannot match.
    //
    // Measured with CspCI (33 bp): comparing a genome against its own reverse
    // complement — identical sequence, so the answer must be 100% — returned
    // 98.81% and lost 91% of that enzyme's anchors. Dropping it returns 99.9999%.
    // Supporting longer tags needs the packing widened past one u64.
    let (ok, too_long): (Vec<_>, Vec<_>) = requested
        .into_iter()
        .partition(|e| e.tag_length <= MAX_PACKABLE_TAG);
    for e in &too_long {
        eprintln!(
            "warning: skipping {} — its {} bp tag exceeds the {} bp packing limit, \
             and truncating it breaks reverse-complement symmetry",
            e.name, e.tag_length, MAX_PACKABLE_TAG
        );
    }
    Ok(ok)
}

/// A digested genome: tags, total length, per-contig lengths, and its id.
struct Digest {
    tags: Vec<GenomeTag>,
    total_length: usize,
    contig_lens: Vec<usize>,
    genome_id: String,
}

/// Digest one genome with the whole enzyme panel into a single tag list.
fn digest_all(path: &Path, enzymes: &[EnzymeConfig]) -> Result<Digest> {
    let mut tags: Vec<GenomeTag> = Vec::new();
    let mut total_length = 0usize;
    let mut contig_lens: Vec<usize> = Vec::new();
    let mut genome_id = String::new();
    for (i, enzyme) in enzymes.iter().enumerate() {
        let set = TagExtractor::extract_from_fasta(path, enzyme)
            .with_context(|| format!("digesting {} with {}", path.display(), enzyme.name))?;
        if i == 0 {
            total_length = set.total_length;
            genome_id = set.genome_id.clone();
            // Contig lengths let chain spans be clamped to the contig they came
            // from, which matters on fragmented assemblies.
            contig_lens = set.sequences.iter().map(|c| c.len()).collect();
        }
        tags.extend(set.tags);
    }
    // Sort by (contig, position) so chaining and window lookups see genome order
    // regardless of which enzyme produced each tag.
    tags.sort_by_key(|t| (t.contig_id, t.position));
    Ok(Digest {
        tags,
        total_length,
        contig_lens,
        genome_id,
    })
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
    let digested: Vec<Digest> = pool.install(|| {
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

    write!(out, "query\treference\tani\tani_uniform\taf_query\taf_reference\tstd_err")?;
    if verbose {
        write!(
            out,
            "\thet_shape\tretention\tani_from_loss\tani_from_hist\tn_anchors\tn_chains\tn_tags\tflag"
        )?;
    }
    writeln!(out)?;

    // Parallelise over the whole query x reference product, not over references
    // within each query. The common shape is many queries against one reference,
    // and an inner-only split leaves that case serial.
    let pairs: Vec<(usize, usize)> = (0..q_sets.len())
        .flat_map(|qi| (0..r_sets.len()).map(move |ri| (qi, ri)))
        .collect();

    let rows: Vec<String> = pool.install(|| {
        pairs
            .par_iter()
            .map(|&(qi, ri)| {
                let q = &q_sets[qi];
                let r = &r_sets[ri];
                let res = chain_ani::compute(
                    &q.tags,
                    &r.tags,
                    &geometry,
                    q.total_length,
                    r.total_length,
                    &q.contig_lens,
                    &r.contig_lens,
                    &cfg,
                );
                let mut line = format!(
                    "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.5}",
                    q.genome_id,
                    r.genome_id,
                    res.ani_het * 100.0,
                    res.ani * 100.0,
                    res.af_query,
                    res.af_reference,
                    res.std_err * 100.0
                );
                if verbose {
                    line.push_str(&format!(
                        "\t{:.3}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
                        res.het_shape,
                        res.retention,
                        res.ani_from_loss * 100.0,
                        res.ani_from_hist * 100.0,
                        res.n_anchors,
                        res.n_chains,
                        res.n_tags_in_chains,
                        if res.below_detection {
                            "BELOW_DETECTION"
                        } else if res.inconsistent {
                            "INCONSISTENT"
                        } else {
                            "ok"
                        }
                    ));
                }
                line
            })
            .collect()
    });
    for line in rows {
        writeln!(out, "{line}")?;
    }
    out.flush()?;
    Ok(())
}
