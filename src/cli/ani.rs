//! `syn2bani ani` — chain-restricted ANI by maximum likelihood.
//!
//! Deliberately separate from `dist`: no GBRT model, no empirical offset, no
//! `--mash-ani` variants. The estimate comes from fitting
//! [`crate::core::mle`]'s likelihood to chain-restricted tag outcomes, so there
//! is nothing to calibrate.
//!
//! `ani` is the rate-heterogeneous estimate; `ani_uniform` is the single-rate
//! fit, kept alongside because the difference between them measures how mosaic
//! the two genomes are. `ani_gated` is the recommended raw estimate: the
//! heterogeneous fit, falling back to the single-rate fit when the two partial
//! estimators (loss rate vs mismatch histogram) disagree by more than 5 ANI
//! points — the regime where the gamma shape is not identifiable and the
//! heterogeneous fit overshoots. The `flag` column reports `BELOW_DETECTION`
//! (retention too low for any estimate), `INCONSISTENT` (gate fallback or
//! disrupted chain structure: >0.5 unconserved within-contig anchor
//! adjacencies per anchor), or `ok`. The `breakpoint_count` column itself is
//! the strict chain-transition count (chain breaks between chained anchors
//! only); the flag's statistic additionally tallies anchors that chaining
//! rejected, which is what its threshold was calibrated on.

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::core::chain_ani::{self, ChainAniConfig, Geometry};
use crate::core::{GenomeTag, TagExtractor, LinearCalModel, load_embedded_cal_model};
use crate::enzyme::{EnzymeConfig, EnzymeRegistry};
use crate::io::{parse_fasta, read_sketch};

/// Longest tag the 2-bit packing can hold. Tags longer than this are truncated,
/// and truncation is not reverse-complement symmetric.
const MAX_PACKABLE_TAG: usize = 32;

pub(crate) fn resolve_enzymes(registry: &EnzymeRegistry, spec: &str) -> Result<Vec<EnzymeConfig>> {
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
pub(crate) struct Digest {
    pub(crate) tags: Vec<GenomeTag>,
    pub(crate) total_length: usize,
    pub(crate) contig_lens: Vec<usize>,
    /// Contig names in `contig_id` order; `struct` needs them for PAF/TSV
    /// output, while `ani` only ever reports the genome id.
    pub(crate) contig_names: Vec<String>,
    pub(crate) genome_id: String,
}

/// Load a genome from a `.s2ba` sketch instead of digesting it.
///
/// The sketch carries tag positions, packed sequences, per-contig lengths and —
/// since format v2 — the enzyme panel itself, which is everything the pairwise
/// comparison needs. Digestion is by far the dominant cost when the same genome
/// is compared many times, so sketching once and reusing is the difference
/// between re-reading every FASTA per invocation and not.
pub(crate) fn load_sketch(path: &Path) -> Result<Digest> {
    let sk = read_sketch(path).with_context(|| format!("reading sketch {}", path.display()))?;
    if sk.enzymes.is_empty() {
        anyhow::bail!(
            "{} is a v1 sketch with no enzyme table, so its enzyme ids cannot be \
             interpreted; rebuild it with `syn2bani sketch`",
            path.display()
        );
    }

    let mut tags: Vec<GenomeTag> = Vec::with_capacity(sk.metadata.tag_count as usize);
    let mut contig_lens = Vec::with_capacity(sk.chromosomes.len());
    let mut contig_names = Vec::with_capacity(sk.chromosomes.len());
    for (cid, chrom) in sk.chromosomes.iter().enumerate() {
        contig_lens.push(chrom.length as usize);
        contig_names.push(chrom.name.clone());
        for st in &chrom.tags {
            let Some(e) = sk.enzymes.get(st.enzyme_id as usize) else {
                anyhow::bail!(
                    "{}: enzyme id {} is outside the sketch's {}-entry panel",
                    path.display(),
                    st.enzyme_id,
                    sk.enzymes.len()
                );
            };
            let seq_len = (e.tag_length as usize).min(32) as u8;
            // Unpack in place; a Vec per tag showed up in the profile.
            let mut sequence = [0u8; 32];
            for i in 0..seq_len as usize {
                sequence[i] = match (st.seq >> (2 * i)) & 0b11 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
            }
            tags.push(GenomeTag {
                position: st.position as usize,
                contig_id: cid,
                sequence,
                packed_sequence: st.seq,
                seq_len,
                direction: if st.direction == 0 { '+' } else { '-' },
                enzyme: e.name.clone(),
            });
        }
    }
    tags.sort_by_key(|t| (t.contig_id, t.position));
    Ok(Digest {
        tags,
        total_length: sk.metadata.total_length as usize,
        contig_lens,
        contig_names,
        genome_id: sk.genome_id,
    })
}

/// Geometry taken from a sketch's own enzyme table, so a sketched run does not
/// depend on `--enzymes` matching what the sketch was built with.
fn geometry_from_sketches(digests: &[Digest], base: &Geometry) -> Geometry {
    let mut g = base.clone();
    for d in digests {
        for t in &d.tags {
            g.entry(t.enzyme.clone()).or_insert((32, 6));
        }
    }
    g
}

/// Digest one genome with the whole enzyme panel into a single tag list.
pub(crate) fn digest_all(path: &Path, enzymes: &[EnzymeConfig]) -> Result<Digest> {
    // The FASTA is parsed ONCE and every enzyme digests the same in-memory
    // records. Going through `extract_from_fasta` per enzyme re-read and
    // re-cloned the whole genome once per enzyme, so a four-enzyme panel paid
    // four times the I/O and allocation for one usable copy. That dominated both
    // runtime and peak memory on multi-pair runs. Contig sequences are never
    // retained here; only their lengths are needed downstream.
    let records = parse_fasta(path)
        .map_err(|e| anyhow::anyhow!("{e}"))
        .with_context(|| format!("parsing {}", path.display()))?;
    if records.is_empty() {
        anyhow::bail!("{} contains no sequences", path.display());
    }

    let genome_id = records[0]
        .id
        .split_whitespace()
        .next()
        .unwrap_or("unknown")
        .to_string();
    let contig_lens: Vec<usize> = records.iter().map(|r| r.sequence.len()).collect();
    let contig_names: Vec<String> = records
        .iter()
        .map(|r| r.id.split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();
    let total_length: usize = contig_lens.iter().sum();

    let mut tags: Vec<GenomeTag> = Vec::new();
    for enzyme in enzymes {
        for (cid, record) in records.iter().enumerate() {
            tags.extend(TagExtractor::extract_from_sequence(
                &record.sequence,
                enzyme,
                cid,
            ));
        }
    }

    // Sort by (contig, position) so chaining and window lookups see genome order
    // regardless of which enzyme produced each tag.
    tags.sort_by_key(|t| (t.contig_id, t.position));
    Ok(Digest {
        tags,
        total_length,
        contig_lens,
        contig_names,
        genome_id,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn run_ani(
    positional: &[PathBuf],
    ql: Option<&Path>,
    rl: Option<&Path>,
    enzymes_spec: &str,
    mismatch_tolerance: usize,
    min_chain_anchors: usize,
    max_gap: usize,
    threads: usize,
    parallel: bool,
    verbose: bool,
    strata_out: Option<&Path>,
    calibrate: bool,
    calibrate_model: Option<&Path>,
    output: Option<&Path>,
) -> Result<()> {
    // Either both list files, or the positional "queries... reference" form.
    let (query, reference): (Vec<PathBuf>, Vec<PathBuf>) = match (ql, rl) {
        (Some(q), Some(r)) => (
            crate::cli::compare::read_path_list(q)?,
            crate::cli::compare::read_path_list(r)?,
        ),
        (None, None) => {
            if positional.len() < 2 {
                anyhow::bail!(
                    "need at least one query and one reference; pass \
                     `<queries...> <reference>` or use --ql/--rl"
                );
            }
            let (qs, rs) = positional.split_at(positional.len() - 1);
            (qs.to_vec(), rs.to_vec())
        }
        _ => anyhow::bail!("--ql and --rl must be given together"),
    };
    if query.is_empty() || reference.is_empty() {
        anyhow::bail!("empty query or reference list");
    }

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

    // Load calibration model if requested.
    let cal_model: Option<LinearCalModel> = if calibrate {
        if let Some(path) = calibrate_model {
            let json = std::fs::read_to_string(path)
                .with_context(|| format!("reading calibration model {}", path.display()))?;
            Some(LinearCalModel::from_json(&json)?)
        } else {
            Some(load_embedded_cal_model())
        }
    } else {
        None
    };

    // Load every genome once, in parallel, then compare all query x reference.
    // A `.s2ba` input skips digestion entirely.
    let all_paths: Vec<&PathBuf> = query.iter().chain(reference.iter()).collect();
    let digested: Vec<Digest> = pool.install(|| {
        all_paths
            .par_iter()
            .map(|p| {
                if p.extension().is_some_and(|e| e == "s2ba") {
                    load_sketch(p)
                } else {
                    digest_all(p, &enzymes)
                }
            })
            .collect::<Result<Vec<_>>>()
    })?;
    // Sketches are self-describing, so honour their panel rather than the flag.
    let geometry = geometry_from_sketches(&digested, &geometry);
    let (q_sets, r_sets) = digested.split_at(query.len());

    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    writeln!(
        out,
        "{}",
        crate::cli::compare::ani_header(cal_model.is_some(), verbose)
    )?;

    // Parallelise over the whole query x reference product, not over references
    // within each query. The common shape is many queries against one reference,
    // and an inner-only split leaves that case serial.
    let pairs: Vec<(usize, usize)> = (0..q_sets.len())
        .flat_map(|qi| (0..r_sets.len()).map(move |ri| (qi, ri)))
        .collect();

    // Per-enzyme sufficient statistics, optionally dumped so that any sub-panel
    // can be re-scored later without touching a genome again. The likelihood is
    // a plain sum over strata, so `estimate(subset)` is exact arithmetic on this
    // table — see `syn2bani panel`.
    let strata_wanted = strata_out.is_some();
    let cal_model = cal_model.clone();

    let collected: Vec<(String, Vec<String>)> = pool.install(|| {
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
                let cal = cal_model.as_ref().map(|m| m.predict_from_result(&res));
                let line =
                    crate::cli::compare::ani_row(&q.genome_id, &r.genome_id, &res, cal, verbose);
                let sl: Vec<String> = if strata_wanted {
                    res.strata
                        .iter()
                        .map(|st| {
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                q.genome_id,
                                r.genome_id,
                                st.enzyme,
                                st.tag_len,
                                st.body_len,
                                st.n_miss,
                                st.hist
                                    .iter()
                                    .map(|c| c.to_string())
                                    .collect::<Vec<_>>()
                                    .join(",")
                            )
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                (line, sl)
            })
            .collect()
    });
    for (line, _) in &collected {
        writeln!(out, "{line}")?;
    }

    if let Some(path) = strata_out {
        let mut sf = BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        );
        writeln!(sf, "query\treference\tenzyme\ttag_len\tbody_len\tn_miss\thist")?;
        for (_, sl) in &collected {
            for row in sl {
                writeln!(sf, "{row}")?;
            }
        }
        sf.flush()?;
    }
    out.flush()?;
    Ok(())
}
