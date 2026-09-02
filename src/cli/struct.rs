//! `syn2bani struct` — structural variation calls from collinear tag chains.
//!
//! This subcommand runs the same chain-restricted pipeline as `ani`
//! ([`crate::core::chain_ani`]) and derives SV calls from the final
//! adaptive-pass chains ([`crate::core::sv`]). It deliberately does not use
//! the v7 `TagMatcher`/`SyntenyBuilder`/`StructureAnalyzer` path: that
//! chaining re-admitted non-collinear anchors and let single chains span a
//! whole genome, so its rearrangement output was not trustworthy.
//!
//! Input is plain FASTA (no sketch support — SV calling is a one-shot
//! pairwise analysis, so the digestion cost is paid once per genome anyway).

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::ani::{digest_all, resolve_enzymes, Digest};
use crate::core::chain_ani::{self, ChainAniConfig, ChainBlock};
use crate::core::sv::{self, SvCall};
use crate::core::SvType;
use crate::enzyme::EnzymeRegistry;

/// Same panel default as `ani`.
const DEFAULT_ENZYMES: &str = "BcgI,AlfI,AloI,FalI";

/// One PAF line per chain. There is no base-level alignment to draw these
/// from, so the last three columns are approximations:
///
/// - `nmatch` is estimated anchor coverage: `n_anchors * mean tag length`,
///   capped at the span. Tags are fixed 25–33 bp windows, so this bounds the
///   bases direct anchor evidence covers.
/// - `alnlen` is the chain span (the longer of the query/reference spans).
/// - `mapq` is a pseudo mapping quality from anchor density: 10 per
///   anchor/kb, capped at 60. A dense chain saturates at 60; sparse chains
///   score lower. It carries no probabilistic meaning.
fn write_paf_line<W: Write>(
    out: &mut W,
    q: &Digest,
    r: &Digest,
    c: &ChainBlock,
    mean_tag_len: usize,
) -> Result<()> {
    let q_span = c.q_end - c.q_start;
    let r_span = c.r_end - c.r_start;
    let alnlen = q_span.max(r_span);
    let nmatch = (c.n_anchors * mean_tag_len).min(alnlen);
    let density = c.n_anchors as f64 / (alnlen.max(1) as f64 / 1_000.0);
    let mapq = (density * 10.0).min(60.0) as u32;
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        q.contig_names[c.q_contig],
        q.contig_lens[c.q_contig],
        c.q_start,
        c.q_end,
        c.orientation,
        r.contig_names[c.r_contig],
        r.contig_lens[c.r_contig],
        c.r_start,
        c.r_end,
        nmatch,
        alnlen,
        mapq,
    )?;
    Ok(())
}

fn write_sv_line<W: Write>(out: &mut W, q: &Digest, r: &Digest, s: &SvCall) -> Result<()> {
    writeln!(
        out,
        "{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        q.genome_id,
        r.genome_id,
        s.sv_type,
        q.contig_names[s.q_contig],
        s.q_start,
        s.q_end,
        r.contig_names[s.r_contig],
        s.r_start,
        s.r_end,
        s.size,
        s.support_left,
        s.support_right,
    )?;
    Ok(())
}

fn sv_type_abbr(t: SvType) -> &'static str {
    match t {
        SvType::Inversion => "INV",
        SvType::Insertion => "INS",
        SvType::Deletion => "DEL",
        SvType::Translocation => "TRA",
        SvType::Duplication => "DUP",
    }
}

/// BED9 colour by SV type (RGB).
fn sv_type_rgb(t: SvType) -> &'static str {
    match t {
        SvType::Inversion => "255,0,0",
        SvType::Insertion => "0,128,0",
        SvType::Deletion => "0,0,255",
        SvType::Translocation => "255,165,0",
        SvType::Duplication => "128,0,128",
    }
}

/// Infer strand for an inversion from the relative orientation of the query
/// and reference intervals. Reverse-orientation mapping is reported on the
/// `-` strand; `.` is used for non-inversions and ambiguous cases.
fn inversion_strand(s: &SvCall) -> char {
    if !matches!(s.sv_type, SvType::Inversion) {
        return '.';
    }
    let q_dir = (s.q_end as i64).saturating_sub(s.q_start as i64).signum();
    let r_dir = (s.r_end as i64).saturating_sub(s.r_start as i64).signum();
    if q_dir == 0 || r_dir == 0 || q_dir == r_dir {
        // An inversion should have opposite orientation; same sign is ambiguous.
        return '.';
    }
    // Query interval progresses forward while reference interval progresses
    // backward → the inverted segment maps to the reverse reference strand.
    '-'
}

fn write_bed_line<W: Write>(out: &mut W, r: &Digest, s: &SvCall) -> Result<()> {
    let start = s.r_start.min(s.r_end);
    let end = s.r_start.max(s.r_end);
    let name = format!("{}_l{}_r{}", sv_type_abbr(s.sv_type), s.support_left, s.support_right);
    let score = s.support_left + s.support_right;
    let strand = inversion_strand(s);
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        r.contig_names[s.r_contig],
        start,
        end,
        name,
        score,
        strand,
        start,
        end,
        sv_type_rgb(s.sv_type),
    )?;
    Ok(())
}

/// Handler for the `struct` subcommand.
///
/// Pairwise query × reference structural variation detection. With `--paf`,
/// emits one PAF record per chain; otherwise one TSV record per SV. The
/// `--rearrangement` / `--indel` flags filter the SV types reported (both
/// unset reports everything); `indel_min` is the smallest offset jump called
/// as an indel.
#[allow(clippy::too_many_arguments)]
pub fn run_struct(
    query: &[PathBuf],
    reference: &[PathBuf],
    output: Option<&Path>,
    paf: bool,
    bed: bool,
    rearrangement: bool,
    indel: bool,
    multi_enzyme: bool,
    enzymes: Option<&str>,
    indel_min: usize,
) -> Result<()> {
    let registry = EnzymeRegistry::new();
    if paf && bed {
        anyhow::bail!("--paf and --bed are mutually exclusive");
    }

    let spec: String = match (enzymes, multi_enzyme) {
        (Some(e), _) => e.to_string(),
        (None, true) => registry
            .all()
            .iter()
            .map(|e| e.name.clone())
            .collect::<Vec<_>>()
            .join(","),
        (None, false) => DEFAULT_ENZYMES.to_string(),
    };
    let enzyme_list = resolve_enzymes(&registry, &spec)?;
    if enzyme_list.is_empty() {
        anyhow::bail!("no enzymes selected");
    }
    let geometry = chain_ani::geometry_from(&enzyme_list);
    let mean_tag_len = enzyme_list
        .iter()
        .map(|e| e.tag_length.min(32))
        .sum::<usize>()
        / enzyme_list.len();
    let cfg = ChainAniConfig::default();

    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };

    if !paf && !bed {
        writeln!(
            writer,
            "query\treference\tsv_type\tq_contig\tq_start\tq_end\tr_contig\tr_start\tr_end\tsize\tsupport_left\tsupport_right"
        )?;
    }

    let q_digests: Vec<Digest> = query
        .iter()
        .map(|p| digest_all(p, &enzyme_list))
        .collect::<Result<_>>()?;
    let r_digests: Vec<Digest> = reference
        .iter()
        .map(|p| digest_all(p, &enzyme_list))
        .collect::<Result<_>>()?;

    for q in &q_digests {
        for r in &r_digests {
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

            if paf {
                for c in &res.chains {
                    write_paf_line(&mut writer, q, r, c, mean_tag_len)?;
                }
                eprintln!(
                    "{}\t{}\t{} chains",
                    q.genome_id,
                    r.genome_id,
                    res.chains.len()
                );
                continue;
            }

            let svs = sv::detect(&res.chains, indel_min);
            // Both filter flags unset means report everything.
            let keep = |s: &SvCall| {
                if !rearrangement && !indel {
                    return true;
                }
                match s.sv_type {
                    SvType::Inversion | SvType::Translocation => rearrangement,
                    SvType::Insertion | SvType::Deletion => indel,
                    SvType::Duplication => false,
                }
            };
            let mut n_reported = 0usize;
            for s in svs.iter().filter(|s| keep(s)) {
                if bed {
                    write_bed_line(&mut writer, r, s)?;
                } else {
                    write_sv_line(&mut writer, q, r, s)?;
                }
                n_reported += 1;
            }
            eprintln!(
                "{}\t{}\t{} chains\t{} SVs reported\tani {:.4}",
                q.genome_id,
                r.genome_id,
                res.chains.len(),
                n_reported,
                res.ani_het * 100.0
            );
        }
    }

    writer.flush()?;
    Ok(())
}
