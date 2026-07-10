use std::io::{Result as IoResult, Write};
use serde::Serialize;

use crate::core::{AniResult, MatchResult, StructuralVariation, SyntenyBlock};

/// skani-compatible TSV formatter.
pub struct TsvFormatter;

impl TsvFormatter {
    /// Write the TSV header line.
    pub fn write_header<W: Write>(writer: &mut W) -> IoResult<()> {
        writeln!(
            writer,
            "query_file\tref_file\tani\taf_q\taf_r\tquery_name\tref_name\tshared_tags\tsv_count"
        )
    }

    /// Write a single TSV record.
    pub fn write_record<W: Write>(
        writer: &mut W,
        query_file: &str,
        ref_file: &str,
        query_name: &str,
        ref_name: &str,
        ani_result: &AniResult,
        sv_count: usize,
    ) -> IoResult<()> {
        writeln!(
            writer,
            "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}",
            query_file,
            ref_file,
            ani_result.ani,
            ani_result.af_query,
            ani_result.af_reference,
            query_name,
            ref_name,
            ani_result.local_ani_profile.len(),
            sv_count,
        )
    }
}

/// Extended TSV formatter with structural variation columns.
pub struct ExtendedTsvFormatter;

impl ExtendedTsvFormatter {
    /// Write the extended TSV header line.
    pub fn write_header<W: Write>(writer: &mut W) -> IoResult<()> {
        writeln!(
            writer,
            "query_file\tref_file\tani\taf_q\taf_r\tquery_name\tref_name\tshared_tags\tsv_count\trearrangements\tindels\tsynteny_blocks"
        )
    }

    /// Write a single extended TSV record.
    pub fn write_record<W: Write>(
        writer: &mut W,
        query_file: &str,
        ref_file: &str,
        query_name: &str,
        ref_name: &str,
        ani_result: &AniResult,
        sv_count: usize,
        rearrangements: usize,
        indels: usize,
        synteny_blocks: usize,
    ) -> IoResult<()> {
        writeln!(
            writer,
            "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            query_file,
            ref_file,
            ani_result.ani,
            ani_result.af_query,
            ani_result.af_reference,
            query_name,
            ref_name,
            ani_result.local_ani_profile.len(),
            sv_count,
            rearrangements,
            indels,
            synteny_blocks,
        )
    }
}

// --- JSON helpers ---

#[derive(Serialize, Debug, Clone)]
struct JsonSyntenyBlock {
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    matched_tags: usize,
    orientation: String,
    block_ani: f64,
}

impl From<&SyntenyBlock> for JsonSyntenyBlock {
    fn from(b: &SyntenyBlock) -> Self {
        Self {
            query_start: b.query_start,
            query_end: b.query_end,
            ref_start: b.ref_start,
            ref_end: b.ref_end,
            matched_tags: b.matched_tags,
            orientation: b.orientation.to_string(),
            block_ani: b.block_ani,
        }
    }
}

#[derive(Serialize, Debug, Clone)]
struct JsonStructuralVariation {
    sv_type: String,
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    size: isize,
    confidence: f64,
}

impl From<&StructuralVariation> for JsonStructuralVariation {
    fn from(sv: &StructuralVariation) -> Self {
        Self {
            sv_type: format!("{:?}", sv.sv_type),
            query_start: sv.query_start,
            query_end: sv.query_end,
            ref_start: sv.ref_start,
            ref_end: sv.ref_end,
            size: sv.size,
            confidence: sv.confidence,
        }
    }
}

#[derive(Serialize, Debug, Clone)]
struct JsonOutput {
    query: String,
    reference: String,
    ani: f64,
    af_query: f64,
    af_reference: f64,
    weighted_ani: f64,
    confidence: f64,
    shared_tags: usize,
    synteny_blocks: Vec<JsonSyntenyBlock>,
    structural_variations: Vec<JsonStructuralVariation>,
}

/// JSON output formatter.
pub struct JsonFormatter;

impl JsonFormatter {
    /// Format a comparison result as a pretty-printed JSON string.
    pub fn format(
        query: &str,
        reference: &str,
        ani_result: &AniResult,
        match_result: &MatchResult,
        svs: &[StructuralVariation],
    ) -> Result<String, serde_json::Error> {
        let output = JsonOutput {
            query: query.to_string(),
            reference: reference.to_string(),
            ani: ani_result.ani,
            af_query: ani_result.af_query,
            af_reference: ani_result.af_reference,
            weighted_ani: ani_result.weighted_ani,
            confidence: ani_result.confidence,
            shared_tags: match_result.matched_pairs.len(),
            synteny_blocks: match_result
                .synteny_blocks
                .iter()
                .map(JsonSyntenyBlock::from)
                .collect(),
            structural_variations: svs.iter().map(JsonStructuralVariation::from).collect(),
        };
        serde_json::to_string_pretty(&output)
    }
}

/// PAF-like formatter for synteny block visualization.
pub struct PafFormatter;

impl PafFormatter {
    /// Write the standard PAF header.
    pub fn write_header<W: Write>(writer: &mut W) -> IoResult<()> {
        writeln!(
            writer,
            "qname\tqlen\tqstart\tqend\tstrand\trname\trlen\trstart\trend\tnmatch\talen\tmapq"
        )
    }

    /// Write synteny blocks in PAF format.
    pub fn write_synteny_blocks<W: Write>(
        writer: &mut W,
        query_name: &str,
        query_len: usize,
        ref_name: &str,
        ref_len: usize,
        blocks: &[SyntenyBlock],
    ) -> IoResult<()> {
        for block in blocks {
            let strand = if block.orientation == '+' { "+" } else { "-" };
            let alen = block.query_end.saturating_sub(block.query_start);
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                query_name,
                query_len,
                block.query_start,
                block.query_end,
                strand,
                ref_name,
                ref_len,
                block.ref_start,
                block.ref_end,
                block.matched_tags,
                alen,
                60,
            )?;
        }
        Ok(())
    }
}

/// Convenience enum for selecting an output format at runtime.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutputFormat {
    Tsv,
    ExtendedTsv,
    Json,
    Paf,
}
