use crate::core::tag_extractor::GenomeTag;
use crate::core::tag_matcher::MatchResult;
use crate::core::synteny_builder::SyntenyBlock;

/// Type of structural variation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SvType {
    Inversion,
    Insertion,
    Deletion,
    Translocation,
    Duplication,
}

/// A structural variation event between query and reference.
#[derive(Debug, Clone)]
pub struct StructuralVariation {
    pub sv_type: SvType,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub size: isize,
    pub confidence: f64,
}

/// Detects structural variations from synteny blocks and unmatched tags.
pub struct StructureAnalyzer;

impl StructureAnalyzer {
    /// Detect rearrangements (inversions, translocations) by inspecting block boundaries.
    pub fn detect_rearrangements(blocks: &[SyntenyBlock]) -> Vec<StructuralVariation> {
        let mut svs = Vec::new();

        for i in 1..blocks.len() {
            let prev = &blocks[i - 1];
            let curr = &blocks[i];

            if curr.orientation != prev.orientation {
                svs.push(StructuralVariation {
                    sv_type: SvType::Inversion,
                    query_start: prev.query_end,
                    query_end: curr.query_start,
                    ref_start: prev.ref_end,
                    ref_end: curr.ref_start,
                    size: (curr.query_end.saturating_sub(prev.query_start)) as isize,
                    confidence: 0.85,
                });
            }

            if curr.ref_start > prev.ref_end + 100000 {
                svs.push(StructuralVariation {
                    sv_type: SvType::Translocation,
                    query_start: prev.query_end,
                    query_end: curr.query_start,
                    ref_start: prev.ref_end,
                    ref_end: curr.ref_start,
                    size: (curr.ref_start.saturating_sub(prev.ref_end)) as isize,
                    confidence: 0.7,
                });
            }
        }

        svs
    }

    /// Detect indels from unmatched query and reference tags.
    pub fn detect_indels(match_result: &MatchResult) -> Vec<StructuralVariation> {
        let mut svs = Vec::new();

        // Insertions from unmatched query tags
        if !match_result.unmatched_query.is_empty() {
            let clusters = Self::cluster_tags(&match_result.unmatched_query, 1000);
            for cluster in clusters {
                let size = cluster.len() as isize * 32; // Approximate tag size
                svs.push(StructuralVariation {
                    sv_type: SvType::Insertion,
                    query_start: cluster.first().map(|t| t.position).unwrap_or(0),
                    query_end: cluster.last().map(|t| t.position).unwrap_or(0),
                    ref_start: 0,
                    ref_end: 0,
                    size,
                    confidence: 0.6,
                });
            }
        }

        // Deletions from unmatched reference tags
        if !match_result.unmatched_ref.is_empty() {
            let clusters = Self::cluster_tags(&match_result.unmatched_ref, 1000);
            for cluster in clusters {
                let size = cluster.len() as isize * 32;
                svs.push(StructuralVariation {
                    sv_type: SvType::Deletion,
                    query_start: 0,
                    query_end: 0,
                    ref_start: cluster.first().map(|t| t.position).unwrap_or(0),
                    ref_end: cluster.last().map(|t| t.position).unwrap_or(0),
                    size,
                    confidence: 0.6,
                });
            }
        }

        svs
    }

    fn cluster_tags(tags: &[GenomeTag], max_gap: usize) -> Vec<Vec<GenomeTag>> {
        if tags.is_empty() {
            return Vec::new();
        }

        let mut sorted: Vec<GenomeTag> = tags.iter().cloned().collect();
        sorted.sort_by_key(|t| t.position);

        let mut clusters = Vec::new();
        let mut current = vec![sorted[0].clone()];

        for tag in sorted.iter().skip(1) {
            if tag.position - current.last().unwrap().position <= max_gap {
                current.push(tag.clone());
            } else {
                clusters.push(current);
                current = vec![tag.clone()];
            }
        }
        clusters.push(current);

        clusters
    }

    /// Produce a simplified PAF-like representation of SVs.
    pub fn to_paf(svs: &[StructuralVariation]) -> String {
        let mut lines = Vec::new();
        lines.push(
            "qname\tqlen\tqstart\tqend\tstrand\trname\trlen\trstart\trend\tnmatch\talen\tmapq\tsv_type\tsv_size"
                .to_string(),
        );

        for sv in svs {
            let strand = if sv.sv_type == SvType::Inversion {
                '-'
            } else {
                '+'
            };
            let line = format!(
                "query\t0\t{}\t{}\t{}\tref\t0\t{}\t{}\t0\t{}\t0\t{:?}\t{}",
                sv.query_start,
                sv.query_end,
                strand,
                sv.ref_start,
                sv.ref_end,
                sv.size.abs(),
                sv.sv_type,
                sv.size
            );
            lines.push(line);
        }

        lines.join("\n")
    }
}
