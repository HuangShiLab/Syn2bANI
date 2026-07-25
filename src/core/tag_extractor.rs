use crate::utils::fxhash::FastHashMap;
use rayon::prelude::*;
use std::path::Path;
use std::sync::Arc;

use thiserror::Error;
use crate::enzyme::{digest_sequence, EnzymeConfig};
use crate::io::fasta_parser::parse_fasta;

/// A single genome tag extracted after enzyme digestion.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomeTag {
    pub position: usize,
    /// Contig index within the source genome (for multi-contig FASTA/MAGs).
    pub contig_id: usize,
    /// Tag sequence as a fixed 32-byte array (zero-padded, actual length may vary).
    pub sequence: [u8; 32],
    /// 2-bit packed sequence (64 bits = 32 bp), aligned with `sequence`.
    pub packed_sequence: u64,
    /// Actual sequence length (tag may be shorter than 32 bp).
    pub seq_len: u8,
    pub direction: char,
    pub enzyme: String,
}

impl GenomeTag {
    /// Strand-canonical packed sequence, for hashing tags across orientations.
    ///
    /// Two genomes that share a locus store the same tag sequence only when
    /// they read it from the same strand. Inversions and arbitrarily oriented
    /// draft contigs break that, so index on this rather than
    /// `packed_sequence`.
    #[inline]
    pub fn canonical(&self) -> u64 {
        canonical_packed(self.packed_sequence, self.seq_len)
    }

    /// Reverse complement of this tag's packed sequence.
    #[inline]
    pub fn packed_revcomp(&self) -> u64 {
        revcomp_packed(self.packed_sequence, self.seq_len)
    }
}

/// A collection of tags from a single genome/contig.
#[derive(Debug, Clone)]
pub struct TagSet {
    pub genome_id: String,
    pub chromosome: String,
    pub tags: Vec<GenomeTag>,
    pub total_length: usize,
    pub gc_content: f64,
    /// Raw contig sequences, indexed by `contig_id`. May be empty for sketch/database modes.
    pub sequences: Vec<Vec<u8>>,
}

/// Tag sets from multiple enzymes for the same genome.
#[derive(Debug, Clone)]
pub struct MultiEnzymeTagSet {
    pub sets: FastHashMap<String, TagSet>,
}

/// Errors that can occur during tag extraction.
#[derive(Error, Debug)]
pub enum ExtractError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid FASTA format: {0}")]
    InvalidFasta(String),
    #[error("Invalid enzyme configuration: {0}")]
    InvalidEnzyme(String),
}

/// Extracts tags from raw sequences or FASTA files.
pub struct TagExtractor;

impl TagExtractor {
    /// Extract tags from a raw sequence slice using the given enzyme configuration.
    pub fn extract_from_sequence(seq: &[u8], enzyme: &EnzymeConfig, contig_id: usize) -> Vec<GenomeTag> {
        let digested = digest_sequence(seq, enzyme);
        digested.into_iter().map(|tag| {
            let packed = pack_bytes(&tag.sequence, tag.seq_len);
            GenomeTag {
                position: tag.position,
                contig_id,
                sequence: tag.sequence,
                packed_sequence: packed,
                seq_len: tag.seq_len,
                direction: match tag.direction {
                    crate::enzyme::Direction::Forward => '+',
                    crate::enzyme::Direction::Reverse => '-',
                },
                enzyme: enzyme.name.clone(),
            }
        }).collect()
    }

    /// Extract tags from a single FASTA file (multi-contig aware).
    pub fn extract_from_fasta(path: &Path, enzyme: &EnzymeConfig) -> Result<TagSet, ExtractError> {
        let records = parse_fasta(path).map_err(|e| ExtractError::InvalidFasta(e.to_string()))?;
        if records.is_empty() {
            return Err(ExtractError::InvalidFasta("Empty file".to_string()));
        }

        let genome_id = records[0]
            .id
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        let mut total_length = 0usize;
        let mut gc_count_total = 0usize;
        let mut tags = Vec::new();
        let mut chromosome_names = Vec::new();

        for (contig_id, record) in records.iter().enumerate() {
            let len = record.sequence.len();
            total_length += len;
            gc_count_total += record
                .sequence
                .iter()
                .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                .count();
            let mut contig_tags = Self::extract_from_sequence(&record.sequence, enzyme, contig_id);
            tags.append(&mut contig_tags);
            chromosome_names.push(record.id.clone());
        }

        let gc_content = gc_count_total as f64 / total_length.max(1) as f64;
        let chromosome = if chromosome_names.len() == 1 {
            chromosome_names[0].clone()
        } else {
            "multi".to_string()
        };
        let sequences: Vec<Vec<u8>> = records.iter().map(|r| r.sequence.clone()).collect();

        Ok(TagSet {
            genome_id,
            chromosome,
            tags,
            total_length,
            gc_content,
            sequences,
        })
    }

    /// Extract tags using multiple enzymes and return a map keyed by enzyme name.
    ///
    /// This is the sequential fallback. For parallel digestion, use
    /// [`extract_multi_enzyme_par`] instead.
    pub fn extract_multi_enzyme(
        path: &Path,
        enzymes: &[EnzymeConfig],
    ) -> Result<MultiEnzymeTagSet, ExtractError> {
        let mut sets = FastHashMap::default();
        for enzyme in enzymes {
            let tag_set = Self::extract_from_fasta(path, enzyme)?;
            let name = enzyme.name.clone();
            sets.insert(name, tag_set);
        }
        Ok(MultiEnzymeTagSet { sets })
    }

    /// Parallel multi-enzyme tag extraction using Rayon.
    ///
    /// The FASTA file is read once into memory, then each enzyme is digested
    /// in parallel across the available CPU cores. This is significantly faster
    /// than the sequential [`extract_multi_enzyme`] when the enzyme panel
    /// contains more than a few enzymes.
    ///
    /// # Performance
    /// On an 8-core machine, the 16-enzyme 2bRAD-M panel is typically
    /// 4–6× faster than sequential digestion.
    pub fn extract_multi_enzyme_par(
        path: &Path,
        enzymes: &[EnzymeConfig],
    ) -> Result<MultiEnzymeTagSet, ExtractError> {
        // Step 1: read the file once (single-threaded I/O)
        let records = parse_fasta(path).map_err(|e| ExtractError::InvalidFasta(e.to_string()))?;
        if records.is_empty() {
            return Err(ExtractError::InvalidFasta("Empty file".to_string()));
        }

        let genome_id = records[0]
            .id
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        let total_length: usize = records.iter().map(|r| r.sequence.len()).sum();
        let gc_count_total: usize = records
            .iter()
            .map(|r| {
                r.sequence
                    .iter()
                    .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
                    .count()
            })
            .sum();
        let gc_content = gc_count_total as f64 / total_length.max(1) as f64;
        let chromosome = if records.len() == 1 {
            records[0].id.clone()
        } else {
            "multi".to_string()
        };
        let sequences: Vec<Vec<u8>> = records.iter().map(|r| r.sequence.clone()).collect();
        let records_arc = Arc::new(records);

        // Step 2: parallel digestion across enzymes
        let sets: Vec<_> = enzymes
            .par_iter()
            .map(|enzyme| {
                let mut tags = Vec::new();
                for (contig_id, record) in records_arc.iter().enumerate() {
                    let mut contig_tags =
                        Self::extract_from_sequence(&record.sequence, enzyme, contig_id);
                    tags.append(&mut contig_tags);
                }
                let tag_set = TagSet {
                    genome_id: genome_id.clone(),
                    chromosome: chromosome.clone(),
                    tags,
                    total_length,
                    gc_content,
                    sequences: sequences.clone(),
                };
                (enzyme.name.clone(), tag_set)
            })
            .collect();

        let mut map = FastHashMap::default();
        for (name, set) in sets {
            map.insert(name, set);
        }

        Ok(MultiEnzymeTagSet { sets: map })
    }
}

/// Pack a DNA sequence (up to 32 bp) into a 64-bit integer using 2-bit encoding.
///
/// Encoding: A/a=0b00, C/c=0b01, G/g=0b10, T/t=0b11.
/// Only the first `len` bases are packed; remaining bits are zero.
#[inline]
pub fn pack_bytes(seq: &[u8; 32], len: u8) -> u64 {
    let mut packed: u64 = 0;
    let n = (len as usize).min(32);
    for i in 0..n {
        let bits = match seq[i] {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00,
        };
        packed |= (bits as u64) << (i * 2);
    }
    packed
}

/// Reverse complement of a 2-bit packed sequence of `len` bases.
///
/// Complement is bitwise NOT under the A=00/C=01/G=10/T=11 encoding, and the
/// base order is reversed.
pub fn revcomp_packed(packed: u64, len: u8) -> u64 {
    let n = (len as usize).min(32);
    let mut out: u64 = 0;
    for i in 0..n {
        let base = (packed >> (i * 2)) & 0b11;
        let comp = (!base) & 0b11;
        out |= comp << ((n - 1 - i) * 2);
    }
    out
}

/// Strand-canonical packed form: `min(forward, reverse_complement)`.
///
/// Without this, a tag inside an inverted segment — or on a draft-assembly
/// contig that happens to be submitted in the opposite orientation — is stored
/// as the reverse complement of its homolog and can never hash-match it. Since
/// roughly half of a draft assembly's contigs are arbitrarily oriented, that
/// silently discards a large share of the genuinely shared tags.
pub fn canonical_packed(packed: u64, len: u8) -> u64 {
    packed.min(revcomp_packed(packed, len))
}

#[cfg(test)]
mod canonical_tests {
    use super::*;

    fn tag_from(seq: &str) -> GenomeTag {
        let mut buf = [0u8; 32];
        buf[..seq.len()].copy_from_slice(seq.as_bytes());
        let len = seq.len() as u8;
        GenomeTag {
            position: 0,
            contig_id: 0,
            sequence: buf,
            packed_sequence: pack_bytes(&buf, len),
            seq_len: len,
            direction: '+',
            enzyme: "test".to_string(),
        }
    }

    #[test]
    fn revcomp_is_an_involution() {
        for seq in ["ACGT", "AAAACCCCGGGGTTTT", "ACGTACGTACGTACGTACGTACGTACGTACGT"] {
            let t = tag_from(seq);
            let once = revcomp_packed(t.packed_sequence, t.seq_len);
            let twice = revcomp_packed(once, t.seq_len);
            assert_eq!(twice, t.packed_sequence, "seq {seq}");
        }
    }

    #[test]
    fn revcomp_matches_expected_sequence() {
        let t = tag_from("ACGT");
        // revcomp(ACGT) = ACGT
        assert_eq!(t.packed_revcomp(), tag_from("ACGT").packed_sequence);
        let t = tag_from("AAAC");
        // revcomp(AAAC) = GTTT
        assert_eq!(t.packed_revcomp(), tag_from("GTTT").packed_sequence);
    }

    #[test]
    fn canonical_agrees_across_strands() {
        // A tag and its reverse complement must hash to the same key, which is
        // what lets an inverted segment still match.
        let fwd = tag_from("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT");
        let rc_seq: String = fwd
            .sequence
            .iter()
            .take(32)
            .rev()
            .map(|&b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => 'N',
            })
            .collect();
        let rev = tag_from(&rc_seq);
        assert_eq!(fwd.canonical(), rev.canonical());
    }
}
