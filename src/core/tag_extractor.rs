use crate::utils::fxhash::FastHashMap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;
use crate::enzyme::{digest_sequence, EnzymeConfig};

/// A single genome tag extracted after enzyme digestion.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomeTag {
    pub position: usize,
    /// Tag sequence as a fixed 32-byte array (zero-padded, actual length may vary).
    pub sequence: [u8; 32],
    /// 2-bit packed sequence (64 bits = 32 bp), aligned with `sequence`.
    pub packed_sequence: u64,
    /// Actual sequence length (tag may be shorter than 32 bp).
    pub seq_len: u8,
    pub direction: char,
    pub enzyme: String,
}

/// A collection of tags from a single genome/contig.
#[derive(Debug, Clone)]
pub struct TagSet {
    pub genome_id: String,
    pub chromosome: String,
    pub tags: Vec<GenomeTag>,
    pub total_length: usize,
    pub gc_content: f64,
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
    pub fn extract_from_sequence(seq: &[u8], enzyme: &EnzymeConfig) -> Vec<GenomeTag> {
        let digested = digest_sequence(seq, enzyme);
        digested.into_iter().map(|tag| {
            let packed = pack_bytes(&tag.sequence, tag.seq_len);
            GenomeTag {
                position: tag.position,
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

    /// Extract tags from a single FASTA file (first sequence only for simplicity).
    pub fn extract_from_fasta(path: &Path, enzyme: &EnzymeConfig) -> Result<TagSet, ExtractError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let header = lines
            .next()
            .ok_or_else(|| ExtractError::InvalidFasta("Empty file".to_string()))??;

        let genome_id = header
            .trim_start_matches('>')
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        let mut sequence = Vec::new();
        for line in lines {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                break; // For simplicity, handle first contig only
            }
            sequence.extend(trimmed.bytes().filter(|&b| b.is_ascii_alphabetic()));
        }

        let total_length = sequence.len();
        let gc_count = sequence
            .iter()
            .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        let gc_content = gc_count as f64 / total_length.max(1) as f64;

        let tags = Self::extract_from_sequence(&sequence, enzyme);

        Ok(TagSet {
            genome_id,
            chromosome: "chrom1".to_string(), // Placeholder for single-contig mode
            tags,
            total_length,
            gc_content,
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
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        let header = lines
            .next()
            .ok_or_else(|| ExtractError::InvalidFasta("Empty file".to_string()))??;

        let genome_id = header
            .trim_start_matches('>')
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();

        let mut sequence = Vec::new();
        for line in lines {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                break;
            }
            sequence.extend(trimmed.bytes().filter(|&b| b.is_ascii_alphabetic()));
        }

        let total_length = sequence.len();
        let gc_count = sequence
            .iter()
            .filter(|&&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
            .count();
        let gc_content = gc_count as f64 / total_length.max(1) as f64;

        // Step 2: parallel digestion across enzymes
        let sets: Vec<_> = enzymes
            .par_iter()
            .map(|enzyme| {
                let tags = Self::extract_from_sequence(&sequence, enzyme);
                let tag_set = TagSet {
                    genome_id: genome_id.clone(),
                    chromosome: "chrom1".to_string(),
                    tags,
                    total_length,
                    gc_content,
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
