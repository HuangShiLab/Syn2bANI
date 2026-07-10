use std::fs::File;
use std::io::Write;
use std::path::Path;
use byteorder::{LittleEndian, WriteBytesExt};

use super::fasta_parser::IoError;

/// Magic bytes for the Syn2bANI binary sketch format.
pub const S2BA_MAGIC: [u8; 4] = *b"S2BA";
/// Current sketch format version.
pub const S2BA_VERSION: u32 = 1;

/// A single 2bRAD tag stored in a compact binary representation.
#[derive(Debug, Clone)]
pub struct SketchTag {
    pub position: u64,
    /// 2-bit packed sequence (up to 32 bp in a u64).
    pub seq: u64,
    /// 0 = forward, 1 = reverse.
    pub direction: u8,
    pub enzyme_id: u16,
}

/// Per-chromosome sketch data.
#[derive(Debug, Clone)]
pub struct ChromSketch {
    pub name: String,
    pub tags: Vec<SketchTag>,
    pub gc_content: f64,
    pub length: u64,
}

/// Metadata summarizing the whole-genome sketch.
#[derive(Debug, Clone)]
pub struct SketchMetadata {
    pub total_length: u64,
    pub gc_content: f64,
    pub tag_count: u64,
}

/// Top-level binary sketch container for a single genome.
#[derive(Debug, Clone)]
pub struct TgtSketch {
    pub magic: [u8; 4],
    pub version: u32,
    pub genome_id: String,
    pub chromosomes: Vec<ChromSketch>,
    pub metadata: SketchMetadata,
}

impl Default for TgtSketch {
    fn default() -> Self {
        Self {
            magic: S2BA_MAGIC,
            version: S2BA_VERSION,
            genome_id: String::new(),
            chromosomes: Vec::new(),
            metadata: SketchMetadata {
                total_length: 0,
                gc_content: 0.0,
                tag_count: 0,
            },
        }
    }
}

/// Pack a DNA sequence into a 64-bit integer using 2-bit encoding.
///
/// Encoding: A=00, C=01, G=10, T=11. Only the first 32 bases are packed.
pub fn pack_sequence(seq: &[u8]) -> u64 {
    let mut packed: u64 = 0;
    let len = seq.len().min(32);
    for i in 0..len {
        let bits = match seq[i].to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1u64,
            b'G' => 2u64,
            b'T' => 3u64,
            _ => 0u64,
        };
        packed |= bits << (2 * i);
    }
    packed
}

/// Unpack a 64-bit 2-bit encoded sequence back into a DNA byte vector.
pub fn unpack_sequence(packed: u64) -> Vec<u8> {
    let mut seq = Vec::with_capacity(32);
    for i in 0..32 {
        let bits = (packed >> (2 * i)) & 0x3;
        let base = match bits {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        };
        seq.push(base);
    }
    seq
}

/// Write a `TgtSketch` to a binary file using little-endian byte order.
pub fn write_sketch(sketch: &TgtSketch, path: &Path) -> Result<(), IoError> {
    let mut file = File::create(path)?;
    file.write_all(&sketch.magic)?;
    file.write_u32::<LittleEndian>(sketch.version)?;

    let genome_id_bytes = sketch.genome_id.as_bytes();
    file.write_u32::<LittleEndian>(genome_id_bytes.len() as u32)?;
    file.write_all(genome_id_bytes)?;

    file.write_u32::<LittleEndian>(sketch.chromosomes.len() as u32)?;

    for chrom in &sketch.chromosomes {
        let name_bytes = chrom.name.as_bytes();
        file.write_u32::<LittleEndian>(name_bytes.len() as u32)?;
        file.write_all(name_bytes)?;
        file.write_u64::<LittleEndian>(chrom.length)?;
        file.write_f64::<LittleEndian>(chrom.gc_content)?;
        file.write_u32::<LittleEndian>(chrom.tags.len() as u32)?;
        for tag in &chrom.tags {
            file.write_u64::<LittleEndian>(tag.position)?;
            file.write_u64::<LittleEndian>(tag.seq)?;
            file.write_u8(tag.direction)?;
            file.write_u16::<LittleEndian>(tag.enzyme_id)?;
        }
    }

    file.write_u64::<LittleEndian>(sketch.metadata.total_length)?;
    file.write_f64::<LittleEndian>(sketch.metadata.gc_content)?;
    file.write_u64::<LittleEndian>(sketch.metadata.tag_count)?;

    Ok(())
}
