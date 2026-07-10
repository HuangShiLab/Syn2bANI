use std::fs::File;
use std::io::Read;
use std::path::Path;
use byteorder::{LittleEndian, ReadBytesExt};

use crate::core::GenomeTag;
use crate::io::{
    ChromSketch, IoError, S2BA_MAGIC, S2BA_VERSION, SketchMetadata, SketchTag, TgtSketch,
};

/// Read a `TgtSketch` from a binary file, validating magic bytes and version.
pub fn read_sketch(path: &Path) -> Result<TgtSketch, IoError> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 4];
    file.read_exact(&mut magic)?;
    if magic != S2BA_MAGIC {
        return Err(IoError::Parse(format!(
            "Invalid magic bytes: expected {:?}, got {:?}",
            S2BA_MAGIC, magic
        )));
    }

    let version = file.read_u32::<LittleEndian>()?;
    if version != S2BA_VERSION {
        return Err(IoError::Parse(format!(
            "Unsupported version: expected {}, got {}",
            S2BA_VERSION, version
        )));
    }

    let genome_id_len = file.read_u32::<LittleEndian>()? as usize;
    let mut genome_id_bytes = vec![0u8; genome_id_len];
    file.read_exact(&mut genome_id_bytes)?;
    let genome_id = String::from_utf8_lossy(&genome_id_bytes).to_string();

    let chrom_count = file.read_u32::<LittleEndian>()? as usize;
    let mut chromosomes = Vec::with_capacity(chrom_count);

    for _ in 0..chrom_count {
        let name_len = file.read_u32::<LittleEndian>()? as usize;
        let mut name_bytes = vec![0u8; name_len];
        file.read_exact(&mut name_bytes)?;
        let name = String::from_utf8_lossy(&name_bytes).to_string();

        let length = file.read_u64::<LittleEndian>()?;
        let gc_content = file.read_f64::<LittleEndian>()?;
        let tag_count = file.read_u32::<LittleEndian>()? as usize;

        let mut tags = Vec::with_capacity(tag_count);
        for _ in 0..tag_count {
            let position = file.read_u64::<LittleEndian>()?;
            let seq = file.read_u64::<LittleEndian>()?;
            let direction = file.read_u8()?;
            let enzyme_id = file.read_u16::<LittleEndian>()?;
            tags.push(SketchTag {
                position,
                seq,
                direction,
                enzyme_id,
            });
        }

        chromosomes.push(ChromSketch {
            name,
            tags,
            gc_content,
            length,
        });
    }

    let total_length = file.read_u64::<LittleEndian>()?;
    let gc_content = file.read_f64::<LittleEndian>()?;
    let tag_count = file.read_u64::<LittleEndian>()?;

    Ok(TgtSketch {
        magic,
        version,
        genome_id,
        chromosomes,
        metadata: SketchMetadata {
            total_length,
            gc_content,
            tag_count,
        },
    })
}

/// Convert a slice of `SketchTag`s back into `GenomeTag`s.
///
/// `enzyme_name` is used to populate the `enzyme` field of each `GenomeTag`.
pub fn sketch_tags_to_genome_tags(sketch_tags: &[SketchTag], enzyme_name: &str) -> Vec<GenomeTag> {
    sketch_tags
        .iter()
        .map(|st| {
            let unpacked = crate::io::unpack_sequence(st.seq);
            let mut sequence = [0u8; 32];
            let copy_len = unpacked.len().min(32);
            sequence[..copy_len].copy_from_slice(&unpacked[..copy_len]);
            GenomeTag {
                position: st.position as usize,
                sequence,
                packed_sequence: st.seq,
                seq_len: copy_len as u8,
                direction: if st.direction == 0 { '+' } else { '-' },
                enzyme: enzyme_name.to_string(),
            }
        })
        .collect()
}
