use std::path::Path;
use needletail::parse_fastx_file;
use thiserror::Error;

/// Errors that can occur during I/O operations.
#[derive(Error, Debug)]
pub enum IoError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
}

/// A single FASTA record (multi-contig FASTA files yield multiple records).
#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

impl FastaRecord {
    /// Compute GC content as a fraction [0.0, 1.0].
    pub fn gc_content(&self) -> f64 {
        let gc = self
            .sequence
            .iter()
            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
            .count();
        gc as f64 / self.sequence.len().max(1) as f64
    }
}

/// Parse a FASTA file into a vector of `FastaRecord`s.
///
/// Handles multi-contig FASTA files (common for MAGs) by returning one record per contig.
pub fn parse_fasta(path: &Path) -> Result<Vec<FastaRecord>, IoError> {
    let mut reader = parse_fastx_file(path).map_err(|e| IoError::Parse(e.to_string()))?;
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| IoError::Parse(e.to_string()))?;
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let seq = rec.seq().to_vec();
        records.push(FastaRecord { id, sequence: seq });
    }
    Ok(records)
}
