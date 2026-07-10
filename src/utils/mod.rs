/// Common utilities for sequence manipulation and statistics.

pub mod fxhash;

/// Reverse-complement a DNA sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            _ => b'N',
        })
        .collect()
}

/// Compute GC content as a fraction [0.0, 1.0].
pub fn gc_content(seq: &[u8]) -> f64 {
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len().max(1) as f64
}

/// Check whether a DNA base matches an IUPAC code.
pub fn iupac_matches(base: u8, code: u8) -> bool {
    match code.to_ascii_uppercase() {
        b'A' => base == b'A' || base == b'a',
        b'C' => base == b'C' || base == b'c',
        b'G' => base == b'G' || base == b'g',
        b'T' => base == b'T' || base == b't',
        b'N' => true,
        b'R' => matches!(base, b'A' | b'a' | b'G' | b'g'),
        b'Y' => matches!(base, b'C' | b'c' | b'T' | b't'),
        b'S' => matches!(base, b'G' | b'g' | b'C' | b'c'),
        b'W' => matches!(base, b'A' | b'a' | b'T' | b't'),
        b'K' => matches!(base, b'G' | b'g' | b'T' | b't'),
        b'M' => matches!(base, b'A' | b'a' | b'C' | b'c'),
        b'B' => !matches!(base, b'A' | b'a'),
        b'D' => !matches!(base, b'C' | b'c'),
        b'H' => !matches!(base, b'G' | b'g'),
        b'V' => !matches!(base, b'T' | b't'),
        _ => false,
    }
}
