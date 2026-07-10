use super::registry::EnzymeConfig;

// ── Lookup table: only A/T/C/G/a/t/c/g are true ────────────────────────────

const ATCG_TABLE: [bool; 256] = {
    let mut table = [false; 256];
    table[b'A' as usize] = true; table[b'a' as usize] = true;
    table[b'T' as usize] = true; table[b't' as usize] = true;
    table[b'C' as usize] = true; table[b'c' as usize] = true;
    table[b'G' as usize] = true; table[b'g' as usize] = true;
    table
};

#[inline]
/// Check if all bytes in the window are A/T/C/G (no degenerate bases).
pub fn is_pure_atcg(window: &[u8]) -> bool {
    window.iter().all(|&b| ATCG_TABLE[b as usize])
}

/// A 2bRAD tag extracted from a genomic sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct Tag {
    /// Position (0-based) of the tag start on the scanned strand.
    pub position: usize,
    /// Strand orientation of the tag.
    pub direction: Direction,
    /// Tag sequence in 5'→3' orientation (up to 32 bp, zero-padded).
    pub sequence: [u8; 32],
    /// Actual length of the tag sequence (may be < 32 for short-tag enzymes).
    pub seq_len: u8,
}

/// Strand orientation of a recognition site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Direction {
    Forward,
    Reverse,
}

impl Direction {
    pub fn as_str(&self) -> &'static str {
        match self {
            Direction::Forward => "F",
            Direction::Reverse => "R",
        }
    }
}

// ── IUPAC matching via lookup table ────────────────────────────────────────

/// Bitmask for each base: A=1, T=2, C=4, G=8.
const BASE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 1; t[b'a' as usize] = 1;
    t[b'T' as usize] = 2; t[b't' as usize] = 2;
    t[b'C' as usize] = 4; t[b'c' as usize] = 4;
    t[b'G' as usize] = 8; t[b'g' as usize] = 8;
    t
};

#[inline]
fn iupac_matches(base: u8, code: u8) -> bool {
    let mask = BASE_MASK[base as usize];
    match code.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => mask == BASE_MASK[code as usize],
        b'N' => true,
        b'R' => mask & (1 | 8) != 0,          // A or G
        b'Y' => mask & (2 | 4) != 0,          // C or T
        b'S' => mask & (4 | 8) != 0,          // G or C
        b'W' => mask & (1 | 2) != 0,          // A or T
        b'K' => mask & (2 | 8) != 0,          // G or T
        b'M' => mask & (1 | 4) != 0,          // A or C
        b'B' => mask & (2 | 4 | 8) != 0,      // C or G or T (not A)
        b'D' => mask & (1 | 2 | 8) != 0,      // A or G or T (not C)
        b'H' => mask & (1 | 2 | 4) != 0,      // A or C or T (not G)
        b'V' => mask & (1 | 4 | 8) != 0,      // A or C or G (not T)
        _ => false,
    }
}

// ── Pattern matching ───────────────────────────────────────────────────────

/// A fixed anchor at a specific offset within the tag window.
struct Anchor { offset: usize, motif: &'static [u8] }

/// An IUPAC degenerate-base constraint at a specific offset.
struct IupacConstraint { offset: usize, allowed: u8 }

/// A pattern that defines forward or reverse recognition within a tag window.
struct Pattern {
    anchors: &'static [Anchor],
    iupac: &'static [IupacConstraint],
}

impl Pattern {
    #[inline]
    fn matches(&self, window: &[u8]) -> bool {
        // Check fixed anchors first (fast rejection)
        for a in self.anchors {
            let end = a.offset + a.motif.len();
            if end > window.len() {
                return false;
            }
            for (i, &expected) in a.motif.iter().enumerate() {
                if !iupac_matches(window[a.offset + i], expected) {
                    return false;
                }
            }
        }
        // Check IUPAC constraints
        for c in self.iupac {
            if c.offset >= window.len() {
                return false;
            }
            if (BASE_MASK[window[c.offset] as usize] & c.allowed) == 0 {
                return false;
            }
        }
        true
    }
}

// ── Enzyme pattern definitions (static, zero-cost) ─────────────────────────

// 1. BcgI (32)
const BCGI_F1: Anchor = Anchor{offset:10, motif:b"CGA"};
const BCGI_F2: Anchor = Anchor{offset:19, motif:b"TGC"};
const BCGI_R1: Anchor = Anchor{offset:10, motif:b"GCA"};
const BCGI_R2: Anchor = Anchor{offset:19, motif:b"TCG"};
const BCGI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BCGI_F1,BCGI_F2], iupac:&[]},
    Pattern{anchors:&[BCGI_R1,BCGI_R2], iupac:&[]},
];

// 2. AlfI (32, palindrome)
const ALFI_A1: Anchor = Anchor{offset:10, motif:b"GCA"};
const ALFI_A2: Anchor = Anchor{offset:19, motif:b"TGC"};
const ALFI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[ALFI_A1,ALFI_A2], iupac:&[]}];

// 3. AloI (27)
const ALOI_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const ALOI_F2: Anchor = Anchor{offset:17, motif:b"TCC"};
const ALOI_R1: Anchor = Anchor{offset:7, motif:b"GGA"};
const ALOI_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const ALOI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[ALOI_F1,ALOI_F2], iupac:&[]},
    Pattern{anchors:&[ALOI_R1,ALOI_R2], iupac:&[]},
];

// 4. BaeI (28, degenerate)
const BAEI_F1: Anchor = Anchor{offset:10, motif:b"AC"};
const BAEI_F2: Anchor = Anchor{offset:16, motif:b"GTA"};
const BAEI_R1: Anchor = Anchor{offset:7, motif:b"G"};
const BAEI_R2: Anchor = Anchor{offset:9, motif:b"TAC"};
const BAEI_FWD_IUPAC: [IupacConstraint;1] = [IupacConstraint{offset:19, allowed:6}]; // Y=[CT]
const BAEI_REV_IUPAC: [IupacConstraint;1] = [IupacConstraint{offset:8, allowed:9}];  // R=[AG]
const BAEI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BAEI_F1,BAEI_F2], iupac:&BAEI_FWD_IUPAC},
    Pattern{anchors:&[BAEI_R1,BAEI_R2], iupac:&BAEI_REV_IUPAC},
];

// 5. BplI (27, palindrome)
const BPLI_A1: Anchor = Anchor{offset:8, motif:b"GAG"};
const BPLI_A2: Anchor = Anchor{offset:16, motif:b"CTC"};
const BPLI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[BPLI_A1,BPLI_A2], iupac:&[]}];

// 6. BsaXI (27)
const BSAXI_F1: Anchor = Anchor{offset:9, motif:b"AC"};
const BSAXI_F2: Anchor = Anchor{offset:16, motif:b"CTCC"};
const BSAXI_R1: Anchor = Anchor{offset:7, motif:b"GGAG"};
const BSAXI_R2: Anchor = Anchor{offset:16, motif:b"GT"};
const BSAXI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSAXI_F1,BSAXI_F2], iupac:&[]},
    Pattern{anchors:&[BSAXI_R1,BSAXI_R2], iupac:&[]},
];

// 7. BslFI (25)
const BSLFI_F1: Anchor = Anchor{offset:6, motif:b"GGGAC"};
const BSLFI_R1: Anchor = Anchor{offset:14, motif:b"GTCCC"};
const BSLFI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSLFI_F1], iupac:&[]},
    Pattern{anchors:&[BSLFI_R1], iupac:&[]},
];

// 8. Bsp24I (27)
const BSP24I_F1: Anchor = Anchor{offset:8, motif:b"GAC"};
const BSP24I_F2: Anchor = Anchor{offset:17, motif:b"TGG"};
const BSP24I_R1: Anchor = Anchor{offset:7, motif:b"CCA"};
const BSP24I_R2: Anchor = Anchor{offset:16, motif:b"GTC"};
const BSP24I_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSP24I_F1,BSP24I_F2], iupac:&[]},
    Pattern{anchors:&[BSP24I_R1,BSP24I_R2], iupac:&[]},
];

// 9. CjeI (28)
const CJEI_F1: Anchor = Anchor{offset:8, motif:b"CCA"};
const CJEI_F2: Anchor = Anchor{offset:17, motif:b"GT"};
const CJEI_R1: Anchor = Anchor{offset:9, motif:b"AC"};
const CJEI_R2: Anchor = Anchor{offset:17, motif:b"TGG"};
const CJEI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CJEI_F1,CJEI_F2], iupac:&[]},
    Pattern{anchors:&[CJEI_R1,CJEI_R2], iupac:&[]},
];

// 10. CjePI (27)
const CJEPI_F1: Anchor = Anchor{offset:7, motif:b"CCA"};
const CJEPI_F2: Anchor = Anchor{offset:17, motif:b"TC"};
const CJEPI_R1: Anchor = Anchor{offset:8, motif:b"GA"};
const CJEPI_R2: Anchor = Anchor{offset:17, motif:b"TGG"};
const CJEPI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CJEPI_F1,CJEPI_F2], iupac:&[]},
    Pattern{anchors:&[CJEPI_R1,CJEPI_R2], iupac:&[]},
];

// 11. CspCI (33)
const CSPCI_F1: Anchor = Anchor{offset:11, motif:b"CAA"};
const CSPCI_F2: Anchor = Anchor{offset:19, motif:b"GTGG"};
const CSPCI_R1: Anchor = Anchor{offset:10, motif:b"CCAC"};
const CSPCI_R2: Anchor = Anchor{offset:19, motif:b"TTG"};
const CSPCI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CSPCI_F1,CSPCI_F2], iupac:&[]},
    Pattern{anchors:&[CSPCI_R1,CSPCI_R2], iupac:&[]},
];

// 12. FalI (27, palindrome)
const FALI_A1: Anchor = Anchor{offset:8, motif:b"AAG"};
const FALI_A2: Anchor = Anchor{offset:16, motif:b"CTT"};
const FALI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[FALI_A1,FALI_A2], iupac:&[]}];

// 13. HaeIV (27, degenerate)
const HAEIV_F1: Anchor = Anchor{offset:7, motif:b"GA"};
const HAEIV_R1: Anchor = Anchor{offset:9, motif:b"GA"};
const HAEIV_FWD_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:9, allowed:6},   // Y=[CT]
    IupacConstraint{offset:15, allowed:9},  // R=[AG]
];
const HAEIV_REV_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:11, allowed:6},  // Y=[CT]
    IupacConstraint{offset:17, allowed:9},  // R=[AG]
];
const HAEIV_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[HAEIV_F1], iupac:&HAEIV_FWD_IUPAC},
    Pattern{anchors:&[HAEIV_R1], iupac:&HAEIV_REV_IUPAC},
];

// 14. Hin4I (27, degenerate)
const HIN4I_F1: Anchor = Anchor{offset:8, motif:b"GA"};
const HIN4I_R1: Anchor = Anchor{offset:8, motif:b"GA"};
const HIN4I_FWD_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:10, allowed:6},   // Y=[CT]
    IupacConstraint{offset:16, allowed:13},  // [GAC]
];
const HIN4I_REV_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:10, allowed:14},  // [CTG]
    IupacConstraint{offset:16, allowed:9},   // R=[AG]
];
const HIN4I_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[HIN4I_F1], iupac:&HIN4I_FWD_IUPAC},
    Pattern{anchors:&[HIN4I_R1], iupac:&HIN4I_REV_IUPAC},
];

// 15. PpiI (27)
const PPII_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const PPII_F2: Anchor = Anchor{offset:16, motif:b"CTC"};
const PPII_R1: Anchor = Anchor{offset:8, motif:b"GAG"};
const PPII_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const PPII_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[PPII_F1,PPII_F2], iupac:&[]},
    Pattern{anchors:&[PPII_R1,PPII_R2], iupac:&[]},
];

// 16. PsrI (27)
const PSRI_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const PSRI_F2: Anchor = Anchor{offset:17, motif:b"TAC"};
const PSRI_R1: Anchor = Anchor{offset:7, motif:b"GTA"};
const PSRI_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const PSRI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[PSRI_F1,PSRI_F2], iupac:&[]},
    Pattern{anchors:&[PSRI_R1,PSRI_R2], iupac:&[]},
];

/// Map enzyme name to its static pattern set and tag length.
fn enzyme_patterns(name: &str) -> Option<(&'static [Pattern], usize)> {
    match name {
        "BcgI"  => Some((&BCGI_PATTERNS, 32)),
        "AlfI"  => Some((&ALFI_PATTERNS, 32)),
        "AloI"  => Some((&ALOI_PATTERNS, 27)),
        "BaeI"  => Some((&BAEI_PATTERNS, 28)),
        "BplI"  => Some((&BPLI_PATTERNS, 27)),
        "BsaXI" => Some((&BSAXI_PATTERNS, 27)),
        "BslFI" => Some((&BSLFI_PATTERNS, 25)),
        "Bsp24I"=> Some((&BSP24I_PATTERNS, 27)),
        "CjeI"  => Some((&CJEI_PATTERNS, 28)),
        "CjePI" => Some((&CJEPI_PATTERNS, 27)),
        "CspCI" => Some((&CSPCI_PATTERNS, 33)),
        "FalI"  => Some((&FALI_PATTERNS, 27)),
        "HaeIV" => Some((&HAEIV_PATTERNS, 27)),
        "Hin4I" => Some((&HIN4I_PATTERNS, 27)),
        "PpiI"  => Some((&PPII_PATTERNS, 27)),
        "PsrI"  => Some((&PSRI_PATTERNS, 27)),
        _ => None,
    }
}

// ── Public API ─────────────────────────────────────────────────────────────

/// In-silico digestion of a DNA sequence using a Type IIB restriction enzyme.
///
/// Fast2bRAD-M style: slides a `tag_length` window across the sequence.
/// If any pattern matches within the window and the window contains only
/// A/T/C/G bases, the window is emitted as a `Tag`.
/// Tags are sorted by position and deduplicated.
pub fn digest_sequence(seq: &[u8], enzyme: &EnzymeConfig) -> Vec<Tag> {
    let mut tags = Vec::new();

    let (patterns, tag_len) = match enzyme_patterns(&enzyme.name) {
        Some(p) => p,
        None => {
            // Fallback: use the old margin-based approach for unknown enzymes
            return digest_sequence_legacy(seq, enzyme);
        }
    };

    if seq.len() < tag_len {
        return tags;
    }

    let max_offset = seq.len() - tag_len;
    let tag_len_u8 = tag_len as u8;

    for pattern in patterns {
        for offset in 0..=max_offset {
            let window = &seq[offset..offset + tag_len];
            if pattern.matches(window) && is_pure_atcg(window) {
                let mut tag_seq = [0u8; 32];
                let copy_len = tag_len.min(32);
                tag_seq[..copy_len].copy_from_slice(&window[..copy_len]);
                tags.push(Tag {
                    position: offset,
                    direction: Direction::Forward,
                    sequence: tag_seq,
                    seq_len: tag_len_u8,
                });
            }
        }
    }

    // Sort by position and deduplicate (same position → same tag)
    tags.sort_by_key(|t| t.position);
    tags.dedup_by_key(|t| t.position);
    tags
}

// ── Legacy: original margin-based digestion ────────────────────────────────

/// Legacy margin-based digestion — kept for benchmark comparison.
/// Scans for recognition sites using left/right anchors + margins,
/// then does a reverse-complement scan separately.
pub fn digest_sequence_legacy(seq: &[u8], enzyme: &EnzymeConfig) -> Vec<Tag> {
    let mut tags = Vec::new();
    let pattern_len = enzyme.pattern_length();

    if seq.len() < enzyme.tag_length || pattern_len == 0 {
        return tags;
    }

    let max_i = seq.len() - pattern_len;
    let tag_len_u8 = enzyme.tag_length as u8;

    // Forward scan
    for i in 0..=max_i {
        if matches_pattern(seq, i, &enzyme.left_anchor, &enzyme.right_anchor, enzyme.spacer_length)
        {
            let start = i.saturating_sub(enzyme.left_margin);
            let end = i + pattern_len + enzyme.right_margin;
            if end <= seq.len() && is_pure_atcg(&seq[start..end]) {
                let mut tag_seq = [0u8; 32];
                let copy_len = (end - start).min(32);
                tag_seq[..copy_len].copy_from_slice(&seq[start..start + copy_len]);
                tags.push(Tag {
                    position: i,
                    direction: Direction::Forward,
                    sequence: tag_seq,
                    seq_len: tag_len_u8,
                });
            }
        }
    }

    // Reverse-complement scan
    let rc_left = reverse_complement(&enzyme.right_anchor);
    let rc_right = reverse_complement(&enzyme.left_anchor);

    for i in 0..=max_i {
        if matches_pattern(seq, i, &rc_left, &rc_right, enzyme.spacer_length) {
            let start = i.saturating_sub(enzyme.right_margin);
            let end = i + pattern_len + enzyme.left_margin;
            if end <= seq.len() && is_pure_atcg(&seq[start..end]) {
                let mut tag_seq = [0u8; 32];
                let copy_len = (end - start).min(32);
                tag_seq[..copy_len].copy_from_slice(&seq[start..start + copy_len]);
                reverse_complement_in_place(&mut tag_seq[..copy_len]);
                tags.push(Tag {
                    position: i,
                    direction: Direction::Reverse,
                    sequence: tag_seq,
                    seq_len: tag_len_u8,
                });
            }
        }
    }

    tags.sort_by_key(|t| t.position);
    tags.dedup_by_key(|t| t.position);
    tags
}

fn matches_pattern(seq: &[u8], pos: usize, left: &str, right: &str, spacer: usize) -> bool {
    let left_b = left.as_bytes();
    let right_b = right.as_bytes();
    let pat_len = left_b.len() + spacer + right_b.len();
    if pos + pat_len > seq.len() { return false; }
    for (j, &expected) in left_b.iter().enumerate() {
        if !iupac_matches(seq[pos + j], expected) { return false; }
    }
    if !right_b.is_empty() {
        for (j, &expected) in right_b.iter().enumerate() {
            if !iupac_matches(seq[pos + left_b.len() + spacer + j], expected) { return false; }
        }
    }
    true
}

fn reverse_complement_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'N' => b'N', b'R' => b'Y', b'Y' => b'R', b'S' => b'S',
        b'W' => b'W', b'K' => b'M', b'M' => b'K',
        b'B' => b'V', b'D' => b'H', b'H' => b'D', b'V' => b'B',
        _ => base,
    }
}

fn reverse_complement(s: &str) -> String {
    s.as_bytes().iter().rev().map(|&b| reverse_complement_base(b) as char).collect()
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    for b in seq.iter_mut() { *b = reverse_complement_base(*b); }
    seq.reverse();
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::registry::EnzymeConfig;

    #[test]
    fn test_is_pure_atcg() {
        assert!(is_pure_atcg(b"ATCG"));
        assert!(!is_pure_atcg(b"ATNG"));
        assert!(!is_pure_atcg(b"ATCGN"));
    }

    #[test]
    fn test_digest_empty_sequence() {
        let tags = digest_sequence(b"", &EnzymeConfig::bcg_i());
        assert!(tags.is_empty());
    }

    #[test]
    fn test_digest_no_sites() {
        let tags = digest_sequence(b"ATATATATATATATATATATATATATATAT", &EnzymeConfig::bcg_i());
        assert!(tags.is_empty());
    }

    #[test]
    fn test_digest_bcgI_real_site() {
        // BcgI: tag_length=32, anchors at offset 10=CGA, offset 19=TGC
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAAAA";
        assert_eq!(seq.len(), 32);
        let tags = digest_sequence(seq, &EnzymeConfig::bcg_i());
        assert_eq!(tags.len(), 1);
        assert_eq!(tags[0].position, 0);
        assert_eq!(tags[0].seq_len, 32);
        assert_eq!(&tags[0].sequence[..32], seq.as_slice());
    }

    #[test]
    fn test_digest_bcgI_with_n_rejected() {
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAANAAA";
        let tags = digest_sequence(seq, &EnzymeConfig::bcg_i());
        assert!(tags.is_empty(), "N-containing window should be rejected");
    }

    #[test]
    fn test_digest_cjeI_correct_site() {
        let seq = b"AAAAAAAACCAAAAAAAGTAAAAAAAAAAAAAAAAAA";
        assert_eq!(seq.len(), 37);
        let tags = digest_sequence(seq, &EnzymeConfig::cje_i());
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn test_digest_baeI_iupac() {
        // BaeI fwd: offset 19 must be C or T (Y=bitmask 6)
        let seq = b"AAAAAAAAAAACAAAAGTACCAAAAAAA"; // 28 bp, C@19
        assert_eq!(seq.len(), 28);
        let tags = digest_sequence(seq, &EnzymeConfig::bae_i());
        assert_eq!(tags.len(), 1);

        let seq_t = b"AAAAAAAAAAACAAAAGTACTAAAAAAA"; // T@19
        let tags_t = digest_sequence(seq_t, &EnzymeConfig::bae_i());
        assert_eq!(tags_t.len(), 1);

        let seq_a = b"AAAAAAAAAAACAAAAGTAACAAAAAAA"; // A@19 (bad)
        let tags_a = digest_sequence(seq_a, &EnzymeConfig::bae_i());
        assert!(tags_a.is_empty());
    }

    #[test]
    fn test_iupac_matching() {
        assert!(iupac_matches(b'A', b'R'));
        assert!(iupac_matches(b'G', b'R'));
        assert!(!iupac_matches(b'C', b'R'));
        assert!(iupac_matches(b'C', b'Y'));
        assert!(iupac_matches(b'T', b'Y'));
        assert!(iupac_matches(b'C', b'N'));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("CGA"), "TCG");
        assert_eq!(reverse_complement("TGC"), "GCA");
    }

    #[test]
    fn test_all_16_enzymes_have_patterns() {
        let registry = crate::enzyme::EnzymeRegistry::new();
        for enzyme in registry.all() {
            assert!(
                enzyme_patterns(&enzyme.name).is_some(),
                "{} missing static patterns",
                enzyme.name
            );
        }
    }
}
