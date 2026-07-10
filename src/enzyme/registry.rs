use std::fmt;

/// Configuration for a Type IIB restriction enzyme.
///
/// Type IIB enzymes recognize asymmetric sequences and cleave on both sides,
/// producing short, iso-length tags. The `left_anchor` and `right_anchor` are the
/// constant regions flanking the degenerate spacer.
#[derive(Debug, Clone, PartialEq)]
pub struct EnzymeConfig {
    pub name: String,
    pub left_anchor: String,
    pub right_anchor: String,
    pub spacer_length: usize,
    pub tag_length: usize,
    pub left_margin: usize,
    pub right_margin: usize,
}

impl EnzymeConfig {
    /// Create a new enzyme configuration.
    pub fn new(
        name: &str,
        left: &str,
        right: &str,
        spacer: usize,
        tag_len: usize,
        left_m: usize,
        right_m: usize,
    ) -> Self {
        Self {
            name: name.to_string(),
            left_anchor: left.to_string(),
            right_anchor: right.to_string(),
            spacer_length: spacer,
            tag_length: tag_len,
            left_margin: left_m,
            right_margin: right_m,
        }
    }

    /// Total length of the recognition sequence (anchors + spacer).
    pub fn pattern_length(&self) -> usize {
        self.left_anchor.len() + self.spacer_length + self.right_anchor.len()
    }

    /// Human-readable recognition pattern, e.g. `CGA-N6-TGC`.
    pub fn recognition_pattern(&self) -> String {
        if self.right_anchor.is_empty() {
            self.left_anchor.clone()
        } else {
            format!(
                "{}-N{}-{}",
                self.left_anchor, self.spacer_length, self.right_anchor
            )
        }
    }

    // --- Predefined Type IIB enzymes (2bRAD-M panel) ---

    pub fn bcg_i() -> Self {
        Self::new("BcgI", "CGA", "TGC", 6, 32, 10, 10)
    }

    pub fn alf_i() -> Self {
        Self::new("AlfI", "GCA", "TGC", 6, 32, 10, 10)
    }

    pub fn alo_i() -> Self {
        Self::new("AloI", "GAAC", "TCC", 6, 27, 7, 7)
    }

    pub fn bae_i() -> Self {
        Self::new("BaeI", "AC", "GTAYC", 4, 28, 10, 7)
    }

    pub fn bpl_i() -> Self {
        Self::new("BplI", "GAG", "CTC", 5, 27, 8, 8)
    }

    pub fn bsa_xi() -> Self {
        Self::new("BsaXI", "AC", "CTCC", 5, 27, 9, 7)
    }

    pub fn bsl_fi() -> Self {
        Self::new("BslFI", "GGGAC", "", 0, 21, 6, 10)
    }

    pub fn bsp24_i() -> Self {
        Self::new("Bsp24I", "GAC", "TGG", 6, 27, 8, 7)
    }

    pub fn cje_i() -> Self {
        Self::new("CjeI", "CCA", "GT", 6, 28, 8, 9)
    }

    pub fn cje_pi() -> Self {
        Self::new("CjePI", "CCA", "TC", 7, 27, 7, 8)
    }

    pub fn csp_ci() -> Self {
        Self::new("CspCI", "CAA", "GTGG", 5, 33, 11, 10)
    }

    pub fn fal_i() -> Self {
        Self::new("FalI", "AAG", "CTT", 5, 27, 8, 8)
    }

    pub fn hae_iv() -> Self {
        Self::new("HaeIV", "GAY", "RTC", 5, 27, 7, 9)
    }

    pub fn hin4_i() -> Self {
        Self::new("Hin4I", "GAY", "VTC", 5, 27, 8, 8)
    }

    pub fn ppi_i() -> Self {
        Self::new("PpiI", "GAAC", "CTC", 5, 28, 7, 9)
    }

    pub fn psr_i() -> Self {
        Self::new("PsrI", "GAAC", "TAC", 6, 27, 7, 7)
    }
}

impl fmt::Display for EnzymeConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.right_anchor.is_empty() {
            write!(
                f,
                "{}: {} ({} bp tag)",
                self.name, self.left_anchor, self.tag_length
            )
        } else {
            write!(
                f,
                "{}: {}-N{}-{} ({} bp tag)",
                self.name,
                self.left_anchor,
                self.spacer_length,
                self.right_anchor,
                self.tag_length
            )
        }
    }
}

/// Registry holding all 16 Type IIB enzymes used in the 2bRAD-M panel.
#[derive(Debug, Clone, Default)]
pub struct EnzymeRegistry {
    enzymes: Vec<EnzymeConfig>,
}

impl EnzymeRegistry {
    /// Create a registry populated with the full 2bRAD-M enzyme panel.
    pub fn new() -> Self {
        Self {
            enzymes: vec![
                EnzymeConfig::bcg_i(),
                EnzymeConfig::alf_i(),
                EnzymeConfig::alo_i(),
                EnzymeConfig::bae_i(),
                EnzymeConfig::bpl_i(),
                EnzymeConfig::bsa_xi(),
                EnzymeConfig::bsl_fi(),
                EnzymeConfig::bsp24_i(),
                EnzymeConfig::cje_i(),
                EnzymeConfig::cje_pi(),
                EnzymeConfig::csp_ci(),
                EnzymeConfig::fal_i(),
                EnzymeConfig::hae_iv(),
                EnzymeConfig::hin4_i(),
                EnzymeConfig::ppi_i(),
                EnzymeConfig::psr_i(),
            ],
        }
    }

    /// Look up an enzyme by name (case-insensitive).
    pub fn get(&self, name: &str) -> Option<&EnzymeConfig> {
        self.enzymes
            .iter()
            .find(|e| e.name.eq_ignore_ascii_case(name))
    }

    /// Return a slice of all registered enzymes.
    pub fn all(&self) -> &[EnzymeConfig] {
        &self.enzymes
    }

    pub fn len(&self) -> usize {
        self.enzymes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.enzymes.is_empty()
    }
}

impl fmt::Display for EnzymeRegistry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "EnzymeRegistry ({} enzymes):", self.len())?;
        for enzyme in &self.enzymes {
            writeln!(f, "  - {}", enzyme)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bcg_i_params() {
        let e = EnzymeConfig::bcg_i();
        assert_eq!(e.name, "BcgI");
        assert_eq!(e.left_anchor, "CGA");
        assert_eq!(e.right_anchor, "TGC");
        assert_eq!(e.spacer_length, 6);
        assert_eq!(e.tag_length, 32);
        assert_eq!(e.left_margin, 10);
        assert_eq!(e.right_margin, 10);
        assert_eq!(e.pattern_length(), 12);
    }

    #[test]
    fn test_registry_len() {
        let reg = EnzymeRegistry::new();
        assert_eq!(reg.len(), 16);
    }

    #[test]
    fn test_registry_lookup() {
        let reg = EnzymeRegistry::new();
        assert!(reg.get("BcgI").is_some());
        assert!(reg.get("bcgi").is_some());
        assert!(reg.get("EcoRI").is_none());
    }

    #[test]
    fn test_recognition_pattern() {
        assert_eq!(EnzymeConfig::bcg_i().recognition_pattern(), "CGA-N6-TGC");
        assert_eq!(EnzymeConfig::bsl_fi().recognition_pattern(), "GGGAC");
    }
}
