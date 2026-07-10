use std::io::Write;
use syn2bani::utils::fxhash::FastHashMap;
use tempfile::NamedTempFile;

use syn2bani::core::{TagExtractor, TagMatcher, AniCalculator, AniConfig, MatchConfig, WeightStrategy, TagSet, MultiEnzymeTagSet};
use syn2bani::enzyme::{EnzymeRegistry, EnzymeConfig};
use syn2bani::io::{parse_fasta, write_sketch, read_sketch, TsvFormatter};
use syn2bani::parallel::parallel_compare;

/// Create a simple synthetic FASTA file for testing.
fn create_test_fasta(seq: &[u8], id: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    writeln!(file, ">{}", id).unwrap();
    writeln!(file, "{}", std::str::from_utf8(seq).unwrap()).unwrap();
    file
}

/// A 5 kb synthetic sequence with one BcgI site (CGA-N6-TGC).
fn synthetic_genome_5kb() -> Vec<u8> {
    let mut seq = vec![b'A'; 5000];
    // Embed a BcgI site at position 1000
    let site = b"CGAAAAAAATGC";
    seq[1000..1000 + site.len()].copy_from_slice(site);
    seq
}

/// A closely related genome with ~2% divergence (one SNP in the BcgI tag region).
fn synthetic_genome_5kb_diverged() -> Vec<u8> {
    let mut seq = synthetic_genome_5kb();
    // Change one base in the tag region
    seq[1005] = b'G';
    seq
}

fn extract_multi_enzyme_from_sequence(seq: &[u8], enzymes: &[EnzymeConfig]) -> MultiEnzymeTagSet {
    let mut sets = FastHashMap::default();
    for enzyme in enzymes {
        let tags = TagExtractor::extract_from_sequence(seq, enzyme);
        let tag_set = TagSet {
            genome_id: "test".to_string(),
            chromosome: "chrom1".to_string(),
            tags,
            total_length: seq.len(),
            gc_content: 0.0,
        };
        sets.insert(enzyme.name.clone(), tag_set);
    }
    MultiEnzymeTagSet { sets }
}

#[test]
fn test_digest_bcg_i() {
    let seq = synthetic_genome_5kb();
    let enzyme = EnzymeConfig::bcg_i();
    let tags = syn2bani::enzyme::digest_sequence(&seq, &enzyme);
    assert!(!tags.is_empty(), "Should find at least one BcgI tag");
    assert_eq!(tags[0].sequence.len(), 32, "BcgI tag length should be 32 bp");
}

#[test]
fn test_tag_matching_and_ani() {
    let q_seq = synthetic_genome_5kb();
    let r_seq = synthetic_genome_5kb_diverged();

    let enzyme = EnzymeConfig::bcg_i();
    let q_tags = TagExtractor::extract_from_sequence(&q_seq, &enzyme);
    let r_tags = TagExtractor::extract_from_sequence(&r_seq, &enzyme);

    let q_set = TagSet {
        genome_id: "query".to_string(),
        chromosome: "chrom1".to_string(),
        tags: q_tags,
        total_length: q_seq.len(),
        gc_content: 0.0,
    };

    let r_set = TagSet {
        genome_id: "ref".to_string(),
        chromosome: "chrom1".to_string(),
        tags: r_tags,
        total_length: r_seq.len(),
        gc_content: 0.0,
    };

    let match_config = MatchConfig::default();
    let match_result = TagMatcher::match_tag_sets(&q_set, &r_set, &match_config);

    assert!(!match_result.matched_pairs.is_empty(), "Should have matched pairs");

    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 1,
        min_af: 0.0,
        debias: false,
        use_gbrt_debias: false,
    };
    let ani_result = AniCalculator::calculate_ani(&match_result, &ani_config);

    assert!(ani_result.ani > 0.95, "ANI should be high for closely related genomes");
    assert!(ani_result.ani <= 1.0, "ANI should not exceed 1.0");
}

#[test]
fn test_fasta_parser() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let file = create_test_fasta(seq, "test_contig");
    let records = parse_fasta(file.path()).unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].id, "test_contig");
    assert_eq!(records[0].sequence, seq);
}

#[test]
fn test_sketch_roundtrip() {
    use syn2bani::io::{TgtSketch, ChromSketch, SketchTag, SketchMetadata, S2BA_MAGIC};

    let sketch = TgtSketch {
        magic: S2BA_MAGIC,
        version: 1,
        genome_id: "test_genome".to_string(),
        chromosomes: vec![ChromSketch {
            name: "chrom1".to_string(),
            tags: vec![SketchTag {
                position: 1000,
                seq: 0b0101,
                direction: 0,
                enzyme_id: 0,
            }],
            gc_content: 0.5,
            length: 10000,
        }],
        metadata: SketchMetadata {
            total_length: 10000,
            gc_content: 0.5,
            tag_count: 1,
        },
    };

    let file = NamedTempFile::new().unwrap();
    write_sketch(&sketch, file.path()).unwrap();
    let loaded = read_sketch(file.path()).unwrap();

    assert_eq!(loaded.genome_id, "test_genome");
    assert_eq!(loaded.chromosomes.len(), 1);
    assert_eq!(loaded.chromosomes[0].tags.len(), 1);
    assert_eq!(loaded.chromosomes[0].tags[0].position, 1000);
}

#[test]
fn test_tsv_formatter() {
    let mut buf = Vec::new();
    TsvFormatter::write_header(&mut buf).unwrap();
    let header = String::from_utf8(buf).unwrap();
    assert!(header.contains("ani"));
    assert!(header.contains("af_q"));
}

#[test]
fn test_parallel_compare() {
    let q_seq = synthetic_genome_5kb();
    let r_seq = synthetic_genome_5kb_diverged();

    let q_file = create_test_fasta(&q_seq, "query");
    let r_file = create_test_fasta(&r_seq, "ref");

    let pairs = vec![
        (q_file.path().to_path_buf(), r_file.path().to_path_buf()),
    ];

    let ani_config = AniConfig {
        weight_strategy: WeightStrategy::Uniform,
        min_shared_tags: 1,
        min_af: 0.0,
        debias: false,
        use_gbrt_debias: false,
    };

    let results = parallel_compare(&pairs, &ani_config, 2);
    assert_eq!(results.len(), 1);
    assert!(results[0].ani > 0.95);
}

#[test]
fn test_registry_all_enzymes() {
    let registry = EnzymeRegistry::new();
    assert_eq!(registry.len(), 16);
    assert!(registry.get("BcgI").is_some());
    assert!(registry.get("bcgi").is_some()); // case-insensitive
    assert!(registry.get("EcoRI").is_none());
}

#[test]
fn test_multi_enzyme_extraction() {
    let seq = synthetic_genome_5kb();
    let registry = EnzymeRegistry::new();
    let enzymes = registry.all().to_vec();
    let multi_set = extract_multi_enzyme_from_sequence(&seq, &enzymes);
    assert!(!multi_set.sets.is_empty());
}
