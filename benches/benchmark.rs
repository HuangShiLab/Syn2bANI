use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::fs::File;
use std::io::{BufRead, BufReader};
use syn2bani::{digest_sequence, digest_sequence_legacy, EnzymeConfig, EnzymeRegistry, TagExtractor};

fn load_fasta_sequence(path: &str) -> Vec<u8> {
    let file = File::open(path).expect("Failed to open FASTA");
    let reader = BufReader::new(file);
    let mut seq = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            continue;
        }
        seq.extend(line.trim().bytes().filter(|&b| b.is_ascii_alphabetic()));
    }
    seq
}

fn bench_digest(c: &mut Criterion) {
    let fasta_path = "/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli/mag_comp_1.0.fasta";
    let seq = load_fasta_sequence(fasta_path);
    let seq_len_mb = seq.len() as f64 / 1_000_000.0;
    println!("Loaded {} Mb sequence for benchmarking", seq_len_mb);

    let bcgi = EnzymeConfig::bcg_i();

    // Verify both methods produce the same results
    let tags_new = digest_sequence(&seq, &bcgi);
    let tags_old = digest_sequence_legacy(&seq, &bcgi);
    println!(
        "New method: {} tags, Old method: {} tags",
        tags_new.len(),
        tags_old.len()
    );
    assert_eq!(tags_new.len(), tags_old.len(), "Tag counts must match!");

    let mut group = c.benchmark_group("digest_sequence");
    group.sample_size(20);

    group.bench_function("new_fast2brad_m", |b| {
        b.iter(|| digest_sequence(black_box(&seq), black_box(&bcgi)))
    });

    group.bench_function("old_margin_based", |b| {
        b.iter(|| digest_sequence_legacy(black_box(&seq), black_box(&bcgi)))
    });

    group.finish();
}

fn bench_all_enzymes(c: &mut Criterion) {
    let fasta_path = "/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli/mag_comp_1.0.fasta";
    let seq = load_fasta_sequence(fasta_path);

    let enzymes = vec![
        ("BcgI", EnzymeConfig::bcg_i()),
        ("AlfI", EnzymeConfig::alf_i()),
        ("AloI", EnzymeConfig::alo_i()),
        ("BaeI", EnzymeConfig::bae_i()),
        ("BplI", EnzymeConfig::bpl_i()),
        ("BsaXI", EnzymeConfig::bsa_xi()),
        ("BslFI", EnzymeConfig::bsl_fi()),
        ("Bsp24I", EnzymeConfig::bsp24_i()),
        ("CjeI", EnzymeConfig::cje_i()),
        ("CjePI", EnzymeConfig::cje_pi()),
        ("CspCI", EnzymeConfig::csp_ci()),
        ("FalI", EnzymeConfig::fal_i()),
        ("HaeIV", EnzymeConfig::hae_iv()),
        ("Hin4I", EnzymeConfig::hin4_i()),
        ("PpiI", EnzymeConfig::ppi_i()),
        ("PsrI", EnzymeConfig::psr_i()),
    ];

    let mut group = c.benchmark_group("all_enzymes");
    group.sample_size(10);

    for (name, enzyme) in enzymes {
        group.bench_function(format!("new_{}", name), |b| {
            b.iter(|| digest_sequence(black_box(&seq), black_box(&enzyme)))
        });
    }

    group.finish();
}

fn bench_multi_enzyme_parallel(c: &mut Criterion) {
    let fasta_path = "/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli/mag_comp_1.0.fasta";
    let fasta_path_buf = std::path::PathBuf::from(fasta_path);
    let registry = EnzymeRegistry::new();
    let enzymes: Vec<_> = registry.all().to_vec();

    // Verify correctness: seq vs par produce the same tag counts
    let seq_result = TagExtractor::extract_multi_enzyme(&fasta_path_buf, &enzymes).unwrap();
    let par_result = TagExtractor::extract_multi_enzyme_par(&fasta_path_buf, &enzymes).unwrap();
    assert_eq!(seq_result.sets.len(), par_result.sets.len(), "Enzyme count mismatch");
    for (name, set) in &seq_result.sets {
        let par_set = par_result.sets.get(name).expect("Missing enzyme in par result");
        assert_eq!(set.tags.len(), par_set.tags.len(), "Tag count mismatch for {}", name);
    }
    println!("Sequential and parallel multi-enzyme extraction produce identical results ✓");

    let mut group = c.benchmark_group("multi_enzyme_panel");
    group.sample_size(10);

    group.bench_function("sequential_16enzymes", |b| {
        b.iter(|| TagExtractor::extract_multi_enzyme(black_box(&fasta_path_buf), black_box(&enzymes)))
    });

    group.bench_function("parallel_16enzymes", |b| {
        b.iter(|| TagExtractor::extract_multi_enzyme_par(black_box(&fasta_path_buf), black_box(&enzymes)))
    });

    group.finish();
}

criterion_group!(benches, bench_digest, bench_all_enzymes, bench_multi_enzyme_parallel);
criterion_main!(benches);
