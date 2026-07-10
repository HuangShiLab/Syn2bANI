# Syn2bANI

> **Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)

## Overview

**Syn2bANI** is a Rust-based bioinformatics tool for estimating Average Nucleotide Identity (ANI) between closely related genomes using **Type IIB restriction enzyme fixed anchors** (2bRAD tags). Unlike k-mer-based methods like skani/FastANI, Syn2bANI leverages the natural positional correspondence of 2bRAD tags to eliminate the costly chaining step, while simultaneously outputting structural variation (SV) and synteny information.

### Core Innovations

1. **Fixed anchors eliminate chaining**: Type IIB restriction sites act as natural positional anchors, replacing random k-mer chaining with O(1) hash-table matching.
2. **ANI + synteny in one pass**: Simultaneously outputs ANI, aligned fraction (AF), structural variations (inversions, indels), and synteny blocks.
3. **Robust to extreme fragmentation**: 2bRAD tags are naturally dispersed short sequences (~32 bp), making Syn2bANI robust against highly fragmented MAGs (N50 < 10 kb).
4. **Experimentally verifiable**: Predicted tags can be directly validated by 2bRAD-M sequencing.
5. **GBRT debiasing**: An embedded Gradient Boosted Regression Tree model corrects systematic ANI overestimation, achieving <0.3% cross-species MAE.

## Installation

### From source (requires Rust ≥ 1.70)

```bash
git clone https://github.com/HuangShiLab/Syn2bANI.git
cd Syn2bANI
cargo build --release
```

The binary will be available at `target/release/syn2bani`.

### Pre-built binaries

Coming soon — see [Releases](https://github.com/HuangShiLab/Syn2bANI/releases).

## Quick Start

### Pairwise ANI (`dist`)

```bash
syn2bani dist -q query.fasta -r reference.fasta
```

### Search against a pre-sketched database (`search`)

```bash
# Step 1: Build sketch database
syn2bani sketch genomes/*.fasta -o db/

# Step 2: Search
syn2bani search -q query.fasta -d db/ -o results.tsv
```

### All-to-all comparison (`triangle`)

```bash
syn2bani triangle genomes/*.fasta -o matrix.tsv
```

### Structural variation analysis (`struct`)

```bash
syn2bani struct -q query.fasta -r reference.fasta --rearrangement --indel -o sv.tsv
```

## Algorithm

Syn2bANI implements a **two-pass fixed-anchor algorithm**:

### Pass 1: Coarse Screening

1. Extract 2bRAD tags from query and reference genomes via in-silico digestion (Fast2bRAD-M style).
2. Count exact shared tags.
3. Estimate coarse ANI using max-containment similarity.
4. Skip distant pairs (ANI < 80%) to save computation.

### Pass 2: Fine ANI Calculation

1. Build a hash index of reference tags.
2. Match query tags against the index (O(1) per tag).
3. Compute Hamming distance for each matched tag pair (32 bp → local ANI).
4. Build synteny blocks from consecutive matched tags.
5. Detect orientation flips (inversions) and gap differences (indels).
6. Compute weighted ANI with optional GBRT debiasing correction.
7. Output ANI, AF, synteny blocks, and structural variations.

## CLI Reference

| Subcommand | Description | skani equivalent |
|-----------|-------------|------------------|
| `dist` | Pairwise ANI between query and reference | `skani dist` |
| `search` | Search query against sketch database | `skani search` |
| `sketch` | Build binary sketch database | `skani sketch` |
| `triangle` | All-to-all pairwise matrix | `skani triangle` |
| `db` | Database management (build, add, remove, list, merge) | — |
| `struct` | Structural variation analysis | **Syn2bANI unique** |

### Common Options

| Flag | Description | Default |
|------|-------------|---------|
| `--enzyme` | Type IIB enzyme to use | `BcgI` |
| `--multi-enzyme` | Use all 16 enzymes for higher coverage | `false` |
| `--threads` | Number of parallel threads | `1` |
| `--structural` | Output structural variation info | `false` |
| `--min-ani` | Minimum ANI threshold for reporting | `80.0` |
| `--gbrt-debias` | Use GBRT model for ANI correction | `true` |

## Output Formats

### TSV (default, skani-compatible)

```
query_file	ref_file	ani	af_q	af_r	query_name	ref_name	shared_tags	sv_count
```

### Extended TSV (with `--structural`)

Adds rearrangement, indel, and synteny block counts.

### JSON (machine-readable)

Full output including per-tag ANI profiles, synteny blocks, and structural variations.

## Architecture

```
syn2bani/
├── src/
│   ├── main.rs              # CLI entry point
│   ├── lib.rs               # Library exports
│   ├── cli/                 # Command handlers (dist, search, sketch, ...)
│   ├── core/                # Core engine
│   │   ├── tag_extractor.rs # In-silico Type IIB digestion (Fast2bRAD-M aligned)
│   │   ├── tag_matcher.rs   # Fixed-anchor hash matching
│   │   ├── ani_calculator.rs# Weighted ANI + GBRT debiasing
│   │   ├── synteny_builder.rs# Synteny block construction
│   │   ├── structure_analyzer.rs # SV detection
│   │   ├── gbrt.rs          # Embedded GBRT model inference
│   │   └── debias.rs        # Simple ANI correction
│   ├── enzyme/              # Enzyme registry & digestion
│   ├── io/                  # FASTA parser, sketch format, output formatters
│   ├── parallel/            # Rayon-based parallelism
│   └── utils/               # Sequence utilities
├── tests/                   # Integration tests
└── benches/                 # Criterion performance benchmarks
```

## Supported Enzymes

All 16 Type IIB enzymes from the 2bRAD-M panel:

| Enzyme | Recognition Pattern | Tag Length |
|--------|-------------------|------------|
| BcgI | CGA-N6-TGC | 32 bp |
| AlfI | GCA-N6-TGC | 32 bp |
| AloI | GAAC-N6-TCC | 27 bp |
| BaeI | AC-N4-GTAYC | 28 bp |
| BplI | GAG-N5-CTC | 27 bp |
| BsaXI | AC-N5-CTCC | 27 bp |
| BslFI | GGGAC | 21 bp |
| Bsp24I | GAC-N6-TGG | 27 bp |
| CjeI | CCA-N6-GT | 28 bp |
| CjePI | CCA-N7-TC | 27 bp |
| CspCI | CAA-N5-GTGG | 33 bp |
| FalI | AAG-N5-CTT | 27 bp |
| HaeIV | GAY-N5-RTC | 27 bp |
| Hin4I | GAY-N5-VTC | 27 bp |
| PpiI | GAAC-N5-CTC | 28 bp |
| PsrI | GAAC-N6-TAC | 27 bp |

## Performance

| Metric | Value | Notes |
|--------|-------|-------|
| Digestion speed | ~107 Mb/s (single thread) | Fast2bRAD-M optimized, BcgI |
| 16-enzyme panel | ~600 ms (4.65 Mb genome) | Single thread |
| Sketch size | ~48 KB per genome | 2-bit packed sequences |
| Memory | ~1 GB for 65k genomes | Compact binary format |

See [`BENCHMARK_REPORT.md`](BENCHMARK_REPORT.md) for detailed benchmark results.

## Citation

If you use Syn2bANI in your research, please cite:

> **Syn2b-ANI: Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes**

## Related Projects

- [Syn2b](https://github.com/HuangShiLab/Syn2b) — Synteny analysis using 2bRAD tags
- [Fast2bRAD-M](https://github.com/HuangShiLab/Fast2bRAD-M) — Fast 2bRAD tag extraction
- [skani](https://github.com/bluenote-1577/skani) — Reference k-mer chaining ANI tool

## License

MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

This project builds upon the 2bRAD-M framework (HuangShiLab) and is inspired by the skani algorithm (Shaw & Yu, *Nature Methods* 2023).
# Syn2bANI
