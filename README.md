# Syn2bANI

> **Strain-level ANI estimation via fixed restriction-site anchors for fragmented metagenome-assembled genomes**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)

## Overview

**Syn2bANI** is a Rust-based bioinformatics tool for estimating Average Nucleotide Identity (ANI) between closely related genomes using **Type IIB restriction enzyme fixed anchors** (2bRAD tags). Unlike k-mer-based methods like skani/FastANI, Syn2bANI leverages the natural positional correspondence of 2bRAD tags to eliminate the costly chaining step, while simultaneously outputting structural variation (SV) and synteny information.

### Core Innovations

1. **Fixed anchors eliminate chaining**: Type IIB restriction sites act as natural positional anchors, replacing random k-mer chaining with O(1) hash-table matching.
2. **ANI + synteny in one pass**: Simultaneously outputs ANI, aligned fraction (AF), structural variations (inversions, indels), and synteny blocks.
3. **Tolerant of fragmentation, with limits**: 2bRAD tags are naturally dispersed short sequences (~32 bp). Measured drift is small down to ~20 kb N50; below that a chain can no longer form inside a single contig and the estimate degrades (ALGORITHM_MLE.md §4.6).
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

### Pairwise ANI, maximum likelihood (`ani`) — recommended

```bash
syn2bani ani query.fasta reference.fasta --verbose -p
```

Chain-restricted, per-enzyme-stratified maximum likelihood. No calibration
model and no empirical offset: ANI comes from fitting a truncated-binomial
likelihood to tag outcomes inside collinear chains, and AF is reported
separately instead of being folded into the ANI estimate.

On simulated ground truth: **MAE 0.074%** over 85–99.9% ANI, and flat as
accessory content varies from 0 to 50%. On 13 real Enterobacteriaceae
chromosomes: **MAE 0.094 vs skani, 0.207 vs FastANI** on the 8 pairs skani also
reports, with the other 5 flagged `BELOW_DETECTION` — the same 5 skani declines
to report.

Note on fragmentation: below ~20 kb N50 the estimate drifts upward more than
skani or FastANI do (+0.65 at 5 kb N50 vs their 0.35–0.38). AF reports the
coverage loss honestly. See ALGORITHM_MLE.md §4.6.

Divergence is modelled with gamma-distributed regional rates, because real
genome pairs are mosaics and a single-rate fit reads systematically high. See
[ALGORITHM_MLE.md](ALGORITHM_MLE.md) for the model, the benchmark tables, and
what is still untested.

### Pairwise ANI over genome sets (`dist`)

```bash
syn2bani dist --ql queries.txt --rl references.txt -o out.tsv
# or the positional form: every path but the last is a query, the last is the reference
syn2bani dist q1.fasta q2.fasta ref.fasta
```

Two stages, the same architecture as skani: a recall-first containment screen
on per-enzyme tag windows rejects pairs that are certainly below the detection
floor, then survivors are refined with the validated chain-restricted MLE
estimator — the identical code path as `ani`, so numbers are byte-identical to
an `ani` run on the same pairs. Output is the `ani` TSV plus a trailing `flag`
column. (This replaced the legacy GBRT-calibrated path, whose screen rejected
94% of true ≥80% ANI pairs and whose debias was miscalibrated near identity.)

### Search against a pre-sketched database (`search`)

```bash
# Step 1: Build sketch database (default panel is now the 4-enzyme ani panel)
syn2bani sketch genomes/*.fasta -o db/

# Step 2: Search (queries may be FASTA or .s2ba)
syn2bani search query.fasta db/ -o results.tsv
```

Screen-then-refine per query against the whole DB; hits are reported best-first
per query with the same estimator and columns as `ani`, filtered by `--min-ani`
(fraction; default 0.8) on the gated estimate.

### All-to-all comparison (`triangle`)

```bash
syn2bani triangle genomes/*.fasta -o matrix.tsv          # matrix, NaN = below floor
syn2bani triangle --ql genomes.txt --edge-list -o edges.tsv
```

All-vs-all lower triangle through the same screen → MLE refine pipeline. The
edge list carries `ani`-style rows for refined pairs; pairs certainly below
the detection floor are omitted (edge list) or written as `NaN` (matrix) —
never `0.0000`.

### Structural variation analysis (`struct`)

```bash
syn2bani struct -q query.fasta -r reference.fasta --rearrangement --indel -o sv.tsv
```

## Algorithm

The estimator is described in [ALGORITHM_MLE.md](ALGORITHM_MLE.md). At
database scale it runs as a two-stage pipeline:

### Stage 1: Containment screen (recall-first)

1. Extract the 4-enzyme panel of Type IIB tags via in-silico digestion (or
   read them from `.s2ba` sketches).
2. Per enzyme, reduce each tag to one strand-canonical key: its centred 18 bp
   window. At 18 bp, calibrated true ≥80% ANI pairs share ≥6 keys (≥29 in the
   80–85% band) while ~83% of random GTDB pairs share 0–2.
3. Merge-intersect the sorted key sets; pass pairs above a calibrated joint
   floor (≥3 shared keys AND ≥0.1% containment of the smaller key set).

### Stage 2: Chain-restricted MLE refine

1. Screen survivors go through `chain_ani::compute` — chaining with an
   adaptive gap, per-enzyme stratified likelihood, gamma rate heterogeneity.
2. Output ANI (heterogeneous, uniform, and gated), AF, synteny blocks,
   breakpoints, and a reliability flag, identical to `ani`.

## CLI Reference

| Subcommand | Description | skani equivalent |
|-----------|-------------|------------------|
| `ani` | Pairwise ANI by chain-restricted maximum likelihood | `skani dist` |
| `dist` | Query × reference ANI (screen + same MLE estimator as `ani`) | `skani dist` |
| `search` | Search queries against sketch database (screen + MLE) | `skani search` |
| `sketch` | Build binary sketch database (default: 4-enzyme panel) | `skani sketch` |
| `triangle` | All-to-all pairwise matrix (screen + MLE) | `skani triangle` |
| `db` | Database management (build, add, remove, list, search, merge) | — |
| `struct` | Structural variation analysis | **Syn2bANI unique** |

### `ani` Options

| Flag | Description | Default |
|------|-------------|---------|
| `--enzymes` | Comma-separated enzyme panel (tags must be ≤ 32 bp) | `BcgI,AlfI,AloI,FalI` |
| `--mismatch-tolerance` | Mismatch budget per tag (`0` = exact only) | `2` |
| `--min-chain-anchors` | Minimum anchors for a trusted chain | `4` |
| `--max-gap` | Max bp between chained anchors | `50000` |
| `--verbose` | Add shape, retention, both partial estimators, chain diagnostics | `false` |

Raising `--mismatch-tolerance` is what extends usable range downward: at 90%
ANI a 32 bp tag matches exactly only 3.4% of the time, leaving too few anchors
to chain.

### Common Options (`dist`, `search`, `triangle`)

| Flag | Description | Default |
|------|-------------|---------|
| `--enzymes` | Comma-separated enzyme panel | `BcgI,AlfI,AloI,FalI` |
| `--threads` / `--parallel` | Thread pool control | auto / off |
| `--min-ani` | Only report pairs with gated ANI ≥ this fraction | `0.0` (`dist`), `0.8` (`search`) |
| `--screen-min-shared` | Screen: minimum shared tag-window keys (AND) | `3` |
| `--screen-min-containment` | Screen: minimum containment of the smaller key set (AND) | `0.001` |
| `--screen-window` | Screen: tag window width (bp) | `18` |
| `--refine-min-approx` | Optional second-tier gate bounding MLE calls | `0.0` (off) |
| `--verbose` | Add shape, retention, partial estimators, chain diagnostics | `false` |

`sketch`/`db build` share `--enzyme` (now accepts a panel; default changed
from BcgI-only to the 4-enzyme panel — a deliberate breaking change; `.s2ba`
files record their enzyme table and stay readable), `--enzymes`, and
`--multi-enzyme`. `db add` digests new genomes with the panel recorded in the
existing database.

## Output Formats

### TSV (`ani`-style; used by `ani`, `dist`, `search`, `triangle --edge-list`, `db search`)

```
query	reference	ani	ani_uniform	af_query	af_reference	std_err	synteny_blocks	synteny_score	breakpoint_count	ani_gated	gate	[flag]
```

`ani_gated` is the recommended estimate (heterogeneous fit with the
disagreement fallback); `gate` records which fit it came from. The database
subcommands append a `flag` column (`ok` / `INCONSISTENT` / `BELOW_DETECTION`).
`triangle` matrix mode writes a full symmetric matrix with `NaN` for pairs
below the detection floor.

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

14x14 all-vs-all (196 pairs), 8 threads, five repeats (`prototype/perf_bench.sh`):

| | min | median | max | peak RSS |
|---|---|---|---|---|
| syn2bani, FASTA input | 0.46 s | 2.52 s | 25.87 s | 343 MB |
| syn2bani, `.s2ba` sketches | **0.34 s** | **0.35 s** | **0.38 s** | **168 MB** |
| skani `dist` | 0.22 s | 2.49 s | 9.19 s | 380 MB |
| skani `triangle` | 0.08 s | 0.43 s | 12.56 s | 260 MB |

Sketch once and reuse — results are bit-identical to the FASTA path:

```bash
syn2bani sketch genomes/*.fasta -o sk --enzymes BcgI,AlfI,AloI,FalI -t 8 -p
syn2bani ani --ql sk_queries.txt --rl sk_refs.txt -t 8 -p
```

That takes the median from 2.52 s to 0.35 s at half the memory, and makes the
runtime reproducible: every FASTA-reading row above, skani's included, spreads by
20-50x between its best and worst run on a machine that is doing anything else,
so treat those medians as indicative only. On best-case runs, sketch-input
syn2bani is within ~1.5x of `skani dist` at under half its peak memory.

Digestion itself runs at ~107 Mb/s single-threaded, and a `.s2ba` sketch is
~120 KB per genome (1.7 MB for 14).

See [`BENCHMARK_REPORT.md`](BENCHMARK_REPORT.md) for the older `dist`-path
benchmarks.

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
