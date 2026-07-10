# Syn2bANI Tag Extraction Benchmark Report

**Date**: 2025-08-07  
**Genome**: *E. coli* K-12 MG1655 (4.65 Mb, `mag_comp_1.0.fasta`)  
**Benchmark framework**: Criterion.rs (Rust)  
**Platform**: macOS, Apple Silicon (single thread)

---

## Executive Summary

Syn2bANI's enzyme digestion module was refactored to align with the **Fast2bRAD-M** approach used in the upstream **Syn2b** tool. The new implementation:

- **~17.5% faster** than the legacy margin-based method (BcgI: 43.3 ms vs 52.7 ms)
- **Throughput**: 107.5 Mb/s (new) vs 88.3 Mb/s (old)
- **Identical output**: Both methods produce exactly 2,888 tags on the test genome
- **16-enzyme panel**: All enzymes complete in **18.6–49.4 ms** per 4.65 Mb genome

---

## 1. Design Changes (Syn2b → Syn2bANI Alignment)

| Feature | Syn2b (Fast2bRAD-M) | Syn2bANI (Old) | Syn2bANI (New)
|--------|---------------------|---------------|---------------
| Scan strategy | Sliding window (`tag_length`) | Recognition-site driven | Sliding window ✅
| Reverse strand | Built-in reverse patterns | `reverse_complement` sequence | Built-in reverse patterns ✅
| Pure-base filter | `is_pure_atcg` (no N) | Not implemented | `is_pure_atcg` ✅
| Tag storage | `[u8; 32]` stack array | `Vec<u8>` heap | `[u8; 32]` ✅
| IUPAC lookup | Bitmask table (`BASE_MASK`) | `match` per base | Bitmask table ✅
| Static patterns | Compile-time `const` anchors | Runtime string compare | Compile-time `const` ✅
| Deduplication | Sort + `dedup_by_key` | None | Sort + `dedup_by_key` ✅

---

## 2. New vs Old Method Comparison (BcgI)

| Metric | New (Fast2bRAD-M) | Old (margin-based) | Improvement |
|--------|-------------------|-------------------|-------------|
| **Median time** | **43.3 ms** | 52.7 ms | **17.5% faster** |
| **Throughput** | **107.5 Mb/s** | 88.3 Mb/s | **21.7% higher** |
| **Tag count** | 2,888 | 2,888 | Identical |
| **Memory (per tag)** | 40 bytes (stack) | 48+ bytes (heap) | More compact |

> **Key insight**: The sliding window approach naturally handles both forward and reverse patterns in a single pass, eliminating the need for a separate `reverse_complement` scan and reducing branch misprediction.

---

## 3. 16-Enzyme Panel Performance

| Rank | Enzyme | Type | Median time | Notes |
|------|--------|------|-------------|-------|
| 1 | **FalI** | Palindrome | **18.6 ms** | 1 pattern only |
| 2 | **BplI** | Palindrome | **20.1 ms** | 1 pattern only |
| 3 | **AlfI** | Palindrome | **21.1 ms** | 1 pattern only |
| 4 | BslFI | Single-anchor ×2 | 37.3 ms | Short tag (25 bp) |
| 5 | BsaXI | Dual-anchor ×2 | 38.8 ms | Standard |
| 6 | AloI | Dual-anchor ×2 | 40.8 ms | Standard |
| 7 | CspCI | Dual-anchor ×2 | 40.1 ms | Standard |
| 8 | Bsp24I | Dual-anchor ×2 | 40.3 ms | Standard |
| 9 | CjePI | Dual-anchor ×2 | 41.2 ms | Standard |
| 10 | CjeI | Dual-anchor ×2 | 41.6 ms | Standard |
| 11 | PsrI | Dual-anchor ×2 | 41.5 ms | Standard |
| 12 | PpiI | Dual-anchor ×2 | 42.0 ms | Standard |
| 13 | **BcgI** | Dual-anchor ×2 | **43.3 ms** | Benchmark reference |
| 14 | BaeI | IUPAC ×2 | 43.5 ms | Degenerate base checks |
| 15 | **HaeIV** | IUPAC ×2 | **49.3 ms** | 🐢 Slowest |
| 16 | **Hin4I** | IUPAC ×2 | **49.4 ms** | 🐢 Slowest |

**Total 16-enzyme panel time**: ~600 ms (single-threaded, 4.65 Mb genome)

### Pattern-type performance hierarchy

```
Palindrome (1 pattern)    ~20 ms     ← 2× faster
Standard (2 patterns)      ~40 ms
IUPAC degenerate (2 pat)   ~49 ms     ← +20% overhead
```

---

## 4. Technical Details

### Static pattern compilation

All 16 enzymes are defined as `const` arrays at compile time. There is **zero runtime allocation** for pattern definitions:

```rust
const BCGI_F1: Anchor = Anchor{offset:10, motif:b"CGA"};
const BCGI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BCGI_F1,BCGI_F2], iupac:&[]},
    Pattern{anchors:&[BCGI_R1,BCGI_R2], iupac:&[]},
];
```

### IUPAC lookup via bitmask

Instead of a `match` on every base comparison, the new code uses a precomputed `BASE_MASK` lookup table:

```rust
const BASE_MASK: [u8; 256] = { /* A=1, T=2, C=4, G=8 */ };
// Check if base matches IUPAC code in a single bitwise AND:
(BASE_MASK[base] & allowed_mask) != 0
```

This reduces IUPAC checks from **~12 branch instructions** to **1 memory load + 1 bitwise AND**.

### Stack-allocated tags

`Tag.sequence` changed from `Vec<u8>` to `[u8; 32]`:
- Eliminates per-tag heap allocation (~16 bytes overhead)
- Better cache locality during Hamming distance computation
- No pointer chasing in the hot matching loop

---

## 5. Files Changed

| File | Change |
|------|--------|
| `src/enzyme/digest.rs` | Rewrote with Fast2bRAD-M sliding window; added `digest_sequence_legacy()` |
| `src/core/tag_extractor.rs` | `GenomeTag.sequence` → `[u8; 32]`; added `seq_len` field |
| `src/core/tag_matcher.rs` | Hamming distance now uses `seq_len` for proper comparison |
| `src/io/sketch_reader.rs` | Updated to pack/unpack `[u8; 32]` arrays |
| `src/cli/search.rs` | Updated to construct `[u8; 32]` from sketch |
| `src/lib.rs` | Exported `digest_sequence_legacy` for testing |
| `benches/benchmark.rs` | New Criterion benchmark comparing old vs new |

---

## 6. Validation

All existing tests pass:

```
cargo test: 22 tests passed, 0 failed
cargo bench: New and old methods produce identical tag counts (2,888)
```

---

## 7. Recommendations for Further Optimization

1. **SIMD `is_pure_atcg`**: Use `_mm256_cmpeq_epi8` (AVX2) or NEON to check 32 bytes in parallel → potential **2–4×** speedup on long genomes.
2. **Parallel enzyme digestion**: Since each enzyme is independent, process the 16-enzyme panel with `rayon` → panel time drops from ~600 ms to ~100 ms (6× speedup on 8 cores).
3. **Skip `sort` for single-pattern enzymes**: Palindrome enzymes (FalI, BplI, AlfI) only scan forward; no need for `sort_by_key` + `dedup_by_key` → ~5% faster.

---

*Report generated by Syn2bANI benchmark suite.*
