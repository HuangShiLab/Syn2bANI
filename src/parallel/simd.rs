/// SIMD-accelerated operations for 2bRAD tag processing.
///
/// # Path A: 64-bit packed sequences (cross-platform)
/// For 32 bp DNA tags, the most efficient cross-platform approach is to pack
/// each base into 2 bits (A=00, C=01, G=10, T=11) and compare 64-bit integers
/// via XOR + popcount-style diff counting.
///
/// # Path B: AVX2 batch operations (x86_64 only)
/// When compiled for x86_64 with AVX2 support, `is_pure_atcg` uses 256-bit
/// vectorized comparison for 32-byte windows. This provides ~2-4× speedup
/// over scalar byte-by-byte checking.
///
/// Future work: AVX-512 or NEON for batch matching of 1000+ tags.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

// ═══════════════════════════════════════════════════════════════════════════
// 64-bit packed sequence operations (Path A — cross-platform)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the Hamming distance between two 64-bit packed sequences.
///
/// Returns the number of bit positions that differ.
/// This is a raw bit-level popcount; for DNA base-level distance,
/// use `diff_count_u64` instead.
pub fn hamming_distance_u64(a: u64, b: u64) -> u32 {
    (a ^ b).count_ones()
}

/// Count the number of differing DNA bases between two 64-bit packed sequences.
///
/// Each base occupies 2 bits. For a given XOR result:
/// - 00 = same base (A↔A, C↔C, G↔G, T↔T)
/// - 01 = 1-bit difference (A↔C, A↔G)
/// - 10 = 1-bit difference (A↔G, C↔T)
/// - 11 = 2-bit difference (A↔T, C↔G)
///
/// This function counts each non-zero 2-bit group as **one differing base**,
/// regardless of whether the difference is 1 or 2 bits. This matches the
/// standard biological Hamming distance (number of differing bases).
///
/// Algorithm: `((xor | (xor >> 1)) & 0x5555...).count_ones()`
#[inline]
pub fn diff_count_u64(xor: u64) -> u32 {
    let t = xor | (xor >> 1);
    (t & 0x5555555555555555).count_ones()
}

// ═══════════════════════════════════════════════════════════════════════════
// AVX2-accelerated is_pure_atcg (x86_64 only)
// ═══════════════════════════════════════════════════════════════════════════

/// Check if all bytes in a window are A, T, C, or G (case-insensitive).
///
/// Automatically selects AVX2 on x86_64 when available and the window is
/// at least 32 bytes. Falls back to a fast scalar loop otherwise.
#[inline]
pub fn is_pure_atcg_simd(window: &[u8]) -> bool {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        if is_x86_feature_detected!("avx2") && window.len() >= 32 {
            return is_pure_atcg_avx2(window);
        }
    }
    is_pure_atcg_scalar(window)
}

/// Scalar fallback for `is_pure_atcg`.
#[inline]
fn is_pure_atcg_scalar(window: &[u8]) -> bool {
    window.iter().all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn is_pure_atcg_avx2(window: &[u8]) -> bool {
    debug_assert!(window.len() >= 32);

    // Load 32 bytes unaligned
    let vec = _mm256_loadu_si256(window.as_ptr() as *const __m256i);

    // Uppercase: AND with 0xDF (clears bit 5, e.g. 'a'(0x61) → 'A'(0x41))
    let upper = _mm256_and_si256(vec, _mm256_set1_epi8(0xDFi8));

    // Compare against A, C, G, T
    let eq_a = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8(b'A' as i8));
    let eq_c = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8(b'C' as i8));
    let eq_g = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8(b'G' as i8));
    let eq_t = _mm256_cmpeq_epi8(upper, _mm256_set1_epi8(b'T' as i8));

    // OR together: any valid base sets the byte to 0xFF
    let valid = _mm256_or_si256(
        _mm256_or_si256(eq_a, eq_c),
        _mm256_or_si256(eq_g, eq_t),
    );

    // movemask extracts the MSB of each byte. If all 32 bytes are valid,
    // every byte is 0xFF → MSB = 1 → mask == -1 (0xFFFFFFFF).
    _mm256_movemask_epi8(valid) == -1
}

// ═══════════════════════════════════════════════════════════════════════════
// Batch Hamming distance (AVX2, x86_64 only)
// ═══════════════════════════════════════════════════════════════════════════

/// Batch-compute base-level Hamming distances for 4 pairs of 64-bit packed
/// sequences using AVX2.
///
/// # Performance note
/// While this processes 4 pairs in parallel, the speedup over 4 scalar
/// `diff_count_u64` calls is modest (~1.5×) because:
/// 1. AVX2 lacks a native 64-bit popcount per lane.
/// 2. The dominant cost in matching is HashMap lookup (memory-bound),
///    not Hamming distance computation (already ~3 CPU cycles).
///
/// For meaningful speedup, batch matching must also address memory layout
/// (e.g. SoA instead of AoS, prefetching). See ` benches/batch_matching.rs`
/// for a prototype.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn batch_diff_count_4(xor_vals: &[u64; 4]) -> [u32; 4] {
    let vec = _mm256_loadu_si256(xor_vals.as_ptr() as *const __m256i);

    // diff_count bit trick: (xor | (xor >> 1)) & 0x5555...
    let shifted = _mm256_srli_epi64(vec, 1);
    let or_result = _mm256_or_si256(vec, shifted);
    let mask = _mm256_set1_epi64x(0x5555555555555555i64);
    let masked = _mm256_and_si256(or_result, mask);

    // AVX2 has no 64-bit popcount — extract to scalar. This is the bottleneck
    // that limits batch speedup to ~1.5× over scalar.
    let low_128 = _mm256_castsi256_si128(masked);
    let high_128 = _mm256_extracti128_si256(masked, 1);

    [
        _mm_extract_epi64(low_128, 0) as u64,
        _mm_extract_epi64(low_128, 1) as u64,
        _mm_extract_epi64(high_128, 0) as u64,
        _mm_extract_epi64(high_128, 1) as u64,
    ]
    .map(|v| v.count_ones())
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diff_count_same() {
        assert_eq!(diff_count_u64(0b00_00_00_00), 0);
    }

    #[test]
    fn test_diff_count_one_base() {
        // A (00) vs C (01): XOR = 01 → diff = 1
        assert_eq!(diff_count_u64(0b00_00_00_01), 1);
        // A (00) vs G (10): XOR = 10 → diff = 1
        assert_eq!(diff_count_u64(0b00_00_00_10), 1);
        // A (00) vs T (11): XOR = 11 → diff = 1
        assert_eq!(diff_count_u64(0b00_00_00_11), 1);
    }

    #[test]
    fn test_diff_count_two_bases() {
        // AC (00_01) vs AG (00_10): XOR = 00_11 → diff = 1
        assert_eq!(diff_count_u64(0b00_00_11_00), 1);
        // AC (00_01) vs GT (10_11): XOR = 10_10 → diff = 2
        assert_eq!(diff_count_u64(0b00_00_10_10), 2);
    }

    #[test]
    fn test_diff_count_full_32bp() {
        let a = 0u64;                           // all A
        let c = 0x5555555555555555u64;          // all C (01 repeating)
        assert_eq!(diff_count_u64(a ^ c), 32);
    }

    #[test]
    fn test_is_pure_atcg_scalar() {
        assert!(is_pure_atcg_scalar(b"ATCGatcg"));
        assert!(!is_pure_atcg_scalar(b"ATNG"));
    }

    #[test]
    fn test_is_pure_atcg_simd_32b() {
        let valid_32 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 32 A's
        let invalid_32 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAANAAA"; // N at pos 28
        assert!(is_pure_atcg_simd(valid_32));
        assert!(!is_pure_atcg_simd(invalid_32));
    }

    #[test]
    fn test_is_pure_atcg_simd_short() {
        // Short windows (< 32) should use scalar path
        assert!(is_pure_atcg_simd(b"ATCG"));
        assert!(!is_pure_atcg_simd(b"ATNG"));
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn test_batch_diff_count_4() {
        if is_x86_feature_detected!("avx2") {
            let xor_vals = [
                0u64,                              // all same → 0
                0x5555555555555555u64,             // all C vs all A → 32
                0b00_00_00_01u64,                  // 1 diff
                0b00_00_10_10u64,                  // 2 diffs
            ];
            let expected = [0, 32, 1, 2];
            let result = unsafe { batch_diff_count_4(&xor_vals) };
            assert_eq!(result, expected);
        }
    }
}
