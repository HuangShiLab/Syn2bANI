//! Fast hash utilities using FxHasher (the same hasher used by rustc).
//!
//! FxHasher is a simple multiplicative hash that trades cryptographic
//! security for raw speed. It is ideal for `HashMap<u64, _>` where the
//! keys are trusted / not attacker-controlled.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// A `HashMap` using `FxHasher`.
pub type FastHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// The FxHasher used internally by the Rust compiler.
///
/// Reference: https://github.com/rust-lang/rustc-hash
#[derive(Default, Clone, Copy, Debug)]
pub struct FxHasher {
    hash: usize,
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, mut bytes: &[u8]) {
        #[cfg(target_pointer_width = "64")]
        const BYTES_WORD: usize = 8;
        #[cfg(target_pointer_width = "32")]
        const BYTES_WORD: usize = 4;

        while bytes.len() >= BYTES_WORD {
            self.add_to_hash(read_usize(bytes));
            bytes = &bytes[BYTES_WORD..];
        }
        if BYTES_WORD == 8 && bytes.len() >= 4 {
            self.add_to_hash(u32::from_le_bytes([
                bytes[0], bytes[1], bytes[2], bytes[3],
            ]) as usize);
            bytes = &bytes[4..];
        }
        for &byte in bytes {
            self.add_to_hash(byte as usize);
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add_to_hash(i as usize);
    }

    #[inline]
    fn write_u16(&mut self, i: u16) {
        self.add_to_hash(i as usize);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add_to_hash(i as usize);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add_to_hash(i as usize);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add_to_hash(i);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl FxHasher {
    #[inline]
    fn add_to_hash(&mut self, i: usize) {
        #[cfg(target_pointer_width = "64")]
        const K: usize = 0x517cc1b727220a95;
        #[cfg(target_pointer_width = "32")]
        const K: usize = 0x9e3779b9;

        self.hash = (self.hash.rotate_left(5) ^ i).wrapping_mul(K);
    }
}

#[inline]
fn read_usize(bytes: &[u8]) -> usize {
    assert!(bytes.len() >= std::mem::size_of::<usize>());
    let mut buf = [0u8; std::mem::size_of::<usize>()];
    buf.copy_from_slice(&bytes[..std::mem::size_of::<usize>()]);
    usize::from_le_bytes(buf)
}
