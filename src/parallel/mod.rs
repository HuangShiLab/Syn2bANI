//! Thread-pool helpers shared by the batch subcommands.
//!
//! The pairwise estimator itself lives in [`crate::core::chain_ani`]; the
//! batch subcommands (`dist`, `search`, `triangle`, `db search`) drive it
//! through [`crate::cli::compare`]. This module only owns the SIMD helpers.

pub mod simd;
