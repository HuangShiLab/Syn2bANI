//! Syn2bANI — Strain-level ANI estimation via Type IIB restriction-site anchors.
//!
//! This library implements the Syn2b algorithm for in-silico digestion of genomic
//! sequences with Type IIB restriction enzymes, producing fixed-length 2bRAD tags
//! for downstream genome comparison and taxonomic profiling.

pub mod core;
pub mod enzyme;
pub mod io;
pub mod cli;
pub mod parallel;
pub mod utils;

pub use enzyme::registry::{EnzymeConfig, EnzymeRegistry};
pub use enzyme::digest::{Tag, Direction, digest_sequence, digest_sequence_legacy};
pub use enzyme::EnzymeError;
pub use core::TagExtractor;
