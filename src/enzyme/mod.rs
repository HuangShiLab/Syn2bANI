//! Enzyme module — Type IIB restriction enzyme registry and in-silico digestion.

pub mod registry;
pub mod digest;

pub use registry::{EnzymeConfig, EnzymeRegistry};
pub use digest::{Tag, Direction, digest_sequence};

use thiserror::Error;

/// Errors that can occur during enzyme operations.
#[derive(Error, Debug)]
pub enum EnzymeError {
    #[error("Unknown enzyme: {0}")]
    UnknownEnzyme(String),
    #[error("Invalid enzyme configuration: {0}")]
    InvalidConfig(String),
    #[error("Sequence too short for digestion: length {0}")]
    SequenceTooShort(usize),
}
