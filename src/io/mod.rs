pub mod fasta_parser;
pub mod sketch;
pub mod sketch_reader;
pub mod formatters;

pub use fasta_parser::{FastaRecord, parse_fasta, IoError};
pub use sketch::{
    TgtSketch, ChromSketch, SketchTag, SketchMetadata,
    pack_sequence, unpack_sequence, write_sketch,
    S2BA_MAGIC, S2BA_VERSION,
};
pub use sketch_reader::{read_sketch, sketch_tags_to_genome_tags};
pub use formatters::{TsvFormatter, ExtendedTsvFormatter, JsonFormatter, PafFormatter, OutputFormat};
