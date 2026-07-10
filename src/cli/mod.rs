use clap::{Parser, Subcommand};
use std::path::PathBuf;

pub mod dist;
pub mod sketch;
pub mod search;
pub mod triangle;
pub mod db;
pub mod r#struct;

/// Build a rayon thread pool according to CLI parallel / threads flags.
///
/// - `parallel=false` → single-thread pool (no parallelism)
/// - `parallel=true, threads=0` → use rayon default (all logical cores)
/// - `parallel=true, threads=N` → pool with exactly N threads
pub fn build_pool(parallel: bool, threads: usize) -> Result<rayon::ThreadPool, rayon::ThreadPoolBuildError> {
    let mut builder = rayon::ThreadPoolBuilder::new();
    if !parallel {
        builder = builder.num_threads(1);
    } else if threads > 0 {
        builder = builder.num_threads(threads);
    }
    builder.build()
}

#[derive(Parser)]
#[command(name = "syn2bani")]
#[command(about = "Strain-level ANI estimation via Type IIB restriction-site anchors")]
#[command(version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    Dist {
        #[arg(required = true)]
        query: Vec<PathBuf>,
        #[arg(required = true)]
        reference: Vec<PathBuf>,
        #[arg(short, long, default_value = "BcgI")]
        enzyme: String,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long)]
        multi_enzyme: bool,
        #[arg(long)]
        structural: bool,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    Search {
        #[arg(required = true)]
        query: Vec<PathBuf>,
        #[arg(required = true)]
        database: PathBuf,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(short, long, default_value = "0.8")]
        min_ani: f64,
    },
    Sketch {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
        #[arg(short, long, default_value = "BcgI")]
        enzyme: String,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long)]
        multi_enzyme: bool,
    },
    Triangle {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long)]
        edge_list: bool,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
    },
    Db {
        #[command(subcommand)]
        command: DbCommands,
    },
    Struct {
        #[arg(required = true)]
        query: Vec<PathBuf>,
        #[arg(required = true)]
        reference: Vec<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long)]
        paf: bool,
        #[arg(long)]
        rearrangement: bool,
        #[arg(long)]
        indel: bool,
    },
}

#[derive(Subcommand)]
pub enum DbCommands {
    Build {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
        #[arg(short, long, default_value = "BcgI")]
        enzyme: String,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long)]
        multi_enzyme: bool,
    },
    Add {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        database: PathBuf,
    },
    Remove {
        #[arg(required = true)]
        genome_ids: Vec<String>,
        #[arg(short, long, required = true)]
        database: PathBuf,
    },
    List {
        #[arg(short, long, required = true)]
        database: PathBuf,
    },
    Merge {
        #[arg(required = true)]
        databases: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
    },
}
