use clap::{Parser, Subcommand};
use std::path::PathBuf;

pub mod dist;
pub mod sketch;
pub mod search;
pub mod triangle;
pub mod db;
pub mod r#struct;

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
        #[arg(short, long, default_value = "1")]
        threads: usize,
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
        #[arg(short, long, default_value = "1")]
        threads: usize,
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
        #[arg(short, long, default_value = "1")]
        threads: usize,
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
        #[arg(short, long, default_value = "1")]
        threads: usize,
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
