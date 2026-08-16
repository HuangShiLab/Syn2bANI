use clap::{Parser, Subcommand};
use std::path::PathBuf;

pub mod ani;
pub mod compare;
pub mod panel;
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
    /// Chain-restricted ANI by maximum likelihood (no calibration model).
    Ani {
        /// Positional form: every path but the last is a query, the last is the
        /// reference. Two greedy lists cannot be split any other way, so use
        /// --ql/--rl when you need more than one reference.
        query: Vec<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing query paths, one per line")]
        ql: Option<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing reference paths, one per line")]
        rl: Option<PathBuf>,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Comma-separated enzyme panel (tags must be <= 32 bp)")]
        enzymes: String,
        #[arg(long, default_value = "2",
               help = "Mismatch budget per tag; 0 = exact match only")]
        mismatch_tolerance: usize,
        #[arg(long, default_value = "4", help = "Minimum anchors for a trusted chain")]
        min_chain_anchors: usize,
        #[arg(long, default_value = "50000", help = "Max bp between chained anchors")]
        max_gap: usize,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long, help = "Also report the two partial estimators and chain diagnostics")]
        verbose: bool,
        #[arg(long, value_name = "FILE",
               help = "Dump per-enzyme sufficient statistics for `panel` re-scoring")]
        strata_out: Option<PathBuf>,
        #[arg(long, help = "Apply embedded linear calibration model; adds an ani_cal column")]
        calibrate: bool,
        #[arg(long, value_name = "FILE",
               help = "Path to a custom calibration model JSON (requires --calibrate)")]
        calibrate_model: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Choose an enzyme panel from statistics already written by `ani --strata-out`.
    Panel {
        #[arg(long, value_name = "FILE", required = true,
               help = "Per-enzyme statistics from `ani --strata-out`")]
        strata: PathBuf,
        #[arg(long, value_name = "FILE", required = true,
               help = "query<TAB>reference<TAB>ani reference values")]
        truth: PathBuf,
        #[arg(long, help = "Run greedy forward selection")]
        greedy: bool,
        #[arg(long, value_name = "LIST",
               help = "Semicolon-separated panels to score, e.g. 'BcgI,AloI;BcgI,AlfI,AloI'")]
        panels: Option<String>,
    },
    /// Pairwise ANI over query x reference sets (screen + chain-restricted MLE).
    Dist {
        /// Positional form: every path but the last is a query, the last is the
        /// reference (same convention as `ani`). Use --ql/--rl for lists.
        query: Vec<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing query paths, one per line")]
        ql: Option<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing reference paths, one per line")]
        rl: Option<PathBuf>,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Comma-separated enzyme panel (tags must be <= 32 bp)")]
        enzymes: String,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long, help = "Also report the partial estimators and chain diagnostics")]
        verbose: bool,
        #[arg(long, default_value = "0.0",
               help = "Only report pairs with gated ANI >= this (0 = report all refined pairs)")]
        min_ani: f64,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_SHARED,
               help = "Screen: pass pairs with at least this many shared tag-window keys")]
        screen_min_shared: usize,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_CONTAINMENT,
               help = "Screen: also pass pairs reaching this shared-key containment")]
        screen_min_containment: f64,
        #[arg(long, default_value_t = crate::core::screen::SCREEN_WINDOW,
               help = "Screen: centred tag window width (bp)")]
        screen_window: usize,
        #[arg(long, default_value = "0.0",
               help = "Second-tier gate: only refine screen survivors whose crude \
                       containment-ANI reaches this (0 = refine all survivors)")]
        refine_min_approx: f64,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Search query genomes against a sketch database (screen + chain-restricted MLE).
    Search {
        /// Query genomes (FASTA or .s2ba), or use --ql.
        query: Vec<PathBuf>,
        /// Sketch database directory (or use --rl with a file of sketch paths).
        database: Option<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing query paths, one per line")]
        ql: Option<PathBuf>,
        #[arg(long, value_name = "FILE",
               help = "File listing database sketch paths (alternative to a directory)")]
        rl: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(short, long, default_value = "0.8",
               help = "Report hits with gated ANI >= this fraction (0.8 = 80%)")]
        min_ani: f64,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Comma-separated enzyme panel for digesting FASTA queries")]
        enzymes: String,
        #[arg(long, help = "Also report the partial estimators and chain diagnostics")]
        verbose: bool,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_SHARED,
               help = "Screen: pass pairs with at least this many shared tag-window keys")]
        screen_min_shared: usize,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_CONTAINMENT,
               help = "Screen: also pass pairs reaching this shared-key containment")]
        screen_min_containment: f64,
        #[arg(long, default_value_t = crate::core::screen::SCREEN_WINDOW,
               help = "Screen: centred tag window width (bp)")]
        screen_window: usize,
        #[arg(long, default_value = "0.0",
               help = "Second-tier gate: only refine screen survivors whose crude \
                       containment-ANI reaches this (0 = refine all survivors)")]
        refine_min_approx: f64,
    },
    /// Build binary sketch files (.s2ba) from genomes.
    ///
    /// BREAKING CHANGE: the default panel is now the validated 4-enzyme panel
    /// (BcgI,AlfI,AloI,FalI), matching `ani`; it used to be BcgI-only.
    /// Sketches record their enzyme table, so old sketches stay readable.
    Sketch {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Enzyme (or comma-separated panel); default changed from BcgI-only")]
        enzyme: String,
        #[arg(long, help = "Comma-separated panel; use the same list as `ani`")]
        enzymes: Option<String>,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(long)]
        multi_enzyme: bool,
    },
    /// All-vs-all pairwise ANI (lower triangle), screen + chain-restricted MLE.
    Triangle {
        /// Genomes (FASTA or .s2ba), or use --ql.
        genomes: Vec<PathBuf>,
        #[arg(long, value_name = "FILE", help = "File listing genome paths, one per line")]
        ql: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long, help = "Emit ani-style rows for refined pairs instead of a matrix")]
        edge_list: bool,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Comma-separated enzyme panel (tags must be <= 32 bp)")]
        enzymes: String,
        #[arg(long, help = "Also report the partial estimators and chain diagnostics")]
        verbose: bool,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_SHARED,
               help = "Screen: pass pairs with at least this many shared tag-window keys")]
        screen_min_shared: usize,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_CONTAINMENT,
               help = "Screen: also pass pairs reaching this shared-key containment")]
        screen_min_containment: f64,
        #[arg(long, default_value_t = crate::core::screen::SCREEN_WINDOW,
               help = "Screen: centred tag window width (bp)")]
        screen_window: usize,
        #[arg(long, default_value = "0.0",
               help = "Second-tier gate: only refine screen survivors whose crude \
                       containment-ANI reaches this (0 = refine all survivors)")]
        refine_min_approx: f64,
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
        #[arg(long)]
        multi_enzyme: bool,
        #[arg(long, help = "Comma-separated enzyme list (overrides --multi-enzyme)")]
        enzymes: Option<String>,
        #[arg(long, default_value = "1000", help = "Minimum offset jump (bp) reported as an indel")]
        indel_min: usize,
    },
}

#[derive(Subcommand)]
pub enum DbCommands {
    Build {
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
        #[arg(short, long, default_value = "BcgI,AlfI,AloI,FalI",
               help = "Enzyme (or comma-separated panel); default changed from BcgI-only")]
        enzyme: String,
        #[arg(long, help = "Comma-separated panel; use the same list as `ani`")]
        enzymes: Option<String>,
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
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
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
    Search {
        #[arg(short, long, required = true)]
        queries: PathBuf,
        #[arg(short, long, required = true)]
        database: PathBuf,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(short, long, default_value = "0", help = "Number of threads (0 = auto)")]
        threads: usize,
        #[arg(short, long, help = "Enable parallel processing")]
        parallel: bool,
        #[arg(short, long, default_value = "0.8",
               help = "Report hits with gated ANI >= this fraction (0.8 = 80%)")]
        min_ani: f64,
        #[arg(long, help = "Also report the partial estimators and chain diagnostics")]
        verbose: bool,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_SHARED,
               help = "Screen: pass pairs with at least this many shared tag-window keys")]
        screen_min_shared: usize,
        #[arg(long, default_value_t = crate::core::screen::DEFAULT_MIN_CONTAINMENT,
               help = "Screen: also pass pairs reaching this shared-key containment")]
        screen_min_containment: f64,
        #[arg(long, default_value_t = crate::core::screen::SCREEN_WINDOW,
               help = "Screen: centred tag window width (bp)")]
        screen_window: usize,
        #[arg(long, default_value = "0.0",
               help = "Second-tier gate: only refine screen survivors whose crude \
                       containment-ANI reaches this (0 = refine all survivors)")]
        refine_min_approx: f64,
    },
    Merge {
        #[arg(required = true)]
        databases: Vec<PathBuf>,
        #[arg(short, long, required = true)]
        output: PathBuf,
    },
}
