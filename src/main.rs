use clap::Parser;
use anyhow::Result;
use env_logger;
use log::info;

use syn2bani::cli::{Cli, Commands, DbCommands};
use syn2bani::cli::dist::run_dist;
use syn2bani::cli::sketch::run_sketch;
use syn2bani::cli::search::run_search;
use syn2bani::cli::triangle::run_triangle;
use syn2bani::cli::db;
use syn2bani::cli::r#struct::run_struct;

fn main() -> Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Dist { query, reference, enzyme, threads, parallel, multi_enzyme, structural, output } => {
            info!("Running dist with enzyme: {}", enzyme);
            run_dist(&query, &reference, &enzyme, threads, parallel, multi_enzyme, structural, output.as_deref())?;
        }
        Commands::Search { query, database, output, threads, parallel, min_ani } => {
            info!("Running search against database: {}", database.display());
            run_search(&query, &database, output.as_deref(), threads, parallel, min_ani)?;
        }
        Commands::Sketch { genomes, output, enzyme, threads, parallel, multi_enzyme } => {
            info!("Running sketch with enzyme: {}", enzyme);
            run_sketch(&genomes, &output, &enzyme, threads, parallel, multi_enzyme)?;
        }
        Commands::Triangle { genomes, output, edge_list, threads, parallel } => {
            info!("Running triangle comparison on {} genomes", genomes.len());
            run_triangle(&genomes, output.as_deref(), edge_list, threads, parallel)?;
        }
        Commands::Db { command } => {
            match command {
                DbCommands::Build { genomes, output, enzyme, threads, parallel, multi_enzyme } => {
                    db::run_db_build(&genomes, &output, &enzyme, threads, parallel, multi_enzyme)?;
                }
                DbCommands::Add { genomes, database } => {
                    db::run_db_add(&genomes, &database)?;
                }
                DbCommands::Remove { genome_ids, database } => {
                    db::run_db_remove(&genome_ids, &database)?;
                }
                DbCommands::List { database } => {
                    db::run_db_list(&database)?;
                }
                DbCommands::Merge { databases, output } => {
                    db::run_db_merge(&databases, &output)?;
                }
            }
        }
        Commands::Struct { query, reference, output, paf, rearrangement, indel } => {
            run_struct(&query, &reference, output.as_deref(), paf, rearrangement, indel)?;
        }
    }

    Ok(())
}
