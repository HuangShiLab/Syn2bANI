use clap::Parser;
use anyhow::Result;
use env_logger;
use log::info;

use syn2bani::cli::{Cli, Commands, DbCommands};
use syn2bani::cli::ani::run_ani;
use syn2bani::cli::panel::run_panel;
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
        Commands::Ani { query, ql, rl, enzymes, mismatch_tolerance, min_chain_anchors, max_gap, threads, parallel, verbose, strata_out, calibrate, calibrate_model, output } => {
            info!("Running chain-restricted MLE ani with enzymes: {}", enzymes);
            run_ani(&query, ql.as_deref(), rl.as_deref(), &enzymes, mismatch_tolerance, min_chain_anchors, max_gap, threads, parallel, verbose, strata_out.as_deref(), calibrate, calibrate_model.as_deref(), output.as_deref())?;
        }
        Commands::Panel { strata, truth, greedy, panels } => {
            run_panel(&strata, &truth, greedy, panels.as_deref())?;
        }
        Commands::Dist { query, ql, rl, enzymes, threads, parallel, verbose, min_ani, screen_min_shared, screen_min_containment, screen_window, refine_min_approx, output } => {
            info!("Running dist with enzymes: {}", enzymes);
            let screen = syn2bani::core::screen::ScreenConfig { min_shared: screen_min_shared, min_containment: screen_min_containment, window: screen_window };
            run_dist(&query, ql.as_deref(), rl.as_deref(), &enzymes, threads, parallel, verbose, min_ani, screen, refine_min_approx, output.as_deref())?;
        }
        Commands::Search { query, database, ql, rl, output, threads, parallel, min_ani, enzymes, verbose, screen_min_shared, screen_min_containment, screen_window, refine_min_approx } => {
            info!("Running search");
            let screen = syn2bani::core::screen::ScreenConfig { min_shared: screen_min_shared, min_containment: screen_min_containment, window: screen_window };
            run_search(&query, ql.as_deref(), database.as_deref(), rl.as_deref(), output.as_deref(), threads, parallel, min_ani, &enzymes, screen, refine_min_approx, verbose)?;
        }
        Commands::Sketch { genomes, output, enzyme, enzymes, threads, parallel, multi_enzyme } => {
            info!("Running sketch with enzyme: {}", enzymes.as_deref().unwrap_or(&enzyme));
            run_sketch(&genomes, &output, &enzyme, threads, parallel, multi_enzyme, enzymes.as_deref())?;
        }
        Commands::Triangle { genomes, ql, output, edge_list, threads, parallel, enzymes, verbose, screen_min_shared, screen_min_containment, screen_window, refine_min_approx } => {
            info!("Running triangle comparison on {} genomes", genomes.len());
            let screen = syn2bani::core::screen::ScreenConfig { min_shared: screen_min_shared, min_containment: screen_min_containment, window: screen_window };
            run_triangle(&genomes, ql.as_deref(), output.as_deref(), edge_list, threads, parallel, &enzymes, screen, refine_min_approx, verbose)?;
        }
        Commands::Db { command } => {
            match command {
                DbCommands::Build { genomes, output, enzyme, enzymes, threads, parallel, multi_enzyme } => {
                    db::run_db_build(&genomes, &output, &enzyme, threads, parallel, multi_enzyme, enzymes.as_deref())?;
                }
                DbCommands::Add { genomes, database, threads, parallel } => {
                    db::run_db_add(&genomes, &database, threads, parallel)?;
                }
                DbCommands::Remove { genome_ids, database } => {
                    db::run_db_remove(&genome_ids, &database)?;
                }
                DbCommands::List { database } => {
                    db::run_db_list(&database)?;
                }
                DbCommands::Search { queries, database, output, threads, parallel, min_ani, verbose, screen_min_shared, screen_min_containment, screen_window, refine_min_approx } => {
                    let screen = syn2bani::core::screen::ScreenConfig { min_shared: screen_min_shared, min_containment: screen_min_containment, window: screen_window };
                    db::run_db_search(&queries, &database, output.as_deref(), threads, parallel, min_ani, screen, refine_min_approx, verbose)?;
                }
                DbCommands::Merge { databases, output } => {
                    db::run_db_merge(&databases, &output)?;
                }
            }
        }
        Commands::Struct { query, reference, output, paf, rearrangement, indel, multi_enzyme, enzymes, indel_min } => {
            run_struct(&query, &reference, output.as_deref(), paf, rearrangement, indel, multi_enzyme, enzymes.as_deref(), indel_min)?;
        }
    }

    Ok(())
}
