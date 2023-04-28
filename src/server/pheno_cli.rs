//! Code for ranking genes on the command line.

use indicatif::{ParallelProgressIterator, ProgressState};
use prost::Message;
use rayon::prelude::*;
use rocksdb::{DBWithThreadMode, MultiThreaded};
use serde::Deserialize;
use std::fmt::Write;
use std::time::Instant;

use clap::Parser;
use hpo::{annotations::AnnotationId, term::HpoGroup, HpoTermId, Ontology};
use tracing::info;

use crate::server::data::SimulationResults;
use crate::{common::trace_rss_now, server::pheno::phenomizer};

/// The version of `varfish-server-worker` package.
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Command line arguments for `server prepare-pheno` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Prepare values for `server pheno`", long_about = None)]
pub struct Args {
    /// Path to the directory with the HPO files.
    #[arg(long, required = true)]
    pub path_hpo_dir: String,

    /// Path to JSON file with the genes to rank.
    #[arg(long)]
    pub path_genes_json: String,
    /// Path to JSON file with HPO IDs of patient.
    #[arg(long)]
    pub path_terms_json: String,
}

/// Struct for loading a gene from JSON.
#[derive(Deserialize, Debug, Clone)]
pub struct Gene {
    /// The gene symbol.
    pub gene_symbol: String,
}

/// Struct for loading an HPO term from JSON.
#[derive(Deserialize, Debug, Clone)]
pub struct HpoTerm {
    /// The term ID.
    pub term_id: String,
}

/// Main entry point for `server pheno-cli` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    if let Some(level) = args_common.verbose.log_level() {
        match level {
            log::Level::Trace | log::Level::Debug => {
                std::env::set_var("RUST_LOG", "debug");
                env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
            }
            _ => (),
        }
    }

    info!("Loading HPO...");
    let before_loading = Instant::now();
    let hpo = Ontology::from_standard(&args.path_hpo_dir)?;
    info!("...done loading HPO in {:?}", before_loading.elapsed());

    info!("Opening RocksDB for reading...");
    let before_rocksdb = Instant::now();
    let path_rocksdb = format!("{}/resnik", args.path_hpo_dir);
    let db = rocksdb::DB::open_for_read_only(&rocksdb::Options::default(), &path_rocksdb, false)?;
    let cf_resnik = db.cf_handle("resnik_pvalues").unwrap();
    info!("...done opening RocksDB in {:?}", before_rocksdb.elapsed());

    info!("Loading genes...");
    let before_load_genes = Instant::now();
    let genes_json = std::fs::read_to_string(&args.path_genes_json)?;
    let genes: Vec<Gene> = serde_json::from_str(&genes_json)?;
    let genes = genes.iter().map(|g| {
        hpo.gene_by_name(&g.gene_symbol)
            .expect(&format!("gene {} not found", &g.gene_symbol))
    });
    info!("... done loadin genes in {:?}", before_load_genes.elapsed());

    info!("Loading genes...");
    let before_load_genes = Instant::now();
    let genes_json = std::fs::read_to_string(&args.path_genes_json)?;
    let genes: Vec<Gene> = serde_json::from_str(&genes_json)?;
    let genes = genes
        .iter()
        .map(|g| {
            hpo.gene_by_name(&g.gene_symbol)
                .expect(&format!("gene {} not found", &g.gene_symbol))
        })
        .collect::<Vec<_>>();
    info!("... done loadin genes in {:?}", before_load_genes.elapsed());

    info!("Loading (patient/query) HPO term ids...");
    let before_load_genes = Instant::now();
    let query_json = std::fs::read_to_string(&args.path_terms_json)?;
    let query: Vec<HpoTerm> = serde_json::from_str(&query_json)?;
    let query = query
        .iter()
        .map(|t| {
            HpoTermId::try_from(t.term_id.as_str())
                .expect(&format!("term {} no valid HPO term ID", &t.term_id))
        })
        .collect::<Vec<_>>();
    let query = {
        let mut group = HpoGroup::new();
        for term in query {
            group.insert(term);
        }
        group
    };
    info!(
        "... done loading HPO IDs in {:?}",
        before_load_genes.elapsed()
    );

    trace_rss_now();

    info!("Starting priorization...");
    let before_priorization = Instant::now();
    for gene in genes {
        let score = phenomizer::score(&query, gene.hpo_terms(), &hpo);
        info!("  {:?} -> {}", gene, score);
    }
    info!(
        "... done with prioritization in {:?}",
        before_priorization.elapsed()
    );

    info!("All done. Have a nice day!");
    Ok(())
}
