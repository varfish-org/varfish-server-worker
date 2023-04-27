//! Code for preparing the empiric P-value distributions for the phenotypes.

use indicatif::ParallelProgressIterator;
use prost::Message;
use rayon::prelude::*;
use rocksdb::{DBWithThreadMode, MultiThreaded};
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
    /// Path to output RocksDB.
    #[arg(long, required = true)]
    pub path_out_rocksdb: String,

    /// Number of simulations to perform for each gene and term set size.
    #[arg(long, default_value_t = 100_000, value_parser = clap::value_parser!(u64).range(2..))]
    pub num_simulations: u64,
    /// Run simulations for `min_terms..=max_terms` terms.
    #[arg(long, default_value_t = 1)]
    pub min_terms: usize,
    /// Run simulations for `min_terms..=max_terms` terms.
    #[arg(long, default_value_t = 10)]
    pub max_terms: usize,
    /// Will sample eCDF to this number of points.
    #[arg(long, default_value_t = 10, value_parser = clap::value_parser!(u64).range(100..))]
    pub ecdf_samples: u64,

    /// Number of threads to use for simulation (default is 1 thread per core).
    #[arg(long)]
    pub num_threads: Option<usize>,
    /// Seed for the random number generator.
    #[arg(long)]
    pub seed: Option<u64>,
}

/// Run simulation using ontology and number of terms.
fn run_simulation(
    db: &DBWithThreadMode<MultiThreaded>,
    ontology: &Ontology,
    args: &Args,
    num_terms: usize,
) -> Result<(), anyhow::Error> {
    info!("  running simulation for {} terms ...", num_terms);
    let before = Instant::now();

    // We want at least two simulations.
    let num_simulations = std::cmp::max(args.num_simulations, 2);

    // Get all HPO terms for phenotypic abnormalities.
    let hpo_abnormality = ontology
        .hpo(HpoTermId::from(String::from("HP:0000118")))
        .ok_or(anyhow::anyhow!(
            "could not find HP:0000118 (phenotypic abnormality)"
        ))?;
    let term_ids = ontology
        .hpos()
        .filter(|t| t.child_of(&hpo_abnormality))
        .map(|t| t.id())
        .collect::<Vec<_>>();
    // Get all genes into a vector so we can use parallel iteration.
    let genes = ontology.genes().collect::<Vec<_>>();

    // Run simulations for each gene in parallel.
    let style = indicatif::ProgressStyle::default_bar();
    genes
        .par_iter()
        .progress_with_style(style)
        .for_each(|gene| {
            // Obtain sorted list of similarity scores from simulations.
            let mut scores = (0..num_simulations)
                .into_par_iter()
                .map(|_| {
                    // Pick `num_terms` random terms with circuit breakers on number of tries.
                    let max_tries = 1000;
                    let mut tries = 0;
                    let mut ts = HpoGroup::new();
                    // let mut ts = Vec::new();
                    while ts.len() < num_terms {
                        tries += 1;
                        if tries > max_tries {
                            panic!("tried too often to pick random terms");
                        }
                        let term_id = term_ids[fastrand::usize(0..term_ids.len())];
                        if !ts.contains(&term_id) {
                            ts.insert(term_id);
                        }
                    }

                    // Compute similarity.
                    phenomizer::score(&ts, &gene.hpo_terms(), &ontology)
                })
                .collect::<prost::alloc::vec::Vec<_>>();
            scores.sort_by(|a, b| a.partial_cmp(b).expect("NaN value"));

            // Copy the scores into the score distribution.
            let ncbi_gene_id = gene.id().as_u32();
            let sim_res = SimulationResults {
                ncbi_gene_id,
                gene_symbol: gene.name().to_string(),
                term_count: num_terms as u32,
                scores,
            };

            // Encode as byte array.
            let sim_res = sim_res.encode_to_vec();

            // Write to RocksDB.
            let cf_resnik = db.cf_handle("resnik_pvalues").unwrap();
            let key = format!("{ncbi_gene_id}:{num_terms}");
            db.put_cf(&cf_resnik, &key.as_bytes(), &sim_res)
                .expect("writing to RocksDB failed");
        });
    info!("  ... done in {:?}", before.elapsed());

    Ok(())
}

/// Main entry point for `server prepare-pheno` sub command.
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
    let ontology = Ontology::from_standard(&args.path_hpo_dir)?;
    info!("...done loading HPO in {:?}", before_loading.elapsed());

    info!("Opening RocksDB for writing...");
    let before_rocksdb = Instant::now();
    let options = {
        let mut options = rocksdb::Options::default();
        options.prepare_for_bulk_load();
        options
    };
    let db = rocksdb::DB::open_cf(&options, &args.path_out_rocksdb, ["meta", "resnik_pvalues"])?;
    // write out metadata
    let cf_meta = db.cf_handle("meta").unwrap();
    db.put_cf(&cf_meta, "hpo-version", ontology.hpo_version())?;
    db.put_cf(&cf_meta, "app-version", VERSION.to_string())?;
    let cf_resnik = db.cf_handle("resnik_pvalues").unwrap();
    info!("...done opening RocksDB in {:?}", before_rocksdb.elapsed());

    info!("Running simulations...");
    let before_simulations = Instant::now();
    if let Some(seed) = args.seed {
        fastrand::seed(seed);
    }
    if let Some(num_threds) = args.num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threds)
            .build_global()?;
    }
    for num_terms in args.min_terms..=args.max_terms {
        run_simulation(&db, &ontology, &args, num_terms)?;
    }
    info!(
        "... done with simulations in {:?}",
        before_simulations.elapsed()
    );

    trace_rss_now();

    tracing::info!("Enforcing manual compaction");
    db.compact_range_cf(&cf_meta, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(&cf_resnik, None::<&[u8]>, None::<&[u8]>);

    let compaction_start = Instant::now();
    let mut last_printed = compaction_start;
    while db
        .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
        .unwrap()
        > 0
        || db
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap()
            > 0
    {
        std::thread::sleep(std::time::Duration::from_millis(100));
        if last_printed.elapsed() > std::time::Duration::from_millis(1000) {
            log::info!(
                "... waiting for compaction for {:?}",
                compaction_start.elapsed()
            );
            last_printed = Instant::now();
        }
    }

    info!("All done. Have a nice day!");
    Ok(())
}
