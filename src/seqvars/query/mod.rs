//! Code implementing the "seqvars query" sub command.

mod schema;

use std::time::Instant;

use clap::{command, Parser};

use rand_core::SeedableRng;

use crate::{common::trace_rss_now, common::GenomeRelease};

/// Command line arguments for `seqvars query` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for seqvars", long_about = None)]
pub struct Args {
    /// Genome release to assume.
    #[arg(long, value_enum)]
    pub genome_release: GenomeRelease,
    /// Path to worker database to use for querying.
    #[arg(long, required = true)]
    pub path_db: String,
    /// Path to query JSON file.
    #[arg(long, required = true)]
    pub path_query_json: String,
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub path_input: String,
    /// Path to the output TSV file.
    #[arg(long, required = true)]
    pub path_output: String,

    /// Optional maximal number of total records to write out.
    #[arg(long)]
    pub max_results: Option<usize>,
    /// Radius around BND sites used when building the database.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: i32,
    /// Radius around INS sites used when building the database.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: i32,
    /// Minimal reciprocal overlap for SVs of the same type, used when building
    /// the database.
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Maximal distance to TAD to consider.
    #[arg(long, default_value_t = 10_000)]
    pub max_tad_distance: i32,
    /// Optional seed for RNG.
    #[arg(long)]
    pub rng_seed: Option<u64>,
}

/// Main entry point for `seqvars query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = Instant::now();
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    // Initialize the random number generator from command line seed if given or local entropy
    // source.
    let _rng = if let Some(rng_seed) = args.rng_seed {
        rand::rngs::StdRng::seed_from_u64(rng_seed)
    } else {
        rand::rngs::StdRng::from_entropy()
    };

    // tracing::info!("Loading query...");
    // let query: CaseQuery = serde_json::from_reader(File::open(&args.path_query_json)?)?;
    // tracing::info!(
    //     "... done loading query = {}",
    //     &serde_json::to_string(&query)?
    // );

    // tracing::info!("Loading worker databases...");
    // let before_loading = Instant::now();
    // let path_worker_db = format!("{}/worker", &args.path_db);
    // let dbs = load_databases(&path_worker_db, args.genome_release, args.max_tad_distance)?;
    // tracing::info!(
    //     "...done loading databases in {:?}",
    //     before_loading.elapsed()
    // );

    trace_rss_now();

    // tracing::info!("Loading mehari tx database...");
    // let before_loading = Instant::now();
    // let path_mehari_tx_db = format!(
    //     "{}/mehari/{}/txs.bin.zst",
    //     &args.path_db,
    //     &args.genome_release.to_string()
    // );
    // tracing::debug!("  path = {}", &path_mehari_tx_db);
    // let mehari_tx_db = mehari::annotate::seqvars::load_tx_db(&path_mehari_tx_db)?;
    // tracing::info!(
    //     "...done loading mehari tx database in {:?}",
    //     before_loading.elapsed()
    // );
    // tracing::info!("Building mehari index data structures...");
    // let before_building = Instant::now();
    // let mehari_tx_idx = TxIntervalTrees::new(&mehari_tx_db, args.genome_release.into());
    // let chrom_to_acc = ASSEMBLY_INFOS[args.genome_release.into()]
    //     .sequences
    //     .iter()
    //     .map(|record| {
    //         (
    //             annonars::common::cli::canonicalize(&record.name),
    //             record.refseq_ac.clone(),
    //         )
    //     })
    //     .collect::<HashMap<_, _>>();
    // tracing::info!(
    //     "...done building mehari index data structures in {:?}",
    //     before_building.elapsed()
    // );

    // trace_rss_now();

    // tracing::info!("Translating gene allow list...");
    // let hgvs_allowlist = if let Some(gene_allowlist) = &query.gene_allowlist {
    //     if gene_allowlist.is_empty() {
    //         None
    //     } else {
    //         Some(translate_gene_allowlist(gene_allowlist, &dbs))
    //     }
    // } else {
    //     None
    // };

    // tracing::info!("Running queries...");
    // let before_query = Instant::now();
    // let query_stats = run_query(
    //     &QueryInterpreter::new(query, hgvs_allowlist),
    //     args,
    //     &dbs,
    //     &mehari_tx_db,
    //     &mehari_tx_idx,
    //     &chrom_to_acc,
    //     &mut rng,
    // )?;
    // tracing::info!("... done running query in {:?}", before_query.elapsed());
    // tracing::info!(
    //     "summary: {} records passed out of {}",
    //     query_stats.count_passed.separate_with_commas(),
    //     query_stats.count_total.separate_with_commas()
    // );
    // tracing::info!("passing records by SV type");
    // for (sv_type, count) in query_stats.by_sv_type.iter() {
    //     tracing::info!("{:?} -- {}", sv_type, count);
    // }

    trace_rss_now();

    tracing::info!(
        "All of `seqvars query` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}
