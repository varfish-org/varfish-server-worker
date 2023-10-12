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
    /// Optional seed for RNG.
    #[arg(long)]
    pub rng_seed: Option<u64>,
    /// Maximal distance to TAD to consider (unused, but required when loading database).
    #[arg(long, default_value_t = 10_000)]
    pub max_tad_distance: i32,
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

    tracing::info!("Loading query...");
    let query: schema::CaseQuery =
        serde_json::from_reader(std::fs::File::open(&args.path_query_json)?)?;
    tracing::info!(
        "... done loading query = {}",
        &serde_json::to_string(&query)?
    );

    tracing::info!("Loading worker databases...");
    let before_loading = Instant::now();
    let path_worker_db = format!("{}/worker", &args.path_db);
    let dbs = crate::strucvars::query::load_databases(
        &path_worker_db,
        args.genome_release,
        args.max_tad_distance,
    )?;
    tracing::info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

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

    let mut file = std::fs::File::create(&args.path_output)?;
    std::io::Write::write_all(&mut file, b"Hello, world!")?;

    trace_rss_now();

    tracing::info!(
        "All of `seqvars query` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    #[tracing_test::traced_test]
    #[test]
    fn smoke_test() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();
        let path_output = format!("{}/out.tsv", tmpdir.to_string_lossy());

        let args_common = Default::default();
        let args = super::Args {
            genome_release: crate::common::GenomeRelease::Grch37,
            path_db: "tests/seqvars/query/db".into(),
            path_query_json: "tests/seqvars/query/Case_1.query.json".into(),
            path_input: "tests/seqvars/query/Case_1.ingested.vcf".into(),
            path_output: path_output,
            max_results: None,
            rng_seed: Some(42),
            max_tad_distance: 10_000,
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(args.path_output.as_str())?);

        Ok(())
    }
}
