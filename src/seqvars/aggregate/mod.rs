//! Implementation of `seqvars aggregate` subcommand.

use std::sync::Arc;

use crate::common;

/// Command line arguments for `seqvars aggregate` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest sequence variant VCF", long_about = None)]
pub struct Args {
    /// The assumed genome build.
    #[clap(long)]
    pub genomebuild: crate::common::GenomeRelease,
    /// Path to the output RocksDB.
    #[clap(long)]
    pub path_out_rocksdb: String,
    /// Path to input VCF file(s).
    #[clap(long)]
    pub path_input: Vec<String>,

    /// Column family name for the count data.
    #[clap(long, default_value = "counts")]
    pub cf_counts: String,
    /// Column family name for the carrier UUID data.
    #[clap(long, default_value = "carriers")]
    pub cf_carriers: String,

    /// Optional path to RocksDB WAL directory.
    #[arg(long)]
    pub path_wal_dir: Option<String>,
}

/// Perform import of VCF files.
fn vcf_import(
    db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
    path_input: &[&str],
    cf_counts: &str,
    cf_carriers: &str,
) -> Result<(), anyhow::Error> {
    Ok(())
}

/// Main entry point for `seqvars aggregate` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    common::trace_rss_now();

    // Build path of all input files to read, read through files given by `@path`.
    let path_input = args
        .path_input
        .iter()
        .map(|path| {
            if path.starts_with("@") {
                std::fs::read_to_string(path.trim_start_matches('@'))
                    .expect("checked above")
                    .lines()
                    .map(|line| line.trim())
                    .filter(|line| !line.is_empty())
                    .map(|line| line.to_string())
                    .collect::<Vec<_>>()
            } else {
                vec![path.clone()]
            }
        })
        .flatten()
        .collect::<Vec<_>>();

    tracing::info!("Opening RocksDB...");
    let options = rocksdb_utils_lookup::tune_options(
        rocksdb::Options::default(),
        args.path_wal_dir.as_ref().map(|s| s.as_ref()),
    );
    let cf_names = &["meta", &args.cf_counts, &args.cf_carriers];
    let db = Arc::new(rocksdb::DB::open_cf_with_opts(
        &options,
        &args.path_out_rocksdb,
        cf_names
            .iter()
            .map(|name| (name.to_string(), options.clone()))
            .collect::<Vec<_>>(),
    )?);
    tracing::info!("  writing meta information");
    let cf_meta = db.cf_handle("meta").unwrap();
    db.put_cf(&cf_meta, "varfish-worker-version", common::worker_version())?;
    db.put_cf(&cf_meta, "db-name", "seqvars-aggregation")?;
    tracing::info!("... done opening RocksDB");

    tracing::info!("Importing VCF files ...");
    let before_import = std::time::Instant::now();
    let paths = path_input.iter().map(|s| s.as_ref()).collect::<Vec<_>>();
    vcf_import(&db, &paths, &args.cf_counts, &args.cf_carriers)?;
    tracing::info!(
        "... done importing VCF files in {:?}",
        before_import.elapsed()
    );

    tracing::info!("Running RocksDB compaction ...");
    let before_compaction = std::time::Instant::now();
    rocksdb_utils_lookup::force_compaction_cf(&db, cf_names, Some("  "), true)?;
    tracing::info!(
        "... done compacting RocksDB in {:?}",
        before_compaction.elapsed()
    );

    tracing::info!(
        "All of `seqvars aggregate` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {}
