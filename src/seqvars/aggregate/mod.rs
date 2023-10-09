//! Implementation of `seqvars aggregate` subcommand.

pub mod ds;

use mehari::common::open_read_maybe_gz;
use noodles_vcf as vcf;
use rayon::prelude::*;
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
    /// Set the number of threads to use, defaults to number of cores.
    #[clap(long)]
    pub num_threads: Option<usize>,

    /// Optional path to RocksDB WAL directory.
    #[arg(long)]
    pub path_wal_dir: Option<String>,
}

/// Extract a PedigreeByName from the VCF header.
fn extract_pedigree(header: &vcf::Header) -> Result<mehari::ped::PedigreeByName, anyhow::Error> {
    todo!()
}

/// Extract counts and carrier data from a single VCF record.
fn handle_record(
    input_record: &vcf::Record,
    pedigree: &mehari::ped::PedigreeByName,
) -> Result<(ds::Counts, ds::CarrierList), anyhow::Error> {
    todo!()
}

/// Import one VCF file into the database.
fn import_vcf(
    db: &Arc<rocksdb::TransactionDB<rocksdb::MultiThreaded>>,
    path_input: &str,
    cf_counts: &str,
    cf_carriers: &str,
) -> Result<(), anyhow::Error> {
    let mut input_reader = open_read_maybe_gz(path_input)
        .map_err(|e| anyhow::anyhow!("could not open file {} for reading: {}", path_input, e))?;
    let mut input_reader = vcf::Reader::new(&mut input_reader);
    let input_header = input_reader.read_header()?;

    let cf_counts = db.cf_handle(cf_counts).expect("checked earlier");
    let cf_carriers = db.cf_handle(cf_carriers).expect("checked earlier");

    let pedigree = extract_pedigree(&input_header)?;
    let mut prev = std::time::Instant::now();
    let mut records = input_reader.records(&input_header);
    loop {
        if let Some(input_record) = records.next() {
            let input_record = input_record?;

            // TODO: we need to lock the database for counts

            // Obtain counts from the current variant.
            let (this_counts_data, this_carrier_data) = handle_record(&input_record, &pedigree)?;
            // Obtain annonars variant key from current allele for RocksDB lookup.
            let vcf_var = annonars::common::keys::Var::from_vcf_allele(&input_record, 0);
            let key: Vec<u8> = vcf_var.clone().into();

            let max_retries = 10;
            let mut retries = 0;
            while retries < max_retries {
                let this_counts_data = this_counts_data.clone();
                let this_carrier_data = this_carrier_data.clone();

                let transaction = db.transaction();

                // Read data for variant from database.
                let mut db_counts_data = transaction.get_cf(&cf_counts, key.clone()).map_err(|e| {
                    anyhow::anyhow!(
                        "problem acessing counts data for variant {:?}: {} (non-existing would be fine)",
                        &vcf_var,
                        e
                    )
                })?.map(|buffer| ds::Counts::from_vec(&buffer)).unwrap_or_default();
                let mut db_carrier_data = transaction.get_cf(&cf_carriers, key.clone()).map_err(|e| {
                    anyhow::anyhow!(
                        "problem acessing carrier data for variant {:?}: {} (non-existing would be fine)",
                        &vcf_var,
                        e
                    )
                })?.map(|buffer| ds::CarrierList::from_vec(&buffer)).unwrap_or_default();

                // Aggregate the data.
                db_counts_data.aggregate(this_counts_data);
                db_carrier_data.aggregate(this_carrier_data);

                // Write data for variant back to database.
                transaction
                    .put_cf(&cf_counts, key.clone(), &db_counts_data.to_vec())
                    .map_err(|e| {
                        anyhow::anyhow!(
                            "problem writing counts data for variant {:?}: {}",
                            &vcf_var,
                            e
                        )
                    })?;
                transaction
                    .put_cf(&cf_carriers, key.clone(), &db_carrier_data.to_vec())
                    .map_err(|e| {
                        anyhow::anyhow!(
                            "problem writing carrier data for variant {:?}: {}",
                            &vcf_var,
                            e
                        )
                    })?;

                let res = transaction.commit();
                match res {
                    Ok(_) => break,
                    Err(e) => {
                        retries += 1;
                        if retries > 5 {
                            tracing::warn!(
                                "problem committing transaction for variant {:?}: {} (retry #{})",
                                &vcf_var,
                                e,
                                retries
                            );
                        }
                    }
                }
            }
            if retries >= max_retries {
                return Err(anyhow::anyhow!(
                    "problem committing transaction for variant {:?}: {} (max retries exceeded)",
                    &vcf_var,
                    retries
                ));
            }

            // Write out progress indicator every 60 seconds.
            if prev.elapsed().as_secs() >= 60 {
                tracing::info!("at {:?}", &vcf_var);
                prev = std::time::Instant::now();
            }
        } else {
            break; // all done
        }
    }

    todo!()
}

/// Perform the parallel import of VCF files.
fn vcf_import(
    db: &Arc<rocksdb::TransactionDB<rocksdb::MultiThreaded>>,
    path_input: &[&str],
    cf_counts: &str,
    cf_carriers: &str,
) -> Result<(), anyhow::Error> {
    path_input
        .par_iter()
        .map(|path_input| {
            import_vcf(db, path_input, cf_counts, cf_carriers)
                .map_err(|e| anyhow::anyhow!("processing VCF file {} failed: {}", path_input, e))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(|_| ())
}

/// Main entry point for `seqvars aggregate` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    if let Some(num_threads) = args.num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .map_err(|e| anyhow::anyhow!("building global Rayon thread pool failed: {}", e))?;
    }

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
    let tx_options = rocksdb::TransactionDBOptions::default();
    let cf_names = &["meta", &args.cf_counts, &args.cf_carriers];
    let cf_descriptors = cf_names
        .iter()
        .map(|name| rocksdb::ColumnFamilyDescriptor::new(*name, options.clone()))
        .collect::<Vec<_>>();

    // scope for the transaction database
    {
        let db: Arc<rocksdb::TransactionDB<rocksdb::MultiThreaded>> =
            Arc::new(rocksdb::TransactionDB::open_cf_descriptors(
                &options,
                &tx_options,
                &args.path_out_rocksdb,
                cf_descriptors,
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
    }

    // scope for cleanup
    {
        let db = Arc::new(rocksdb::DB::open_cf_with_opts(
            &options,
            &args.path_out_rocksdb,
            cf_names
                .iter()
                .map(|name| (name.to_string(), options.clone()))
                .collect::<Vec<_>>(),
        )?);
        tracing::info!("Running RocksDB compaction ...");
        let before_compaction = std::time::Instant::now();
        rocksdb_utils_lookup::force_compaction_cf(&db, cf_names, Some("  "), true)?;
        tracing::info!(
            "... done compacting RocksDB in {:?}",
            before_compaction.elapsed()
        );
    }

    tracing::info!(
        "All of `seqvars aggregate` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {}
