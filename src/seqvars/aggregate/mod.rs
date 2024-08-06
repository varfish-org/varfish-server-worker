//! Implementation of `seqvars aggregate` subcommand.

pub mod ds;

use mehari::common::noodles::open_vcf_reader;
use noodles::vcf;
use rayon::prelude::*;
use std::{str::FromStr as _, sync::Arc};

use crate::common::{self, genotype_to_string, Chrom, Genotype};

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

/// Returns whether the given coordinate is in PAR for `chrom`, `pos` (1-based) and `genombuild`.
fn is_par(chrom: Chrom, pos: usize, genomebuild: crate::common::GenomeRelease) -> bool {
    match (chrom, genomebuild) {
        (Chrom::X, crate::common::GenomeRelease::Grch37) => {
            (60001..=2699520).contains(&pos) || (154931044..=155260560).contains(&pos)
        }
        (Chrom::X, crate::common::GenomeRelease::Grch38) => {
            (10001..=2781479).contains(&pos) || (155701383..=156030895).contains(&pos)
        }
        (Chrom::Y, crate::common::GenomeRelease::Grch37) => {
            (10001..=2649520).contains(&pos) || (59034050..=59363566).contains(&pos)
        }
        (Chrom::Y, crate::common::GenomeRelease::Grch38) => {
            (10001..=2781479).contains(&pos) || (56887903..=57217415).contains(&pos)
        }
        _ => false,
    }
}

/// Extract counts and carrier data from a single VCF record.
fn handle_record(
    input_record: &vcf::variant::RecordBuf,
    input_header: &vcf::Header,
    pedigree: &mehari::ped::PedigreeByName,
    case_uuid: &uuid::Uuid,
    genomebuild: crate::common::GenomeRelease,
) -> Result<(ds::Counts, ds::CarrierList), anyhow::Error> {
    let chrom: Chrom = annonars::common::cli::canonicalize(
        input_record.reference_sequence_name().to_string().as_str(),
    )
    .as_str()
    .parse()?;

    let mut res_counts = ds::Counts::default();
    let mut res_carriers = ds::CarrierList::default();

    let start: usize = input_record
        .variant_start()
        .ok_or_else(|| anyhow::anyhow!("missing variant start in record {:?}", &input_record))?
        .into();

    for (name, sample) in input_header
        .sample_names()
        .iter()
        .zip(input_record.samples().values())
    {
        let individual = pedigree
            .individuals
            .get(name)
            .ok_or_else(|| anyhow::anyhow!("individual {} not found in pedigree", name))?;

        use noodles::vcf::variant::record::samples::keys::key;
        let genotype = if let Some(Some(
            vcf::variant::record_buf::samples::sample::value::Value::Genotype(gt),
        )) = sample.get(key::GENOTYPE)
        {
            Genotype::from_str(&genotype_to_string(&gt)?)?
        } else {
            anyhow::bail!("invalid genotype value in {:?}", &sample)
        };

        // Ac-hoc enum for readable PAR status.
        #[derive(Debug, PartialEq, Eq)]
        enum _IsPar {
            IsPar,
            NoPar,
        }
        use _IsPar::*;
        let is_par = if is_par(chrom, start, genomebuild) {
            IsPar
        } else {
            NoPar
        };

        let carrier_genotype = match (chrom, is_par, individual.sex, genotype) {
            (_, _, _, Genotype::WithNoCall) => continue,
            // On the autosomes, male/female are handled the same.
            (Chrom::Auto, _, _, Genotype::HomRef) => {
                res_counts.count_homref += 1;
                ds::Genotype::HomRef
            }
            (Chrom::Auto, _, _, Genotype::Het) => {
                res_counts.count_het += 1;
                ds::Genotype::Het
            }
            (Chrom::Auto, _, _, Genotype::HomAlt) => {
                res_counts.count_homalt += 2;
                ds::Genotype::HomAlt
            }
            // On the gonomosomes, we handle call male variant calls as hemizygous outside PAR.
            (Chrom::X, NoPar, mehari::ped::Sex::Male, Genotype::HomRef)
            | (Chrom::Y, NoPar, mehari::ped::Sex::Male, Genotype::HomRef) => {
                res_counts.count_hemiref += 1;
                ds::Genotype::HemiRef
            }
            (Chrom::X, NoPar, mehari::ped::Sex::Male, Genotype::Het)
            | (Chrom::Y, NoPar, mehari::ped::Sex::Male, Genotype::Het)
            | (Chrom::X, NoPar, mehari::ped::Sex::Male, Genotype::HomAlt)
            | (Chrom::Y, NoPar, mehari::ped::Sex::Male, Genotype::HomAlt) => {
                res_counts.count_hemialt += 1;
                ds::Genotype::HemiAlt
            }
            // For female samples, we handle chrX as biallelic, same as male inside PAR on
            // chrX/chrY.  Note that read mapping pipelines N-mask the PAR on chrY, so we
            // don't count twice.
            (Chrom::X, IsPar, mehari::ped::Sex::Male, Genotype::HomRef)
            | (Chrom::Y, IsPar, mehari::ped::Sex::Male, Genotype::HomRef)
            | (Chrom::X, _, mehari::ped::Sex::Female, Genotype::HomRef) => {
                res_counts.count_homref += 1;
                ds::Genotype::HomRef
            }
            (Chrom::X, IsPar, mehari::ped::Sex::Male, Genotype::Het)
            | (Chrom::Y, IsPar, mehari::ped::Sex::Male, Genotype::Het)
            | (Chrom::X, _, mehari::ped::Sex::Female, Genotype::Het) => {
                res_counts.count_het += 1;
                ds::Genotype::Het
            }
            (Chrom::X, IsPar, mehari::ped::Sex::Male, Genotype::HomAlt)
            | (Chrom::Y, IsPar, mehari::ped::Sex::Male, Genotype::HomAlt)
            | (Chrom::X, _, mehari::ped::Sex::Female, Genotype::HomAlt) => {
                res_counts.count_homalt += 2;
                ds::Genotype::HomAlt
            }
            // We ignore calls to chrY for female samples.
            (Chrom::Y, _, mehari::ped::Sex::Female, Genotype::HomRef)
            | (Chrom::Y, _, mehari::ped::Sex::Female, Genotype::Het)
            | (Chrom::Y, _, mehari::ped::Sex::Female, Genotype::HomAlt) => ds::Genotype::HomRef,
            // Do not count samples with unknown sex on gonomosomes.
            (Chrom::X, _, mehari::ped::Sex::Unknown, _)
            | (Chrom::Y, _, mehari::ped::Sex::Unknown, _) => ds::Genotype::HomRef,
        };

        if carrier_genotype != ds::Genotype::HomRef {
            res_carriers.carriers.push(ds::Carrier {
                uuid: *case_uuid,
                index: pedigree
                    .individuals
                    .get_index_of(name.as_str())
                    .ok_or_else(|| anyhow::anyhow!("individual {} not found in pedigree", &name))?
                    as u8,
                genotype: carrier_genotype,
            });
        }
    }

    Ok((res_counts, res_carriers))
}

/// Import one VCF file into the database.
///
/// This function is `async` because we potentially need to read from S3.
async fn import_vcf(
    db: &Arc<rocksdb::TransactionDB<rocksdb::MultiThreaded>>,
    path_input: &str,
    cf_counts: &str,
    cf_carriers: &str,
    genomebuild: crate::common::GenomeRelease,
) -> Result<(), anyhow::Error> {
    let mut input_reader = open_vcf_reader(path_input)
        .await
        .map_err(|e| anyhow::anyhow!("could not open file {} for reading: {}", path_input, e))?;
    let input_header = input_reader.read_header().await?;

    let cf_counts = db.cf_handle(cf_counts).expect("checked earlier");
    let cf_carriers = db.cf_handle(cf_carriers).expect("checked earlier");

    let (pedigree, case_uuid) = common::extract_pedigree_and_case_uuid(&input_header)?;
    let mut prev = std::time::Instant::now();
    let mut record_buf = vcf::variant::RecordBuf::default();
    loop {
        let bytes_read = input_reader
            .read_record_buf(&input_header, &mut record_buf)
            .await
            .map_err(|e| anyhow::anyhow!("problem reading VCF file {}: {}", path_input, e))?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Obtain counts from the current variant.
        let (this_counts_data, this_carrier_data) = handle_record(
            &record_buf,
            &input_header,
            &pedigree,
            &case_uuid,
            genomebuild,
        )?;
        // Obtain annonars variant key from current allele for RocksDB lookup.
        let vcf_var = annonars::common::keys::Var::from_vcf_allele(&record_buf, 0);
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
                })?
                .map(|buffer| ds::CarrierList::try_from(buffer.as_slice()))
                .transpose()
                .map_err(|e| {
                    anyhow::anyhow!(
                        "problem decoding carrier data for variant {:?}: {}",
                        &vcf_var,
                        e
                    )
                })?
                .unwrap_or_default();

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
    }

    Ok(())
}

/// Perform the parallel import of VCF files.
fn vcf_import(
    db: &Arc<rocksdb::TransactionDB<rocksdb::MultiThreaded>>,
    path_input: &[&str],
    cf_counts: &str,
    cf_carriers: &str,
    genomebuild: crate::common::GenomeRelease,
) -> Result<(), anyhow::Error> {
    path_input
        .par_iter()
        .map(|path_input| {
            // We create a Tokio scheduler for the current file as we need it
            // to wait / block for the VCF import running in the current Rayon
            // thread.
            tokio::runtime::Builder::new_current_thread()
                .build()
                .map_err(|e| {
                    anyhow::anyhow!(
                        "building Tokio runtime for VCF file {} failed: {}",
                        path_input,
                        e
                    )
                })?
                .block_on(import_vcf(
                    db,
                    path_input,
                    cf_counts,
                    cf_carriers,
                    genomebuild,
                ))
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
        .flat_map(|path| {
            if path.starts_with('@') {
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
        vcf_import(
            &db,
            &paths,
            &args.cf_counts,
            &args.cf_carriers,
            args.genomebuild,
        )?;
        tracing::info!(
            "... done importing VCF files in {:?}",
            before_import.elapsed()
        );
    }

    // scope for compaction
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
mod test {
    use super::*;

    #[test]
    fn test_is_par() {
        assert_eq!(
            super::is_par(super::Chrom::X, 60000, crate::common::GenomeRelease::Grch37),
            false
        );
        assert_eq!(
            super::is_par(super::Chrom::X, 60001, crate::common::GenomeRelease::Grch37),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                2699520,
                crate::common::GenomeRelease::Grch37
            ),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                2699521,
                crate::common::GenomeRelease::Grch37
            ),
            false
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                154931043,
                crate::common::GenomeRelease::Grch37
            ),
            false
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                154931044,
                crate::common::GenomeRelease::Grch37
            ),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                155260560,
                crate::common::GenomeRelease::Grch37
            ),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                155260561,
                crate::common::GenomeRelease::Grch37
            ),
            false
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                155260561,
                crate::common::GenomeRelease::Grch38
            ),
            false
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                155701383,
                crate::common::GenomeRelease::Grch38
            ),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                156030895,
                crate::common::GenomeRelease::Grch38
            ),
            true
        );
        assert_eq!(
            super::is_par(
                super::Chrom::X,
                156030896,
                crate::common::GenomeRelease::Grch38
            ),
            false
        );
        assert_eq!(
            super::is_par(super::Chrom::Y, 10000, crate::common::GenomeRelease::Grch37),
            false
        );
        assert_eq!(
            super::is_par(super::Chrom::Y, 10001, crate::common::GenomeRelease::Grch37),
            true
        );
    }

    #[tracing_test::traced_test]
    #[test]
    fn handle_record_snapshot() -> Result<(), anyhow::Error> {
        let path = "tests/seqvars/aggregate/ingest.vcf";
        let mut vcf_reader = vcf::io::reader::Builder::default()
            .build_from_path(path)
            .unwrap();
        let header = vcf_reader.read_header().unwrap();

        let mut record_buf = vcf::variant::RecordBuf::default();
        loop {
            let bytes_read = vcf_reader
                .read_record_buf(&header, &mut record_buf)
                .map_err(|e| anyhow::anyhow!("problem reading VCF file {}: {}", path, e))?;
            if bytes_read == 0 {
                break; // EOF
            }

            let (pedigree, case_uuid) = common::extract_pedigree_and_case_uuid(&header)?;
            let (counts, carriers) = super::handle_record(
                &record_buf,
                &header,
                &pedigree,
                &case_uuid,
                crate::common::GenomeRelease::Grch37,
            )?;

            insta::assert_debug_snapshot!(counts);
            insta::assert_debug_snapshot!(carriers);
        }

        Ok(())
    }
}
