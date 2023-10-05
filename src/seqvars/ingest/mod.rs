//! Implementation of `seqvars ingest` subcommand.

use std::sync::Arc;

use crate::common::{self, GenomeRelease};
use mehari::annotate::seqvars::provider::MehariProvider;
use noodles_vcf as vcf;
use thousands::Separable;

pub mod header;

/// Command line arguments for `seqvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest sequence variant VCF", long_about = None)]
pub struct Args {
    /// Maximal number of variants to write out; optional.
    #[clap(long)]
    pub max_var_count: Option<usize>,
    /// The path to the mehari database.
    #[clap(long)]
    pub path_mehari_db: String,
    /// The assumed genome build.
    #[clap(long)]
    pub genomebuild: GenomeRelease,
    /// Path to the pedigree file.
    #[clap(long)]
    pub path_ped: String,
    /// Path to input file.
    #[clap(long)]
    pub path_in: String,
    /// Path to output file.
    #[clap(long)]
    pub path_out: String,
}

/// Return the version of the `varfish-server-worker` crate and `x.y.z` in tests.
fn worker_version() -> &'static str {
    if cfg!(test) {
        "x.y.z"
    } else {
        env!("CARGO_PKG_VERSION")
    }
}

/// Return path component fo rth egiven assembly.
pub fn path_component(genomebuild: GenomeRelease) -> &'static str {
    match genomebuild {
        GenomeRelease::Grch37 => "grch37",
        GenomeRelease::Grch38 => "grch38",
    }
}

/// Process the variants from `input_reader` to `output_writer`.
fn process_variants(
    output_writer: &mut vcf::Writer<std::io::BufWriter<std::fs::File>>,
    input_reader: &mut vcf::Reader<std::io::BufReader<std::fs::File>>,
    output_header: &vcf::Header,
    input_header: &vcf::Header,
    args: &Args,
) -> Result<(), anyhow::Error> {
    // Open the frequency RocksDB database in read only mode.
    tracing::info!("Opening frequency database");
    let rocksdb_path = format!(
        "{}/{}/seqvars/freqs",
        &args.path_mehari_db,
        path_component(args.genomebuild)
    );
    tracing::debug!("RocksDB path = {}", &rocksdb_path);
    let options = rocksdb::Options::default();
    let db_freq = rocksdb::DB::open_cf_for_read_only(
        &options,
        &rocksdb_path,
        ["meta", "autosomal", "gonosomal", "mitochondrial"],
        false,
    )?;

    let cf_autosomal = db_freq.cf_handle("autosomal").unwrap();
    let cf_gonosomal = db_freq.cf_handle("gonosomal").unwrap();
    let cf_mtdna = db_freq.cf_handle("mitochondrial").unwrap();

    // Open the ClinVar RocksDB database in read only mode.
    tracing::info!("Opening ClinVar database");
    let rocksdb_path = format!(
        "{}/{}/seqvars/clinvar/rocksdb",
        &args.path_mehari_db,
        path_component(args.genomebuild)
    );
    tracing::debug!("RocksDB path = {}", &rocksdb_path);
    let options = rocksdb::Options::default();
    let db_clinvar =
        rocksdb::DB::open_cf_for_read_only(&options, &rocksdb_path, ["meta", "clinvar"], false)?;

    let cf_clinvar = db_clinvar.cf_handle("clinvar").unwrap();

    // Open the serialized transcripts.
    tracing::info!("Opening transcript database");
    let tx_db = mehari::annotate::seqvars::load_tx_db(&format!(
        "{}/{}/txs.bin.zst",
        &args.path_mehari_db,
        path_component(args.genomebuild)
    ))?;
    tracing::info!("Building transcript interval trees ...");
    let assembly = if args.genomebuild == GenomeRelease::Grch37 {
        hgvs::static_data::Assembly::Grch37p10
    } else {
        hgvs::static_data::Assembly::Grch38
    };
    let provider = Arc::new(MehariProvider::new(tx_db, assembly));
    let predictor = mehari::annotate::seqvars::csq::ConsequencePredictor::new(provider, assembly);
    tracing::info!("... done building transcript interval trees");

    let start = std::time::Instant::now();
    let mut prev = std::time::Instant::now();
    let mut total_written = 0usize;
    let mut records = input_reader.records(&input_header);
    loop {
        if let Some(record) = records.next() {
            let mut vcf_record = record?;

            // TODO: ignores all but the first alternative allele!
            let vcf_var = annonars::common::keys::Var::from_vcf_allele(&vcf_record, 0);

            // Skip records with a deletion as alternative allele.
            if vcf_var.alternative == "*" {
                continue;
            }

            if prev.elapsed().as_secs() >= 60 {
                tracing::info!("at {:?}", &vcf_var);
                prev = std::time::Instant::now();
            }

            // Only attempt lookups into RocksDB for canonical contigs.
            if annonars::common::cli::is_canonical(vcf_var.chrom.as_str()) {
                // Build key for RocksDB database from `vcf_var`.
                let key: Vec<u8> = vcf_var.clone().into();

                // Annotate with frequency.
                if mehari::annotate::seqvars::CHROM_AUTO.contains(vcf_var.chrom.as_str()) {
                    mehari::annotate::seqvars::annotate_record_auto(
                        &db_freq,
                        &cf_autosomal,
                        &key,
                        &mut vcf_record,
                    )?;
                } else if mehari::annotate::seqvars::CHROM_XY.contains(vcf_var.chrom.as_str()) {
                    mehari::annotate::seqvars::annotate_record_xy(
                        &db_freq,
                        &cf_gonosomal,
                        &key,
                        &mut vcf_record,
                    )?;
                } else if mehari::annotate::seqvars::CHROM_MT.contains(vcf_var.chrom.as_str()) {
                    mehari::annotate::seqvars::annotate_record_mt(
                        &db_freq,
                        &cf_mtdna,
                        &key,
                        &mut vcf_record,
                    )?;
                } else {
                    tracing::trace!(
                        "Record @{:?} on non-canonical chromosome, skipping.",
                        &vcf_var
                    );
                }

                // Annotate with ClinVar information.
                mehari::annotate::seqvars::annotate_record_clinvar(
                    &db_clinvar,
                    &cf_clinvar,
                    &key,
                    &mut vcf_record,
                )?;
            }

            let annonars::common::keys::Var {
                chrom,
                pos,
                reference,
                alternative,
            } = vcf_var;

            // Annotate with variant effect.
            if let Some(ann_fields) =
                predictor.predict(&mehari::annotate::seqvars::csq::VcfVariant {
                    chromosome: chrom,
                    position: pos,
                    reference,
                    alternative,
                })?
            {
                if !ann_fields.is_empty() {
                    vcf_record.info_mut().insert(
                        "ANN".parse()?,
                        Some(vcf::record::info::field::Value::Array(
                            vcf::record::info::field::value::Array::String(
                                ann_fields.iter().map(|ann| Some(ann.to_string())).collect(),
                            ),
                        )),
                    );
                }
            }

            // Write out the record.
            output_writer.write_record(&output_header, &vcf_record)?;
        } else {
            break; // all done
        }

        total_written += 1;
        if let Some(max_var_count) = args.max_var_count {
            if total_written >= max_var_count {
                tracing::warn!(
                    "Stopping after {} records as requested by --max-var-count",
                    total_written
                );
                break;
            }
        }
    }
    tracing::info!(
        "... annotated {} records in {:?}",
        total_written.separate_with_commas(),
        start.elapsed()
    );

    todo!()
}

/// Main entry point for `seqvars ingest` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    common::trace_rss_now();

    tracing::info!("loading pedigree...");
    let pedigree = mehari::ped::PedigreeByName::from_path(&args.path_ped)
        .map_err(|e| anyhow::anyhow!("problem parsing PED file: {}", e))?;
    tracing::info!("pedigre = {:#?}", &pedigree);

    tracing::info!("opening input file...");
    let mut input_reader = {
        let file = std::fs::File::open(&args.path_in)
            .map_err(|e| anyhow::anyhow!("could not open input file {}: {}", &args.path_in, e))
            .map(std::io::BufReader::new)?;
        vcf::reader::Builder
            .build_from_reader(file)
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?
    };

    tracing::info!("processing header...");
    let input_header = input_reader
        .read_header()
        .map_err(|e| anyhow::anyhow!("problem reading VCF header: {}", e))?;
    let output_header = header::build_output_header(
        &input_header,
        &Some(pedigree),
        args.genomebuild,
        worker_version(),
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    let mut output_writer = {
        let writer = std::fs::File::create(&args.path_out).map_err(|e| {
            anyhow::anyhow!(
                "could not output file for writing {}: {}",
                &args.path_out,
                e
            )
        })?;
        let writer = std::io::BufWriter::new(writer);
        vcf::writer::Writer::new(writer)
    };
    output_writer
        .write_header(&output_header)
        .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

    process_variants(
        &mut output_writer,
        &mut input_reader,
        &output_header,
        &input_header,
        args,
    )?;

    tracing::info!(
        "All of `seqvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::common::GenomeRelease;

    macro_rules! set_snapshot_suffix {
        ($($expr:expr),*) => {
            let mut settings = insta::Settings::clone_current();
            settings.set_snapshot_suffix(format!($($expr,)*));
            let _guard = settings.bind_to_scope();
        }
    }

    #[rstest]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn smoke_test_run(#[case] path: &str) -> Result<(), anyhow::Error> {
        set_snapshot_suffix!("{:?}", path.split('/').last().unwrap().replace('.', "_"));

        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_mehari_db: Default::default(),
            path_ped: path.replace(".vcf", ".ped"),
            genomebuild: GenomeRelease::Grch37,
            path_in: path.into(),
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
        };
        super::run(&args_common, &args)
    }
}
