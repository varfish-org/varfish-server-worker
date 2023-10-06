//! Implementation of `seqvars ingest` subcommand.

use std::sync::{Arc, OnceLock};

use crate::common::{self, open_read_maybe_gz, open_write_maybe_gz, GenomeRelease};
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

/// Known keys information and logic for `FORMAT`.
#[derive(Debug)]
struct KnownFormatKeys {
    /// The keys that will be written to the output.
    output_keys: Vec<vcf::record::genotypes::keys::Key>,
    /// The keys that are known from the input keys.
    known_keys: Vec<vcf::record::genotypes::keys::Key>,
    /// Mapping from known to output keys where it is not identity
    known_to_output_map: std::collections::HashMap<
        vcf::record::genotypes::keys::Key,
        vcf::record::genotypes::keys::Key,
    >,
}

impl Default for KnownFormatKeys {
    /// Constructor.
    fn default() -> Self {
        Self {
            output_keys: vec![
                vcf::record::genotypes::keys::key::GENOTYPE, // GT
                vcf::record::genotypes::keys::key::CONDITIONAL_GENOTYPE_QUALITY, // GQ
                vcf::record::genotypes::keys::key::READ_DEPTH, // DP
                vcf::record::genotypes::keys::key::READ_DEPTHS, // AD
                vcf::record::genotypes::keys::key::PHASE_SET, // PS
            ],
            known_keys: vec![
                vcf::record::genotypes::keys::key::GENOTYPE,
                vcf::record::genotypes::keys::key::CONDITIONAL_GENOTYPE_QUALITY,
                vcf::record::genotypes::keys::key::READ_DEPTH,
                vcf::record::genotypes::keys::key::READ_DEPTHS,
                vcf::record::genotypes::keys::key::PHASE_SET, // PS
                "SQ".parse().expect("invalid key: SQ"),       // written as AD
            ],
            known_to_output_map: vec![(
                "SQ".parse().expect("invalid key: SQ"),
                vcf::record::genotypes::keys::key::CONDITIONAL_GENOTYPE_QUALITY,
            )]
            .into_iter()
            .collect(),
        }
    }
}

impl KnownFormatKeys {
    /// Map from known to output key.
    pub fn known_to_output(
        &self,
        key: &vcf::record::genotypes::keys::Key,
    ) -> vcf::record::genotypes::keys::Key {
        self.known_to_output_map.get(key).unwrap_or(key).clone()
    }
}

/// The known `FORMAT` keys.
static KNOWN_FORMAT_KEYS: OnceLock<KnownFormatKeys> = OnceLock::new();

/// Transform the ``FORMAT`` key if known.
fn transform_format_value(
    value: &Option<&vcf::record::genotypes::sample::Value>,
    key: &vcf::record::genotypes::keys::Key,
    allele_no: usize,
) -> Option<Option<vcf::record::genotypes::sample::Value>> {
    if let Some(value) = value {
        Some(Some(match key.as_ref() {
            "AD" => {
                // Only write out current allele as AD.
                match *value {
                    vcf::record::genotypes::sample::Value::Array(
                        vcf::record::genotypes::sample::value::Array::Integer(ad_values),
                    ) => vcf::record::genotypes::sample::Value::Integer(
                        ad_values[allele_no].expect("SQ should be integer value"),
                    ),
                    _ => return None, // unreachable!("FORMAT/AD must be array of integer"),
                }
            }
            "SQ" => {
                // SQ is written as AD.
                match *value {
                    vcf::record::genotypes::sample::Value::Float(sq_value) => {
                        vcf::record::genotypes::sample::Value::Float(*sq_value)
                    }
                    vcf::record::genotypes::sample::Value::Array(
                        vcf::record::genotypes::sample::value::Array::Float(sq_values),
                    ) => vcf::record::genotypes::sample::Value::Float(
                        sq_values[allele_no].expect("SQ should be float value"),
                    ),
                    _ => return None, // unreachable!("FORMAT/PS must be integer"),
                }
            }
            _ => return None, // unreachable!("unknown key: {:?}", key),
        }))
    } else {
        Some(None)
    }
}

/// Copy the `FORMAT/GQ` fields for all samples.
///
/// The implementation assumes that there are no duplicates in the output keys when mapped
/// from input keys.
fn copy_format(
    input_record: &vcf::Record,
    builder: vcf::record::Builder,
    idx_output_to_input: &[usize],
    allele_no: usize,
    known_format_keys: &KnownFormatKeys,
) -> Result<vcf::record::Builder, anyhow::Error> {
    let keys_from_input_known = input_record
        .genotypes()
        .keys()
        .iter()
        .filter(|k| known_format_keys.known_keys.contains(*k))
        .cloned()
        .collect::<Vec<_>>();
    let output_keys = keys_from_input_known
        .iter()
        .map(|k| known_format_keys.known_to_output(k).clone())
        .collect::<Vec<_>>();

    let values = idx_output_to_input
        .iter()
        .copied()
        .map(|input_idx| {
            let sample = input_record
                .genotypes()
                .get_index(input_idx)
                .expect("input_idx must be valid here");
            keys_from_input_known
                .iter()
                .map(|key| {
                    let input_value = sample.get(key).expect("key must be valid");
                    if let Some(value) = transform_format_value(&input_value, key, allele_no) {
                        value
                    } else if known_format_keys.output_keys.contains(key) {
                        input_value.cloned()
                    } else {
                        unreachable!("don't know how to handle key: {:?}", key)
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let genotypes = vcf::record::Genotypes::new(
        vcf::record::genotypes::Keys::try_from(output_keys).expect("invalid keys"),
        values,
    );

    Ok(builder.set_genotypes(genotypes))
}

/// Process the variants from `input_reader` to `output_writer`.
fn process_variants<R, W>(
    output_writer: &mut vcf::Writer<W>,
    input_reader: &mut vcf::Reader<R>,
    output_header: &vcf::Header,
    input_header: &vcf::Header,
    args: &Args,
) -> Result<(), anyhow::Error>
where
    R: std::io::BufRead,
    W: std::io::Write,
{
    // Open the frequency RocksDB database in read only mode.
    tracing::info!("Opening frequency database");
    let rocksdb_path = format!(
        "{}/{}/seqvars/freqs/rocksdb",
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

    // Build mapping from output sample index to input sample index.
    let idx_output_to_input = {
        let output_sample_to_idx = output_header
            .sample_names()
            .iter()
            .enumerate()
            .map(|(idx, name)| (name, idx))
            .collect::<std::collections::HashMap<_, _>>();
        let mut res = vec![usize::MAX; output_header.sample_names().len()];
        for (input_idx, sample) in input_header.sample_names().iter().enumerate() {
            res[output_sample_to_idx[sample]] = input_idx;
        }
        res
    };

    // Read through input file, construct output records, and annotate these.
    let start = std::time::Instant::now();
    let mut prev = std::time::Instant::now();
    let mut total_written = 0usize;
    let mut records = input_reader.records(input_header);
    let known_format_keys = KNOWN_FORMAT_KEYS.get_or_init(Default::default);
    loop {
        if let Some(input_record) = records.next() {
            let input_record = input_record?;

            for (allele_no, alt_allele) in input_record.alternate_bases().iter().enumerate() {
                // Construct record with first few fields describing one variant allele.
                let builder = vcf::Record::builder()
                    .set_chromosome(input_record.chromosome().clone())
                    .set_position(input_record.position())
                    .set_reference_bases(input_record.reference_bases().clone())
                    .set_alternate_bases(vcf::record::AlternateBases::from(vec![
                        alt_allele.clone()
                    ]));

                // Copy over the well-known FORMAT fields and construct output record.
                let builder = copy_format(
                    &input_record,
                    builder,
                    &idx_output_to_input,
                    allele_no,
                    known_format_keys,
                )?;

                let mut output_record = builder.build()?;

                // Obtain annonars variant key from current allele for RocksDB lookup.
                let vcf_var = annonars::common::keys::Var::from_vcf_allele(&output_record, 0);

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
                            &mut output_record,
                        )?;
                    } else if mehari::annotate::seqvars::CHROM_XY.contains(vcf_var.chrom.as_str()) {
                        mehari::annotate::seqvars::annotate_record_xy(
                            &db_freq,
                            &cf_gonosomal,
                            &key,
                            &mut output_record,
                        )?;
                    } else if mehari::annotate::seqvars::CHROM_MT.contains(vcf_var.chrom.as_str()) {
                        mehari::annotate::seqvars::annotate_record_mt(
                            &db_freq,
                            &cf_mtdna,
                            &key,
                            &mut output_record,
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
                        &mut output_record,
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
                        output_record.info_mut().insert(
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
                output_writer.write_record(output_header, &output_record)?;
                total_written += 1;
            }
        } else {
            break; // all done
        }

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

    Ok(())
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
        vcf::reader::Builder
            .build_from_reader(open_read_maybe_gz(&args.path_in)?)
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

    let mut output_writer = { vcf::writer::Writer::new(open_write_maybe_gz(&args.path_out)?) };
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
    use std::io::Read;

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
    #[case("tests/seqvars/ingest/NA12878_dragen.vcf")]
    #[case("tests/seqvars/ingest/Case_1.vcf")]
    fn result_snapshot_test(#[case] path: &str) -> Result<(), anyhow::Error> {
        set_snapshot_suffix!("{}", path.split('/').last().unwrap().replace('.', "_"));

        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_mehari_db: "tests/seqvars/ingest/db".into(),
            path_ped: path.replace(".vcf", ".ped"),
            genomebuild: GenomeRelease::Grch37,
            path_in: path.into(),
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }

    fn read_to_bytes<P>(path: P) -> Result<Vec<u8>, anyhow::Error>
    where
        P: AsRef<std::path::Path>,
    {
        let mut f = std::fs::File::open(&path).expect("no file found");
        let metadata = std::fs::metadata(&path).expect("unable to read metadata");
        let mut buffer = vec![0; metadata.len() as usize];
        f.read_exact(&mut buffer).expect("buffer overflow");
        Ok(buffer)
    }

    #[test]
    fn result_snapshot_test_gz() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let path_in: String = "tests/seqvars/ingest/NA12878_dragen.vcf.gz".into();
        let path_ped = path_in.replace(".vcf.gz", ".ped");
        let path_out = tmpdir
            .join("out.vcf.gz")
            .to_str()
            .expect("invalid path")
            .into();
        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_mehari_db: "tests/seqvars/ingest/db".into(),
            path_ped,
            genomebuild: GenomeRelease::Grch37,
            path_in,
            path_out,
        };
        super::run(&args_common, &args)?;

        let mut buffer = Vec::new();
        hxdmp::hexdump(&read_to_bytes(&args.path_out)?, &mut buffer)?;
        insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }
}
