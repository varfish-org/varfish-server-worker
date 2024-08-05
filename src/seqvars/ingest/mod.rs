//! Implementation of `seqvars ingest` subcommand.

use std::sync::{Arc, OnceLock};

use crate::{
    common::{self, genotype_to_string, strip_gt_leading_slash, worker_version, GenomeRelease},
    flush_and_shutdown,
};
use mehari::annotate::seqvars::provider::Provider as MehariProvider;
use mehari::common::noodles::{open_vcf_reader, open_vcf_writer, AsyncVcfReader, AsyncVcfWriter};
use noodles::vcf;
use thousands::Separable;
use tokio::io::AsyncWriteExt;

pub mod header;

/// Command line arguments for `seqvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest sequence variant VCF", long_about = None)]
pub struct Args {
    /// Value to write to `##fileDate`.
    #[arg(long)]
    pub file_date: String,
    /// The case UUID to write out.
    #[clap(long)]
    pub case_uuid: uuid::Uuid,
    /// The assumed genome build.
    #[clap(long)]
    pub genomebuild: GenomeRelease,

    /// The path to the mehari database.
    #[clap(long)]
    pub path_mehari_db: String,
    /// Path to the pedigree file.
    #[clap(long)]
    pub path_ped: String,
    /// Path to input file.
    #[clap(long)]
    pub path_in: String,
    /// Path to output file.
    #[clap(long)]
    pub path_out: String,

    /// Maximal number of variants to write out; optional.
    #[clap(long)]
    pub max_var_count: Option<usize>,
    /// Per-file identifier mapping, either a JSON or @-prefixed path to JSON.
    #[clap(long)]
    pub id_mapping: Option<String>,
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
    output_keys: Vec<String>,
    /// The keys that are known from the input keys.
    known_keys: Vec<String>,
    /// Mapping from known to output keys where it is not identity
    known_to_output_map: std::collections::HashMap<String, String>,
}

impl Default for KnownFormatKeys {
    /// Constructor.
    fn default() -> Self {
        use noodles::vcf::variant::record::samples::keys::key;
        Self {
            output_keys: vec![
                key::GENOTYPE.to_string(),                     // GT
                key::CONDITIONAL_GENOTYPE_QUALITY.to_string(), // GQ
                key::READ_DEPTH.to_string(),                   // DP
                key::READ_DEPTHS.to_string(),                  // AD
                key::PHASE_SET.to_string(),                    // PS
            ],
            known_keys: vec![
                key::GENOTYPE.to_string(),
                key::CONDITIONAL_GENOTYPE_QUALITY.to_string(),
                key::READ_DEPTH.to_string(),
                key::READ_DEPTHS.to_string(),
                key::PHASE_SET.to_string(), // PS
                "SQ".to_string(),           // written as AD
            ],
            known_to_output_map: vec![(
                "SQ".to_string(),
                key::CONDITIONAL_GENOTYPE_QUALITY.to_string(),
            )]
            .into_iter()
            .collect(),
        }
    }
}

impl KnownFormatKeys {
    /// Map from known to output key.
    pub fn known_to_output(&self, key: &str) -> String {
        self.known_to_output_map
            .get(key)
            .cloned()
            .unwrap_or(key.to_string())
    }
}

/// The known `FORMAT` keys.
static KNOWN_FORMAT_KEYS: OnceLock<KnownFormatKeys> = OnceLock::new();

/// Regular expression for parsing `GT` values.
static GT_RE: OnceLock<regex::Regex> = OnceLock::new();

/// Transform the ``FORMAT`` key if known.
fn transform_format_value(
    value: &Option<&vcf::variant::record_buf::samples::sample::value::Value>,
    key: &str,
    allele_no: usize,
    sample: &vcf::variant::record_buf::samples::Sample<'_>,
) -> Option<Option<vcf::variant::record_buf::samples::sample::value::Value>> {
    let gt_re = GT_RE
        .get_or_init(|| regex::Regex::new(r"([^\|]+)([/|])([^\|]+)").expect("could not parse RE"));

    let curr_allele = format!("{}", allele_no);

    fn transform_allele(allele_to_transform: &str, curr_allele: &str) -> &'static str {
        if allele_to_transform == curr_allele {
            "1"
        } else {
            "0"
        }
    }

    if let Some(value) = value {
        Some(Some(match key {
            "GT" => {
                let gt = if let Some(Some(
                    vcf::variant::record_buf::samples::sample::value::Value::Genotype(gt),
                )) =
                    sample.get(noodles::vcf::variant::record::samples::keys::key::GENOTYPE)
                {
                    strip_gt_leading_slash(
                        &genotype_to_string(&gt)
                            .unwrap_or_else(|_| panic!("invalid genotype: {:?}", &gt)),
                    )
                    .to_string()
                } else {
                    unreachable!("FORMAT/GT must be string")
                };
                if ["./.", ".|.", "."].contains(&gt.as_str()) {
                    // no need to transform no-call
                    vcf::variant::record_buf::samples::sample::value::Value::String(gt)
                } else {
                    // transform all others
                    let gt_captures = gt_re
                        .captures(&gt)
                        .unwrap_or_else(|| panic!("FORMAT/GT cannot be parsed: {}", &gt));
                    let gt_1 = gt_captures.get(1).expect("must be captured").as_str();
                    let gt_2 = gt_captures.get(2).expect("must be captured").as_str();
                    let gt_3 = gt_captures.get(3).expect("must be captured").as_str();

                    let new_gt = format!(
                        "{}{}{}",
                        transform_allele(gt_1, &curr_allele),
                        gt_2,
                        transform_allele(gt_3, &curr_allele),
                    );

                    vcf::variant::record_buf::samples::sample::value::Value::String(new_gt)
                }
            }
            "AD" => {
                let dp = match sample
                    .get(noodles::vcf::variant::record::samples::keys::key::READ_DEPTH)
                    .expect("FORMAT/DP must be present")
                    .cloned()
                    .unwrap_or_else(|| unreachable!("FORMAT/DP must be present and not None"))
                {
                    vcf::variant::record_buf::samples::sample::value::Value::Integer(dp) => dp,
                    _ => unreachable!("FORMAT/DP must be integer"),
                };

                // Only write out reference and current allele as AD.
                match *value {
                    vcf::variant::record_buf::samples::sample::value::Value::Array(
                        vcf::variant::record_buf::samples::sample::value::Array::Integer(ad_values),
                    ) => {
                        let ad = ad_values[allele_no].expect("AD should be integer value");
                        vcf::variant::record_buf::samples::sample::value::Value::Array(
                            vcf::variant::record_buf::samples::sample::value::Array::Integer(vec![
                                Some(dp - ad),
                                Some(ad),
                            ]),
                        )
                    }
                    _ => return None, // unreachable!("FORMAT/AD must be array of integer"),
                }
            }
            "SQ" => {
                // SQ is written as AD.
                match *value {
                    vcf::variant::record_buf::samples::sample::value::Value::Float(sq_value) => {
                        vcf::variant::record_buf::samples::sample::value::Value::Float(*sq_value)
                    }
                    vcf::variant::record_buf::samples::sample::value::Value::Array(
                        vcf::variant::record_buf::samples::sample::value::Array::Float(sq_values),
                    ) => vcf::variant::record_buf::samples::sample::value::Value::Integer(
                        sq_values[allele_no - 1]
                            .expect("SQ should be float value")
                            .round() as i32,
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
    record_buf: &vcf::variant::RecordBuf,
    builder: vcf::variant::record_buf::builder::Builder,
    idx_output_to_input: &[usize],
    allele_no: usize,
    known_format_keys: &KnownFormatKeys,
) -> Result<vcf::variant::record_buf::builder::Builder, anyhow::Error> {
    let keys_from_input_known = record_buf
        .samples()
        .keys()
        .as_ref()
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
            let sample = record_buf
                .samples()
                .get_index(input_idx)
                .expect("input_idx must be valid here");
            keys_from_input_known
                .iter()
                .map(|key| {
                    let input_value = sample.get(key).expect("key must be valid");
                    if let Some(value) =
                        transform_format_value(&input_value, key, allele_no, &sample)
                    {
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

    let genotypes =
        vcf::variant::record_buf::samples::Samples::new(output_keys.into_iter().collect(), values);

    Ok(builder.set_samples(genotypes))
}

/// Process the variants from `input_reader` to `output_writer`.
async fn process_variants(
    output_writer: &mut AsyncVcfWriter,
    input_reader: &mut AsyncVcfReader,
    output_header: &vcf::Header,
    input_header: &vcf::Header,
    id_mapping: &Option<indexmap::IndexMap<String, String>>,
    args: &Args,
) -> Result<(), anyhow::Error> {
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
        biocommons_bioutils::assemblies::Assembly::Grch37p10
    } else {
        biocommons_bioutils::assemblies::Assembly::Grch38
    };
    let provider = Arc::new(MehariProvider::new(tx_db, assembly, Default::default()));
    let predictor = mehari::annotate::seqvars::csq::ConsequencePredictor::new(
        provider,
        assembly,
        Default::default(),
    );
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
            let sample = if let Some(id_mapping) = id_mapping {
                id_mapping.get(sample).expect("checked earlier")
            } else {
                sample
            };
            res[output_sample_to_idx[sample]] = input_idx;
        }
        res
    };

    // Read through input file, construct output records, and annotate these.
    let start = std::time::Instant::now();
    let mut prev = std::time::Instant::now();
    let mut total_written = 0usize;
    let mut input_record = vcf::variant::RecordBuf::default();
    let known_format_keys = KNOWN_FORMAT_KEYS.get_or_init(Default::default);
    loop {
        let bytes_read = input_reader
            .read_record_buf(input_header, &mut input_record)
            .await
            .map_err(|e| anyhow::anyhow!("problem reading VCF file: {}", e))?;
        if bytes_read == 0 {
            break; // EOF
        }

        for (allele_no, alt_allele) in input_record.alternate_bases().as_ref().iter().enumerate() {
            let allele_no = allele_no + 1;
            // Construct record with first few fields describing one variant allele.
            let builder = noodles::vcf::variant::RecordBuf::builder()
                .set_reference_sequence_name(input_record.reference_sequence_name())
                .set_variant_start(
                    input_record
                        .variant_start()
                        .ok_or_else(|| anyhow::anyhow!("missing start position"))?,
                )
                .set_reference_bases(input_record.reference_bases())
                .set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                    alt_allele.clone(),
                ]));

            // Copy over the well-known FORMAT fields and construct output record.
            let builder = copy_format(
                &input_record,
                builder,
                &idx_output_to_input,
                allele_no,
                known_format_keys,
            )?;

            // Build the output `RecordBuf`.
            let mut output_record = builder.build();

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
                    tracing::debug!(
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
                        Some(vcf::variant::record_buf::info::field::Value::Array(
                            vcf::variant::record_buf::info::field::value::Array::String(
                                ann_fields.iter().map(|ann| Some(ann.to_string())).collect(),
                            ),
                        )),
                    );
                }
            }

            // Write out the record.
            output_writer
                .write_variant_record(output_header, &output_record)
                .await?;
            total_written += 1;
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
pub async fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    common::trace_rss_now();

    tracing::info!("loading pedigree...");
    let pedigree = mehari::ped::PedigreeByName::from_path(&args.path_ped)
        .map_err(|e| anyhow::anyhow!("problem parsing PED file: {}", e))?;
    tracing::info!("pedigre = {:#?}", &pedigree);

    tracing::info!("opening input file...");
    let mut input_reader = open_vcf_reader(&args.path_in)
        .await
        .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?;

    tracing::info!("loading file identifier mappings...");
    let id_mapping = args
        .id_mapping
        .as_ref()
        .map(
            |id_mapping| -> Result<Option<indexmap::IndexMap<_, _>>, anyhow::Error> {
                let id_mappings = if id_mapping.starts_with('@') {
                    crate::common::id_mapping::FileIdentifierMappings::load_from_path(
                        id_mapping.trim_start_matches('@'),
                    )
                    .map_err(|e| {
                        anyhow::anyhow!(
                            "could not load ID mapping from file {:?}: {}",
                            id_mapping.trim_start_matches('@'),
                            e
                        )
                    })
                } else {
                    crate::common::id_mapping::FileIdentifierMappings::load_from_json(id_mapping)
                        .map_err(|e| {
                            anyhow::anyhow!(
                                "could not load ID mapping from JSON string {:?}: {}",
                                &id_mapping,
                                &e,
                            )
                        })
                }?;

                if id_mappings.file_names().contains(&args.path_in) {
                    tracing::debug!("- we have an ID mapping for {}", &args.path_in);
                    Ok(Some(
                        id_mappings
                            .mapping_for_file(&args.path_in)
                            .cloned()
                            .expect("checked above"),
                    ))
                } else {
                    tracing::debug!("- we have no ID mapping for {}", &args.path_in);
                    Ok(None)
                }
            },
        )
        .transpose()?
        .flatten();

    tracing::info!("processing header...");
    let mut input_header = input_reader
        .read_header()
        .await
        .map_err(|e| anyhow::anyhow!("problem reading VCF header: {}", e))?;
    let output_header = header::build_output_header(
        &input_header,
        &Some(pedigree),
        &id_mapping,
        args.genomebuild,
        &args.file_date,
        &args.case_uuid,
        worker_version(),
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    // Work around glnexus issue with RNC.
    if let Some(format) = input_header.formats_mut().get_mut("RNC") {
        *format.number_mut() = vcf::header::record::value::map::format::Number::Count(1);
        *format.type_mut() = vcf::header::record::value::map::format::Type::String;
    }

    // Use output file helper.
    let out_path_helper = crate::common::s3::OutputPathHelper::new(&args.path_out)?;

    {
        let mut output_writer = open_vcf_writer(out_path_helper.path_out()).await?;
        output_writer
            .write_header(&output_header)
            .await
            .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

        process_variants(
            &mut output_writer,
            &mut input_reader,
            &output_header,
            &input_header,
            &id_mapping,
            args,
        )
        .await?;

        flush_and_shutdown!(output_writer);
    }

    out_path_helper.create_tbi_for_bgzf().await?;
    out_path_helper.upload_for_s3().await?;

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

    #[rstest]
    #[case::clair3_glnexus("tests/seqvars/ingest/clair3_glnexus.vcf")]
    #[case::dragen_07_021_624_3_10_4("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case::dragen_07_021_624_3_10_9("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case::gatk_hc_3_7("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case::gatk_hc_4_4_0_0("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    #[case::dragen_na12787("tests/seqvars/ingest/NA12878_dragen.vcf")]
    #[case::gatk_hc_case_1("tests/seqvars/ingest/Case_1.vcf")]
    #[tokio::test]
    async fn result_snapshot_test(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!(
            "{}",
            path.split('/').last().unwrap().replace('.', "_")
        );

        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            file_date: String::from("20230421"),
            case_uuid: uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
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
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }

    #[tokio::test]
    async fn result_snapshot_test_gz() -> Result<(), anyhow::Error> {
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
            file_date: String::from("20230421"),
            case_uuid: uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
            max_var_count: None,
            path_mehari_db: "tests/seqvars/ingest/db".into(),
            path_ped,
            genomebuild: GenomeRelease::Grch37,
            path_in,
            path_out,
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;
        insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }

    #[rstest]
    #[case::dragen_na12787("tests/seqvars/ingest/NA12878_dragen.vcf")]
    #[case::gatk_hc_case_1("tests/seqvars/ingest/Case_1.vcf")]
    #[tokio::test]
    async fn result_snapshot_test_with_id_map(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!(
            "{}",
            path.split('/').last().unwrap().replace('.', "_")
        );

        let tmpdir = temp_testdir::TempDir::default();

        let path_ped = path.replace(".vcf", ".custom_id.ped");
        let path_out = tmpdir
            .join("out.vcf")
            .to_str()
            .expect("invalid path")
            .into();
        let args_common = Default::default();
        let args = super::Args {
            file_date: String::from("20230421"),
            case_uuid: uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
            max_var_count: None,
            path_mehari_db: "tests/seqvars/ingest/db".into(),
            path_ped,
            genomebuild: GenomeRelease::Grch37,
            path_in: path.into(),
            path_out,
            id_mapping: Some(
                r#"
                {
                    "mappings": [
                        {
                            "path": "tests/seqvars/ingest/NA12878_dragen.vcf",
                            "entries": [
                                {
                                    "src": "NA12878",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/seqvars/ingest/Case_1.vcf",
                            "entries": [
                                {
                                    "src": "Case_1_index-N1-DNA1-WGS1",
                                    "dst": "Case_1_index"
                                },
                                {
                                    "src": "Case_1_mother-N1-DNA1-WGS1",
                                    "dst": "Case_1_mother"
                                },
                                {
                                    "src": "Case_1_father-N1-DNA1-WGS1",
                                    "dst": "Case_1_father"
                                }
                            ]
                        }
                    ]
                }
                "#
                .to_string(),
            ),
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
}
