//! Implementation of `strucvars ingest` subcommand.

use crate::common::noodles::open_vcf_readers;
use crate::common::{self, worker_version, GenomeRelease};
use crate::flush_and_shutdown;
use futures::future::join_all;
use mehari::annotate::strucvars::bnd::Breakend;
use mehari::annotate::strucvars::guess_sv_caller;
use mehari::common::noodles::{
    open_vcf_writer, AsyncVcfWriter, NoodlesVariantReader as _, VariantReader,
};
use noodles::vcf;
use rand_core::SeedableRng;
use tokio::io::AsyncWriteExt;

pub mod header;

/// Command line arguments for `strucvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest structural variant VCF", long_about = None)]
pub struct Args {
    /// Value to write to `##fileDate`.
    #[arg(long)]
    pub file_date: String,
    /// Value to write out for `##x-varfish-case-uuid`.
    #[arg(long)]
    pub case_uuid: String,
    /// The assumed genome build.
    #[clap(long)]
    pub genomebuild: GenomeRelease,

    /// Path to the pedigree file.
    #[clap(long)]
    pub path_ped: String,
    /// Path to input files.
    #[clap(long, required = true)]
    pub path_in: Vec<String>,
    /// Path to coverage VCF files from maelstrom; optional.
    #[clap(long)]
    pub path_cov_vcf: Vec<String>,
    /// Path to output file.
    #[clap(long)]
    pub path_out: String,

    /// Minimal reciprocal overlap to require.
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Slack to use around break-ends.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: i32,
    /// Slack to use around insertions.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: i32,

    /// Seed for random number generator (UUIDs), if any.
    #[arg(long)]
    pub rng_seed: Option<u64>,
    /// Maximal number of variants to write out; optional.
    #[clap(long)]
    pub max_var_count: Option<usize>,
    /// Per-file identifier mapping, either a JSON or @-prefixed path to JSON.
    #[clap(long)]
    pub id_mapping: Option<String>,
}

async fn write_ingest_record(
    output_header: &vcf::Header,
    writer: &mut AsyncVcfWriter,
    input_record: &vcf::variant::RecordBuf,
) -> Result<(), anyhow::Error> {
    // copy over CHROM, POS, REF
    let builder = vcf::variant::record_buf::builder::Builder::default()
        .set_reference_sequence_name(input_record.reference_sequence_name())
        .set_variant_start(input_record.variant_start().expect("no variant_start?"))
        .set_reference_bases(input_record.reference_bases());

    // copy over first ALT allele, remove any SV sub types
    if input_record.alternate_bases().as_ref().len() != 1 {
        anyhow::bail!(
            "unexpected number of ALT alleles (should be ==1) in: {:?}",
            input_record.alternate_bases()
        );
    }
    let alt_0 = &input_record
        .alternate_bases()
        .as_ref()
        .iter()
        .next()
        .expect("alternate_bases cannot be empty");
    let (sv_type, bnd, mut builder) = if alt_0.contains('[') || alt_0.contains(']') {
        (
            "BND".to_string(),
            Some(Breakend::from_ref_alt_str(
                input_record.reference_bases(),
                &input_record
                    .alternate_bases()
                    .as_ref()
                    .iter()
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("no alternate allele?"))?
                    .to_string(),
            )?),
            builder.set_alternate_bases(input_record.alternate_bases().clone()),
        )
    } else if alt_0.contains('<') && alt_0.contains('>') {
        let sv_type = alt_0
            .split('<')
            .nth(1)
            .ok_or_else(|| anyhow::anyhow!("no < in SV type"))?
            .split('>')
            .next()
            .ok_or_else(|| anyhow::anyhow!("no > in SV type"))?
            .split(':')
            .next()
            .expect("empty SVTYPE?");
        (
            sv_type.to_string(),
            None,
            builder.set_alternate_bases(vcf::variant::record_buf::AlternateBases::from(vec![
                format!("<{}>", sv_type),
            ])),
        )
    } else {
        anyhow::bail!("unexpected alternate base type: {:?}", &alt_0)
    };

    // copy over FORMAT tags, all except FT
    let mut keys_with_value = std::collections::HashSet::<String>::new();
    let output_format_values = input_record
        .samples()
        .values()
        .map(|g| {
            g.keys()
                .as_ref()
                .iter()
                .zip(g.values().iter())
                .filter(|(k, _)| k.as_str() != "FT")
                .map(|(k, v)| {
                    if v.is_some() {
                        keys_with_value.insert(k.clone());
                    }

                    v.clone()
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let output_keys = input_record
        .samples()
        .keys()
        .as_ref()
        .iter()
        .filter(|k| k.as_str() != "FT")
        .cloned()
        .map(|k| {
            if k.as_str() == "CN" {
                "cn".parse().expect("invalid key: cn")
            } else {
                k
            }
        })
        .collect();
    builder = builder.set_samples(vcf::variant::record_buf::samples::Samples::new(
        output_keys,
        output_format_values,
    ));

    // copy over INFO tags
    // Note: annsv will be added only in "strucvars query"
    let mut info: vcf::variant::record_buf::Info = Default::default();
    match sv_type.as_str() {
        "DEL" | "DUP" | "CNV" => {
            let claim = if keys_with_value.contains("pev") || keys_with_value.contains("srv") {
                "DJ"
            } else {
                "D"
            };
            info.insert(
                vcf::variant::record::info::field::key::SV_CLAIM.to_string(),
                Some(vcf::variant::record_buf::info::field::Value::Array(
                    vcf::variant::record_buf::info::field::value::Array::String(vec![Some(
                        claim.to_string(),
                    )]),
                )),
            );
        }
        "INS" | "INV" | "BND" => {
            info.insert(
                vcf::variant::record::info::field::key::SV_CLAIM.to_string(),
                Some(vcf::variant::record_buf::info::field::Value::Array(
                    vcf::variant::record_buf::info::field::value::Array::String(vec![Some(
                        "J".to_string(),
                    )]),
                )),
            );
        }
        _ => anyhow::bail!("unexpected SV type: {}", sv_type),
    }
    info.insert(
        vcf::variant::record::info::field::key::SV_TYPE.to_string(),
        Some(vcf::variant::record_buf::info::field::Value::String(
            sv_type.to_string(),
        )),
    );
    if let Some(Some(vcf::variant::record_buf::info::field::Value::Integer(end))) = input_record
        .info()
        .get(vcf::variant::record::info::field::key::END_POSITION)
    {
        info.insert(
            vcf::variant::record::info::field::key::END_POSITION.to_string(),
            Some(vcf::variant::record_buf::info::field::Value::Integer(*end)),
        );

        if sv_type == "BND" {
            info.insert(
                "chr2".to_string(),
                Some(vcf::variant::record_buf::info::field::Value::String(
                    bnd.expect("must be set here").chrom.clone(),
                )),
            );
        } else {
            let pos: usize = input_record
                .variant_start()
                .expect("no variant_start?")
                .into();
            let sv_len: usize = *end as usize - pos + 1;
            info.insert(
                vcf::variant::record::info::field::key::SV_LENGTHS.to_string(),
                Some(vcf::variant::record_buf::info::field::Value::Array(
                    vcf::variant::record_buf::info::field::value::Array::Integer(vec![Some(
                        sv_len as i32,
                    )]),
                )),
            );
        }
    }

    fn map_caller(caller: &str) -> Result<Option<String>, anyhow::Error> {
        if caller.starts_with("DELLYv") {
            Ok(Some("Delly".to_string()))
        } else if caller.starts_with("DRAGEN_CNVv") {
            Ok(Some("DragenCnv".to_string()))
        } else if caller.starts_with("DRAGEN_SVv") {
            Ok(Some("DragenSv".to_string()))
        } else if caller.starts_with("GATK_GCNVv") {
            Ok(Some("Gcnv".to_string()))
        } else if caller.starts_with("MANTAv") {
            Ok(Some("Manta".to_string()))
        } else if caller.starts_with("POPDELv") {
            Ok(Some("Popdel".to_string()))
        } else if caller.starts_with("MELTv") {
            Ok(Some("Melt".to_string()))
        } else if caller.starts_with("SNIFFLESv") {
            Ok(Some("Sniffles".to_string()))
        } else {
            anyhow::bail!("unknown caller: {}", caller)
        }
    }

    if let Some(Some(callers)) = input_record.info().get("callers") {
        if let vcf::variant::record_buf::info::field::Value::Array(
            vcf::variant::record_buf::info::field::value::Array::String(callers),
        ) = callers
        {
            let output_callers = callers
                .iter()
                .flatten()
                .map(|caller| map_caller(caller))
                .collect::<Result<Vec<_>, _>>()?;
            info.insert(
                "callers".to_string(),
                Some(vcf::variant::record_buf::info::field::Value::Array(
                    vcf::variant::record_buf::info::field::value::Array::String(output_callers),
                )),
            );
        } else if let vcf::variant::record_buf::info::field::Value::String(caller) = callers {
            let output_callers = vec![map_caller(caller)?];
            info.insert(
                "callers".to_string(),
                Some(vcf::variant::record_buf::info::field::Value::Array(
                    vcf::variant::record_buf::info::field::value::Array::String(output_callers),
                )),
            );
        }
    } else {
        anyhow::bail!("no callers INFO tag found");
    }

    builder = builder.set_info(info);

    let record = builder.build();

    writer
        .write_variant_record(output_header, &record)
        .await
        .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
}

/// Write out variants from input files.
async fn process_variants(
    pedigree: &mehari::ped::PedigreeByName,
    output_header: &vcf::Header,
    output_writer: &mut AsyncVcfWriter,
    input_readers: Vec<VariantReader>,
    input_header: &[vcf::Header],
    input_sv_callers: &[mehari::annotate::strucvars::SvCaller],
    args: &Args,
) -> Result<(), anyhow::Error> {
    // Initialize the random number generator from command line seed if given or local entropy
    // source.
    let mut rng = if let Some(rng_seed) = args.rng_seed {
        rand::rngs::StdRng::seed_from_u64(rng_seed)
    } else {
        rand::rngs::StdRng::from_entropy()
    };

    // Create temporary directory.  We will create one temporary file (containing `jsonl`
    // seriealized `VarFishStrucvarTsvRecord`s) for each SV type and contig.
    let tmp_dir = tempfile::TempDir::new()?;

    // Read through input VCF files and write out to temporary files.
    tracing::info!("converting input VCF files to temporary files...");
    let mut input_readers = input_readers;
    for (mut reader, sv_caller, header) in itertools::izip!(
        input_readers.drain(..),
        input_sv_callers.iter(),
        input_header.iter()
    ) {
        mehari::annotate::strucvars::run_vcf_to_jsonl(
            pedigree,
            &mut reader,
            header,
            sv_caller,
            &tmp_dir,
            &mut std::collections::HashMap::new(),
            &mut rng,
        )
        .await?;
    }
    tracing::info!("... done converting input files");

    tracing::info!("clustering SVs to output...");
    // Read through temporary files by contig, cluster by overlap as configured, and write to `writer`.
    for contig_no in 1..=25 {
        tracing::info!(
            "  contig: {}",
            annonars::common::cli::CANONICAL[contig_no - 1]
        );
        let clusters = mehari::annotate::strucvars::read_and_cluster_for_contig(
            &tmp_dir,
            contig_no,
            args.slack_ins,
            args.slack_bnd,
            args.min_overlap,
        )?;
        for record in clusters {
            write_ingest_record(output_header, output_writer, &record.try_into()?).await?;
        }
    }
    tracing::info!("... done clustering SVs to output");

    Ok(())
}

/// Main entry point for `strucvars ingest` sub command.
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
    let mut input_readers = open_vcf_readers(&args.path_in).await?;

    tracing::info!("loading file identifier mappings...");
    let id_mappings = args
        .id_mapping
        .as_ref()
        .map(
            |id_mapping| -> Result<crate::common::id_mapping::FileIdentifierMappings, anyhow::Error> {
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

                let paths_in = args.path_in.iter().cloned().collect::<indexmap::IndexSet<_>>();
                let paths_mapped = id_mappings.file_names().iter().cloned().collect::<indexmap::IndexSet<_>>();

                if paths_in != paths_mapped {
                    return Err(anyhow::anyhow!(
                        "input files and ID mappings do not match: {:?} != {:?}",
                        paths_in,
                        paths_mapped
                    ));
                }

                Ok(id_mappings)
            },
        )
        .transpose()?;

    tracing::info!("guessing SV callers...");
    let input_sv_callers = {
        let mut sv_callers = Vec::new();
        for mut reader in open_vcf_readers(&args.path_in).await? {
            sv_callers.push(guess_sv_caller(&mut reader).await?);
        }
        sv_callers
    };

    tracing::info!("processing header...");
    let input_headers = join_all(
        input_readers
            .iter_mut()
            .map(|input_reader| input_reader.read_header()),
    )
    .await
    .into_iter()
    .collect::<Result<Vec<_>, _>>()
    .map_err(|e| anyhow::anyhow!("problem reading header: {}", e))?;
    let orig_sample_names = input_headers
        .first()
        .expect("must have at least one input file")
        .sample_names();
    let sample_names = orig_sample_names
        .iter()
        .map(|name| {
            if let Some(id_mappings) = &id_mappings {
                let mapping = id_mappings
                    .mapping_for_file(args.path_in.first().expect("count checked above"))
                    .expect("checked above");
                mapping
                    .get(name)
                    .cloned()
                    .ok_or_else(|| anyhow::anyhow!("no mapping for sample name: {}", name))
            } else {
                Ok(name.clone())
            }
        })
        .collect::<Result<indexmap::IndexSet<_>, _>>()?;
    for (indexno, other_input_header) in input_headers.iter().enumerate().skip(1) {
        let other_sample_names = if let Some(id_mappings) = &id_mappings {
            let mapping = id_mappings
                .mapping_for_file(&args.path_in[indexno])
                .expect("checked above");
            other_input_header
                .sample_names()
                .iter()
                .map(|s| {
                    mapping.get(s).cloned().ok_or_else(|| {
                        anyhow::anyhow!(
                            "no mapping for sample name: {} in file: {}",
                            s,
                            &args.path_in[indexno]
                        )
                    })
                })
                .collect::<Result<indexmap::IndexSet<_>, _>>()?
        } else {
            other_input_header
                .sample_names()
                .iter()
                .cloned()
                .collect::<indexmap::IndexSet<_>>()
        };
        if other_sample_names != sample_names {
            return Err(anyhow::anyhow!(
                "input file #{} has different sample names than first one: {}",
                indexno,
                &args.path_in[indexno]
            ));
        }
    }
    let output_header = header::build_output_header(
        orig_sample_names,
        &input_sv_callers.iter().collect::<Vec<_>>(),
        id_mappings.as_ref().map(|id_mappings| {
            id_mappings
                .mapping_for_file(args.path_in.first().expect("count checked above"))
                .expect("checked above")
        }),
        Some(&pedigree),
        args.genomebuild,
        &args.file_date,
        worker_version(),
        &args.case_uuid,
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    // Use output file helper.
    let out_path_helper = crate::common::s3::OutputPathHelper::new(&args.path_out)?;

    {
        // Map sample names in input headers.
        let mapped_input_headers = if let Some(id_mappings) = &id_mappings {
            let mut mapped_input_headers = Vec::new();
            for (i, input_header) in input_headers.iter().enumerate() {
                let mapping = id_mappings
                    .mapping_for_file(&args.path_in[i])
                    .expect("checked above");
                let mut input_header = input_header.clone();
                let orig_sample_names = input_header.sample_names().clone();
                input_header.sample_names_mut().clear();
                for sample_name in orig_sample_names {
                    let mapped_sample_name =
                        mapping.get(&sample_name).cloned().ok_or_else(|| {
                            anyhow::anyhow!("no mapping for sample name: {}", sample_name)
                        })?;
                    input_header.sample_names_mut().insert(mapped_sample_name);
                }
                mapped_input_headers.push(input_header);
            }
            mapped_input_headers
        } else {
            input_headers.clone()
        };

        // Perform actual writing
        let mut output_writer = open_vcf_writer(out_path_helper.path_out()).await?;
        output_writer
            .write_header(&output_header)
            .await
            .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

        process_variants(
            &pedigree,
            &output_header,
            &mut output_writer,
            input_readers,
            &mapped_input_headers,
            &input_sv_callers,
            args,
        )
        .await?;

        flush_and_shutdown!(output_writer);
    }

    out_path_helper.create_tbi_for_bgzf().await?;
    out_path_helper.upload_for_s3().await?;

    tracing::info!(
        "All of `strucvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    use crate::common::GenomeRelease;

    #[tracing_test::traced_test]
    #[tokio::test]
    async fn smoke_test_trio() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default().permanent();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_in: vec![
                String::from("tests/strucvars/ingest/delly2-min.vcf"),
                String::from("tests/strucvars/ingest/popdel-min.vcf"),
            ],
            path_cov_vcf: vec![],
            path_ped: "tests/strucvars/ingest/delly2-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
            rng_seed: Some(42),
            file_date: String::from("20230421"),
            case_uuid: String::from("d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c"),
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
    #[tracing_test::traced_test]
    #[tokio::test]
    async fn smoke_test_singleton() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_in: vec![
                String::from("tests/strucvars/ingest/dragen-cnv-min.vcf"),
                String::from("tests/strucvars/ingest/dragen-sv-min.vcf"),
                String::from("tests/strucvars/ingest/gcnv-min.vcf"),
                String::from("tests/strucvars/ingest/manta-min.vcf"),
                String::from("tests/strucvars/ingest/melt-min.vcf"),
                String::from("tests/strucvars/ingest/sniffles2-min.vcf"),
            ],
            path_cov_vcf: vec![],
            path_ped: "tests/strucvars/ingest/dragen-cnv-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
            rng_seed: Some(42),
            file_date: String::from("20230421"),
            case_uuid: String::from("d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c"),
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[tokio::test]
    async fn smoke_test_trio_gz() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_in: vec![
                String::from("tests/strucvars/ingest/delly2-min.vcf.gz"),
                String::from("tests/strucvars/ingest/popdel-min.vcf.gz"),
            ],
            path_cov_vcf: vec![],
            path_ped: "tests/strucvars/ingest/delly2-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf.gz")
                .to_str()
                .expect("invalid path")
                .into(),
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
            rng_seed: Some(42),
            file_date: String::from("20230421"),
            case_uuid: String::from("d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c"),
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;

        Ok(())
    }

    #[tracing_test::traced_test]
    #[tokio::test]
    async fn smoke_test_singleton_gz() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_in: vec![
                String::from("tests/strucvars/ingest/dragen-cnv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/dragen-sv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/gcnv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/manta-min.vcf.gz"),
                String::from("tests/strucvars/ingest/melt-min.vcf.gz"),
                String::from("tests/strucvars/ingest/sniffles2-min.vcf.gz"),
            ],
            path_cov_vcf: vec![],
            path_ped: "tests/strucvars/ingest/dragen-cnv-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf.gz")
                .to_str()
                .expect("invalid path")
                .into(),
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
            rng_seed: Some(42),
            file_date: String::from("20230421"),
            case_uuid: String::from("d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c"),
            id_mapping: None,
        };
        super::run(&args_common, &args).await?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;

        Ok(())
    }

    #[tracing_test::traced_test]
    #[tokio::test]
    async fn smoke_test_singleton_with_id_mapping() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            max_var_count: None,
            path_in: vec![
                String::from("tests/strucvars/ingest/dragen-cnv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/dragen-sv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/gcnv-min.vcf.gz"),
                String::from("tests/strucvars/ingest/manta-min.vcf.gz"),
                String::from("tests/strucvars/ingest/melt-min.vcf.gz"),
                String::from("tests/strucvars/ingest/sniffles2-min.vcf.gz"),
            ],
            path_cov_vcf: vec![],
            path_ped: "tests/strucvars/ingest/dragen-cnv-min.custom_id.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
            rng_seed: Some(42),
            file_date: String::from("20230421"),
            case_uuid: String::from("d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c"),
            id_mapping: Some(
                r#"
                {
                    "mappings": [
                        {
                            "path": "tests/strucvars/ingest/dragen-cnv-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/strucvars/ingest/dragen-sv-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/strucvars/ingest/gcnv-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/strucvars/ingest/manta-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/strucvars/ingest/melt-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        },
                        {
                            "path": "tests/strucvars/ingest/sniffles2-min.vcf.gz",
                            "entries": [
                                {
                                    "src": "SAMPLE",
                                    "dst": "my-custom-id"
                                }
                            ]
                        }
                    ]
                }
                "#
                .into(),
            ),
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
}
