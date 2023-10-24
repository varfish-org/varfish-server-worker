//! Implementation of `strucvars ingest` subcommand.

use crate::common::{self, worker_version, GenomeRelease};
use crate::flush_and_shutdown;
use futures::future::join_all;
use mehari::annotate::strucvars::guess_sv_caller;
use mehari::common::io::std::is_gz;
use mehari::common::noodles::{open_vcf_readers, open_vcf_writer, AsyncVcfReader, AsyncVcfWriter};
use noodles_vcf as vcf;
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
}

async fn write_ingest_record(
    writer: &mut AsyncVcfWriter,
    input_record: &vcf::Record,
) -> Result<(), anyhow::Error> {
    // copy over CHROM, POS, REF
    let mut builder = vcf::Record::builder()
        .set_chromosome(input_record.chromosome().clone())
        .set_position(input_record.position())
        .set_reference_bases(input_record.reference_bases().clone());

    // copy over first ALT allele, remove any SV sub types
    if input_record.alternate_bases().len() != 1 {
        anyhow::bail!(
            "unexpected number of ALT alleles (should be ==1) in: {:?}",
            input_record.alternate_bases()
        );
    }
    let alt_0 = &input_record.alternate_bases()[0];
    let sv_type;
    let bnd;
    match alt_0 {
        vcf::record::alternate_bases::Allele::Breakend(bnd_string) => {
            sv_type = "BND".parse()?;
            builder =
                builder.set_alternate_bases(vcf::record::AlternateBases::from(vec![alt_0.clone()]));
            bnd = Some(
                mehari::annotate::strucvars::bnd::Breakend::from_ref_alt_str(
                    &format!("{}", input_record.reference_bases()),
                    bnd_string,
                )?,
            );
        }
        vcf::record::alternate_bases::Allele::Symbol(symbol) => match symbol {
            vcf::record::alternate_bases::allele::Symbol::StructuralVariant(sv) => {
                sv_type = sv.ty();
                builder = builder
                .set_alternate_bases(vcf::record::AlternateBases::from(
                    vec![vcf::record::alternate_bases::Allele::Symbol(
                        vcf::record::alternate_bases::allele::Symbol::StructuralVariant(
                            vcf::record::alternate_bases::allele::symbol::structural_variant::StructuralVariant::from(
                                sv.ty()
                            )
                        )
                    )]
                ));
                bnd = None;
            }
            _ => anyhow::bail!("unexpected symbolic allele: {:?}", &symbol),
        },
        _ => anyhow::bail!("unexpected alternate base type: {:?}", &alt_0),
    }

    // copy over FORMAT tags, all except FT
    let mut keys_with_value = std::collections::HashSet::<String>::new();
    let output_format_values = input_record
        .genotypes()
        .values()
        .map(|g| {
            g.keys()
                .iter()
                .zip(g.values().iter())
                .filter(|(k, _)| k.as_ref() != "FT")
                .map(|(k, v)| {
                    if v.is_some() {
                        keys_with_value.insert(k.as_ref().to_string());
                    }

                    v.clone()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let output_keys = vcf::record::genotypes::Keys::try_from(
        input_record
            .genotypes()
            .keys()
            .iter()
            .filter(|k| k.as_ref() != "FT")
            .cloned()
            .map(|k| {
                if k.as_ref() == "CN" {
                    "cn".parse().expect("invalid key: cn")
                } else {
                    k
                }
            })
            .collect::<Vec<_>>(),
    )?;
    builder = builder.set_genotypes(vcf::record::Genotypes::new(
        output_keys,
        output_format_values,
    ));

    // copy over INFO tags
    // Note: annsv will be added only in "strucvars query"
    let mut info: noodles_vcf::record::Info = Default::default();
    match sv_type {
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Deletion |
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Duplication |
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::CopyNumberVariation => {
            let claim = if keys_with_value.contains("pev") || keys_with_value.contains("srv") {
                "DJ"
            } else {
                "D"
            };
            info.insert(
                vcf::record::info::field::key::SV_CLAIM,
                Some(vcf::record::info::field::Value::Array(
                    vcf::record::info::field::value::Array::String(vec![Some(claim.to_string())]),
                )
            ));

        }
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Insertion |
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Inversion |
        vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Breakend => {
            info.insert(
                vcf::record::info::field::key::SV_CLAIM,
                Some(vcf::record::info::field::Value::Array(
                    vcf::record::info::field::value::Array::String(vec![Some("J".to_string())]),
                )
            ));
        },
    }
    info.insert(
        vcf::record::info::field::key::SV_TYPE,
        Some(vcf::record::info::field::Value::String(sv_type.to_string())),
    );
    if let Some(Some(vcf::record::info::field::Value::Integer(end))) = input_record
        .info()
        .get(&vcf::record::info::field::key::END_POSITION)
    {
        info.insert(
            vcf::record::info::field::key::END_POSITION,
            Some(vcf::record::info::field::Value::Integer(*end)),
        );

        if sv_type
            == vcf::record::alternate_bases::allele::symbol::structural_variant::Type::Breakend
        {
            info.insert(
                "chr2".parse()?,
                Some(vcf::record::info::field::value::Value::String(
                    bnd.expect("must be set here").chrom.clone(),
                )),
            );
        } else {
            let pos: usize = input_record.position().into();
            let sv_len: usize = *end as usize - pos + 1;
            info.insert(
                vcf::record::info::field::key::SV_LENGTHS,
                Some(vcf::record::info::field::Value::Array(
                    vcf::record::info::field::value::Array::Integer(vec![Some(sv_len as i32)]),
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
        } else {
            anyhow::bail!("unknown caller: {}", caller)
        }
    }

    let key_callers: vcf::record::info::field::Key = "callers".parse()?;
    if let Some(Some(callers)) = input_record.info().get(&key_callers) {
        if let vcf::record::info::field::Value::Array(
            vcf::record::info::field::value::Array::String(callers),
        ) = callers
        {
            let output_callers = callers
                .iter()
                .flatten()
                .map(|caller| map_caller(caller))
                .collect::<Result<Vec<_>, _>>()?;
            info.insert(
                key_callers,
                Some(vcf::record::info::field::Value::Array(
                    vcf::record::info::field::value::Array::String(output_callers),
                )),
            );
        } else if let vcf::record::info::field::Value::String(caller) = callers {
            let output_callers = vec![map_caller(caller)?];
            info.insert(
                key_callers,
                Some(vcf::record::info::field::Value::Array(
                    vcf::record::info::field::value::Array::String(output_callers),
                )),
            );
        }
    } else {
        anyhow::bail!("no callers INFO tag found");
    }

    builder = builder.set_info(info);

    let record = builder.build()?;

    writer
        .write_record(&record)
        .await
        .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
}

/// Write out variants from input files.
async fn process_variants(
    pedigree: &mehari::ped::PedigreeByName,
    output_writer: &mut AsyncVcfWriter,
    input_readers: &mut [AsyncVcfReader],
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
    for (reader, sv_caller, header) in itertools::izip!(
        input_readers.iter_mut(),
        input_sv_callers.iter(),
        input_header.iter()
    ) {
        mehari::annotate::strucvars::run_vcf_to_jsonl(
            pedigree,
            reader,
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
            write_ingest_record(output_writer, &record.try_into()?).await?;
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
    let sample_names = input_headers
        .first()
        .expect("must have at least one input file")
        .sample_names();
    for (indexno, other_input_header) in input_headers.iter().enumerate().skip(1) {
        if other_input_header.sample_names() != sample_names {
            return Err(anyhow::anyhow!(
                "input file #{} has different sample names than first one: {}",
                indexno,
                &args.path_in[indexno]
            ));
        }
    }
    let output_header = header::build_output_header(
        sample_names,
        &input_sv_callers.iter().collect::<Vec<_>>(),
        Some(&pedigree),
        args.genomebuild,
        &args.file_date,
        worker_version(),
        &args.case_uuid,
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    {
        let mut output_writer = open_vcf_writer(&args.path_out).await?;
        output_writer
            .write_header(&output_header)
            .await
            .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

        process_variants(
            &pedigree,
            &mut output_writer,
            &mut input_readers,
            &input_headers,
            &input_sv_callers,
            args,
        )
        .await?;

        flush_and_shutdown!(output_writer);
    }

    if is_gz(&args.path_out) {
        tracing::info!("Creating TBI index for BGZF VCF file...");
        crate::common::noodles::build_tbi(&args.path_out, &format!("{}.tbi", &args.path_out))
            .await
            .map_err(|e| anyhow::anyhow!("problem building TBI: {}", e))?;
        tracing::info!("... done writing TBI index");
    } else {
        tracing::info!("(not building TBI index for plain text VCF file");
    }

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
        let tmpdir = temp_testdir::TempDir::default();

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
        };
        super::run(&args_common, &args).await?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;

        Ok(())
    }
}
