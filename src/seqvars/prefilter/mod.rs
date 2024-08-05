//! Implementation of `seqvars prefilter` subcommand.

use std::io::BufRead;

use mehari::annotate::seqvars::ann::AnnField;
use mehari::common::noodles::{open_vcf_reader, open_vcf_writer, AsyncVcfReader, AsyncVcfWriter};
use noodles::vcf;
use thousands::Separable;
use tokio::io::AsyncWriteExt;

use crate::{common, flush_and_shutdown};

/// Arguments for the `seqvars prefilter` subcommand.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
struct PrefilterParams {
    /// Path to output file.
    pub prefilter_path: String,
    /// Maximal allele population frequency.
    pub max_freq: f64,
    /// Maximal distance to exon.
    pub max_exon_dist: i32,
}

/// Command line arguments for `seqvars prefilter` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "prefilter an ingested variant VCF", long_about = None)]
pub struct Args {
    /// Path to input file.
    #[clap(long)]
    pub path_in: String,
    /// Prefilter parameters or @ with path to JSONL file.
    #[clap(long)]
    pub params: Vec<String>,
}

/// Load prefilter params from strings or files with such strings.
fn load_params(params: &[String]) -> Result<Vec<PrefilterParams>, anyhow::Error> {
    let mut result = Vec::new();

    for param in params {
        if param.starts_with('@') {
            let path = param.trim_start_matches('@');
            let file = std::fs::File::open(path)
                .map_err(|e| anyhow::anyhow!("failed to open prefilter params file: {}", e))?;
            let reader = std::io::BufReader::new(file);
            for (lineno, line) in reader.lines().enumerate() {
                let params: PrefilterParams = serde_json::from_str(&line?).map_err(|e| {
                    anyhow::anyhow!(
                        "failed to parse prefilter params from line {}: {}",
                        lineno,
                        e
                    )
                })?;
                result.push(params);
            }
        } else {
            let params: PrefilterParams = serde_json::from_str(param)
                .map_err(|e| anyhow::anyhow!("failed to parse prefilter params: {}", e))?;
            result.push(params);
        }
    }

    Ok(result)
}

/// Extract an `i32` from a VCF record's `INFO`.
fn get_info_i32(input_record: &vcf::variant::RecordBuf, key: &str) -> Result<i32, anyhow::Error> {
    if let Some(Some(vcf::variant::record_buf::info::field::Value::Integer(gnomad_exomes_an))) =
        input_record.info().get(key)
    {
        Ok(*gnomad_exomes_an)
    } else {
        Ok(0)
    }
}

/// Extract largest population frequency and exon distance from input_record.
///
/// Note that all variants on chrMT will be returned.
fn get_freq_and_distance(
    input_record: &vcf::variant::RecordBuf,
) -> Result<(f64, Option<i32>), anyhow::Error> {
    if annonars::common::cli::canonicalize(input_record.reference_sequence_name()) == "MT" {
        return Ok((0.0, Some(0))); // all variants on chrMT are returned
    }

    let gnomad_exomes_an = get_info_i32(input_record, "gnomad_exomes_an")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_exomes_an: {}", e))?;
    let gnomad_exomes_hom = get_info_i32(input_record, "gnomad_exomes_hom")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_exomes_hom: {}", e))?;
    let gnomad_exomes_het = get_info_i32(input_record, "gnomad_exomes_het")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_exomes_het: {}", e))?;
    let gnomad_exomes_hemi = get_info_i32(input_record, "gnomad_exomes_hemi")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_exomes_hemi: {}", e))?;
    let gnomad_exomes_an_pos = gnomad_exomes_hom * 2 + gnomad_exomes_het + gnomad_exomes_hemi;
    let gnomad_exomes_freq = if gnomad_exomes_an > 0 {
        gnomad_exomes_an_pos as f64 / gnomad_exomes_an as f64
    } else {
        0f64
    };

    let gnomad_genomes_an = get_info_i32(input_record, "gnomad_genomes_an")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_genomes_an: {}", e))?;
    let gnomad_genomes_hom = get_info_i32(input_record, "gnomad_genomes_hom")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_genomes_hom: {}", e))?;
    let gnomad_genomes_het = get_info_i32(input_record, "gnomad_genomes_het")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_genomes_het: {}", e))?;
    let gnomad_genomes_hemi = get_info_i32(input_record, "gnomad_genomes_hemi")
        .map_err(|e| anyhow::anyhow!("failed to get gnomad_genomes_hemi: {}", e))?;
    let gnomad_genomes_an_pos = gnomad_genomes_hom * 2 + gnomad_genomes_het + gnomad_genomes_hemi;
    let gnomad_genomes_freq = if gnomad_genomes_an > 0 {
        gnomad_genomes_an_pos as f64 / gnomad_genomes_an as f64
    } else {
        0f64
    };

    let exon_dist = if let Some(Some(ann)) = input_record.info().get("ANN") {
        if let vcf::variant::record_buf::info::field::Value::Array(
            vcf::variant::record_buf::info::field::value::Array::String(anns),
        ) = ann
        {
            let ann = anns
                .first()
                .expect("ANN field is empty")
                .as_ref()
                .expect("ANN field entry is empty");
            let record: AnnField = ann
                .parse()
                .map_err(|e| anyhow::anyhow!("failed to parse ANN field from {}: {}", ann, e))?;
            record.distance
        } else {
            anyhow::bail!("wrong data type for ANN in VCF record: {:?}", &input_record)
        }
    } else {
        None
    };

    Ok((f64::max(gnomad_exomes_freq, gnomad_genomes_freq), exon_dist))
}

/// Perform the actual prefiltration.
async fn run_filtration(
    input_reader: &mut AsyncVcfReader,
    input_header: &vcf::Header,
    output_writers: &mut [AsyncVcfWriter],
    output_headers: &[vcf::Header],
    params: &[PrefilterParams],
) -> Result<(), anyhow::Error> {
    let start = std::time::Instant::now();
    let mut prev = std::time::Instant::now();
    let mut total_written = 0usize;
    let mut input_record = vcf::variant::RecordBuf::default();
    loop {
        let bytes_read = input_reader
            .read_record_buf(input_header, &mut input_record)
            .await
            .map_err(|e| anyhow::anyhow!("problem reading VCF file: {}", e))?;
        if bytes_read == 0 {
            break; // EOF
        }

        let (frequency, exon_distance) = get_freq_and_distance(&input_record)?;
        if let Some(exon_distance) = exon_distance {
            for ((writer_params, output_writer), output_header) in params
                .iter()
                .zip(output_writers.iter_mut())
                .zip(output_headers.iter())
            {
                if frequency <= writer_params.max_freq
                    && exon_distance <= writer_params.max_exon_dist
                {
                    output_writer
                        .write_variant_record(output_header, &input_record)
                        .await
                        .map_err(|e| anyhow::anyhow!("failed to write record: {}", e))?;
                }
            }
        }

        let vcf_var = annonars::common::keys::Var::from_vcf_allele(&input_record, 0);
        if prev.elapsed().as_secs() >= 60 {
            tracing::info!("at {:?}", &vcf_var);
            prev = std::time::Instant::now();
        }

        total_written += 1;
    }

    tracing::info!(
        "... annotated {} records in {:?}",
        total_written.separate_with_commas(),
        start.elapsed()
    );

    Ok(())
}

/// Main entry point for `seqvars prefilter` sub command.
pub async fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    tracing::info!("loading prefilter params...");
    let params_list = load_params(&args.params)?;
    tracing::info!("opening input file...");
    let mut reader = open_vcf_reader(&args.path_in)
        .await
        .map_err(|e| anyhow::anyhow!("could not open input file: {}", e))?;
    let header = reader
        .read_header()
        .await
        .map_err(|e| anyhow::anyhow!("problem reading header: {}", e))?;

    {
        tracing::info!("opening output files...");
        let mut output_writers = Vec::new();
        let mut out_path_helpers = Vec::new();
        let output_headers = {
            let mut output_headers = Vec::new();
            for params in params_list.iter() {
                let header_params = PrefilterParams {
                    prefilter_path: "<stripped>".into(),
                    ..params.clone()
                };
                let mut header = header.clone();
                header.insert(
                    "x-varfish-prefilter-params"
                        .parse()
                        .map_err(|e| anyhow::anyhow!("{}", e))?,
                    vcf::header::record::Value::from(
                        serde_json::to_string(&header_params)
                            .map_err(|e| anyhow::anyhow!("failed to serialize params: {}", e))?,
                    ),
                )?;

                out_path_helpers.push(crate::common::s3::OutputPathHelper::new(
                    &params.prefilter_path,
                )?);
                let mut writer =
                    open_vcf_writer(out_path_helpers.last().expect("just pushed").path_out())
                        .await?;
                writer.write_header(&header).await.map_err(|e| {
                    anyhow::anyhow!(
                        "could not write header to {}: {}",
                        &params.prefilter_path,
                        e
                    )
                })?;
                output_writers.push(writer);
                output_headers.push(header);
            }
            output_headers
        };

        common::trace_rss_now();

        tracing::info!("starting filtration...");
        run_filtration(
            &mut reader,
            &header,
            &mut output_writers,
            &output_headers,
            &params_list,
        )
        .await?;
        tracing::info!("... done with filtration");

        for output_writer in output_writers.drain(..) {
            flush_and_shutdown!(output_writer);
        }
        for out_path_helper in out_path_helpers.drain(..) {
            out_path_helper.create_tbi_for_bgzf().await?;
            out_path_helper.upload_for_s3().await?;
        }
    };

    tracing::info!(
        "All of `seqvars prefilter` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    #[tokio::test]
    async fn single_output_arg() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args = super::Args {
            path_in: "tests/seqvars/prefilter/ingest.vcf".into(),
            params: vec![format!(
                r#"{{
                    "prefilter_path": "{}/out-1.vcf",
                    "max_freq": 0.01,
                    "max_exon_dist": 200
                }}"#,
                tmpdir.to_path_buf().to_str().unwrap()
            )],
        };

        super::run(&crate::common::Args::default(), &args).await?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }

    #[tokio::test]
    async fn single_output_file() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let params_json = format!(
            r#"{{"prefilter_path": "{}/out-1.vcf", "max_freq": 0.01, "max_exon_dist": 200}}"#,
            tmpdir.to_path_buf().to_str().unwrap()
        );

        let params_file = tmpdir.to_path_buf().join("params.json");
        std::fs::write(&params_file, &params_json)?;

        let args = super::Args {
            path_in: "tests/seqvars/prefilter/ingest.vcf".into(),
            params: vec![format!("@{}", params_file.to_str().unwrap())],
        };

        super::run(&crate::common::Args::default(), &args).await?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }

    #[tokio::test]
    async fn two_output_arg() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args = super::Args {
            path_in: "tests/seqvars/prefilter/ingest.vcf".into(),
            params: vec![
                format!(
                    r#"{{
                        "prefilter_path": "{}/out-1.vcf",
                        "max_freq": 0.01,
                        "max_exon_dist": 200
                    }}"#,
                    tmpdir.to_path_buf().to_str().unwrap()
                ),
                format!(
                    r#"{{
                        "prefilter_path": "{}/out-2.vcf",
                        "max_freq": 0,
                        "max_exon_dist": 20
                    }}"#,
                    tmpdir.to_path_buf().to_str().unwrap()
                ),
            ],
        };

        super::run(&crate::common::Args::default(), &args).await?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());
        assert!(std::path::Path::new(&format!(
            "{}/out-2.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(format!(
            "{}/out-1.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);
        insta::assert_snapshot!(std::fs::read_to_string(format!(
            "{}/out-2.vcf",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }
}
