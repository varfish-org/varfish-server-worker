//! Implementation of `seqvars prefilter` subcommand.

use std::io::BufRead;

use noodles_vcf as vcf;
use thousands::Separable;

use crate::common::{self, open_read_maybe_gz};

/// Arguments for the `seqvars prefilter` subcommand.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
struct PrefilterParams {
    /// Path to output file.
    pub path_out: String,
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
        if param.starts_with("@") {
            let path = param.trim_start_matches("@");
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
            let params: PrefilterParams = serde_json::from_str(&param)
                .map_err(|e| anyhow::anyhow!("failed to parse prefilter params: {}", e))?;
            result.push(params);
        }
    }

    Ok(result)
}

/// Extract largest population frequency and exon distance from input_record.
///
/// Note that all variants on chrMT will be returned.
fn get_freq_and_distance(input_record: &vcf::Record) -> Result<(f64, i32), anyhow::Error> {
    let value = input_record.info().get("gnomad_exomes_an".parse()?);
}

// ##INFO=<ID=gnomad_exomes_an,Number=1,Type=Integer,Description="Number of samples in gnomAD exomes">
// ##INFO=<ID=gnomad_exomes_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in gnomAD exomes">
// ##INFO=<ID=gnomad_exomes_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in gnomAD exomes">
// ##INFO=<ID=gnomad_exomes_hemi,Number=1,Type=Integer,Description="Number of hemi. alt. carriers in gnomAD exomes">
// ##INFO=<ID=gnomad_genomes_an,Number=1,Type=Integer,Description="Number of samples in gnomAD genomes">
// ##INFO=<ID=gnomad_genomes_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in gnomAD genomes">
// ##INFO=<ID=gnomad_genomes_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in gnomAD genomes">
// ##INFO=<ID=gnomad_genomes_hemi,Number=1,Type=Integer,Description="Number of hemi. alt. carriers in gnomAD genomes">
// ##INFO=<ID=helix_an,Number=1,Type=Integer,Description="Number of samples in HelixMtDb">
// ##INFO=<ID=helix_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in HelixMtDb">
// ##INFO=<ID=helix_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in HelixMtDb">


/// Perform the actual prefiltration.
fn run_filtration(
    input_reader: &mut vcf::Reader<Box<dyn std::io::BufRead>>,
    input_header: &vcf::Header,
    output_writers: &mut [vcf::Writer<Box<dyn std::io::Write>>],
    params: &[PrefilterParams],
) -> Result<(), anyhow::Error> {
    let start = std::time::Instant::now();
    let mut prev = std::time::Instant::now();
    let mut records = input_reader.records(input_header);
    let mut total_written = 0usize;
    loop {
        if let Some(input_record) = records.next() {
            let input_record = input_record?;

            let (frequency, exon_distance) = get_freq_and_distance(&input_record)?;
            for (writer_params, output_writer) in params.iter().zip(output_writers.iter_mut()) {
                if frequency <= writer_params.max_freq && exon_distance <= writer_params.max_exon_dist {
                    output_writer
                        .write_record(&input_header, &input_record)
                        .map_err(|e| anyhow::anyhow!("failed to write record: {}", e))?;
                }
            }

            let vcf_var = annonars::common::keys::Var::from_vcf_allele(&input_record, 0);
            if prev.elapsed().as_secs() >= 60 {
                tracing::info!("at {:?}", &vcf_var);
                prev = std::time::Instant::now();
            }

            total_written += 1;
        } else {
            break; // all done
        }
    }
    tracing::info!(
        "... annotated {} records in {:?}",
        total_written.separate_with_commas(),
        start.elapsed()
    );

    Ok(())
}

/// Main entry point for `seqvars prefilter` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    tracing::info!("loading prefilter params...");
    let params = load_params(&args.params)?;
    tracing::info!("opening input file...");
    let reader = open_read_maybe_gz(&args.path_in)
        .map_err(|e| anyhow::anyhow!("could not open input file: {}", e))?;
    let mut reader = vcf::Reader::new(reader);
    let header = reader
        .read_header()
        .map_err(|e| anyhow::anyhow!("problem reading header: {}", e))?;

    tracing::info!("opening output files...");
    let mut output_files = params
        .iter()
        .map(|params| {
            let mut header = header.clone();
            header.insert(
                "x-varfish-prefilter-params"
                    .parse()
                    .map_err(|e| anyhow::anyhow!("{}", e))?,
                vcf::header::record::Value::from(""),
            )?;

            let mut writer =
                vcf::Writer::new(common::open_write_maybe_gz(&params.path_out).map_err(|e| {
                    anyhow::anyhow!("could not open output file {}: {}", &params.path_out, e)
                })?);
            writer.write_header(&header).map_err(|e| {
                anyhow::anyhow!("could not write header to {}: {}", &params.path_out, e)
            })?;

            Ok(writer)
        })
        .collect::<Result<Vec<_>, anyhow::Error>>()?;

    common::trace_rss_now();

    tracing::info!("starting filtration...");
    run_filtration(&mut reader, &header, &mut output_files, &params)?;
    tracing::info!("... done with filtration");

    tracing::info!(
        "All of `seqvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    #[test]
    fn single_output_arg() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args = super::Args {
            path_in: "tests/data/seqvars/ingest/ingested.vcf.gz".into(),
            params: vec![format!(
                r#"{{
                    "path_out": "{}/out-1.vcf.gz",
                    "max_freq": 0.01,
                    "max_exon_dist": 200
                }}"#,
                tmpdir.to_path_buf().to_str().unwrap()
            )],
        };

        super::run(&crate::common::Args::default(), &args)?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }

    #[test]
    fn single_output_file() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let params_json = format!(
            r#"{{
                "path_out": "{}/out-1.vcf.gz",
                "max_freq": 0.01,
                "max_exon_dist": 200
            }}"#,
            tmpdir.to_path_buf().to_str().unwrap()
        );

        let params_file = tmpdir.to_path_buf().join("params.json");
        std::fs::write(&params_file, &params_json)?;

        let args = super::Args {
            path_in: "tests/data/seqvars/ingest/ingested.vcf.gz".into(),
            params: vec![format!("@{}", params_file.to_str().unwrap())],
        };

        super::run(&crate::common::Args::default(), &args)?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }

    #[test]
    fn two_output_arg() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let args = super::Args {
            path_in: "tests/data/seqvars/ingest/ingested.vcf.gz".into(),
            params: vec![
                format!(
                    r#"{{
                        "path_out": "{}/out-1.vcf.gz",
                        "max_freq": 0.01,
                        "max_exon_dist": 200
                    }}"#,
                    tmpdir.to_path_buf().to_str().unwrap()
                ),
                format!(
                    r#"{{
                        "path_out": "{}/out-2.vcf.gz",
                        "max_freq": 0.005,
                        "max_exon_dist": 20
                    }}"#,
                    tmpdir.to_path_buf().to_str().unwrap()
                ),
            ],
        };

        super::run(&crate::common::Args::default(), &args)?;

        assert!(std::path::Path::new(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());
        assert!(std::path::Path::new(&format!(
            "{}/out-2.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))
        .exists());

        insta::assert_snapshot!(std::fs::read_to_string(&format!(
            "{}/out-1.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);
        insta::assert_snapshot!(std::fs::read_to_string(&format!(
            "{}/out-2.vcf.gz",
            tmpdir.to_path_buf().to_str().unwrap()
        ))?);

        Ok(())
    }
}
