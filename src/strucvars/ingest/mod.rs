//! Implementation of `strucvars ingest` subcommand.

use std::io::{BufRead, Write};

use crate::common::{self, open_write_maybe_gz, worker_version, GenomeRelease};
use mehari::{annotate::seqvars::AnnotatedVcfWriter, common::open_read_maybe_gz};
use rand_core::SeedableRng;

use mehari::annotate::strucvars::guess_sv_caller;
use noodles_vcf as vcf;

pub mod header;

/// Command line arguments for `strucvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest structural variant VCF", long_about = None)]
pub struct Args {
    /// Maximal number of variants to write out; optional.
    #[clap(long)]
    pub max_var_count: Option<usize>,
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
    /// Value to write to `##fileDate`.
    #[arg(long)]
    pub file_date: String,
}

/// Wrapper around noodle's VCF writer that adjusts the record for the worker.
pub struct WriterWrapper {
    inner: vcf::Writer<Box<dyn Write>>,
}

impl WriterWrapper {
    pub fn new(inner: vcf::Writer<Box<dyn Write>>) -> Self {
        Self { inner }
    }
}

impl mehari::annotate::seqvars::AnnotatedVcfWriter for WriterWrapper {
    fn write_header(&mut self, header: &vcf::Header) -> Result<(), anyhow::Error> {
        self.inner
            .write_header(header)
            .map_err(|e| anyhow::anyhow!("Error writing VCF header: {}", e))
    }

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::Record,
    ) -> Result<(), anyhow::Error> {
        eprintln!("foo");
        self.inner
            .write_record(header, record)
            .map_err(|e| anyhow::anyhow!("Error writing VCF record: {}", e))
    }
}

/// Write out variants from input files.
fn process_variants(
    pedigree: &mehari::ped::PedigreeByName,
    output_writer: &mut dyn mehari::annotate::seqvars::AnnotatedVcfWriter,
    input_readers: &mut [vcf::Reader<Box<dyn BufRead>>],
    output_header: &vcf::Header,
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
    let tmp_dir = tempdir::TempDir::new("mehari")?;

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
        )?;
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
            output_writer.write_record(output_header, &record.try_into()?)?;
        }
    }
    tracing::info!("... done clustering SVs to output");

    Ok(())
}

/// Main entry point for `strucvars ingest` sub command.
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
    let mut input_readers = args
        .path_in
        .iter()
        .map(|path_in| {
            vcf::reader::Builder
                .build_from_reader(open_read_maybe_gz(path_in)?)
                .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))
        })
        .collect::<Result<Vec<_>, _>>()?;

    tracing::info!("guessing SV callers...");
    let input_sv_callers = args
        .path_in
        .iter()
        .map(|path| {
            let reader = open_read_maybe_gz(path)?;
            guess_sv_caller(reader)
        })
        .collect::<Result<Vec<_>, _>>()?;

    tracing::info!("processing header...");
    let input_headers = input_readers
        .iter_mut()
        .map(|input_reader| {
            input_reader
                .read_header()
                .map_err(|e| anyhow::anyhow!("problem reading VCF header: {}", e))
        })
        .collect::<Result<Vec<_>, _>>()?;
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
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    let mut output_writer = WriterWrapper::new(vcf::writer::Writer::new(open_write_maybe_gz(
        &args.path_out,
    )?));
    output_writer
        .write_header(&output_header)
        .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

    process_variants(
        &pedigree,
        &mut output_writer,
        &mut input_readers,
        &output_header,
        &input_headers,
        &input_sv_callers,
        args,
    )?;

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
    #[test]
    fn smoke_test_trio() -> Result<(), anyhow::Error> {
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
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
    #[tracing_test::traced_test]
    #[test]
    fn smoke_test_singleton() -> Result<(), anyhow::Error> {
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
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn smoke_test_trio_gz() -> Result<(), anyhow::Error> {
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
        };
        super::run(&args_common, &args)?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;

        Ok(())
    }

    #[tracing_test::traced_test]
    #[test]
    fn smoke_test_singleton_gz() -> Result<(), anyhow::Error> {
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
        };
        super::run(&args_common, &args)?;

        let mut buffer: Vec<u8> = Vec::new();
        hxdmp::hexdump(&crate::common::read_to_bytes(&args.path_out)?, &mut buffer)?;

        Ok(())
    }
}
