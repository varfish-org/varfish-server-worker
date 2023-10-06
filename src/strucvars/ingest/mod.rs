//! Implementation of `strucvars ingest` subcommand.

use crate::common::{self, open_write_maybe_gz, worker_version, GenomeRelease};
use mehari::common::open_read_maybe_gz;

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
    /// Path to output file.
    #[clap(long)]
    pub path_out: String,
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
        .map(guess_sv_caller)
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
        &Some(pedigree),
        args.genomebuild,
        worker_version(),
    )
    .map_err(|e| anyhow::anyhow!("problem building output header: {}", e))?;

    let mut output_writer = { vcf::writer::Writer::new(open_write_maybe_gz(&args.path_out)?) };
    output_writer
        .write_header(&output_header)
        .map_err(|e| anyhow::anyhow!("problem writing header: {}", e))?;

    // process_variants(
    //     &mut output_writer,
    //     &mut input_readers,
    //     &output_header,
    //     &input_header,
    //     args,
    // )?;

    tracing::info!(
        "All of `strucvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    use crate::common::GenomeRelease;

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
            path_ped: "tests/strucvars/ingest/delly2-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
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
            path_ped: "tests/strucvars/ingest/dragen-cnv-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
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
            path_ped: "tests/strucvars/ingest/delly2-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf.gz")
                .to_str()
                .expect("invalid path")
                .into(),
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
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
            path_ped: "tests/strucvars/ingest/dragen-cnv-min.ped".into(),
            genomebuild: GenomeRelease::Grch37,
            path_out: tmpdir
                .join("out.vcf.gz")
                .to_str()
                .expect("invalid path")
                .into(),
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(&args.path_out)?);

        Ok(())
    }
}
