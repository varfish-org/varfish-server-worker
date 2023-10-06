//! Implementation of `strucvars ingest` subcommand.

use crate::common::{self, open_read_maybe_gz, open_write_maybe_gz, worker_version, GenomeRelease};

use mehari::annotate::seqvars::provider::MehariProvider;
use noodles_vcf as vcf;
use thousands::Separable;

pub mod header;

/// Command line arguments for `strucvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest structural variant VCF", long_about = None)]
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
        input_header.sample_names(),
        &vec![],
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
    //     &mut input_reader,
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
mod test {}
