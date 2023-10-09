//! Implementation of `seqvars prefilter` subcommand.

use crate::common;

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

/// Main entry point for `seqvars prefilter` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    common::trace_rss_now();


    tracing::info!(
        "All of `seqvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}
