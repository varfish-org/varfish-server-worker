//! Code for converting from background DB TSV file to binary background DB format.

use clap::{command, Parser};
use tracing::info;

/// Command line arguments for `sv build-inhouse-db` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Convert background TSV to binary file", long_about = None)]
pub struct Args {
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub input_tsv: String,
    /// Path to output binary file.
    #[arg(long, required = true)]
    pub output_bin: String,
}

/// Main entry point for the `sv convert-bgdb` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting sv convert-bgdb");
    info!("common_args = {:?}", &common_args);
    info!("args = {:?}", &args);

    Ok(())
}
