use clap::{command, Parser};
use tracing::{info, instrument};

#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for SVs", long_about = None)]
pub struct Args {}

#[instrument]
pub(crate) fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);
    Ok(())
}
