use clap::Parser;

use crate::common::Args as CommonArgs;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Path to input file
    #[arg(short, long)]
    pub path_input: String,
}

pub(crate) fn run(
    term: &console::Term,
    common: &CommonArgs,
    args: &Args,
) -> Result<(), anyhow::Error> {
    todo!()
}
