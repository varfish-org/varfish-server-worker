mod common;
mod err;
mod sv_query;

use clap::{Parser, Subcommand};
use console::{Emoji, Term};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Commonly used arguments
    #[command(flatten)]
    common: common::Args,

    /// The sub command to run
    #[command(subcommand)]
    command: Commands,
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug)]
enum Commands {
    /// Create contigs with synthetic sequence
    SvQuery(sv_query::Args),
}

fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    let term = Term::stderr();
    match &cli.command {
        Commands::SvQuery(args) => {
            sv_query::run(&term, &cli.common, args)?;
        }
    }
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
