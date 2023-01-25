pub mod common;
pub mod err;
pub mod sv_build_bgdb;
pub mod sv_query;

use clap::{Args, Parser, Subcommand};
use console::{Emoji, Term};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Varfish heavy lifting",
    long_about = "This tool performs the heavy lifting for varfish-server"
)]
struct Cli {
    /// Commonly used arguments
    #[command(flatten)]
    common: common::Args,

    /// The sub command to run
    #[command(subcommand)]
    command: Commands,
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum Commands {
    /// SV related commands
    Sv(Sv),
}

#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Sv {
    /// The sub command to run
    #[command(subcommand)]
    command: SvCommands,
}

#[derive(Debug, Subcommand)]
enum SvCommands {
    Query(sv_query::Args),
    BuildBgDb(sv_build_bgdb::Args),
}

fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    let term = Term::stderr();
    match &cli.command {
        Commands::Sv(sv) => match &sv.command {
            SvCommands::Query(args) => {
                sv_query::run(&term, &cli.common, args)?;
            }
            SvCommands::BuildBgDb(args) => {
                sv_build_bgdb::run(&term, &cli.common, args)?;
            }
        },
    }
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
