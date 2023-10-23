//! VarFish Server Worker main executable

pub mod common;
pub mod seqvars;
pub mod strucvars;

use clap::{Args, Parser, Subcommand};
use console::{Emoji, Term};

/// CLI parser based on clap.
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Varfish Server heavy lifting",
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

/// Enum supporting the parsing of top-level commands.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum Commands {
    /// Structural variant related commands.
    Strucvars(Strucvars),
    /// Sequence variant related commands.
    Seqvars(Seqvars),
}

/// Parsing of "strucvars *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Strucvars {
    /// The sub command to run
    #[command(subcommand)]
    command: StrucvarsCommands,
}

/// Enum supporting the parsing of "strucvars *" sub commands.
#[derive(Debug, Subcommand)]
enum StrucvarsCommands {
    Aggregate(strucvars::aggregate::cli::Args),
    Ingest(strucvars::ingest::Args),
    Query(strucvars::query::Args),
    TxtToBin(strucvars::txt_to_bin::cli::Args),
}

/// Parsing of "seqvars *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Seqvars {
    /// The sub command to run
    #[command(subcommand)]
    command: SeqvarsCommands,
}

/// Enum supporting the parsing of "strucvars *" sub commands.
#[derive(Debug, Subcommand)]
enum SeqvarsCommands {
    Aggregate(seqvars::aggregate::Args),
    Ingest(seqvars::ingest::Args),
    Prefilter(seqvars::prefilter::Args),
    Query(seqvars::query::Args),
}

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    // Build a tracing subscriber according to the configuration in `cli.common`.
    let collector = tracing_subscriber::fmt()
        .with_target(false)
        .with_max_level(match cli.common.verbose.log_level() {
            Some(level) => match level {
                log::Level::Error => tracing::Level::ERROR,
                log::Level::Warn => tracing::Level::WARN,
                log::Level::Info => tracing::Level::INFO,
                log::Level::Debug => tracing::Level::DEBUG,
                log::Level::Trace => tracing::Level::TRACE,
            },
            None => tracing::Level::INFO,
        })
        .compact()
        .finish();
    tracing::subscriber::set_global_default(collector)?;

    // Install collector and go into sub commands.
    let term = Term::stderr();
    match &cli.command {
        Commands::Seqvars(seqvars) => match &seqvars.command {
            SeqvarsCommands::Aggregate(args) => {
                seqvars::aggregate::run(&cli.common, args)?;
            }
            SeqvarsCommands::Ingest(args) => {
                seqvars::ingest::run(&cli.common, args)?;
            }
            SeqvarsCommands::Prefilter(args) => {
                seqvars::prefilter::run(&cli.common, args)?;
            }
            SeqvarsCommands::Query(args) => {
                seqvars::query::run(&cli.common, args)?;
            }
        },
        Commands::Strucvars(strucvars) => match &strucvars.command {
            StrucvarsCommands::Aggregate(args) => {
                strucvars::aggregate::cli::run(&cli.common, args).await?;
            }
            StrucvarsCommands::Ingest(args) => {
                strucvars::ingest::run(&cli.common, args).await?;
            }
            StrucvarsCommands::Query(args) => {
                strucvars::query::run(&cli.common, args).await?;
            }
            StrucvarsCommands::TxtToBin(args) => {
                strucvars::txt_to_bin::cli::run(&cli.common, args)?;
            }
        },
    }
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
