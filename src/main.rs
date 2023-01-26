//! VarFish Server Worker main executable

pub mod common;
pub mod sv;

#[allow(
    non_snake_case,
    unused_imports,
    clippy::extra_unused_lifetimes,
    clippy::missing_safety_doc,
    clippy::derivable_impls,
    clippy::size_of_in_element_count
)]
#[path = "../target/flatbuffers/world_generated.rs"]
pub mod world_flatbuffers;

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
    /// SV related commands
    Sv(Sv),
}

/// Parsing of top-level CLI commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Sv {
    /// The sub command to run
    #[command(subcommand)]
    command: SvCommands,
}

/// Enum supporting the parsing of "sv *" sub commands.
#[derive(Debug, Subcommand)]
enum SvCommands {
    BgDbToBin(sv::bg_db_to_bin::Args),
    InhouseDbBuild(sv::inhouse_db_build::Args),
    Query(sv::query::Args),
}

fn main() -> Result<(), anyhow::Error> {
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

    // Install collector and go into sub commands.
    let term = Term::stderr();
    tracing::subscriber::with_default(collector, || {
        match &cli.command {
            Commands::Sv(sv) => match &sv.command {
                SvCommands::BgDbToBin(args) => {
                    sv::bg_db_to_bin::run(&cli.common, args)?;
                }
                SvCommands::InhouseDbBuild(args) => {
                    sv::inhouse_db_build::run(&cli.common, args)?;
                }
                SvCommands::Query(args) => {
                    sv::query::run(&cli.common, args)?;
                }
            },
        }

        Ok::<(), anyhow::Error>(())
    })?;
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
