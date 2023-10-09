//! VarFish Server Worker main executable

pub mod common;
pub mod db;
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
    /// Database building related commands.
    Db(Db),
    /// Structural variant related commands.
    Strucvars(Strucvars),
    /// Sequence variant related commands.
    Seqvars(Seqvars),
}

/// Parsing of "db *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Db {
    /// The sub command to run
    #[command(subcommand)]
    command: DbCommands,
}

/// Enum supporting the parsing of "db *" sub commands.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum DbCommands {
    ToBin(db::to_bin::cli::Args),
    MkInhouse(db::mk_inhouse::cli::Args),
}

/// Parsing of "sv *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Strucvars {
    /// The sub command to run
    #[command(subcommand)]
    command: StrucvarsCommands,
}

/// Enum supporting the parsing of "sv *" sub commands.
#[derive(Debug, Subcommand)]
enum StrucvarsCommands {
    Ingest(strucvars::ingest::Args),
    Query(strucvars::query::Args),
}

/// Parsing of "seqvars *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Seqvars {
    /// The sub command to run
    #[command(subcommand)]
    command: SeqvarsCommands,
}

/// Enum supporting the parsing of "sv *" sub commands.
#[derive(Debug, Subcommand)]
enum SeqvarsCommands {
    Ingest(seqvars::ingest::Args),
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
            Commands::Db(db) => match &db.command {
                DbCommands::MkInhouse(args) => {
                    db::mk_inhouse::cli::run(&cli.common, args)?;
                }
                DbCommands::ToBin(args) => {
                    db::to_bin::cli::run(&cli.common, args)?;
                }
            },
            Commands::Seqvars(seqvars) => match &seqvars.command {
                SeqvarsCommands::Ingest(args) => {
                    seqvars::ingest::run(&cli.common, args)?;
                }
            },
            Commands::Strucvars(strucvars) => match &strucvars.command {
                StrucvarsCommands::Ingest(args) => {
                    strucvars::ingest::run(&cli.common, args)?;
                }
                StrucvarsCommands::Query(args) => {
                    strucvars::query::run(&cli.common, args)?;
                }
            },
        }

        Ok::<(), anyhow::Error>(())
    })?;
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
