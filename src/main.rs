//! VarFish Server Worker main executable

pub mod common;
pub mod db;
pub mod pheno;
pub mod server;
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
    /// Database-related commands.
    Db(Db),
    /// Phenotype-related commands.
    Pheno(Pheno),
    /// SV related commands.
    Sv(Sv),
    /// Server related commands.
    Server(Server),
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
    Copy(annonars::db_utils::cli::copy::Args),
    Compile(db::compile::Args),
    MkInhouse(db::mk_inhouse::Args),
    Genes(Genes),
    Seqvars(Seqvars),
}

/// Parsing of "db genes *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Genes {
    /// The sub command to run
    #[command(subcommand)]
    command: GenesCommands,
}

/// Enum supporting the parsing of "db genes *" sub commands.
#[derive(Debug, Subcommand)]
enum GenesCommands {
    Build(db::genes::build::Args),
}

/// Parsing of "db seqvars *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Seqvars {
    /// The sub command to run
    #[command(subcommand)]
    command: SeqvarsCommands,
}

/// Enum supporting the parsing of "db seqvars *" sub commands.
#[derive(Debug, Subcommand)]
enum SeqvarsCommands {
    Build(SeqvarsBuild),
}

/// Enum supporting the parsing of "db seqvars *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct SeqvarsBuild {
    /// The sub command to run
    #[command(subcommand)]
    command: SeqvarsBuildCommands,
}

/// Enum supporting the parsing of "db seqvars *" sub commands.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum SeqvarsBuildCommands {
    Dbsnp(annonars::dbsnp::cli::import::Args),
    GnomadNuclear(annonars::gnomad_nuclear::cli::import::Args),
    GnomadMtdna(annonars::gnomad_mtdna::cli::import::Args),
    Helixmtdb(annonars::helixmtdb::cli::import::Args),
    UcscConservation(annonars::cons::cli::import::Args),
    Tsv(annonars::tsv::cli::import::Args),
}

/// Parsing of "pheno *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Pheno {
    /// The sub command to run
    #[command(subcommand)]
    command: PhenoCommands,
}

/// Enum supporting the parsing of "db *" sub commands.
#[derive(Debug, Subcommand)]
enum PhenoCommands {
    Prepare(pheno::prepare::Args),
    Query(pheno::query::Args),
}

/// Parsing of "sv *" sub commands.
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
    Query(sv::query::Args),
}

/// Enum supporting the parsing of "server *" sub commands.
#[derive(Debug, Subcommand)]
enum ServerCommands {
    Annos(server::annos::Args),
    Genes(server::genes::Args),
    Rest(server::rest::Args),
    Pheno(server::pheno::Args),
}
/// Parsing of "sv *" sub commands.
#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
struct Server {
    /// The sub command to run
    #[command(subcommand)]
    command: ServerCommands,
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

    // Common variables for CLI commands from annonars.
    let annonars_common = annonars::common::cli::Args {
        verbose: cli.common.verbose.clone(),
    };

    // Install collector and go into sub commands.
    let term = Term::stderr();
    tracing::subscriber::with_default(collector, || {
        match &cli.command {
            Commands::Db(db) => match &db.command {
                DbCommands::Compile(args) => {
                    db::compile::run(&cli.common, args)?;
                }
                DbCommands::MkInhouse(args) => {
                    db::mk_inhouse::run(&cli.common, args)?;
                }
                DbCommands::Genes(args) => match &args.command {
                    GenesCommands::Build(args) => {
                        db::genes::build::run(&cli.common, args)?;
                    }
                },
                DbCommands::Seqvars(args) => match &args.command {
                    SeqvarsCommands::Build(args) => match &args.command {
                        SeqvarsBuildCommands::Dbsnp(args) => {
                            annonars::dbsnp::cli::import::run(&annonars_common, args)?;
                        }
                        SeqvarsBuildCommands::GnomadNuclear(args) => {
                            annonars::gnomad_nuclear::cli::import::run(&annonars_common, args)?;
                        }
                        SeqvarsBuildCommands::GnomadMtdna(args) => {
                            annonars::gnomad_mtdna::cli::import::run(&annonars_common, args)?;
                        }
                        SeqvarsBuildCommands::Helixmtdb(args) => {
                            annonars::helixmtdb::cli::import::run(&annonars_common, args)?;
                        }
                        SeqvarsBuildCommands::UcscConservation(args) => {
                            annonars::cons::cli::import::run(&annonars_common, args)?;
                        }
                        SeqvarsBuildCommands::Tsv(args) => {
                            annonars::tsv::cli::import::run(&annonars_common, args)?;
                        }
                    },
                },
                DbCommands::Copy(args) => {
                    annonars::db_utils::cli::copy::run(&annonars_common, args)?
                }
            },
            Commands::Pheno(pheno) => match &pheno.command {
                PhenoCommands::Prepare(args) => pheno::prepare::run(&cli.common, args)?,
                PhenoCommands::Query(args) => pheno::query::run(&cli.common, args)?,
            },
            Commands::Sv(sv) => match &sv.command {
                SvCommands::Query(args) => {
                    sv::query::run(&cli.common, args)?;
                }
            },
            Commands::Server(server) => match &server.command {
                ServerCommands::Annos(args) => server::annos::run(&cli.common, args)?,
                ServerCommands::Genes(args) => server::genes::run(&cli.common, args)?,
                ServerCommands::Rest(args) => server::rest::run(&cli.common, args)?,
                ServerCommands::Pheno(args) => server::pheno::run(&cli.common, args)?,
            },
        }

        Ok::<(), anyhow::Error>(())
    })?;
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
