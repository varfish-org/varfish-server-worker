//! Code supporting the `server rest` sub command.

use clap::Parser;
use tracing::info;

/// Command line arguments for `server rest` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run REST API server", long_about = None)]
pub struct Args {
    /// Path to database to use for querying.
    #[arg(long, required = true)]
    pub path_db: String,
    /// Path to configuration file, defaults to `${path_db}/conf.toml`.
    #[arg(long)]
    pub path_conf: Option<String>,
    /// IP to listen on.
    #[arg(long, default_value = "127.0.0.1")]
    pub listen_host: String,
    /// Port to listen on.
    #[arg(long, default_value_t = 8081)]
    pub listen_port: u32,
}

/// Main entry point for `sv query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    info!("Loading databases...");

    info!("Launching server ...");

    info!("All done. Have a nice day!");
    Ok(())
}
