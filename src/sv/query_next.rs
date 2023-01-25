use std::path::{Path, PathBuf};

use anyhow::anyhow;
use clap::{command, Parser};
use tracing::{error, info};

use crate::sv::conf::{sanity_check_db, DbConf};

#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for SVs", long_about = None)]
pub struct Args {
    /// Path to database to use for querying.
    #[arg(long)]
    pub path_db: String,
    /// Path to configuration file, defaults to `${path_db}/conf.toml`.
    #[arg(long)]
    pub path_conf: Option<String>,
    /// Disable checksum verifiation.
    #[arg(long, default_value_t = false)]
    pub disable_checksums: bool,
}

fn load_conf(args: &Args) -> Result<DbConf, anyhow::Error> {
    // Get full path to database configuration.
    let path_config = if let Some(path_conf) = &args.path_conf {
        PathBuf::from(path_conf)
    } else {
        Path::new(&args.path_db).join("conf.toml")
    };

    // Perform sanity checks on database.
    if let Some(error_msgs) = sanity_check_db(
        Path::new(&args.path_db),
        &path_config,
        !args.disable_checksums,
    )? {
        error!("Found {} errors in your database", error_msgs.len());
        for msg in &error_msgs {
            error!("error: {}", &msg);
        }
        return Err(anyhow!("Errors found in database sanity heck"));
    }

    // Load configuration
    let toml_str = std::fs::read_to_string(path_config)?;
    let conf: DbConf = toml::from_str(&toml_str)?;

    Ok(conf)
}

pub(crate) fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    let db_conf = load_conf(args)?;

    Ok(())
}
