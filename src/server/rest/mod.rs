//! Code supporting the `server rest` sub command.

use std::{collections::HashMap, path::PathBuf, str::FromStr, time::Instant};

use actix_web::web::Data;
use clap::Parser;
use enum_map::{enum_map, EnumMap};
use tracing::info;

use crate::{
    common::{build_chrom_map, trace_rss_now},
    db::conf::{GenomeRelease, Top},
    sv::query::{
        bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db, masked::load_masked_dbs,
        pathogenic::load_patho_dbs, tads::load_tads, Databases,
    },
};

pub mod actix_server;

pub struct WebServerData {
    pub chrom_map: HashMap<String, usize>,
    pub dbs: EnumMap<GenomeRelease, Databases>,
}

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
    pub listen_port: u16,
}

/// Main entry point for `server rest` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    if let Some(level) = args_common.verbose.log_level() {
        match level {
            log::Level::Trace | log::Level::Debug => {
                std::env::set_var("RUST_LOG", "debug");
                env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
            }
            _ => (),
        }
    }

    info!("Loading database config...");
    let db_conf: Top = {
        let path_conf = if let Some(path_conf) = &args.path_conf {
            PathBuf::from_str(path_conf)?
        } else {
            PathBuf::from_str(&args.path_db)?.join("conf.toml")
        };
        let toml_str = std::fs::read_to_string(&path_conf)?;
        toml::from_str(&toml_str)?
    };

    info!("Loading databases...");
    let before_loading = Instant::now();
    let dbs = enum_map! {
        GenomeRelease::Grch37 => Databases {
            bg_dbs: load_bg_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            patho_dbs: load_patho_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            tad_sets: load_tads(&args.path_db, &db_conf.features[GenomeRelease::Grch37])?,
            masked: load_masked_dbs(&args.path_db, &db_conf.features[GenomeRelease::Grch37])?,
            genes: load_gene_db(&args.path_db, &db_conf.genes, &db_conf.features[GenomeRelease::Grch37])?,
            clinvar_sv: load_clinvar_sv(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
        },
        GenomeRelease::Grch38 => Default::default(),
    };
    info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    let data = Data::new(WebServerData {
        chrom_map: build_chrom_map(),
        dbs,
    });

    trace_rss_now();

    info!("Launching server ...");
    actix_server::main(args, data)?;

    info!("All done. Have a nice day!");
    Ok(())
}
