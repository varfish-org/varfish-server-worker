//! Code supporting the `server rest` sub command.

use std::{path::PathBuf, str::FromStr, time::Instant};

use clap::Parser;
use enum_map::{enum_map, EnumMap};
use tracing::info;

use crate::{
    common::trace_rss_now,
    db::conf::{GenomeRelease, Top},
    sv::query::{
        bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db,
        pathogenic::load_patho_dbs, tads::load_tads, Databases,
    },
};

/// Implementation of the actix server.
pub mod actix_server {
    use actix_web::{App, HttpServer};
    use enum_map::EnumMap;

    use crate::{db::conf::GenomeRelease, sv::query::Databases};

    use super::Args;

    #[actix_web::main]
    pub async fn main(
        args: &Args,
        _dbs: &EnumMap<GenomeRelease, Databases>,
    ) -> std::io::Result<()> {
        HttpServer::new(|| App::new())
            .bind((args.listen_host.as_str(), args.listen_port))?
            .run()
            .await
    }
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

/// Main entry point for `sv query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    info!("Loading database config...");
    let db_conf: Top = {
        let path_conf = if let Some(path_conf) = &args.path_conf {
            PathBuf::from_str(&path_conf)?
        } else {
            PathBuf::from_str(&args.path_db)?.join("conf.toml")
        };
        let toml_str = std::fs::read_to_string(&path_conf)?;
        toml::from_str(&toml_str)?
    };

    info!("Loading databases...");
    let before_loading = Instant::now();
    let dbs: EnumMap<GenomeRelease, Databases> = enum_map! {
        GenomeRelease::Grch37 => Databases {
            bg_dbs: load_bg_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            patho_dbs: load_patho_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            tad_sets: load_tads(&args.path_db, &db_conf.features[GenomeRelease::Grch37])?,
            genes: load_gene_db(&args.path_db, &db_conf.genes, &db_conf.features[GenomeRelease::Grch37])?,
            clinvar_sv: load_clinvar_sv(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
        },
        GenomeRelease::Grch38 => Default::default(),
    };
    info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    info!("Launching server ...");
    actix_server::main(&args, &dbs)?;

    info!("All done. Have a nice day!");
    Ok(())
}
