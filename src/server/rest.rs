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
        bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db,
        pathogenic::load_patho_dbs, tads::load_tads, Databases,
    },
};

pub struct WebServerData {
    pub chrom_map: HashMap<String, usize>,
    pub dbs: EnumMap<GenomeRelease, Databases>,
}

/// Implementation of the actix server.
pub mod actix_server {
    use std::str::FromStr;

    use actix_web::{
        get,
        web::{self, Data, Json, Path},
        App, HttpServer, Responder, ResponseError,
    };
    use serde::Serialize;
    use tracing::trace;

    use crate::{common::CHROMS, db::conf::TadSet as TadSetChoice};
    use crate::{db::conf::GenomeRelease, sv::query::records::ChromRange};

    use super::{Args, WebServerData};

    #[derive(Debug)]
    struct MyError {
        err: anyhow::Error,
    }

    impl std::fmt::Display for MyError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{:?}", self.err)
        }
    }

    impl MyError {
        fn new(err: anyhow::Error) -> Self {
            MyError { err }
        }
    }

    impl ResponseError for MyError {}

    /// Result type of `fetch_tads`.
    #[derive(Serialize, Debug)]
    struct Tad {
        chromosome: String,
        begin: u32,
        end: u32,
    }

    /// List the TADs of the given TAD set
    #[get("/public/svs/tads/{release}/{tad_set}/")]
    async fn fetch_tads(
        data: Data<WebServerData>,
        path: Path<(String, String)>,
        chrom_range: web::Query<ChromRange>,
    ) -> actix_web::Result<impl Responder, MyError> {
        let genome_release =
            GenomeRelease::from_str(&path.0).map_err(|e| MyError::new(e.into()))?;
        trace!("genome_release = {:?}", genome_release);
        let tad_set = TadSetChoice::from_str(&path.1).map_err(|e| MyError::new(e.into()))?;
        trace!("tad_set = {:?}", tad_set);
        let tads = data.dbs[genome_release]
            .tad_sets
            .fetch_tads(tad_set, &chrom_range, &data.chrom_map)
            .into_iter()
            .map(|record| Tad {
                chromosome: CHROMS[record.chrom_no as usize].to_string(),
                begin: record.begin,
                end: record.end,
            })
            .collect::<Vec<Tad>>();
        trace!("tads = {:?}", &tads);
        Ok(Json(tads))
    }

    #[actix_web::main]
    pub async fn main(args: &Args, dbs: Data<WebServerData>) -> std::io::Result<()> {
        HttpServer::new(move || App::new().app_data(dbs.clone()).service(fetch_tads))
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

    if let Some(level) = args_common.verbose.log_level() {
        match level {
            log::Level::Trace | log::Level::Debug => {
                std::env::set_var("RUST_LOG", "debug");
                env_logger::init();
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
