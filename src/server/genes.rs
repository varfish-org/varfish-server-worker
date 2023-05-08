//! Code supporting the `server genes` sub command.

use std::time::Instant;

use actix_web::web::Data;
use clap::Parser;
use tracing::info;

use crate::common::trace_rss_now;

pub struct WebServerData {
    pub db: rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
}

/// Implementation of the actix server.
pub mod actix_server {
    use actix_web::{middleware::Logger, web::Data, App, HttpServer, ResponseError};
    use serde::{Deserialize, Serialize};

    use super::{Args, WebServerData};

    #[derive(Debug)]
    struct CustomError {
        err: anyhow::Error,
    }

    impl std::fmt::Display for CustomError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{:?}", self.err)
        }
    }

    impl CustomError {
        fn new(err: anyhow::Error) -> Self {
            CustomError { err }
        }
    }

    impl ResponseError for CustomError {}

    /// Specify how to perform query matches in the API calls.
    #[derive(Deserialize, Debug, Clone, Copy, Default, PartialEq, Eq)]
    #[serde(rename_all = "lowercase")]
    enum Match {
        #[default]
        Exact,
        Prefix,
        Suffix,
        Contains,
    }

    /// Representation of a gene.
    #[derive(Serialize, Debug, Clone)]
    struct ResultGene {
        /// The HPO ID.
        pub gene_id: u32,
        /// The description.
        pub gene_symbol: String,
    }

    /// Representation of an HPO term.
    #[derive(Serialize, Debug, Clone)]
    struct ResultHpoTerm {
        /// The HPO ID.
        pub term_id: String,
        /// The description.
        pub name: String,
    }

    /// Code for `/genes/details`.
    pub mod genes_details {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use linked_hash_map::LinkedHashMap;
        use serde::{Deserialize, Serialize};
        use serde_with::{formats::CommaSeparator, StringWithSeparator};

        use crate::{db::genes, pheno::prepare::VERSION, server::genes::WebServerData};

        use super::CustomError;

        /// Parameters for `fetch_genes_details`.
        ///
        /// This allows to query for details on one or more genes.  The first given of the
        /// following is interpreted.
        ///
        /// - `hgnc_id` -- specify HGNC ID(s); comma-separated list
        ///
        /// The following propery defines how matches are performed:
        ///
        /// - `match` -- how to match
        #[serde_with::skip_serializing_none]
        #[serde_with::serde_as]
        #[derive(Deserialize, Debug, Clone)]
        struct Request {
            /// The HGNC IDs to search for.
            #[serde_as(as = "Option<StringWithSeparator::<CommaSeparator, String>>")]
            pub hgnc_id: Option<Vec<String>>,
        }

        /// Result for `handle`.
        #[derive(Serialize, Debug, Clone)]
        struct Result {
            /// Version of the server code.
            pub server_version: String,
            /// Version of the builder code.
            pub builder_version: String,
            // TODO: add data version
            /// The resulting gene information.
            pub genes: LinkedHashMap<String, genes::data::Record>,
        }

        /// Query for genes in the HPO database.
        #[get("/genes/details")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            eprintln!("{:?}", &query);
            let cf_genes = data
                .db
                .cf_handle("genes")
                .expect("no 'genes' column family");
            let mut genes = LinkedHashMap::new();
            if let Some(hgnc_id) = query.hgnc_id.as_ref() {
                for hgnc_id in hgnc_id {
                    let raw_json = data
                        .db
                        .get_cf(&cf_genes, hgnc_id)
                        .map_err(|e| {
                            CustomError::new(anyhow::anyhow!("problem querying database: {}", e))
                        })?
                        .ok_or_else(|| {
                            CustomError::new(anyhow::anyhow!("no such gene: {}", hgnc_id))
                        })?;
                    let json = std::str::from_utf8(&raw_json).map_err(|e| {
                        CustomError::new(anyhow::anyhow!("problem decoding value: {}", e))
                    })?;
                    let record: genes::data::Record = serde_json::from_str(json).map_err(|e| {
                        CustomError::new(anyhow::anyhow!(
                            "problem converting to gene record: {}",
                            e
                        ))
                    })?;
                    genes.insert(hgnc_id.to_string(), record);
                }
            }

            let cf_meta = data.db.cf_handle("meta").expect("no 'meta' column family");
            let raw_builder_version = &data
                .db
                .get_cf(&cf_meta, "builder-version")
                .map_err(|e| CustomError::new(anyhow::anyhow!("problem querying database: {}", e)))?
                .expect("database missing 'builder-version' key?");
            let builder_version = std::str::from_utf8(raw_builder_version)
                .map_err(|e| CustomError::new(anyhow::anyhow!("problem decoding value: {}", e)))?
                .to_string();

            Ok(Json(Result {
                server_version: VERSION.to_string(),
                builder_version,
                genes,
            }))
        }
    }

    #[actix_web::main]
    pub async fn main(args: &Args, dbs: Data<WebServerData>) -> std::io::Result<()> {
        HttpServer::new(move || {
            let app = App::new()
                .app_data(dbs.clone())
                .service(genes_details::handle);
            app.wrap(Logger::default())
        })
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

    info!("Opening database...");
    let before_opening = Instant::now();
    let path_rocksdb = format!("{}/genes/db", &args.path_db);
    let db = rocksdb::DB::open_cf_for_read_only(
        &rocksdb::Options::default(),
        path_rocksdb,
        ["meta", "genes"],
        true,
    )?;
    info!("...done opening database in {:?}", before_opening.elapsed());

    let data = Data::new(WebServerData { db });

    trace_rss_now();

    info!(
        "Launching server main on http://{}:{} ...",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/genes/details?hgnc_id=HGNC:1097,HGNC:11319",
        args.listen_host.as_str(),
        args.listen_port
    );
    actix_server::main(args, data)?;

    info!("All done. Have a nice day!");
    Ok(())
}
