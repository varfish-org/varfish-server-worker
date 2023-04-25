//! Implementation of phenotype-related REST API server.

use std::time::Instant;

use actix_web::web::Data;
use clap::Parser;
use hpo::Ontology;
use tracing::info;

use crate::common::trace_rss_now;

/// Data to keep in the web server.
pub struct WebServerData {
    pub ontology: Ontology,
}

/// Implementation of the actix server.
pub mod actix_server {
    use actix_web::{middleware::Logger, web::Data, App, HttpServer, ResponseError};

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

    // Code for `/hpo/genes`.
    pub mod hpo_genes {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use hpo::{
            annotations::{AnnotationId, Gene, GeneId},
            Ontology,
        };
        use serde::{Deserialize, Serialize};

        use crate::server::pheno::WebServerData;

        use super::CustomError;

        /// Specify how to perform query matches in the API calls, sucha s `fetch_genes`.
        #[derive(Deserialize, Debug, Clone, Copy, Default, PartialEq, Eq)]
        #[serde(rename_all = "lowercase")]
        enum Match {
            #[default]
            Exact,
            Prefix,
            Suffix,
            Contains,
        }

        /// Parameters for `fetch_hpo_genes`.
        ///
        /// This allows to query for genes.  The first given of the following is interpreted.
        ///
        /// - `gene_id` -- specify gene ID
        /// - `gene_symbol` -- specify the gene symbol
        /// - `max_results` -- the maximnum number of records to return
        /// - `hpo_terms` -- whether to include `"hpo_terms"` in result
        ///
        /// The following propery defines how matches are performed:
        ///
        /// - `match` -- how to match
        #[derive(Deserialize, Debug, Clone)]
        struct FetchHpoGenesRequest {
            /// The gene ID to search for.
            pub gene_id: Option<String>,
            /// The gene symbol to search for.
            pub gene_symbol: Option<String>,
            /// The match mode.
            #[serde(alias = "match")]
            pub match_: Option<Match>,
            /// Maximal number of results to return.
            #[serde(default = "_default_max_results")]
            pub max_results: usize,
            /// Whether to include HPO terms.
            #[serde(default = "_default_hpo_terms")]
            pub hpo_terms: bool,
        }

        /// Return default of `FetchHpoGenesRequest::max_results`.
        fn _default_max_results() -> usize {
            100
        }

        /// Return default of `FetchHpoGenesRequest::hpo_terms`.
        fn _default_hpo_terms() -> bool {
            false
        }

        /// Representation of an HPO term.
        #[derive(Serialize, Debug, Clone)]
        struct HpoTerm {
            /// The HPO ID.
            pub term_id: String,
            /// The description.
            pub name: String,
        }

        /// Result entry for `fetch_hpo_genes`.
        #[derive(Serialize, Debug, Clone)]
        struct ResultEntry {
            /// The gene's NCBI ID.
            pub gene_ncbi_id: u32,
            /// The gene's HGNC symbol.
            pub gene_symbol: String,
            /// The gene's associated HPO terms.
            #[serde(default = "Option::default", skip_serializing_if = "Option::is_none")]
            pub hpo_terms: Option<Vec<HpoTerm>>,
        }

        impl ResultEntry {
            pub fn from_gene_with_ontology(
                gene: &Gene,
                ontology: &Ontology,
                hpo_terms: bool,
            ) -> Self {
                let hpo_terms = if hpo_terms {
                    Some(
                        gene.hpo_terms()
                            .into_iter()
                            .map(|term_id| ontology.hpo(term_id))
                            .filter(|term| term.is_some())
                            .map(|term| {
                                let term = term.expect("filtered above");
                                HpoTerm {
                                    term_id: term.id().to_string(),
                                    name: term.name().to_string(),
                                }
                            })
                            .collect(),
                    )
                } else {
                    None
                };
                ResultEntry {
                    gene_ncbi_id: gene.id().as_u32(),
                    gene_symbol: gene.name().to_string(),
                    hpo_terms,
                }
            }
        }

        /// Query for genes in the HPO database.
        #[get("/hpo/genes")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<FetchHpoGenesRequest>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let ontology = &data.ontology;
            let match_ = query.match_.unwrap_or_default();
            let mut result: Vec<ResultEntry> = Vec::new();

            if match_ == Match::Exact {
                let gene = if let Some(gene_ncbi_id) = &query.gene_id {
                    let gene_id = GeneId::from(
                        gene_ncbi_id
                            .parse::<u32>()
                            .map_err(|e| CustomError::new(anyhow::anyhow!(e)).into())?,
                    );
                    ontology.gene(&gene_id)
                } else if let Some(gene_symbol) = &query.gene_symbol {
                    ontology.gene_by_name(&gene_symbol)
                } else {
                    None
                };
                if let Some(gene) = gene {
                    result.push(ResultEntry::from_gene_with_ontology(
                        gene,
                        ontology,
                        query.hpo_terms,
                    ));
                }
            } else {
                if let Some(gene_symbol) = &query.gene_symbol {
                    let mut it = ontology.genes();
                    let mut gene = it.next();
                    while gene.is_some() && result.len() < query.max_results {
                        let symbol = gene.expect("checked above").symbol();
                        let is_match = match query.match_.unwrap_or_default() {
                            Match::Prefix => symbol.starts_with(gene_symbol),
                            Match::Suffix => symbol.ends_with(gene_symbol),
                            Match::Contains => symbol.contains(gene_symbol),
                            Match::Exact => panic!("cannot happen here"),
                        };
                        if is_match {
                            result.push(ResultEntry::from_gene_with_ontology(
                                gene.expect("checked above"),
                                ontology,
                                query.hpo_terms,
                            ));
                        }

                        gene = it.next();
                    }
                }
            }

            Ok(Json(result))
        }
    }

    #[actix_web::main]
    pub async fn main(args: &Args, dbs: Data<WebServerData>) -> std::io::Result<()> {
        HttpServer::new(move || {
            App::new()
                .app_data(dbs.clone())
                .service(hpo_genes::handle)
                .wrap(Logger::default())
        })
        .bind((args.listen_host.as_str(), args.listen_port))?
        .run()
        .await
    }
}

/// Command line arguments for `server pheno` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run pheno REST API server", long_about = None)]
pub struct Args {
    /// Path to the directory with the HPO files.
    #[arg(long, required = true)]
    pub path_hpo_dir: String,

    /// IP to listen on.
    #[arg(long, default_value = "127.0.0.1")]
    pub listen_host: String,
    /// Port to listen on.
    #[arg(long, default_value_t = 8081)]
    pub listen_port: u16,
}

/// Main entry point for `server pheno` sub command.
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

    info!("Loading HPO HPO...");
    let before_loading = Instant::now();
    let ontology = Ontology::from_standard(&args.path_hpo_dir)?;
    let data = Data::new(WebServerData { ontology });
    info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    info!("Launching server ...");
    actix_server::main(args, data)?;

    info!("All done. Have a nice day!");
    Ok(())
}
