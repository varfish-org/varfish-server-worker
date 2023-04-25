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

        use super::{CustomError, Match, ResultHpoTerm};

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
        struct Request {
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

        /// Return default of `Request::max_results`.
        fn _default_max_results() -> usize {
            100
        }

        /// Return default of `Request::hpo_terms`.
        fn _default_hpo_terms() -> bool {
            false
        }

        /// Result entry for `handle`.
        #[derive(Serialize, Debug, Clone)]
        struct ResultEntry {
            /// The gene's NCBI ID.
            pub gene_ncbi_id: u32,
            /// The gene's HGNC symbol.
            pub gene_symbol: String,
            /// The gene's associated HPO terms.
            #[serde(default = "Option::default", skip_serializing_if = "Option::is_none")]
            pub hpo_terms: Option<Vec<ResultHpoTerm>>,
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
                                ResultHpoTerm {
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
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let ontology = &data.ontology;
            let match_ = query.match_.unwrap_or_default();
            let mut result: Vec<ResultEntry> = Vec::new();

            if match_ == Match::Exact {
                let gene = if let Some(gene_ncbi_id) = &query.gene_id {
                    let gene_id = GeneId::from(
                        gene_ncbi_id
                            .parse::<u32>()
                            .map_err(|e| CustomError::new(anyhow::anyhow!(e)))?,
                    );
                    ontology.gene(&gene_id)
                } else if let Some(gene_symbol) = &query.gene_symbol {
                    ontology.gene_by_name(gene_symbol)
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
            } else if let Some(gene_symbol) = &query.gene_symbol {
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

            Ok(Json(result))
        }
    }

    // Code for `/hpo/terms`.
    pub mod hpo_terms {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use hpo::{annotations::AnnotationId, HpoTerm, HpoTermId, Ontology};
        use serde::{Deserialize, Serialize};

        use crate::server::pheno::WebServerData;

        use super::{CustomError, Match, ResultGene};

        /// Parameters for `handle`.
        ///
        /// This allows to query for terms.  The first given of the following is interpreted.
        ///
        /// - `term_id` -- specify term ID
        /// - `gene_symbol` -- specify the gene symbol
        /// - `max_results` -- the maximum number of records to return
        /// - `genes` -- whether to include `"genes"` in result
        ///
        /// The following propery defines how matches are performed:
        ///
        /// - `match` -- how to match
        #[derive(Deserialize, Debug, Clone)]
        struct Request {
            /// The term ID to search for.
            pub term_id: Option<String>,
            /// The term name to search for.
            pub name: Option<String>,
            /// The match mode.
            #[serde(alias = "match")]
            pub match_: Option<Match>,
            /// Maximal number of results to return.
            #[serde(default = "_default_max_results")]
            pub max_results: usize,
            /// Whether to include genes.
            #[serde(default = "_default_genes")]
            pub genes: bool,
        }

        /// Return default of `Request::max_results`.
        fn _default_max_results() -> usize {
            100
        }

        /// Return default of `Request::genes`.
        fn _default_genes() -> bool {
            false
        }

        /// Result entry for `fetch_hpo_genes`.
        #[derive(Serialize, Debug, Clone)]
        struct ResultEntry {
            /// The HPO term's ID.
            pub term_id: String,
            /// The HPO term's name.
            pub name: String,
            /// The gene's associated HPO terms.
            #[serde(default = "Option::default", skip_serializing_if = "Option::is_none")]
            pub genes: Option<Vec<ResultGene>>,
        }

        impl ResultEntry {
            pub fn from_term_with_ontology(
                term: &HpoTerm,
                ontology: &Ontology,
                genes: bool,
            ) -> Self {
                let genes = if genes {
                    Some(
                        term.gene_ids()
                            .iter()
                            .map(|gene_id| ontology.gene(gene_id))
                            .filter(|term| term.is_some())
                            .map(|term| {
                                let gene = term.expect("filtered above");
                                ResultGene {
                                    gene_id: gene.id().as_u32(),
                                    gene_symbol: gene.name().to_string(),
                                }
                            })
                            .collect(),
                    )
                } else {
                    None
                };
                ResultEntry {
                    term_id: term.id().to_string(),
                    name: term.name().to_string(),
                    genes,
                }
            }
        }

        /// Query for terms in the HPO database.
        #[get("/hpo/terms")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let ontology = &data.ontology;
            let match_ = query.match_.unwrap_or_default();
            let mut result: Vec<ResultEntry> = Vec::new();

            if match_ == Match::Exact {
                let term = if let Some(term_id) = &query.term_id {
                    let term_id = HpoTermId::from(term_id.clone());
                    ontology.hpo(term_id)
                } else if let Some(name) = &query.name {
                    let mut term = None;
                    let mut it = ontology.hpos();
                    let mut tmp = it.next();
                    while tmp.is_some() && term.is_none() {
                        if tmp.expect("checked above").name() == name {
                            term = tmp;
                        }
                        tmp = it.next();
                    }
                    term
                } else {
                    None
                };
                if let Some(term) = &term {
                    result.push(ResultEntry::from_term_with_ontology(
                        term,
                        ontology,
                        query.genes,
                    ));
                }
            } else if let Some(name) = &query.name {
                let mut it = ontology.hpos();
                let mut term = it.next();
                while term.is_some() && result.len() < query.max_results {
                    let term_name = term.as_ref().expect("checked above").name();
                    let is_match = match query.match_.unwrap_or_default() {
                        Match::Prefix => term_name.starts_with(name),
                        Match::Suffix => term_name.ends_with(name),
                        Match::Contains => term_name.contains(name),
                        Match::Exact => panic!("cannot happen here"),
                    };
                    if is_match {
                        result.push(ResultEntry::from_term_with_ontology(
                            term.as_ref().expect("checked above"),
                            ontology,
                            query.genes,
                        ));
                    }

                    term = it.next();
                }
            }

            Ok(Json(result))
        }
    }

    // Code for `/hpo/omims`.
    pub mod hpo_omims {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use hpo::{
            annotations::{OmimDisease, OmimDiseaseId},
            Ontology,
        };
        use serde::{Deserialize, Serialize};

        use crate::server::pheno::WebServerData;

        use super::{CustomError, Match, ResultHpoTerm};

        /// Parameters for `handle`.
        ///
        /// This allows to query for diseases.  The first given of the following is interpreted.
        ///
        /// - `omim_id` -- specify disease ID
        /// - `name` -- specify the name to query for
        /// - `max_results` -- the maximum number of records to return
        /// - `hpo_terms` -- whether to include `"hpo_terms"` in result
        ///
        /// The following propery defines how matches are performed:
        ///
        /// - `match` -- how to match
        #[derive(Deserialize, Debug, Clone)]
        struct Request {
            /// The OMIM ID to search for.
            pub omim_id: Option<String>,
            /// The disease name to search for.
            pub name: Option<String>,
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

        /// Return default of `Request::max_results`.
        fn _default_max_results() -> usize {
            100
        }

        /// Return default of `Request::hpo_terms`.
        fn _default_hpo_terms() -> bool {
            false
        }

        /// Result entry for `handle`.
        #[derive(Serialize, Debug, Clone)]
        struct ResultEntry {
            /// The OMIM ID.
            pub omim_id: String,
            /// The OMIM disease name.
            pub name: String,
            /// The gene's associated HPO terms.
            #[serde(default = "Option::default", skip_serializing_if = "Option::is_none")]
            pub hpo_terms: Option<Vec<ResultHpoTerm>>,
        }

        impl ResultEntry {
            pub fn from_omim_disease_with_ontology(
                omim_disease: &OmimDisease,
                ontology: &Ontology,
                hpo_terms: bool,
            ) -> Self {
                let hpo_terms = if hpo_terms {
                    Some(
                        omim_disease
                            .hpo_terms()
                            .into_iter()
                            .map(|term_id| ontology.hpo(term_id))
                            .filter(|term| term.is_some())
                            .map(|term| {
                                let term = term.expect("filtered above");
                                ResultHpoTerm {
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
                    omim_id: omim_disease.id().to_string(),
                    name: omim_disease.name().to_string(),
                    hpo_terms,
                }
            }
        }

        /// Query for OMIM diseases in the HPO database.
        #[get("/hpo/omims")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let ontology = &data.ontology;
            let match_ = query.match_.unwrap_or_default();
            let mut result: Vec<ResultEntry> = Vec::new();

            if match_ == Match::Exact {
                let omim_disease = if let Some(omim_id) = &query.omim_id {
                    let omim_id = OmimDiseaseId::try_from(omim_id.as_ref())
                        .map_err(|e| CustomError::new(anyhow::anyhow!(e)))?;
                    ontology.omim_disease(&omim_id)
                } else if let Some(name) = &query.name {
                    let mut omim_disease = None;
                    let mut it = ontology.omim_diseases();
                    let mut tmp = it.next();
                    while tmp.is_some() && omim_disease.is_none() {
                        if tmp.expect("checked above").name() == name {
                            omim_disease = tmp;
                        }
                        tmp = it.next();
                    }
                    omim_disease
                } else {
                    None
                };
                if let Some(omim_disease) = &omim_disease {
                    result.push(ResultEntry::from_omim_disease_with_ontology(
                        omim_disease,
                        ontology,
                        query.hpo_terms,
                    ));
                }
            } else if let Some(name) = &query.name {
                let mut it = ontology.omim_diseases();
                let mut omim_disease = it.next();
                while omim_disease.is_some() && result.len() < query.max_results {
                    let omim_name = omim_disease.as_ref().expect("checked above").name();
                    let is_match = match query.match_.unwrap_or_default() {
                        Match::Prefix => omim_name.starts_with(name),
                        Match::Suffix => omim_name.ends_with(name),
                        Match::Contains => omim_name.contains(name),
                        Match::Exact => panic!("cannot happen here"),
                    };
                    if is_match {
                        result.push(ResultEntry::from_omim_disease_with_ontology(
                            omim_disease.as_ref().expect("checked above"),
                            ontology,
                            query.hpo_terms,
                        ));
                    }

                    omim_disease = it.next();
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
                .service(hpo_terms::handle)
                .service(hpo_omims::handle)
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
