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

/// Similarity computation using the Phenomizer method.
pub(crate) mod phenomizer {
    use hpo::{
        similarity::{Builtins, Similarity},
        term::{HpoGroup, InformationContentKind},
        Ontology,
    };

    /// Compute symmetric similarity score.
    pub(crate) fn score(q: &HpoGroup, d: &HpoGroup, o: &Ontology) -> f32 {
        let s = Builtins::Resnik(InformationContentKind::Omim);
        0.5 * score_dir(q, d, o, &s) + 0.5 * score_dir(d, q, o, &s)
    }

    /// "Directed" score part of phenomizer score.
    fn score_dir(qs: &HpoGroup, ds: &HpoGroup, o: &Ontology, s: &impl Similarity) -> f32 {
        // Handle case of empty `qs`.
        if qs.is_empty() {
            return 0f32;
        }

        // For each `q in qs` compute max similarity to any `d in ds`.
        let mut tmp: Vec<f32> = Vec::new();
        for q in qs {
            if let Some(q) = o.hpo(q) {
                tmp.push(
                    ds.iter()
                        .filter_map(|d| o.hpo(d).map(|d| q.similarity_score(&d, s)))
                        .max_by(|a, b| a.partial_cmp(b).expect("try to compare NaN"))
                        .unwrap_or_default(),
                );
            }
        }

        tmp.iter().sum::<f32>() / (qs.len() as f32)
    }
}

/// Implementation of the actix server.
pub mod actix_server {
    use actix_web::{middleware::Logger, web::Data, App, HttpServer, ResponseError};
    use serde::{Deserialize, Deserializer, Serialize};

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

    /// Helper to deserialize a comma-separated list of strings.
    fn vec_str_deserialize<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let str_sequence = String::deserialize(deserializer)?;
        Ok(str_sequence
            .split(',')
            .map(|item| item.to_owned())
            .collect())
    }

    /// Helper to deserialize a comma-separated list of strings.
    fn option_vec_str_deserialize<'de, D>(deserializer: D) -> Result<Option<Vec<String>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let str_sequence = String::deserialize(deserializer)?;
        if str_sequence.is_empty() {
            Ok(None)
        } else {
            Ok(Some(
                str_sequence
                    .split(',')
                    .map(|item| item.to_owned())
                    .collect(),
            ))
        }
    }

    /// Code for `/hpo/genes`.
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

    /// Code for `/hpo/terms`.
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

    /// Code for `/hpo/omims`.
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

    /// Code for `/hpo/sim/{term-term,term-gene}` endpoint.
    pub mod hpo_sim {
        /// Entry point `/hpo/sim/term-term` allows the pairwise similary computation between two
        /// sets of HPO terms.
        pub mod term_term {
            use std::str::FromStr;

            use actix_web::{
                get,
                web::{self, Data, Json, Path},
                Responder,
            };
            use hpo::{
                similarity::{Builtins, Similarity},
                term::InformationContentKind,
                HpoTermId, Ontology,
            };
            use itertools::Itertools;
            use serde::{Deserialize, Serialize};

            use crate::server::pheno::{actix_server::CustomError, WebServerData};

            /// Enum for representing similarity method to use.
            #[derive(Default, Debug, Clone, Copy, derive_more::Display)]
            pub enum SimilarityMethod {
                #[default]
                #[display(fmt = "resnik::omim")]
                ResnikOmim,
            }

            impl FromStr for SimilarityMethod {
                type Err = anyhow::Error;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    Ok(match s {
                        "resnik::omim" => Self::ResnikOmim,
                        _ => anyhow::bail!("unknown similarity method: {}", s),
                    })
                }
            }

            impl From<SimilarityMethod> for Builtins {
                fn from(val: SimilarityMethod) -> Self {
                    match val {
                        SimilarityMethod::ResnikOmim => {
                            Builtins::Resnik(InformationContentKind::Omim)
                        }
                    }
                }
            }

            /// Parameters for `handle`.
            ///
            /// This allows to compute differences between
            ///
            /// - `lhs` -- first set of terms to compute similarity for
            /// - `rhs` -- econd set of terms to compute similarity for
            #[derive(Deserialize, Debug, Clone)]
            struct Request {
                /// The one set of HPO terms to compute similarity for.
                #[serde(deserialize_with = "super::super::vec_str_deserialize")]
                pub lhs: Vec<String>,
                /// The second set of HPO terms to compute similarity for.
                #[serde(deserialize_with = "super::super::vec_str_deserialize")]
                pub rhs: Vec<String>,
                /// The similarity method to use.
                #[serde(
                    deserialize_with = "help::similarity_deserialize",
                    default = "help::default_sim"
                )]
                pub sim: SimilarityMethod,
            }

            /// Helpers for deserializing `Request`.
            mod help {
                /// Helper to deserialize a similarity
                pub fn similarity_deserialize<'de, D>(
                    deserializer: D,
                ) -> Result<super::SimilarityMethod, D::Error>
                where
                    D: serde::Deserializer<'de>,
                {
                    let s = <String as serde::Deserialize>::deserialize(deserializer)?;
                    std::str::FromStr::from_str(&s).map_err(serde::de::Error::custom)
                }

                /// Default value for `Request::sim`.
                pub fn default_sim() -> super::SimilarityMethod {
                    super::SimilarityMethod::ResnikOmim
                }
            }

            /// Result entry for `handle`.
            #[derive(Serialize, Debug, Clone)]
            struct ResultEntry {
                /// The lhs entry.
                pub lhs: String,
                /// The rhs entry.
                pub rhs: String,
                /// The similarity score.
                pub score: f32,
                /// The score type that was used to compute the similarity for.
                pub sim: String,
            }

            /// Query for pairwise term similarity.
            ///
            /// In the case of Resnik, this corresponds to `IC(MICA(t_1, t_2))`.
            #[get("/hpo/sim/term-term")]
            async fn handle(
                data: Data<WebServerData>,
                _path: Path<()>,
                query: web::Query<Request>,
            ) -> actix_web::Result<impl Responder, CustomError> {
                let ontology: &Ontology = &data.ontology;
                let mut result = Vec::new();

                let ic: Builtins = query.sim.into();

                // Translate strings from the query into HPO terms.
                let lhs = query
                    .lhs
                    .iter()
                    .filter_map(|lhs| ontology.hpo(HpoTermId::from(lhs.clone())))
                    .collect::<Vec<_>>();
                let rhs = query
                    .rhs
                    .iter()
                    .filter_map(|rhs| ontology.hpo(HpoTermId::from(rhs.clone())))
                    .collect::<Vec<_>>();

                // Compute the similarity for each pair.
                for (lhs, rhs) in lhs.iter().cartesian_product(rhs.iter()) {
                    let similarity = ic.calculate(lhs, rhs);
                    let elem = ResultEntry {
                        lhs: lhs.id().to_string(),
                        rhs: rhs.id().to_string(),
                        score: similarity,
                        sim: query.sim.to_string(),
                    };
                    result.push(elem);
                }

                Ok(Json(result))
            }
        }

        /// Entry point `/hpo/sim/term-gene` that allows the similarity computation between a
        /// set of terms and a gene.
        pub mod term_gene {
            use std::str::FromStr;

            use actix_web::{
                get,
                web::{self, Data, Json, Path},
                Responder,
            };

            use hpo::{
                annotations::{AnnotationId, GeneId},
                term::HpoGroup,
                HpoTermId, Ontology,
            };
            use serde::{Deserialize, Serialize};

            use super::super::CustomError;
            use crate::server::pheno::{phenomizer, WebServerData};

            /// Enum for representing similarity method to use.
            #[derive(Default, Debug, Clone, Copy, derive_more::Display)]
            pub enum SimilarityMethod {
                #[default]
                #[display(fmt = "phenomizer")]
                Phenomizer,
            }

            impl FromStr for SimilarityMethod {
                type Err = anyhow::Error;

                fn from_str(s: &str) -> Result<Self, Self::Err> {
                    Ok(match s {
                        "resnik::omim" => Self::Phenomizer,
                        _ => anyhow::bail!("unknown similarity method: {}", s),
                    })
                }
            }

            /// Parameters for `handle`.
            ///
            /// This allows to compute differences between
            ///
            /// - `terms` -- set of terms to use as query
            /// - `gene_ids` -- set of ids for genes to use as "database"
            /// - `gene_symbols` -- set of symbols for genes to use as "database"
            #[derive(Deserialize, Debug, Clone)]
            struct Request {
                /// Set of terms to use as query.
                #[serde(deserialize_with = "super::super::vec_str_deserialize")]
                pub terms: Vec<String>,
                /// The set of ids for genes to use as "database".
                #[serde(
                    default = "Option::default",
                    skip_serializing_if = "Option::is_none",
                    deserialize_with = "super::super::option_vec_str_deserialize"
                )]
                pub gene_ids: Option<Vec<String>>,
                /// The set of symbols for genes to use as "database".
                #[serde(
                    default = "Option::default",
                    skip_serializing_if = "Option::is_none",
                    deserialize_with = "super::super::option_vec_str_deserialize"
                )]
                pub gene_symbols: Option<Vec<String>>,
                /// The similarity method to use.
                #[serde(
                    deserialize_with = "help::similarity_deserialize",
                    default = "help::default_sim"
                )]
                pub sim: SimilarityMethod,
            }

            /// Helpers for deserializing `Request`.
            mod help {
                /// Helper to deserialize a similarity method.
                pub fn similarity_deserialize<'de, D>(
                    deserializer: D,
                ) -> Result<super::SimilarityMethod, D::Error>
                where
                    D: serde::Deserializer<'de>,
                {
                    let s = <String as serde::Deserialize>::deserialize(deserializer)?;
                    std::str::FromStr::from_str(&s).map_err(serde::de::Error::custom)
                }

                /// Default value for `Request::sim`.
                pub fn default_sim() -> super::SimilarityMethod {
                    super::SimilarityMethod::Phenomizer
                }
            }

            /// Result entry for one gene.
            #[derive(Serialize, Debug, Clone)]
            struct ResultEntry {
                /// The gene ID of the entry.
                pub gene_id: u32,
                /// The gene HGNC symbol of the gene.
                pub gene_symbol: String,
                /// The similarity score.
                pub score: f32,
                /// The score type that was used to compute the similarity for.
                pub sim: String,
            }

            /// Query for similarity between a set of terms to each entry in a list of genes.
            #[get("/hpo/sim/term-gene")]
            async fn handle(
                data: Data<WebServerData>,
                _path: Path<()>,
                query: web::Query<Request>,
            ) -> actix_web::Result<impl Responder, CustomError> {
                let ontology: &Ontology = &data.ontology;

                // Translate strings from the query into an `HpoGroup`.
                let query_terms = {
                    let mut query_terms = HpoGroup::new();
                    for term in &query.terms {
                        if let Some(term) = ontology.hpo(HpoTermId::from(term.clone())) {
                            query_terms.insert(term.id());
                        }
                    }
                    query_terms
                };

                // Translate strings from the query into genes via symbol or gene ID.
                let genes = if let Some(gene_ids) = &query.gene_ids {
                    Ok(gene_ids
                        .iter()
                        .filter_map(|gene_id| {
                            if let Ok(gene_id) = gene_id.parse::<u32>() {
                                ontology.gene(&GeneId::from(gene_id))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>())
                } else if let Some(gene_symbols) = &query.gene_symbols {
                    Ok(gene_symbols
                        .iter()
                        .filter_map(|gene_symbol| ontology.gene_by_name(gene_symbol))
                        .collect::<Vec<_>>())
                } else {
                    Err(CustomError::new(anyhow::anyhow!(
                        "either `gene_ids` or `gene_symbols` must be given"
                    )))
                }?;

                // Compute the similarity for the given term for each.  Create one result record
                // for each gene.
                let mut result = Vec::new();
                for gene in genes {
                    let gene_terms = gene.hpo_terms();
                    let score = phenomizer::score(&query_terms, gene_terms, ontology);

                    result.push(ResultEntry {
                        gene_id: gene.id().as_u32(),
                        gene_symbol: gene.name().to_string(),
                        score,
                        sim: query.sim.to_string(),
                    });
                }

                Ok(Json(result))
            }
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
                .service(hpo_sim::term_term::handle)
                .service(hpo_sim::term_gene::handle)
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

    info!("Loading HPO...");
    let before_loading = Instant::now();
    let ontology = Ontology::from_standard(&args.path_hpo_dir)?;
    let data = Data::new(WebServerData { ontology });
    info!("...done loading HPO in {:?}", before_loading.elapsed());

    trace_rss_now();

    info!(
        "Launching server main on http://{}:{} ...",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/genes?gene_symbol=TGDS",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/genes?gene_id=23483&hpo_terms=true",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/omims?omim_id=616145&hpo_terms=true",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/terms?term_id=HP:0000023&genes=true",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/sim/term-term?lhs=HP:0001166,HP:0040069&rhs=HP:0005918,\
        HP:0004188&sim=resnik::omim",
        args.listen_host.as_str(),
        args.listen_port
    );
    info!(
        "  try: http://{}:{}/hpo/sim/term-gene?terms=HP:0001166,HP:0000098&gene_symbols=FBN1,TGDS,TTN",
        args.listen_host.as_str(),
        args.listen_port
    );
    actix_server::main(args, data)?;

    info!("All done. Have a nice day!");
    Ok(())
}
