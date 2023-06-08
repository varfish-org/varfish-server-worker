//! Code for ranking genes on the command line.

use hpo::similarity::Builtins;
use prost::Message;
use rocksdb::{DBWithThreadMode, MultiThreaded};
use serde::{Deserialize, Serialize};
use std::time::Instant;

use clap::Parser;
use hpo::{annotations::AnnotationId, term::HpoGroup, HpoTermId, Ontology};
use tracing::info;

use crate::common::trace_rss_now;
use crate::pheno::algos::phenomizer;
use crate::pheno::pbs::SimulationResults;
use crate::pheno::prepare::VERSION;
use crate::pheno::query::query_result::TermDetails;
use crate::server::pheno::actix_server::hpo_sim::term_gene::SimilarityMethod;

/// Command line arguments for `server prepare-pheno` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Prepare values for `server pheno`", long_about = None)]
pub struct Args {
    /// Path to the directory with the HPO files.
    #[arg(long, required = true)]
    pub path_hpo_dir: String,

    /// Path to JSON file with the genes to rank.
    #[arg(long)]
    pub path_genes_json: String,
    /// Path to JSON file with HPO IDs of patient.
    #[arg(long)]
    pub path_terms_json: String,
}

/// Struct for loading a gene from JSON.
#[derive(Deserialize, Debug, Clone)]
pub struct Gene {
    /// The gene symbol.
    pub gene_symbol: String,
}

/// Struct for loading an HPO term from JSON.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct HpoTerm {
    /// The term ID.
    pub term_id: String,
    /// The term name (optional).
    #[serde(default = "Option::default")]
    pub term_name: Option<String>,
}

/// Query result records.
pub mod query_result {
    use serde::Serialize;

    use super::HpoTerm;

    /// Result container data structure.
    #[derive(Serialize, Debug, Clone)]
    pub struct Container {
        /// Version of the HPO.
        pub hpo_version: String,
        /// Version of the `varfish-server-worker` package.
        pub varfish_version: String,
        /// The scoring method used.
        pub score_method: String,
        /// The original query records.
        pub query: Vec<HpoTerm>,
        /// The resulting records for the scored genes.
        pub result: Vec<Record>,
    }

    /// Store score for a record with information on individual terms.
    #[derive(Serialize, Debug, Clone)]
    pub struct Record {
        /// The gene symbol.
        pub gene_symbol: String,
        /// The estimate for empirical P-value
        pub p_value: f32,
        /// The score (`-10 * log10(p_value)`).
        pub score: f32,
        /// Details on individual terms.
        #[serde(default = "Option::default")]
        pub terms: Option<Vec<TermDetails>>,
    }

    /// Detailed term scores.
    #[derive(Serialize, Debug, Clone)]
    pub struct TermDetails {
        /// The query HPO term.
        pub term_query: Option<HpoTerm>,
        /// The gene's HPO term.
        pub term_gene: HpoTerm,
        /// The similarity score.
        pub score: f32,
    }
}

/// Run the actual phenotypic similarity query for patient terms and list of
/// genes.
///
/// # Arguments
///
/// * `patient`: The query/patient HPO terms.
/// * `genes`: The list of genes to score.
/// * `hpo`: The HPO ontology.
/// * `db`: The RocksDB instance for the Resnik P-values.
///
/// # Returns
///
/// * `Ok(query_result::Container)` if successful.
pub fn run_query(
    patient: &HpoGroup,
    genes: &Vec<&hpo::annotations::Gene>,
    hpo: &Ontology,
    db: &DBWithThreadMode<MultiThreaded>,
) -> Result<query_result::Container, anyhow::Error> {
    let cf_resnik = db.cf_handle("resnik_pvalues").unwrap();

    let num_terms = std::cmp::min(10, patient.len());
    let mut result = query_result::Container {
        hpo_version: hpo.hpo_version(),
        varfish_version: VERSION.to_string(),
        score_method: SimilarityMethod::Phenomizer.to_string(),
        query: patient
            .iter()
            .map(|t| {
                let term = hpo.hpo(t).expect("could not resolve HPO term");
                HpoTerm {
                    term_id: term.id().to_string(),
                    term_name: Some(term.name().to_string()),
                }
            })
            .collect(),
        result: Vec::new(),
    };
    for gene in genes {
        let ncbi_gene_id = gene.id().as_u32();
        let key = format!("{ncbi_gene_id}:{num_terms}");
        let data = db
            .get_cf(&cf_resnik, key.as_bytes())?
            .expect("key not found");
        let res = SimulationResults::decode(&data[..])?;
        let score = phenomizer::score(
            patient,
            &HpoGroup::from_iter(
                gene.to_hpo_set(hpo)
                    .child_nodes()
                    .without_modifier()
                    .into_iter(),
            ),
            hpo,
        );

        tracing::debug!("gene = {:?}", gene);
        tracing::debug!("score = {:?}", score);

        let lower_bound = res.scores[..].partition_point(|x| *x < score);
        let upper_bound = res.scores[..].partition_point(|x| *x <= score);
        tracing::debug!("scores = {:?}", res.scores);
        let idx = (lower_bound + upper_bound) / 2;
        let idx = std::cmp::min(idx, res.scores.len() - 1);
        let p = 1.0 - (idx as f64) / (res.scores.len() as f64);
        let log_p = -10.0 * p.log10();

        tracing::debug!("gene = {:?}", gene);
        tracing::debug!("score = {:?}", score);
        tracing::debug!("p = {:?}", p);

        // For each term in the gene, provide query term with the highest similarity.
        let mut terms = HpoGroup::from_iter(
            gene.to_hpo_set(hpo)
                .child_nodes()
                .without_modifier()
                .into_iter(),
        )
        .iter()
        .map(|gene_term_id| {
            let gene_term = hpo.hpo(gene_term_id).expect("gene HPO term not found");
            let (best_term, best_score) = patient
                .iter()
                .map(|query_term_id| {
                    let query_term = hpo.hpo(query_term_id).expect("query HPO term not found");
                    let score = gene_term.similarity_score(
                        &query_term,
                        &Builtins::Resnik(hpo::term::InformationContentKind::Gene),
                    );
                    (query_term, score)
                })
                .max_by(|(_, score1), (_, score2)| score1.partial_cmp(score2).unwrap())
                .expect("could not determine best query term");

            let term_query = if best_score > 0.0 {
                Some(HpoTerm {
                    term_id: best_term.id().to_string(),
                    term_name: Some(best_term.name().to_string()),
                })
            } else {
                None
            };

            TermDetails {
                term_query,
                term_gene: HpoTerm {
                    term_id: gene_term.id().to_string(),
                    term_name: Some(gene_term.name().to_string()),
                },
                score: best_score,
            }
        })
        .collect::<Vec<_>>();
        terms.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

        result.result.push(query_result::Record {
            gene_symbol: gene.name().to_string(),
            p_value: p as f32,
            score: log_p as f32,
            terms: Some(terms),
        });
    }

    result
        .result
        .sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    Ok(result)
}

/// Main entry point for `server pheno-cli` sub command.
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
    let hpo = Ontology::from_standard(&args.path_hpo_dir)?;
    info!("...done loading HPO in {:?}", before_loading.elapsed());

    info!("Opening RocksDB for reading...");
    let before_rocksdb = Instant::now();
    let path_rocksdb = format!("{}/resnik", args.path_hpo_dir);
    let db = rocksdb::DB::open_cf_for_read_only(
        &rocksdb::Options::default(),
        &path_rocksdb,
        ["meta", "resnik_pvalues"],
        true,
    )?;
    info!("...done opening RocksDB in {:?}", before_rocksdb.elapsed());

    info!("Loading genes...");
    let before_load_genes = Instant::now();
    let genes_json = std::fs::read_to_string(&args.path_genes_json)?;
    let genes: Vec<Gene> = serde_json::from_str(&genes_json)?;
    let mut missing_genes = Vec::new();
    let genes = genes
        .iter()
        .filter_map(|g| {
            let mapped = hpo.gene_by_name(&g.gene_symbol);
            if mapped.is_none() {
                missing_genes.push(g.clone());
            }
            mapped
        })
        .collect::<Vec<_>>();
    info!("... done loadin genes in {:?}", before_load_genes.elapsed());

    info!("Loading (patient/query) HPO term ids...");
    let before_load_genes = Instant::now();
    let query_json = std::fs::read_to_string(&args.path_terms_json)?;
    let query: Vec<HpoTerm> = serde_json::from_str(&query_json)?;
    let query = query
        .iter()
        .map(|t| {
            HpoTermId::try_from(t.term_id.as_str())
                .unwrap_or_else(|_| panic!("term {} no valid HPO term ID", &t.term_id))
        })
        .collect::<Vec<_>>();
    let query = {
        let mut group = HpoGroup::new();
        for term in query {
            group.insert(term);
        }
        group
    };
    info!(
        "... done loading HPO IDs in {:?}",
        before_load_genes.elapsed()
    );

    trace_rss_now();

    info!("Starting priorization...");
    let before_priorization = Instant::now();
    let result = run_query(&query, &genes, &hpo, &db)?;
    info!(
        "... done with prioritization in {:?}",
        before_priorization.elapsed()
    );

    println!("{:#?}", result);

    info!(
        "{: >4} | {: <10} | {: >10} | {: >10}",
        "rank", "gene", "P-value", "score"
    );
    info!("     |            |            |");
    for (i, gene) in result.result.iter().enumerate() {
        info!(
            "{: >4} | {: <10} | {: >10.5} | {: >10.2}",
            i + 1,
            gene.gene_symbol,
            gene.p_value,
            gene.score
        );
    }

    info!("All done. Have a nice day!");
    Ok(())
}
