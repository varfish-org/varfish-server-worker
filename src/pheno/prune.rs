//! Pruning of `phenotype_to_genes.txt` for recent HPO releases.
//!
//! Also see here: https://github.com/anergictcell/hpo/issues/44

use std::collections::HashSet;
use std::io::{BufRead, Write};
use std::time::Instant;

use clap::Parser;
use hpo::Ontology;

/// Command line arguments for `pheno prune` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Prune phenotype to gene association file", long_about = None)]
pub struct Args {
    /// Path to the directory with the HPO files.
    #[arg(long, required = true)]
    pub path_hpo_dir: String,

    /// Path to output phenotype_to_genes.txt file.
    #[arg(long)]
    pub path_out_phenotype_to_genes: String,
}

/// Build pruning list.
///
/// Returns set of (HPO term ID, gene symbol) pairs.
fn build_pruning_list(hpo: &Ontology) -> Result<HashSet<(String, String)>, anyhow::Error> {
    let mut result = HashSet::new();

    for gene in hpo.genes() {
        let mut to_prune_for_gene = HashSet::new();
        for hpo_term_id in gene.hpo_terms() {
            let hpo_term = hpo.hpo(hpo_term_id).unwrap();
            for parent in hpo_term.parents() {
                if parent.id() != hpo_term_id {
                    to_prune_for_gene.insert(parent.id());
                }
            }
        }

        for entry in to_prune_for_gene {
            result.insert((entry.to_string(), gene.symbol().to_string()));
        }
    }

    Ok(result)
}

/// Main entry point for `server pheno-cli` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    if let Some(level) = args_common.verbose.log_level() {
        match level {
            log::Level::Trace | log::Level::Debug => {
                std::env::set_var("RUST_LOG", "debug");
                env_logger::init_from_env(env_logger::Env::new().default_filter_or("info"));
            }
            _ => (),
        }
    }

    tracing::info!("Loading HPO...");
    let before_loading = Instant::now();
    let hpo = Ontology::from_standard(&args.path_hpo_dir)?;
    tracing::info!("...done loading HPO in {:?}", before_loading.elapsed());

    tracing::info!("Building pruning list...");
    let before_list = Instant::now();
    let pruning_list = build_pruning_list(&hpo)?;
    tracing::info!(
        "...done building pruning list in {:?}",
        before_list.elapsed()
    );

    tracing::info!("Writing pruned phenotype_to_genes.txt...");
    let before_writing = Instant::now();
    let mut writer =
        std::io::BufWriter::new(std::fs::File::create(&args.path_out_phenotype_to_genes)?);
    let reader = std::io::BufReader::new(std::fs::File::open(&format!(
        "{}/phenotype_to_genes.txt",
        args.path_hpo_dir
    ))?);
    for line in reader.lines() {
        let line = line?;
        let fields = line.split('\t').into_iter().collect::<Vec<_>>();
        let hpo_term_id = fields[0];
        let gene_symbol = fields[3];
        if !pruning_list.contains(&(hpo_term_id.to_string(), gene_symbol.to_string())) {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }
    }
    tracing::info!(
        "...done writing pruned phenotype_to_genes.txt in {:?}",
        before_writing.elapsed()
    );

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
