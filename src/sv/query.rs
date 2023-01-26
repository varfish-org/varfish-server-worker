pub mod bgdbs;
pub mod genes;
pub mod interpreter;
pub mod pathogenic;
pub mod records;
pub mod schema;
pub mod tads;

use std::{
    collections::BTreeMap,
    fs::File,
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::anyhow;
use clap::{command, Parser};
use thousands::Separable;
use tracing::{error, info, warn};

use crate::{
    common::{build_chrom_map, open_maybe_gz, trace_rss_now},
    sv::{
        conf::{sanity_check_db, DbConf},
        query::{
            bgdbs::load_bg_dbs, genes::load_gene_regions, interpreter::QueryInterpreter,
            pathogenic::load_patho_dbs, records::StructuralVariant as RecordSv, schema::CaseQuery,
            schema::StructuralVariant as SchemaSv, tads::load_tads,
        },
    },
};

use self::{bgdbs::BgDbBundle, pathogenic::PathoDbBundle, schema::SvType};

/// Command line arguments for `sv query` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for SVs", long_about = None)]
pub struct Args {
    /// Path to database to use for querying.
    #[arg(long, required = true)]
    pub path_db: String,
    /// Path to configuration file, defaults to `${path_db}/conf.toml`.
    #[arg(long)]
    pub path_conf: Option<String>,
    /// Path to query JSON file.
    #[arg(long, required = true)]
    pub path_query_json: String,
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub path_input_svs: String,
    /// Disable checksum verifiation.
    #[arg(long, default_value_t = false)]
    pub disable_checksums: bool,
    /// Slack to use around BND.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: u32,
    /// Slack to use around INS.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: u32,
}

/// Load database configuration and perform sanity checks as configured.
fn load_db_conf(args: &Args) -> Result<DbConf, anyhow::Error> {
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

/// Utility struct to store statistics about counts.
#[derive(Debug, Default)]
struct QueryStats {
    pub count_passed: usize,
    pub count_total: usize,
    pub by_sv_type: BTreeMap<SvType, usize>,
}

/// Open the SV file at `path_sv_tsv` and run through the given `interpreter`.
fn run_query(
    interpreter: &QueryInterpreter,
    bg_dbs: &BgDbBundle,
    patho_dbs: &PathoDbBundle,
    args: &Args,
) -> Result<QueryStats, anyhow::Error> {
    let chrom_map = build_chrom_map();
    let mut stats = QueryStats::default();

    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.path_input_svs)?);
    for record in csv_reader.deserialize() {
        stats.count_total += 1;
        let record_sv: RecordSv = record?;
        let schema_sv: SchemaSv = record_sv.into();

        if interpreter.passes(&schema_sv, |sv: &SchemaSv| {
            bg_dbs.count_overlaps(
                sv,
                &interpreter.query,
                &chrom_map,
                args.slack_ins,
                args.slack_bnd,
            )
        })? {
            stats.count_passed += 1;
            *stats.by_sv_type.entry(schema_sv.sv_type).or_default() += 1;
            if patho_dbs.count_overlaps(&schema_sv, &chrom_map) > 0 {
                warn!("found overlap with pathogenic {:?}", &schema_sv);
            }
        }
    }

    Ok(stats)
}

/// Main entry point for `sv query` sub command.
pub(crate) fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = Instant::now();
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    info!("Loading databases...");
    let db_conf = load_db_conf(args)?;
    let before_loading = Instant::now();
    let bg_dbs = load_bg_dbs(&args.path_db, &db_conf.background_dbs)?;
    let patho_dbs = load_patho_dbs(&args.path_db, &db_conf.known_pathogenic)?;
    let _tad_sets = load_tads(&args.path_db, &db_conf.tads)?;
    let _gene_regions = load_gene_regions(&args.path_db, &db_conf.genes)?;
    info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    info!("Loading query...");
    let query: CaseQuery = serde_json::from_reader(File::open(&args.path_query_json)?)?;
    info!(
        "... done loading query = {}",
        &serde_json::to_string(&query)?
    );

    info!("Running queries...");
    let before_query = Instant::now();
    let query_stats = run_query(&QueryInterpreter::new(query), &bg_dbs, &patho_dbs, args)?;
    info!("... done running query in {:?}", before_query.elapsed());
    info!(
        "summary: {} records passed out of {}",
        query_stats.count_passed.separate_with_commas(),
        query_stats.count_total.separate_with_commas()
    );
    info!("passing records by SV type");
    for (sv_type, count) in query_stats.by_sv_type.iter() {
        info!("{:?} -- {}", sv_type, count);
    }

    trace_rss_now();

    info!(
        "All of `sv query` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}
