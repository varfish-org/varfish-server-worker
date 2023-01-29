pub mod bgdbs;
pub mod clinvar;
pub mod genes;
pub mod interpreter;
pub mod pathogenic;
pub mod records;
pub mod schema;
pub mod tads;

use std::{
    collections::{BTreeMap, HashSet},
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
            bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db,
            interpreter::QueryInterpreter, pathogenic::load_patho_dbs,
            records::StructuralVariant as RecordSv, schema::CaseQuery,
            schema::StructuralVariant as SchemaSv, tads::load_tads,
        },
    },
};

use self::{
    bgdbs::BgDbBundle,
    clinvar::ClinvarSv,
    genes::GeneDb,
    pathogenic::PathoDbBundle,
    schema::{Pathogenicity, SvType},
    tads::{TadSetBundle, TadSetChoice},
};

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
    args: &Args,
    dbs: &Databases,
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
            dbs.bg_dbs.count_overlaps(
                sv,
                &interpreter.query,
                &chrom_map,
                args.slack_ins,
                args.slack_bnd,
            )
        })? {
            stats.count_passed += 1;
            *stats.by_sv_type.entry(schema_sv.sv_type).or_default() += 1;
            // perform some more queries for measuring time
            if dbs
                .patho_dbs
                .count_overlaps(&schema_sv, &chrom_map, Some(0.8))
                > 0
            {
                warn!("found overlap with pathogenic {:?}", &schema_sv);
            }
            let vcvs = dbs.clinvar_sv.overlapping_vcvs(
                &schema_sv,
                &chrom_map,
                Some(Pathogenicity::LikelyPathogenic),
                Some(0.8),
            );
            if !vcvs.is_empty() {
                warn!("found overlap with clinvar {:?}", &schema_sv);
                for vcv in vcvs {
                    warn!(" --> VCV{:09}", vcv);
                }
            }

            let gene_ids = {
                let mut gene_ids = dbs.genes.overlapping_gene_ids(
                    interpreter.query.database,
                    *chrom_map
                        .get(&schema_sv.chrom)
                        .expect("cannot translate chromosome") as u32,
                    schema_sv.pos.saturating_sub(1)..schema_sv.end,
                );
                gene_ids.sort();
                gene_ids
            };
            let tad_gene_ids = {
                let tads =
                    dbs.tad_sets
                        .overlapping_tads(TadSetChoice::Hesc, &schema_sv, &chrom_map);
                let mut tad_gene_ids = Vec::new();
                tads.iter()
                    .map(|tad| {
                        dbs.genes.overlapping_gene_ids(
                            interpreter.query.database,
                            tad.chrom_no,
                            tad.begin..tad.end,
                        )
                    })
                    .for_each(|mut v| tad_gene_ids.append(&mut v));
                let tad_gene_ids: HashSet<_> = HashSet::from_iter(tad_gene_ids.into_iter());
                let mut tad_gene_ids = Vec::from_iter(tad_gene_ids);
                tad_gene_ids.sort();
                tad_gene_ids
            };

            gene_ids.iter().for_each(|gene_id| {
                dbs.genes
                    .xlink
                    .gene_id_to_symbols(interpreter.query.database, *gene_id);
            });
            tad_gene_ids.iter().for_each(|gene_id| {
                dbs.genes
                    .xlink
                    .gene_id_to_symbols(interpreter.query.database, *gene_id);
            });
        }
    }

    Ok(stats)
}

/// Bundle the used database to reduce argument count.
pub struct Databases {
    pub bg_dbs: BgDbBundle,
    pub patho_dbs: PathoDbBundle,
    pub tad_sets: TadSetBundle,
    pub genes: GeneDb,
    pub clinvar_sv: ClinvarSv,
}

/// Main entry point for `sv query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = Instant::now();
    info!("args_common = {:?}", &args_common);
    info!("args = {:?}", &args);

    info!("Loading databases...");
    let before_loading = Instant::now();
    let db_conf = load_db_conf(args)?;
    let dbs = Databases {
        bg_dbs: load_bg_dbs(&args.path_db, &db_conf.background_dbs)?,
        patho_dbs: load_patho_dbs(&args.path_db, &db_conf.known_pathogenic)?,
        tad_sets: load_tads(&args.path_db, &db_conf.tads)?,
        genes: load_gene_db(&args.path_db, &db_conf.genes)?,
        clinvar_sv: load_clinvar_sv(&args.path_db, &db_conf.clinvar)?,
    };
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
    let query_stats = run_query(&QueryInterpreter::new(query), args, &dbs)?;
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
