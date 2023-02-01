//! Code implementing the "sv query" sub command.

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
use csv::QuoteStyle;
use indexmap::IndexMap;
use log::warn;
use serde::Serialize;
use thousands::Separable;
use tracing::{error, info};
use uuid::Uuid;

use crate::{
    common::{build_chrom_map, open_maybe_gz, trace_rss_now},
    sv::{
        conf::{sanity_check_db, DbConf},
        query::{
            bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db,
            interpreter::QueryInterpreter, pathogenic::load_patho_dbs,
            pathogenic::Record as KnownPathogenicRecord, records::StructuralVariant as RecordSv,
            schema::CaseQuery, schema::StructuralVariant as SchemaSv, tads::load_tads,
        },
    },
};

use self::{
    bgdbs::{BgDbBundle, BgDbOverlaps},
    clinvar::ClinvarSv,
    genes::GeneDb,
    pathogenic::PathoDbBundle,
    schema::{CallInfo, Database, StrandOrientation, SvSubType, SvType},
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
    /// Path to the output TSV file.
    #[arg(long, required = true)]
    pub path_output_svs: String,
    /// Value to write to the result set id column.
    #[arg(long, default_value_t = 0)]
    pub result_set_id: u64,
    /// Disable checksum verifiation.
    #[arg(long, default_value_t = false)]
    pub disable_checksums: bool,
    /// Slack to use around BND.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: u32,
    /// Slack to use around INS.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: u32,
    /// Optional maximal number of total records to write out.
    #[arg(long)]
    pub max_results: Option<usize>,
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

/// Gene information.
#[derive(Debug, Default, Serialize)]
struct Gene {
    /// Gene symbol
    symbol: Option<String>,
    /// ENSEMBL gene ID
    ensembl_id: Option<String>,
    /// Entrez gene ID
    entrez_id: Option<u32>,
    /// Whether the gene is in the ACMG list for incidental findings.
    is_acmg: bool,
    /// Whether the gene is linked to an OMIM disease.
    is_disease_gene: bool,
}

/// The structured result information of the result record.
#[derive(Debug, Default, Serialize)]
struct ResultPayload {
    /// The overlapping VCVs
    clinvar_ovl_vcvs: Vec<String>,
    /// The directly overlapping genes.
    ovl_genes: Vec<Gene>,
    /// Genes that are not directly overlapping but contained in overlapping TADs.
    tad_genes: Vec<Gene>,
    /// Overlapping known pathogenic SV records.
    known_pathogenic: Vec<KnownPathogenicRecord>,
    /// Information about the call support from the structural variant.
    call_info: IndexMap<String, CallInfo>,
    /// Whether there is an overlap with a disease gene in the overlap.
    ovl_disease_gene: bool,
    /// Whether there is an overlap with a disease gene in the overlapping TADs.
    tad_disease_gene: bool,
    /// The size of the SV, None for ins and BND
    sv_length: Option<u32>,
    /// Overlap counts with background databases.
    overlap_counts: BgDbOverlaps,
    /// Distance to next TAD boundary.
    tad_boundary_distance: Option<u32>,
}

/// A result record from the query.
#[derive(Debug, Default, Serialize)]
struct ResultRecord {
    sodar_uuid: Uuid,
    svqueryresultset: u64,
    release: String,
    chromosome: String,
    chromosome_no: u32,
    bin: u32,
    chromosome2: String,
    chromosome_no2: u32,
    bin2: u32,
    start: u32,
    end: u32,
    pe_orientation: StrandOrientation,
    sv_type: SvType,
    sv_sub_type: SvSubType,
    payload: String,
}

fn resolve_gene_id(database: Database, gene_db: &GeneDb, gene_id: u32) -> Vec<Gene> {
    let record_idxs = match database {
        Database::Refseq => gene_db.xlink.from_entrez.get_vec(&gene_id),
        Database::Ensembl => gene_db.xlink.from_ensembl.get_vec(&gene_id),
    };
    if let Some(record_idxs) = record_idxs {
        record_idxs
            .iter()
            .map(|record_idx| {
                let record = &gene_db.xlink.records[*record_idx as usize];
                Gene {
                    symbol: Some(record.symbol.clone()),
                    ensembl_id: Some(format!("ENSG{:011}", record.ensembl_gene_id)),
                    entrez_id: Some(record.entrez_id),
                    is_acmg: gene_db.acmg.contains(record.entrez_id),
                    is_disease_gene: gene_db.omim.contains(record.entrez_id),
                }
            })
            .collect()
    } else {
        match database {
            Database::Refseq => vec![Gene {
                entrez_id: Some(gene_id),
                symbol: None,
                ensembl_id: None,
                is_acmg: gene_db.acmg.contains(gene_id),
                is_disease_gene: gene_db.omim.contains(gene_id),
            }],
            Database::Ensembl => vec![Gene {
                ensembl_id: Some(format!("ENSG{gene_id:011}")),
                symbol: None,
                entrez_id: None,
                is_acmg: false,
                is_disease_gene: false,
            }],
        }
    }
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

    // Construct reader and writer for CSV records
    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.path_input_svs)?);
    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .quote_style(QuoteStyle::Never)
        .from_path(&args.path_output_svs)?;

    // Read through input records using the query interpreter as a filter
    for record in csv_reader.deserialize() {
        stats.count_total += 1;
        let record_sv: RecordSv = record?;
        let schema_sv: SchemaSv = record_sv.clone().into();

        let mut result_payload = ResultPayload {
            call_info: schema_sv.call_info.clone(),
            ..ResultPayload::default()
        };

        let is_pass = interpreter.passes(&schema_sv, &mut |sv: &SchemaSv| {
            result_payload.overlap_counts = dbs.bg_dbs.count_overlaps(
                sv,
                &interpreter.query,
                &chrom_map,
                args.slack_ins,
                args.slack_bnd,
            );
            result_payload.overlap_counts.clone()
        })?;

        if is_pass {
            if schema_sv.sv_type != SvType::Ins && schema_sv.sv_type != SvType::Bnd {
                result_payload.sv_length = Some(schema_sv.end - schema_sv.pos + 1);
            }

            // Count passing record in statistics
            stats.count_passed += 1;
            *stats.by_sv_type.entry(schema_sv.sv_type).or_default() += 1;

            // Get overlaps with known pathogenic SVs and ClinVar SVs
            result_payload.known_pathogenic = dbs.patho_dbs.overlapping_records(
                &schema_sv,
                &chrom_map,
                interpreter.query.known_pathogenic_min_overlap,
            );
            result_payload.clinvar_ovl_vcvs = dbs
                .clinvar_sv
                .overlapping_vcvs(
                    &schema_sv,
                    &chrom_map,
                    interpreter.query.clinvar_sv_min_pathogenicity,
                    interpreter.query.clinvar_sv_min_overlap,
                )
                .into_iter()
                .map(|vcv| format!("VCV{vcv:09}"))
                .collect();

            // Get overlapping genes and genes in overlapping TADs
            let ovl_gene_ids = {
                let mut ovl_gene_ids = dbs.genes.overlapping_gene_ids(
                    interpreter.query.database,
                    *chrom_map
                        .get(&schema_sv.chrom)
                        .expect("cannot translate chromosome") as u32,
                    schema_sv.pos.saturating_sub(1)..schema_sv.end,
                );
                ovl_gene_ids.sort();
                ovl_gene_ids
            };
            let tad_gene_ids = {
                let gene_ids: HashSet<_> = HashSet::from_iter(ovl_gene_ids.iter());
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
                tad_gene_ids.retain(|gene_id| !gene_ids.contains(gene_id));
                tad_gene_ids.sort();
                tad_gene_ids
            };
            result_payload.tad_boundary_distance =
                dbs.tad_sets
                    .boundary_dist(TadSetChoice::Hesc, &schema_sv, &chrom_map);

            // Convert the genes into more verbose records and put them into the result
            ovl_gene_ids.iter().for_each(|gene_id| {
                result_payload.ovl_genes.append(&mut resolve_gene_id(
                    interpreter.query.database,
                    &dbs.genes,
                    *gene_id,
                ))
            });
            result_payload.ovl_disease_gene = result_payload
                .ovl_genes
                .iter()
                .any(|gene| gene.is_disease_gene);
            tad_gene_ids.iter().for_each(|gene_id| {
                result_payload.tad_genes.append(&mut resolve_gene_id(
                    interpreter.query.database,
                    &dbs.genes,
                    *gene_id,
                ))
            });
            result_payload.tad_disease_gene = result_payload
                .tad_genes
                .iter()
                .any(|gene| gene.is_disease_gene);

            if let Some(max_results) = args.max_results {
                if stats.count_total > max_results {
                    warn!(
                        "stopping writing {} records but there are more results!",
                        stats.count_total
                    );
                }
            }

            // Finally, write out the record.
            csv_writer.serialize(&ResultRecord {
                sodar_uuid: uuid::Uuid::new_v4(),
                svqueryresultset: args.result_set_id,
                release: record_sv.release.clone(),
                chromosome: record_sv.chromosome.clone(),
                chromosome_no: record_sv.chromosome_no,
                start: record_sv.start,
                bin: record_sv.bin,
                chromosome2: record_sv
                    .chromosome2
                    .unwrap_or(record_sv.chromosome)
                    .clone(),
                chromosome_no2: record_sv.chromosome_no2,
                bin2: record_sv.bin2,
                end: record_sv.end,
                pe_orientation: record_sv.pe_orientation,
                sv_type: record_sv.sv_type,
                sv_sub_type: record_sv.sv_sub_type,
                payload: serde_json::to_string(&result_payload)?,
            })?;
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
