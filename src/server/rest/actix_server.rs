//! Implementation of the actix server.

use std::str::FromStr;

use actix_web::{
    get,
    middleware::Logger,
    web::{self, Data, Json, Path},
    App, HttpServer, Responder, ResponseError,
};
use serde::{Deserialize, Serialize};
use thousands::Separable;

use crate::{
    common::CHROMS,
    db::conf::TadSet as TadSetChoice,
    sv::query::{
        bgdbs::BgDbType,
        schema::{Pathogenicity, SvType, VariationType},
    },
};
use crate::{db::conf::GenomeRelease, sv::query::records::ChromRange};

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

/// Chromosomal range for parsing from REST.
///
/// We allow floats here as igv.js likes to put floats for coordinates.
#[derive(Deserialize, Debug)]
struct QueryChromRange {
    chromosome: String,
    begin: f64,
    end: f64,
}

impl From<QueryChromRange> for ChromRange {
    fn from(val: QueryChromRange) -> Self {
        ChromRange {
            chromosome: val.chromosome,
            begin: val.begin as i32,
            end: val.end as i32,
        }
    }
}

/// Parameters for `fetch_tads`.
#[derive(Deserialize, Debug, Clone)]
struct FetchTadsRequest {
    chromosome: String,
    begin: f64,
    end: f64,
    /// A padding of 1 forces igv.js to display book-ended domains to be
    /// displayed on different vertical positions.
    padding: Option<i32>,
}

impl From<FetchTadsRequest> for ChromRange {
    fn from(val: FetchTadsRequest) -> Self {
        ChromRange {
            chromosome: val.chromosome,
            begin: val.begin as i32,
            end: val.end as i32,
        }
    }
}

/// Result type of ""/public/svs/tads/{release}/{tad_set}/".
#[derive(Serialize, Debug)]
struct Tad {
    chromosome: String,
    begin: i32,
    end: i32,
}

/// List the overlapping TADs of the given TAD set.
#[get("/public/svs/tads/{release}/{tad_set}/")]
async fn fetch_tads(
    data: Data<WebServerData>,
    path: Path<(String, String)>,
    query: web::Query<FetchTadsRequest>,
) -> actix_web::Result<impl Responder, CustomError> {
    let genome_release =
        GenomeRelease::from_str(&path.0).map_err(|e| CustomError::new(e.into()))?;
    let tad_set = TadSetChoice::from_str(&path.1).map_err(|e| CustomError::new(e.into()))?;
    let chrom_range: ChromRange = query.clone().into_inner().into();
    let tads = data.dbs[genome_release]
        .tad_sets
        .fetch_tads(tad_set, &chrom_range, &data.chrom_map)
        .into_iter()
        .map(|record| Tad {
            chromosome: CHROMS[record.chrom_no as usize].to_string(),
            begin: record.begin - query.padding.unwrap_or_default(),
            end: record.end + query.padding.unwrap_or_default(),
        })
        .collect::<Vec<Tad>>();
    Ok(Json(tads))
}

/// Result type of "/public/svs/pathogenic/{release}/".
#[derive(Serialize, Debug)]
struct KnownPathogenic {
    chromosome: String,
    begin: i32,
    end: i32,
    sv_type: SvType,
    id: String,
}

/// List the overlapping TADs of the given TAD set.
#[get("/public/svs/pathogenic/{release}/")]
async fn fetch_pathogenic(
    data: Data<WebServerData>,
    path: Path<(String,)>,
    chrom_range: web::Query<QueryChromRange>,
) -> actix_web::Result<impl Responder, CustomError> {
    let genome_release =
        GenomeRelease::from_str(&path.0).map_err(|e| CustomError::new(e.into()))?;
    let chrom_range: ChromRange = chrom_range.into_inner().into();
    let known_pathogenics = data.dbs[genome_release]
        .patho_dbs
        .fetch_records(&chrom_range, &data.chrom_map)
        .into_iter()
        .map(|record| KnownPathogenic {
            chromosome: chrom_range.chromosome.clone(),
            begin: record.begin,
            end: record.end,
            sv_type: record.sv_type,
            id: record.id,
        })
        .collect::<Vec<KnownPathogenic>>();
    Ok(Json(known_pathogenics))
}

/// Result type of "/public/svs/clinvar/{release}/".
#[derive(Serialize, Debug)]
struct ClinvarEntry {
    chromosome: String,
    begin: i32,
    end: i32,
    variation_type: VariationType,
    pathogenicity: Pathogenicity,
    vcv: String,
    name: String,
}

/// We allow begin and end to be floats as igv.js likes to put floats.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub struct ClinVarSvQuery {
    /// Chromosome name.
    pub chromosome: String,
    /// 0-based begin position.
    pub begin: f64,
    /// 0-based end position.
    pub end: f64,
    /// Minimal pathogenicity.
    pub min_pathogenicity: Option<Pathogenicity>,
}

/// List the overlapping TADs of the given TAD set.
#[get("/public/svs/clinvar/{release}/")]
async fn fetch_clinvar_sv(
    data: Data<WebServerData>,
    path: Path<(String,)>,
    query: web::Query<ClinVarSvQuery>,
) -> actix_web::Result<impl Responder, CustomError> {
    let genome_release =
        GenomeRelease::from_str(&path.0).map_err(|e| CustomError::new(e.into()))?;
    let chrom_range = ChromRange {
        chromosome: query.chromosome.clone(),
        begin: query.begin as i32,
        end: query.end as i32,
    };
    let clinvar_entries = data.dbs[genome_release]
        .clinvar_sv
        .fetch_records(&chrom_range, &data.chrom_map, query.min_pathogenicity)
        .into_iter()
        .map(|record| ClinvarEntry {
            chromosome: chrom_range.chromosome.clone(),
            begin: record.start - 1,
            end: record.stop,
            variation_type: crate::sv::query::clinvar::pbs::VariationType::from_i32(
                record.variation_type,
            )
            .expect("unknown variation type")
            .try_into()
            .expect("problem with variation type"),
            pathogenicity: crate::sv::query::clinvar::pbs::Pathogenicity::from_i32(
                record.pathogenicity,
            )
            .expect("unknown pathogenicity")
            .try_into()
            .expect("problem with pathogenicity"),
            vcv: format!("VCV{:09}", record.vcv),
            name: format!(
                "{:?} @ {}:{}-{} ({})",
                record.variation_type,
                &chrom_range.chromosome,
                record.start.separate_with_commas(),
                record.stop.separate_with_commas(),
                record.pathogenicity
            ),
        })
        .collect::<Vec<ClinvarEntry>>();
    Ok(Json(clinvar_entries))
}

#[derive(Serialize)]
struct BgdbResponseRecord {
    pub chromosome: String,
    pub begin: i32,
    pub end: i32,
    pub sv_type: SvType,
    pub count: u32,
    pub name: String,
}

/// List the overlapping background database entries.
#[get("/public/svs/bgdb/{release}/{database}/")]
async fn fetch_bgdb_records(
    data: Data<WebServerData>,
    path: Path<(String, String)>,
    query: web::Query<QueryChromRange>,
) -> actix_web::Result<impl Responder, CustomError> {
    let genome_release =
        GenomeRelease::from_str(&path.0).map_err(|e| CustomError::new(e.into()))?;
    let database = BgDbType::from_str(&path.1).map_err(|e| CustomError::new(e.into()))?;
    let chrom_range = ChromRange {
        chromosome: query.chromosome.clone(),
        begin: query.begin as i32,
        end: query.end as i32,
    };
    let records: Vec<_> = data.dbs[genome_release]
        .bg_dbs
        .fetch_records(&chrom_range, &data.chrom_map, database)
        .into_iter()
        .map(|record| BgdbResponseRecord {
            chromosome: query.chromosome.clone(),
            begin: record.begin,
            end: record.end,
            sv_type: record.sv_type,
            count: record.count,
            name: format!(
                "{:?} @ {}:{}-{} (carriers: {})",
                record.sv_type,
                &query.chromosome,
                (record.begin + 1).separate_with_commas(),
                record.end.separate_with_commas(),
                record.count
            ),
        })
        .collect();
    Ok(Json(records))
}

#[actix_web::main]
pub async fn main(args: &Args, dbs: Data<WebServerData>) -> std::io::Result<()> {
    let tracks_dir = std::path::Path::new(&args.path_db.clone())
        .join("public")
        .join("tracks");

    HttpServer::new(move || {
        let app = App::new()
            .app_data(dbs.clone())
            .service(fetch_tads)
            .service(fetch_pathogenic)
            .service(fetch_clinvar_sv)
            .service(fetch_bgdb_records);
        let app = app.service(actix_files::Files::new("/public/tracks", &tracks_dir));
        app.wrap(Logger::default())
    })
    .bind((args.listen_host.as_str(), args.listen_port))?
    .run()
    .await
}
