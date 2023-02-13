//! Code supporting the `server rest` sub command.

use std::{collections::HashMap, path::PathBuf, str::FromStr, time::Instant};

use actix_web::web::Data;
use clap::Parser;
use enum_map::{enum_map, EnumMap};
use tracing::info;

use crate::{
    common::{build_chrom_map, trace_rss_now},
    db::conf::{GenomeRelease, Top},
    sv::query::{
        bgdbs::load_bg_dbs, clinvar::load_clinvar_sv, genes::load_gene_db, masked::load_masked_dbs,
        pathogenic::load_patho_dbs, tads::load_tads, Databases,
    },
};

pub struct WebServerData {
    pub chrom_map: HashMap<String, usize>,
    pub dbs: EnumMap<GenomeRelease, Databases>,
}

/// Implementation of the actix server.
pub mod actix_server {
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
                begin: val.begin as u32,
                end: val.end as u32,
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
        padding: Option<u32>,
    }

    impl From<FetchTadsRequest> for ChromRange {
        fn from(val: FetchTadsRequest) -> Self {
            ChromRange {
                chromosome: val.chromosome,
                begin: val.begin as u32,
                end: val.end as u32,
            }
        }
    }

    /// Result type of ""/public/svs/tads/{release}/{tad_set}/".
    #[derive(Serialize, Debug)]
    struct Tad {
        chromosome: String,
        begin: u32,
        end: u32,
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
                begin: record
                    .begin
                    .saturating_sub(query.padding.unwrap_or_default()),
                end: record.end.saturating_add(query.padding.unwrap_or_default()),
            })
            .collect::<Vec<Tad>>();
        Ok(Json(tads))
    }

    /// Result type of "/public/svs/pathogenic/{release}/".
    #[derive(Serialize, Debug)]
    struct KnownPathogenic {
        chromosome: String,
        begin: u32,
        end: u32,
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
        begin: u32,
        end: u32,
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
            begin: query.begin as u32,
            end: query.end as u32,
        };
        let clinvar_entries = data.dbs[genome_release]
            .clinvar_sv
            .fetch_records(&chrom_range, &data.chrom_map, query.min_pathogenicity)
            .into_iter()
            .map(|record| ClinvarEntry {
                chromosome: chrom_range.chromosome.clone(),
                begin: record.begin,
                end: record.end,
                variation_type: record.variation_type,
                pathogenicity: record.pathogenicity,
                vcv: format!("VCV{:09}", record.vcv),
                name: format!(
                    "{:?} @ {}:{}-{} ({})",
                    record.variation_type,
                    &chrom_range.chromosome,
                    (record.begin + 1).separate_with_commas(),
                    record.end.separate_with_commas(),
                    record.pathogenicity
                ),
            })
            .collect::<Vec<ClinvarEntry>>();
        Ok(Json(clinvar_entries))
    }

    #[derive(Serialize)]
    struct BgdbResponseRecord {
        pub chromosome: String,
        pub begin: u32,
        pub end: u32,
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
            begin: query.begin as u32,
            end: query.end as u32,
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
}

/// Command line arguments for `server rest` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run REST API server", long_about = None)]
pub struct Args {
    /// Path to database to use for querying.
    #[arg(long, required = true)]
    pub path_db: String,
    /// Path to configuration file, defaults to `${path_db}/conf.toml`.
    #[arg(long)]
    pub path_conf: Option<String>,
    /// IP to listen on.
    #[arg(long, default_value = "127.0.0.1")]
    pub listen_host: String,
    /// Port to listen on.
    #[arg(long, default_value_t = 8081)]
    pub listen_port: u16,
}

/// Main entry point for `sv query` sub command.
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

    info!("Loading database config...");
    let db_conf: Top = {
        let path_conf = if let Some(path_conf) = &args.path_conf {
            PathBuf::from_str(path_conf)?
        } else {
            PathBuf::from_str(&args.path_db)?.join("conf.toml")
        };
        let toml_str = std::fs::read_to_string(&path_conf)?;
        toml::from_str(&toml_str)?
    };

    info!("Loading databases...");
    let before_loading = Instant::now();
    let dbs = enum_map! {
        GenomeRelease::Grch37 => Databases {
            bg_dbs: load_bg_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            patho_dbs: load_patho_dbs(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
            tad_sets: load_tads(&args.path_db, &db_conf.features[GenomeRelease::Grch37])?,
            masked: load_masked_dbs(&args.path_db, &db_conf.features[GenomeRelease::Grch37])?,
            genes: load_gene_db(&args.path_db, &db_conf.genes, &db_conf.features[GenomeRelease::Grch37])?,
            clinvar_sv: load_clinvar_sv(&args.path_db, &db_conf.vardbs[GenomeRelease::Grch37].strucvar)?,
        },
        GenomeRelease::Grch38 => Default::default(),
    };
    info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    let data = Data::new(WebServerData {
        chrom_map: build_chrom_map(),
        dbs,
    });

    trace_rss_now();

    info!("Launching server ...");
    actix_server::main(args, data)?;

    info!("All done. Have a nice day!");
    Ok(())
}
