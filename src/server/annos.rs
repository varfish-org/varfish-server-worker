//! Code supporting the `server annos` sub command.

use std::time::Instant;

use clap::Parser;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use crate::common::{rocksdb_utils::fetch_meta, trace_rss_now, GenomeRelease};

/// Encode annotation database.
#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    strum::Display,
    strum::EnumString,
    strum::EnumIter,
    serde::Serialize,
    enum_map::Enum,
)]
#[strum(serialize_all = "kebab-case")]
#[serde(rename_all = "kebab-case")]
pub enum AnnoDb {
    /// Other database.
    #[default]
    Other,
    /// CADD annotations.
    Cadd,
    /// dbSNP annotations.
    Dbsnp,
    /// dbNSFP annotations.
    Dbnsfp,
    /// dbscSNV annotations.
    Dbscsnv,
    /// gnomAD mtDNA annotations.
    GnomadMtdna,
    /// gnomAD exomes annotations.
    GnomadExomes,
    /// gnomAD genomes annotations.
    GnomadGenomes,
    /// HelixMtDb annotations.
    Helixmtdb,
    /// UCSC conservation annotations.
    UcscConservation,
}

impl AnnoDb {
    /// Return the expected column family name of the database.
    pub fn cf_name(self) -> &'static str {
        match self {
            AnnoDb::Cadd => "tsv_data",
            AnnoDb::Dbsnp => "dbsnp_data",
            AnnoDb::Dbnsfp => "tsv_data",
            AnnoDb::Dbscsnv => "tsv_data",
            AnnoDb::GnomadMtdna => "gnomad_mtdna_data",
            AnnoDb::GnomadExomes => "gnomad_nuclear_data",
            AnnoDb::GnomadGenomes => "gnomad_nuclear_data",
            AnnoDb::Helixmtdb => "helixmtdb_data",
            AnnoDb::UcscConservation => "ucsc_conservation",
            AnnoDb::Other => panic!("cannot get CF name for 'Other'"),
        }
    }

    /// Return the key for the database version.
    fn db_version_meta(&self) -> Option<&'static str> {
        match self {
            AnnoDb::Cadd => Some("db-version"),
            AnnoDb::Dbsnp => Some("db-version"),
            AnnoDb::Dbnsfp => Some("db-version"),
            AnnoDb::Dbscsnv => Some("db-version"),
            AnnoDb::GnomadMtdna => Some("gnomad-version"),
            AnnoDb::GnomadExomes => Some("gnomad-version"),
            AnnoDb::GnomadGenomes => Some("gnomad-version"),
            AnnoDb::Helixmtdb => None,
            AnnoDb::UcscConservation => None,
            AnnoDb::Other => panic!("cannot get meta version name name for 'Other'"),
        }
    }
}

/// Genome-release specific annotation for each database.
pub type ReleaseAnnos =
    enum_map::EnumMap<AnnoDb, Option<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>>;

/// Database information
#[derive(serde::Serialize, Debug, Clone, Default)]
pub struct DbInfo {
    /// Identifier of the database.
    pub name: AnnoDb,
    /// Version of the database.
    pub db_version: Option<String>,
    /// Version of the builder code.
    pub builder_version: String,
}

/// Fetch database information from the given RocksDB.
fn fetch_db_info(
    db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
    name: AnnoDb,
) -> Result<(GenomeRelease, DbInfo), anyhow::Error> {
    let genome_release: GenomeRelease = fetch_meta(db, "genome-release")?
        .ok_or(anyhow::anyhow!("meta:genome-release not found in data"))?
        .as_str()
        .parse()?;
    let db_version = name
        .db_version_meta()
        .map(|db_version_meta| {
            fetch_meta(db, db_version_meta)?.ok_or(anyhow::anyhow!(
                "meta:{} not found in database",
                db_version_meta
            ))
        })
        .transpose()?;
    let builder_version = fetch_meta(db, "annonars-version")?.ok_or(anyhow::anyhow!(
        "meta:annonars-version not found in database {}",
        db.path().display()
    ))?;
    let db_info = DbInfo {
        name,
        db_version,
        builder_version,
    };
    Ok((genome_release, db_info))
}

/// Data for the web server.
#[derive(Default)]
pub struct WebServerData {
    /// Release-specific annotations for each `GenomeRelease`.
    pub annos: enum_map::EnumMap<GenomeRelease, ReleaseAnnos>,
    /// Version information for each database.
    pub db_infos: enum_map::EnumMap<GenomeRelease, enum_map::EnumMap<AnnoDb, Option<DbInfo>>>,
}

/// Implementation of the actix server.
pub mod actix_server {
    use actix_web::{middleware::Logger, web::Data, App, HttpServer, ResponseError};

    use crate::common::{db_keys, rocksdb_utils::fetch_meta};

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

    /// Function to fetch prost Message from a variant database.
    fn fetch_var_protobuf<T>(
        db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
        cf_name: &str,
        key: db_keys::Var,
    ) -> Result<Option<serde_json::Value>, CustomError>
    where
        T: prost::Message + serde::Serialize + Default,
    {
        let cf_data = db
            .cf_handle(cf_name)
            .unwrap_or_else(|| panic!("unknown column family: {}", cf_name));
        let key: Vec<u8> = key.into();

        let raw_data = db
            .get_cf(&cf_data, key)
            .map_err(|e| CustomError::new(anyhow::anyhow!("problem querying database: {}", e)))?;
        raw_data
            .map(|raw_data| {
                let msg: T = prost::Message::decode(&raw_data[..]).map_err(|e| {
                    CustomError::new(anyhow::anyhow!(
                        "problem decoding protobuf from database (cf_name={}): {}",
                        cf_name,
                        e
                    ))
                })?;
                serde_json::to_value(msg).map_err(|e| {
                    CustomError::new(anyhow::anyhow!("problem decoding JSON from database: {e}",))
                })
            })
            .transpose()
    }

    /// Function to fetch prost Message from a position database.
    fn fetch_pos_protobuf<T>(
        db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
        cf_name: &str,
        start: db_keys::Pos,
        stop: db_keys::Pos,
    ) -> Result<Option<serde_json::Value>, CustomError>
    where
        T: prost::Message + serde::Serialize + Default,
    {
        let stop = annonars::common::keys::Pos {
            chrom: stop.chrom.clone(),
            pos: stop.pos as i32,
        };

        let cf_data = db.cf_handle(cf_name).unwrap();
        let mut iter = db.raw_iterator_cf(&cf_data);
        let start: Vec<u8> = start.into();
        iter.seek(&start);

        let mut result = Vec::new();
        while iter.valid() {
            if let Some(raw_value) = iter.value() {
                let iter_key = iter.key().unwrap();
                let iter_pos: annonars::common::keys::Pos = iter_key.into();

                if iter_pos > stop {
                    break;
                }

                let msg: T = prost::Message::decode(raw_value).map_err(|e| {
                    CustomError::new(anyhow::anyhow!(
                        "problem decoding protobuf from database (cf_name={}): {}",
                        cf_name,
                        e
                    ))
                })?;
                result.push(serde_json::to_value(msg).map_err(|e| {
                    CustomError::new(anyhow::anyhow!("problem decoding JSON from database: {e}",))
                })?);

                iter.next();
            }
        }

        Ok(Some(serde_json::Value::Array(result)))
    }

    /// Function to fetch a annonars::tsv record from a database by variant.
    fn fetch_tsv_json(
        db: &rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>,
        cf_name: &str,
        key: db_keys::Var,
    ) -> Result<Option<serde_json::Value>, CustomError> {
        let db_schema: annonars::tsv::schema::FileSchema = fetch_meta(db, "db-schema")
            .map_err(CustomError::new)?
            .map(|s| {
                serde_json::from_str(&s).map_err(|e| {
                    CustomError::new(anyhow::anyhow!(
                        "problem loading schema from JSON (cf_name={}): {}",
                        cf_name,
                        e
                    ))
                })
            })
            .transpose()?
            .ok_or(CustomError::new(anyhow::anyhow!(
                "db-schema not found in TSV data"
            )))?;
        let infer_config: annonars::tsv::schema::infer::Config = fetch_meta(db, "db-infer-config")
            .map_err(CustomError::new)?
            .map(|s| {
                serde_json::from_str(&s).map_err(|e| {
                    CustomError::new(anyhow::anyhow!(
                        "problem loading inference from JSON: {}",
                        e
                    ))
                })
            })
            .transpose()?
            .ok_or(CustomError::new(anyhow::anyhow!(
                "db-schema not found in TSV data"
            )))?;
        let cf_data = db
            .cf_handle(cf_name)
            .ok_or(CustomError::new(anyhow::anyhow!(
                "TSV data does not have a column family named {}",
                cf_name
            )))?;
        let ctx = annonars::tsv::coding::Context::new(infer_config, db_schema.clone());
        let key: Vec<u8> = key.into();
        let raw_values = db.get_cf(&cf_data, key).map_err(|e| {
            CustomError::new(anyhow::anyhow!(
                "problem querying database (cf_name={}): {}",
                cf_name,
                e
            ))
        })?;
        let values = raw_values
            .map(|v| {
                ctx.decode_values(&v).map_err(|e| {
                    CustomError::new(anyhow::anyhow!(
                        "problem decoding data (cf_name={}): {}",
                        cf_name,
                        e
                    ))
                })
            })
            .transpose()?;

        Ok(values.as_ref().map(|values| {
            let mut m = serde_json::Map::new();
            for (col, value) in db_schema.columns.iter().zip(values.iter()) {
                m.insert(col.name.clone(), value.clone());
            }
            serde_json::Value::Object(m)
        }))
    }

    /// Code for `/annos/variant`.
    pub mod annos_variant {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use serde::{Deserialize, Serialize};
        use strum::IntoEnumIterator;

        use crate::{common::db_keys, pheno::prepare::VERSION, server::annos::AnnoDb};

        use super::{
            fetch_pos_protobuf, fetch_tsv_json, fetch_var_protobuf, CustomError, WebServerData,
        };

        /// Parameters for `variant_annos::handle`.
        ///
        /// Defines a variant in SPDI format with a genome release specification.
        #[serde_with::skip_serializing_none]
        #[serde_with::serde_as]
        #[derive(Serialize, Deserialize, Debug, Clone)]
        #[serde(rename_all = "kebab-case")]
        struct Request {
            /// Genome release specification.
            #[allow(dead_code)]
            pub genome_release: String,
            /// Chromosome name.
            pub chromosome: String,
            /// 1-based position for SPDI.
            pub pos: u32,
            /// Reference allele bases.
            pub reference: String,
            /// Alterantive allele bases.
            pub alternative: String,
        }

        impl From<Request> for db_keys::Var {
            fn from(value: Request) -> Self {
                db_keys::Var {
                    chrom: value.chromosome,
                    pos: value.pos,
                    reference: value.reference,
                    alternative: value.alternative,
                }
            }
        }

        impl From<Request> for db_keys::Pos {
            fn from(value: Request) -> Self {
                db_keys::Pos {
                    chrom: value.chromosome,
                    pos: value.pos,
                }
            }
        }

        /// Result for `handle`.
        #[derive(Serialize, Debug, Clone)]
        #[serde_with::skip_serializing_none]
        struct Container {
            /// Version of the server code.
            pub server_version: String,
            /// The query parameters.
            pub query: Request,
            /// Annotations for the variant from each database.
            pub result: std::collections::BTreeMap<AnnoDb, Option<serde_json::Value>>,
        }

        /// Query for annotations for one variant.
        #[allow(clippy::option_map_unit_fn)]
        #[get("/annos/variant")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let genome_release = query
                .clone()
                .into_inner()
                .genome_release
                .parse()
                .map_err(CustomError::new)?;

            let mut annotations = std::collections::BTreeMap::default();
            for anno_db in AnnoDb::iter() {
                match anno_db {
                    AnnoDb::Other => (),
                    AnnoDb::Cadd | AnnoDb::Dbnsfp | AnnoDb::Dbscsnv => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                fetch_tsv_json(
                                    db,
                                    anno_db.cf_name(),
                                    query.clone().into_inner().into(),
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::Dbsnp => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                fetch_var_protobuf::<annonars::dbsnp::pbs::Record>(
                                    db,
                                    anno_db.cf_name(),
                                    query.clone().into_inner().into(),
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::Helixmtdb => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                fetch_var_protobuf::<annonars::helixmtdb::pbs::Record>(
                                    db,
                                    anno_db.cf_name(),
                                    query.clone().into_inner().into(),
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::GnomadMtdna => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                fetch_var_protobuf::<annonars::gnomad_pbs::mtdna::Record>(
                                    db,
                                    anno_db.cf_name(),
                                    query.clone().into_inner().into(),
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::GnomadExomes => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                fetch_var_protobuf::<annonars::gnomad_pbs::gnomad2::Record>(
                                    db,
                                    anno_db.cf_name(),
                                    query.clone().into_inner().into(),
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::GnomadGenomes => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                let db_version = data.db_infos[genome_release][anno_db]
                                    .as_ref()
                                    .expect("must have db info here")
                                    .db_version
                                    .as_ref()
                                    .expect("gnomAD must have db version");
                                if db_version.starts_with("2.") {
                                    fetch_var_protobuf::<annonars::gnomad_pbs::gnomad2::Record>(
                                        db,
                                        anno_db.cf_name(),
                                        query.clone().into_inner().into(),
                                    )
                                } else if db_version.starts_with("3.") {
                                    fetch_var_protobuf::<annonars::gnomad_pbs::gnomad3::Record>(
                                        db,
                                        anno_db.cf_name(),
                                        query.clone().into_inner().into(),
                                    )
                                } else {
                                    Err(CustomError::new(anyhow::anyhow!(
                                        "don't know how to tread gnomAD version {}",
                                        db_version
                                    )))
                                }
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                    AnnoDb::UcscConservation => {
                        data.annos[genome_release][anno_db]
                            .as_ref()
                            .map(|db| {
                                let start: db_keys::Pos = query.clone().into_inner().into();
                                let start = db_keys::Pos {
                                    chrom: start.chrom,
                                    pos: start.pos - 2,
                                };
                                let stop = query.clone().into_inner().into();
                                fetch_pos_protobuf::<annonars::cons::pbs::RecordList>(
                                    db,
                                    anno_db.cf_name(),
                                    start,
                                    stop,
                                )
                            })
                            .transpose()?
                            .map(|v| annotations.insert(anno_db, v));
                    }
                }
            }

            let result = Container {
                server_version: VERSION.to_string(),
                query: query.into_inner(),
                result: annotations,
            };

            Ok(Json(result))
        }
    }

    /// Code for `/annos/range`.
    pub mod annos_range {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use serde::{Deserialize, Serialize};

        use super::{CustomError, WebServerData};

        /// Parameters for `variant_annos::handle`.
        #[serde_with::skip_serializing_none]
        #[serde_with::serde_as]
        #[derive(Deserialize, Debug, Clone)]
        #[serde(rename_all = "kebab-case")]
        struct Request {
            /// Genome release version.
            pub genome_release: String,
            /// Chromosome name.
            pub chromosome: String,
            /// 1-based start position.
            pub pos: u32,
            /// Reference allele.
            pub reference: String,
            /// Alternative allele.
            pub alternative: String,
        }

        /// Result for `handle`.
        #[derive(Serialize, Debug, Clone)]
        struct Result {
            /// Version of the server code.
            pub server_version: String,
            /// Version of the builder code.
            pub builder_version: String,
        }

        /// Query for annotations for one variant.
        #[get("/annos/range")]
        async fn handle(
            _data: Data<WebServerData>,
            _path: Path<()>,
            _query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            Ok(Json(1))
        }
    }

    /// Code for `/annos/db-info`.
    pub mod annos_db_info {
        use actix_web::{
            get,
            web::{self, Data, Json, Path},
            Responder,
        };
        use serde::{Deserialize, Serialize};

        use crate::server::annos::DbInfo;

        use super::{CustomError, WebServerData};

        /// Parameters for `variant_annos::handle`.
        #[serde_with::skip_serializing_none]
        #[serde_with::serde_as]
        #[derive(Deserialize, Debug, Clone)]
        #[serde(rename_all = "kebab-case")]
        struct Request {
            pub genome_release: String,
        }

        /// Result for `handle`.
        #[derive(Serialize, Debug, Clone)]
        struct ResultEntry {
            /// Information for each database.
            pub db_info: linked_hash_map::LinkedHashMap<String, DbInfo>,
        }

        /// Query for annotations for one variant.
        #[get("/annos/db-info")]
        async fn handle(
            data: Data<WebServerData>,
            _path: Path<()>,
            query: web::Query<Request>,
        ) -> actix_web::Result<impl Responder, CustomError> {
            let genome_release = query
                .into_inner()
                .genome_release
                .parse()
                .map_err(CustomError::new)?;
            Ok(Json(data.db_infos[genome_release].clone()))
        }
    }

    #[actix_web::main]
    pub async fn main(args: &Args, dbs: Data<WebServerData>) -> std::io::Result<()> {
        HttpServer::new(move || {
            let app = App::new()
                .app_data(dbs.clone())
                .service(annos_variant::handle)
                .service(annos_range::handle)
                .service(annos_db_info::handle);
            app.wrap(Logger::default())
        })
        .bind((args.listen_host.as_str(), args.listen_port))?
        .run()
        .await
    }
}

/// Command line arguments for `server rest` sub command.
///
/// Each path can be given more than one time to support multiple releases.  When the server
/// is started, it needs to be given a file for each database with each release.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run variant annotation REST API", long_about = None)]
pub struct Args {
    /// CADD database(s), one for each release.
    #[arg(long)]
    pub path_cadd: Vec<String>,
    /// dbSNP database(s), one for each release.
    #[arg(long)]
    pub path_dbsnp: Vec<String>,
    /// dbNSFP database(s), one for each release.
    #[arg(long)]
    pub path_dbnsfp: Vec<String>,
    /// PdbscSNV database(s), one for each release.
    #[arg(long)]
    pub path_dbscsnv: Vec<String>,
    /// gnomAD mtDNA database(s), one for each release.
    #[arg(long)]
    pub path_gnomad_mtdna: Vec<String>,
    /// gnomAD-exomes database(s), one for each release.
    #[arg(long)]
    pub path_gnomad_exomes: Vec<String>,
    /// gnomAD-genomes database(s), one for each release.
    #[arg(long)]
    pub path_gnomad_genomes: Vec<String>,
    /// HelixMtDB database(s), one for each release.
    #[arg(long)]
    pub path_helixmtdb: Vec<String>,
    /// UCSC conservation database(s), one for each release.
    #[arg(long)]
    pub path_ucsc_conservation: Vec<String>,

    /// IP to listen on.
    #[arg(long, default_value = "127.0.0.1")]
    pub listen_host: String,
    /// Port to listen on.
    #[arg(long, default_value_t = 8081)]
    pub listen_port: u16,
}

/// Open a RocksDB database.
///
/// # Arguments
///
/// * `path` - Path to the database.
/// * `cf_name` - Name of the column family to open (besides the mandatory `meta` column family).
fn open_db(
    path: &str,
    cf_name: &str,
) -> Result<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>, anyhow::Error> {
    tracing::info!("Opening database {}...", path);
    let before_open = Instant::now();
    let res = rocksdb::DB::open_cf_for_read_only(
        &rocksdb::Options::default(),
        path,
        ["meta", cf_name],
        true,
    )
    .map_err(|e| anyhow::anyhow!("problem opening database: {}", e));
    tracing::info!("...done opening database in {:?}", before_open.elapsed());
    res
}

/// Main entry point for `server rest` sub command.
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

    tracing::info!("Opening databases...");
    let before_opening = Instant::now();
    // Argument lists from the command line with the corresponding database enum value.
    let paths_db_pairs = vec![
        (&args.path_cadd, AnnoDb::Cadd),
        (&args.path_dbnsfp, AnnoDb::Dbnsfp),
        (&args.path_dbsnp, AnnoDb::Dbsnp),
        (&args.path_dbscsnv, AnnoDb::Dbscsnv),
        (&args.path_gnomad_mtdna, AnnoDb::GnomadMtdna),
        (&args.path_gnomad_exomes, AnnoDb::GnomadExomes),
        (&args.path_gnomad_genomes, AnnoDb::GnomadGenomes),
        (&args.path_helixmtdb, AnnoDb::Helixmtdb),
        (&args.path_ucsc_conservation, AnnoDb::UcscConservation),
    ];
    // "Unpack" the list of paths to single paths.
    let path_db_pairs = paths_db_pairs
        .iter()
        .map(|(paths, anno_db)| {
            paths
                .iter()
                .map(|path| (path.clone(), *anno_db))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    // Open the corresponding databases in parallel and extract database infos.  Store the
    // resulting database infos in `data`.
    let mut data = WebServerData::default();
    path_db_pairs
        .par_iter()
        .progress_with(annonars::common::cli::progress_bar(path_db_pairs.len()))
        .map(|(path, anno_db)| -> Result<_, anyhow::Error> {
            let db = open_db(path, anno_db.cf_name())?;
            let (genome_release, db_info) = fetch_db_info(&db, *anno_db)?;

            Ok((db_info, genome_release, db))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .for_each(|(db_info, genome_release, db)| {
            let name = db_info.name;
            data.db_infos[genome_release][name] = Some(db_info);
            data.annos[genome_release][name] = Some(db);
        });
    tracing::info!(
        "...done opening databases in {:?}",
        before_opening.elapsed()
    );

    trace_rss_now();

    tracing::info!(
        "Launching server main on http://{}:{} ...",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/annos/db-info?genome-release=grch37",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/annos/variant?genome-release=grch37&chromosome=1&pos=55505599&reference=C&alternative=G",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/annos/variant?genome-release=grch37&chromosome=1&pos=10001&reference=T&alternative=A",
        args.listen_host.as_str(),
        args.listen_port
    );
    tracing::info!(
        "  try: http://{}:{}/annos/range?genome-release=grch37&chromosome=1&start=1&stop=55516888",
        args.listen_host.as_str(),
        args.listen_port
    );
    actix_server::main(args, actix_web::web::Data::new(data))?;

    tracing::info!("All done. Have a nice day!");
    Ok(())
}
