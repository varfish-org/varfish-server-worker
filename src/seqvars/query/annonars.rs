//! Code connecting to annonars RocksDB databases for CADD and dbNSFP.

use std::{path::Path, sync::Arc};

use crate::{common::GenomeRelease, seqvars::ingest::path_component};

pub struct AnnonarsDbs {
    /// ClinVar database as annonars RocksDB.
    pub clinvar_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// ClinVar metadata from annonars.
    pub clinvar_meta: annonars::clinvar_minimal::cli::query::Meta,
    /// CADD database as annonars RocksDB.
    pub cadd_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// CADD metadata from annonars.
    pub cadd_meta: annonars::tsv::cli::query::Meta,
    /// dbNSFP database as annonars RocksDB.
    pub dbnsfp_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// dbNSFP metadata from annonars.
    pub dbnsfp_meta: annonars::tsv::cli::query::Meta,
}

impl AnnonarsDbs {
    /// Initialize from path that contains the annonars databases.
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        genome_release: GenomeRelease,
    ) -> Result<Self, anyhow::Error> {
        let path_base = path
            .as_ref()
            .join("annonars")
            .join(path_component(genome_release));
        let (clinvar_db, clinvar_meta) = {
            let path: std::path::PathBuf = path_base.join("clinvar-minimal").join("rocksdb");
            annonars::clinvar_minimal::cli::query::open_rocksdb(&path, "clinvar", "meta").map_err(
                |e| {
                    anyhow::anyhow!(
                        "problem opening ClinVar database at {}: {}",
                        path.as_os_str().to_string_lossy(),
                        e
                    )
                },
            )?
        };
        let (cadd_db, cadd_meta) = {
            let path: std::path::PathBuf = path_base.join("cadd").join("rocksdb");
            annonars::tsv::cli::query::open_rocksdb(&path, "tsv_data", "meta").map_err(|e| {
                anyhow::anyhow!(
                    "problem opening CADD database at {}: {}",
                    path.as_os_str().to_string_lossy(),
                    e
                )
            })?
        };
        let (dbnsfp_db, dbnsfp_meta) = {
            let path = path_base.join("dbnsfp").join("rocksdb");
            annonars::tsv::cli::query::open_rocksdb(&path, "tsv_data", "meta").map_err(|e| {
                anyhow::anyhow!(
                    "problem opening dbNSFP database at {}: {}",
                    path.as_os_str().to_string_lossy(),
                    e
                )
            })?
        };

        Ok(Self {
            clinvar_db,
            clinvar_meta,
            cadd_db,
            cadd_meta,
            dbnsfp_db,
            dbnsfp_meta,
        })
    }
}
