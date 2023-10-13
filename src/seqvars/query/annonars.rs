//! Code connecting to annonars RocksDB databases for CADD and dbNSFP.

use std::{path::Path, sync::Arc};

use crate::{common::GenomeRelease, seqvars::ingest::path_component};

pub struct AnnonarsDbs {
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
            cadd_db,
            cadd_meta,
            dbnsfp_db,
            dbnsfp_meta,
        })
    }
}
