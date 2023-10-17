//! Code connecting to annonars RocksDB databases for CADD and dbNSFP.

use std::{path::Path, sync::Arc};

use prost::Message;

use crate::{common::GenomeRelease, seqvars::ingest::path_component};

use super::schema::SequenceVariant;

/// Bundle the types needed for databases.
struct AnnonarsDbs {
    /// ClinVar database as annonars RocksDB.
    clinvar_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// ClinVar metadata from annonars.
    clinvar_meta: annonars::clinvar_minimal::cli::query::Meta,
    /// dbSNP database as annonars RocksDB.
    dbsnp_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// dbSNP metadata from annonars.
    dbsnp_meta: annonars::dbsnp::cli::query::Meta,
    /// CADD database as annonars RocksDB.
    #[allow(dead_code)]
    cadd_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// CADD metadata from annonars.
    #[allow(dead_code)]
    cadd_meta: annonars::tsv::cli::query::Meta,
    /// dbNSFP database as annonars RocksDB.
    #[allow(dead_code)]
    dbnsfp_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
    /// dbNSFP metadata from annonars.
    #[allow(dead_code)]
    dbnsfp_meta: annonars::tsv::cli::query::Meta,
    /// annonaars gene RocksDB.
    #[allow(dead_code)]
    genes_db: Arc<rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>>,
}

impl AnnonarsDbs {
    /// Initialize from path that contains the annonars databases.
    fn with_path<P: AsRef<Path>>(
        path: P,
        genome_release: GenomeRelease,
    ) -> Result<Self, anyhow::Error> {
        let path_annonars = path.as_ref().join("annonars");
        let path_genome_release = path_annonars.join(path_component(genome_release));

        macro_rules! open_rocksdb {
            ($path_token:expr, $module:ident, $db_name:expr, $meta_name:expr) => {{
                let path: std::path::PathBuf =
                    path_genome_release.join($path_token).join("rocksdb");
                annonars::$module::cli::query::open_rocksdb(&path, $db_name, $meta_name).map_err(
                    |e| {
                        anyhow::anyhow!(
                            "problem opening {} metadata at {}: {}",
                            $db_name,
                            path.as_os_str().to_string_lossy(),
                            e
                        )
                    },
                )?
            }};
        }

        let (clinvar_db, clinvar_meta) =
            open_rocksdb!("clinvar-minimal", clinvar_minimal, "clinvar", "meta");
        let (cadd_db, cadd_meta) = open_rocksdb!("cadd", tsv, "tsv_data", "meta");
        let (dbnsfp_db, dbnsfp_meta) = open_rocksdb!("dbnsfp", tsv, "tsv_data", "meta");
        let (dbsnp_db, dbsnp_meta) = open_rocksdb!("dbsnp", dbsnp, "dbsnp_data", "meta");

        let path_rocksdb = path_annonars.join("genes").join("rocksdb");
        let genes_db = annonars::genes::cli::query::open_rocksdb(&path_rocksdb, "genes", "meta")
            .map_err(|e| {
                anyhow::anyhow!(
                    "problem opening genes metadata at {}: {}",
                    path_rocksdb.as_os_str().to_string_lossy(),
                    e
                )
            })?;

        Ok(Self {
            clinvar_db,
            clinvar_meta,
            dbsnp_db,
            dbsnp_meta,
            cadd_db,
            cadd_meta,
            dbnsfp_db,
            dbnsfp_meta,
            genes_db,
        })
    }
}

/// Utility for sequence variant annotation with annonars.
pub struct Annotator {
    /// Annonars database bundles.
    annonars_dbs: AnnonarsDbs,
}

impl Annotator {
    /// Construct with path to annonars databases.
    ///
    /// # Errors
    ///
    /// If there is a problem opening the databases.
    pub fn with_path<P: AsRef<Path>>(
        path: P,
        genome_release: GenomeRelease,
    ) -> Result<Self, anyhow::Error> {
        let annonars_dbs = AnnonarsDbs::with_path(path.as_ref(), genome_release).map_err(|e| {
            anyhow::anyhow!(
                "problem opening annonars databases at {}: {}",
                path.as_ref().as_os_str().to_string_lossy(),
                e
            )
        })?;
        Ok(Self { annonars_dbs })
    }

    /// Query `clinvar-minimal` database for a given variant.
    ///
    /// # Errors
    ///
    /// If there is a problem querying the database.
    pub fn query_clinvar_minimal(
        &self,
        seqvar: &SequenceVariant,
    ) -> Result<annonars::clinvar_minimal::pbs::Record, anyhow::Error> {
        let cf_data = self
            .annonars_dbs
            .clinvar_db
            .cf_handle("clinvar")
            .ok_or_else(|| anyhow::anyhow!("could not get clinvar column family"))?;

        let variant = annonars::common::spdi::Var::new(
            annonars::common::cli::canonicalize(&seqvar.chrom),
            seqvar.pos,
            seqvar.reference.clone(),
            seqvar.alternative.clone(),
        );

        annonars::clinvar_minimal::cli::query::query_for_variant(
            &variant,
            &self.annonars_dbs.clinvar_meta,
            &self.annonars_dbs.clinvar_db,
            &cf_data,
        )
        .map_err(|e| anyhow::anyhow!("problem querying clinvar-minimal database: {}", e))
    }

    /// Query `dbsnp` database for a given variant.
    ///
    /// # Errors
    ///
    /// If there is a problem querying the database.
    pub fn query_dbsnp(
        &self,
        seqvar: &SequenceVariant,
    ) -> Result<annonars::dbsnp::pbs::Record, anyhow::Error> {
        let cf_data = self
            .annonars_dbs
            .dbsnp_db
            .cf_handle("dbsnp")
            .ok_or_else(|| anyhow::anyhow!("could not get dbsnp column family"))?;

        let variant = annonars::common::spdi::Var::new(
            annonars::common::cli::canonicalize(&seqvar.chrom),
            seqvar.pos,
            seqvar.reference.clone(),
            seqvar.alternative.clone(),
        );

        annonars::dbsnp::cli::query::query_for_variant(
            &variant,
            &self.annonars_dbs.dbsnp_meta,
            &self.annonars_dbs.dbsnp_db,
            &cf_data,
        )
        .map_err(|e| anyhow::anyhow!("problem querying dbsnp database: {}", e))
    }

    /// Query `genes` database for a given HGNC ID.
    ///
    /// # Errors
    ///
    /// If there is a problem querying the database.
    pub fn query_genes(
        &self,
        hgnc_id: &str,
    ) -> Result<Option<annonars::genes::pbs::Record>, anyhow::Error> {
        let cf_data = self
            .annonars_dbs
            .genes_db
            .cf_handle("genes")
            .ok_or_else(|| anyhow::anyhow!("could not get genes column family"))?;

        let raw_value = self
            .annonars_dbs
            .genes_db
            .get_cf(&cf_data, hgnc_id.as_bytes())
            .map_err(|e| {
                anyhow::anyhow!(
                    "problem querying genes database for HGNC ID {}: {}",
                    hgnc_id,
                    e
                )
            })?;

        raw_value
            .map(|raw_value| {
                annonars::genes::pbs::Record::decode(&mut std::io::Cursor::new(&raw_value)).map_err(
                    |e| {
                        anyhow::anyhow!(
                            "problem decoding record from genes database for HGNC ID {}: {}",
                            hgnc_id,
                            e
                        )
                    },
                )
            })
            .transpose()
    }
}
