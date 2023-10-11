//! Code for supporting annotation with overlapping genes.

use std::{collections::HashSet, path::Path, time::Instant};

use mehari::common::open_read_maybe_gz;
use prost::Message;
use serde::Deserialize;
use tracing::info;

use crate::{common::GenomeRelease, strucvars::pbs};

/// Information to store for the interlink table.
#[derive(Default, Debug)]
pub struct XlinkDbRecord {
    pub entrez_id: u32,
    pub ensembl_gene_id: u32,
    pub symbol: String,
    pub hgnc_id: String,
}

/// The interlink DB.
#[derive(Debug, Default)]
pub struct XlinkDb {
    /// The interlink database with symbols.
    pub records: Vec<XlinkDbRecord>,
    /// Link from entrez ID to indices in records.
    pub from_entrez: multimap::MultiMap<u32, u32>,
    /// Link from ensembl ID to indices in records.
    pub from_ensembl: multimap::MultiMap<u32, u32>,
    /// Link from HGNC ID to indices in records.
    pub from_hgnc: multimap::MultiMap<String, u32>,
}

#[tracing::instrument]
fn load_xlink_db(path: &Path) -> Result<XlinkDb, anyhow::Error> {
    tracing::debug!("loading binary xlink records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = XlinkDb::default();

    let fcontents =
        std::fs::read(path).map_err(|e| anyhow::anyhow!("error reading {:?}: {}", &path, e))?;
    let xlink_db = pbs::XlinkDatabase::decode(std::io::Cursor::new(fcontents))
        .map_err(|e| anyhow::anyhow!("error decoding {:?}: {}", &path, e))?;

    let mut total_count = 0;
    for record in xlink_db.records.into_iter() {
        result
            .from_entrez
            .insert(record.entrez_id, result.records.len() as u32);
        result
            .from_ensembl
            .insert(record.ensembl_id, result.records.len() as u32);
        result
            .from_hgnc
            .insert(record.hgnc_id.clone(), result.records.len() as u32);
        result.records.push(XlinkDbRecord {
            entrez_id: record.entrez_id,
            ensembl_gene_id: record.ensembl_id,
            symbol: record.symbol,
            hgnc_id: record.hgnc_id,
        });
        total_count += 1;
    }
    tracing::debug!(
        "... done loading {} records and building maps in {:?}",
        total_count,
        before_loading.elapsed(),
    );

    Ok(result)
}

/// Information from ACMG file.
#[derive(Deserialize, Default, Clone, Debug)]
pub struct AcmgRecord {
    /// HGNC gene symbol.
    pub gene_symbol: String,
    /// ENSEMBL gene ID.
    pub ensembl_gene_id: String,
    /// Entrez / NCBI gene ID.
    #[serde(alias = "ncbi_gene_id")]
    pub entrez_id: u32,
}

/// Container for set of genes in ACMG
#[derive(Default, Clone, Debug)]
pub struct AcmgDb {
    pub entrez_ids: HashSet<u32>,
}

#[tracing::instrument]
fn load_acmg_db(path: &Path) -> Result<AcmgDb, anyhow::Error> {
    tracing::debug!("loading ACMG TSV records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = AcmgDb::default();

    // Setup CSV reader for ordinary TSV file - header is written on top and used;
    // no comment.
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path.to_str().unwrap())?);

    let mut total_count = 0;
    for record in reader.deserialize() {
        let record: AcmgRecord = record?;
        result.entrez_ids.insert(record.entrez_id);
        total_count += 1;
    }
    tracing::debug!(
        "... done loading {} records in {:?}",
        total_count,
        before_loading.elapsed(),
    );

    Ok(result)
}

impl AcmgDb {
    pub fn contains(&self, entrez_id: u32) -> bool {
        self.entrez_ids.contains(&entrez_id)
    }
}

/// Information to store for an OMIM entry.
#[derive(Deserialize, Default, Clone, Debug)]
pub struct OmimRecord {
    /// OMIM ID
    pub omim_id: u32,
    /// Entrez gene ID
    pub entrez_id: u32,
}

/// Container for set of genes in OMIM
#[derive(Default, Clone, Debug)]
pub struct OmimDb {
    pub entrez_ids: HashSet<u32>,
}

impl OmimDb {
    pub fn contains(&self, entrez_id: u32) -> bool {
        self.entrez_ids.contains(&entrez_id)
    }
}

#[tracing::instrument]
fn load_mim2gene_db(path: &Path) -> Result<OmimDb, anyhow::Error> {
    tracing::debug!("loading OMIM TSV records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = OmimDb::default();

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path.to_str().unwrap())?);

    let mut total_count = 0;
    for record in reader.deserialize() {
        let record: OmimRecord = record?;
        result.entrez_ids.insert(record.entrez_id);
        total_count += 1;
    }
    tracing::debug!(
        "... done loading {} records in {:?}",
        total_count,
        before_loading.elapsed(),
    );

    Ok(result)
}

/// Bundle of gene region DBs and the xlink info packaged with VarFish.
#[derive(Default, Debug)]
pub struct GeneDb {
    pub xlink: XlinkDb,
    pub acmg: AcmgDb,
    pub mim2gene: OmimDb,
}

// Load all gene information, such as region, id mapping and symbols.
#[tracing::instrument]
pub fn load_gene_db(path_db: &str, genome_release: GenomeRelease) -> Result<GeneDb, anyhow::Error> {
    info!("Loading gene dbs");

    let result = GeneDb {
        xlink: load_xlink_db(Path::new(path_db).join("noref/genes/xlink.bin").as_path())?,
        acmg: load_acmg_db(Path::new(path_db).join("noref/genes/acmg.tsv").as_path())?,
        mim2gene: load_mim2gene_db(
            Path::new(path_db)
                .join("noref/genes/mim2gene.tsv")
                .as_path(),
        )?,
    };

    Ok(result)
}
