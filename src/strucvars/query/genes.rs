//! Code for supporting annotation with overlapping genes.

use std::{collections::HashSet, ops::Range, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use mehari::common::open_read_maybe_gz;
use prost::Message;
use serde::Deserialize;
use tracing::info;

use crate::{
    common::CHROMS,
    db::{
        conf::{Database, GenomeRelease},
        pbs,
    },
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

/// Information to store for a TAD set entry.
#[derive(Default, Clone, Debug)]
pub struct GeneRegionRecord {
    /// 0-based begin position.
    pub begin: i32,
    /// End position.
    pub end: i32,
    /// gene ID
    pub gene_id: u32,
}

/// Gene region overlapping information.
#[derive(Default, Debug)]
pub struct GeneRegionDb {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<GeneRegionRecord>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl GeneRegionDb {
    /// Return IDs of overlapping genes.
    pub fn overlapping_gene_ids(&self, chrom_no: u32, query: Range<i32>) -> Vec<u32> {
        self.trees[chrom_no as usize]
            .find(query)
            .iter()
            .map(|cursor| self.records[chrom_no as usize][*cursor.data() as usize].gene_id)
            .collect()
    }
}

#[tracing::instrument]
fn load_gene_regions_db(path: &Path) -> Result<GeneRegionDb, anyhow::Error> {
    tracing::debug!("loading binary gene region records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = GeneRegionDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let fcontents =
        std::fs::read(path).map_err(|e| anyhow::anyhow!("error reading {:?}: {}", &path, e))?;
    let gene_region_db = pbs::GeneRegionDatabase::decode(std::io::Cursor::new(fcontents))
        .map_err(|e| anyhow::anyhow!("error decoding {:?}: {}", &path, e))?;

    let mut total_count = 0;
    for record in gene_region_db.records.into_iter() {
        let chrom_no = record.chrom_no as usize;
        let key = (record.start - 1)..record.stop;
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(GeneRegionRecord {
            begin: record.start - 1,
            end: record.stop,
            gene_id: record.gene_id,
        });

        total_count += 1;
    }
    result.trees.iter_mut().for_each(|tree| tree.index());
    tracing::debug!(
        "... done loading {} records and building trees in {:?}",
        total_count,
        before_loading.elapsed(),
    );

    Ok(result)
}

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

impl XlinkDb {
    fn entrez_id_to_symbols(&self, entrez_id: u32) -> Vec<String> {
        let mut tmp = self
            .from_entrez
            .get_vec(&entrez_id)
            .map_or(Vec::new(), |idxs| {
                idxs.iter()
                    .map(|idx| self.records[*idx as usize].symbol.clone())
                    .collect::<Vec<String>>()
            });
        tmp.sort();
        tmp.dedup();
        tmp
    }

    fn ensembl_to_symbols(&self, ensembl_id: u32) -> Vec<String> {
        let mut tmp = self
            .from_ensembl
            .get_vec(&ensembl_id)
            .map_or(Vec::new(), |idxs| {
                idxs.iter()
                    .map(|idx| self.records[*idx as usize].symbol.clone())
                    .collect::<Vec<String>>()
            });
        tmp.sort();
        tmp.dedup();
        tmp
    }

    pub fn gene_id_to_symbols(&self, database: Database, gene_id: u32) -> Vec<String> {
        match database {
            Database::RefSeq => self.entrez_id_to_symbols(gene_id),
            Database::Ensembl => self.ensembl_to_symbols(gene_id),
        }
    }
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
fn load_omim_db(path: &Path) -> Result<OmimDb, anyhow::Error> {
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
    pub refseq: GeneRegionDb,
    pub ensembl: GeneRegionDb,
    pub xlink: XlinkDb,
    pub acmg: AcmgDb,
    pub omim: OmimDb,
}

impl GeneDb {
    /// Return IDs of overlapping genes.
    pub fn overlapping_gene_ids(
        &self,
        database: Database,
        chrom_no: u32,
        query: Range<i32>,
    ) -> Vec<u32> {
        match database {
            Database::RefSeq => self.refseq.overlapping_gene_ids(chrom_no, query),
            Database::Ensembl => self.ensembl.overlapping_gene_ids(chrom_no, query),
        }
    }
}

// Load all gene information, such as region, id mapping and symbols.
#[tracing::instrument]
pub fn load_gene_db(path_db: &str, genome_release: GenomeRelease) -> Result<GeneDb, anyhow::Error> {
    info!("Loading gene dbs");

    let result = GeneDb {
        refseq: load_gene_regions_db(
            Path::new(path_db)
                .join(format!("{}/genes/refseq_regions.bin", genome_release))
                .as_path(),
        )?,
        ensembl: load_gene_regions_db(
            Path::new(path_db)
                .join(format!("{}/genes/ensembl_regions.bin", genome_release))
                .as_path(),
        )?,
        xlink: load_xlink_db(Path::new(path_db).join("noref/genes/xlink.bin").as_path())?,
        acmg: load_acmg_db(Path::new(path_db).join("noref/genes/acmg.tsv").as_path())?,
        omim: load_omim_db(Path::new(path_db).join("noref/genes/omim.tsv").as_path())?,
    };

    Ok(result)
}
