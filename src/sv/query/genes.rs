//! Code for supporting annotation with overlapping genes.

use std::{fs::File, ops::Range, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use tracing::{debug, info};

use crate::{
    common::CHROMS,
    sv::conf::GenesConf,
    world_flatbuffers::var_fish_server_worker::{GeneRegionDatabase, XlinkDatabase},
};

use super::schema::Database;

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Information to store for a TAD set entry.
#[derive(Default, Clone, Debug)]
pub struct GeneRegionRecord {
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
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
    pub fn overlapping_gene_ids(&self, chrom_no: u32, query: Range<u32>) -> Vec<u32> {
        self.trees[chrom_no as usize]
            .find(query)
            .iter()
            .map(|cursor| self.records[chrom_no as usize][*cursor.data() as usize].gene_id)
            .collect()
    }
}

#[tracing::instrument]
fn load_gene_regions_db(path: &Path) -> Result<GeneRegionDb, anyhow::Error> {
    debug!("loading binary gene region records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = GeneRegionDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let gene_region_db = flatbuffers::root::<GeneRegionDatabase>(&mmap)?;
    let records = gene_region_db
        .records()
        .expect("no records in gene region db");

    let mut total_count = 0;
    for record in &records {
        let chrom_no = record.chrom_no() as usize;
        let key = record.begin()..record.end();
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(GeneRegionRecord {
            begin: record.begin(),
            end: record.end(),
            gene_id: record.gene_id(),
        });

        total_count += 1;
    }
    result.trees.iter_mut().for_each(|tree| tree.index());
    debug!(
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
            Database::Refseq => self.entrez_id_to_symbols(gene_id),
            Database::Ensembl => self.ensembl_to_symbols(gene_id),
        }
    }
}

#[tracing::instrument]
fn load_xlink_db(path: &Path) -> Result<XlinkDb, anyhow::Error> {
    debug!("loading binary xlink records from {:?}...", path);

    let before_loading = Instant::now();
    let mut result = XlinkDb::default();

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let xlink_db = flatbuffers::root::<XlinkDatabase>(&mmap)?;
    let records = xlink_db.records().expect("no records in xlink db");
    let strings = xlink_db.symbols().expect("no sybmols in xlink db");

    let mut total_count = 0;
    for (idx, record) in records.iter().enumerate() {
        result
            .from_entrez
            .insert(record.entrez_id(), result.records.len() as u32);
        result
            .from_ensembl
            .insert(record.ensembl_id(), result.records.len() as u32);
        result.records.push(XlinkDbRecord {
            entrez_id: record.entrez_id(),
            ensembl_gene_id: record.ensembl_id(),
            symbol: strings.get(idx).to_owned(),
        });
        total_count += 1;
    }
    debug!(
        "... done loading {} records and building maps in {:?}",
        total_count,
        before_loading.elapsed(),
    );

    Ok(result)
}

/// Bundle of gene region DBs and the xlink info packaged with VarFish.
pub struct GeneDb {
    pub refseq: GeneRegionDb,
    pub ensembl: GeneRegionDb,
    pub xlink: XlinkDb,
}

impl GeneDb {
    /// Return IDs of overlapping genes.
    pub fn overlapping_gene_ids(
        &self,
        database: Database,
        chrom_no: u32,
        query: Range<u32>,
    ) -> Vec<u32> {
        match database {
            Database::Refseq => self.refseq.overlapping_gene_ids(chrom_no, query),
            Database::Ensembl => self.ensembl.overlapping_gene_ids(chrom_no, query),
        }
    }
}

// Load all gene information, such as region, id mapping and symbols.
#[tracing::instrument]
pub fn load_gene_db(path_db: &str, conf: &GenesConf) -> Result<GeneDb, anyhow::Error> {
    info!("Loading gene dbs");
    let result = GeneDb {
        refseq: load_gene_regions_db(Path::new(path_db).join(&conf.refseq.regions.path).as_path())?,
        ensembl: load_gene_regions_db(
            Path::new(path_db)
                .join(&conf.ensembl.regions.path)
                .as_path(),
        )?,
        xlink: load_xlink_db(Path::new(path_db).join(&conf.xlink.path).as_path())?,
    };

    Ok(result)
}