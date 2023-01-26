//! Code for supporting annotation with overlapping genes.

use std::{fs::File, ops::Range, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use tracing::{debug, info};

use crate::{
    common::CHROMS, sv::conf::GenesConf,
    world_flatbuffers::var_fish_server_worker::GeneRegionDatabase,
};

use super::schema::Database;

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Information to store for a TAD set entry.
#[derive(Default, Clone, Debug)]
pub struct Record {
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
    pub records: Vec<Vec<Record>>,
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

/// Bundle of gene DBs packaged with VarFish.
pub struct GeneRegionDbBundle {
    pub refseq: GeneRegionDb,
    pub ensembl: GeneRegionDb,
}

impl GeneRegionDbBundle {
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
        result.records[chrom_no].push(Record {
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

// Load all pathogenic SV databases from database given the configuration.
#[tracing::instrument]
pub fn load_gene_regions(
    path_db: &str,
    conf: &GenesConf,
) -> Result<GeneRegionDbBundle, anyhow::Error> {
    info!("Loading gene region dbs");
    let result = GeneRegionDbBundle {
        refseq: load_gene_regions_db(Path::new(path_db).join(&conf.refseq.regions.path).as_path())?,
        ensembl: load_gene_regions_db(
            Path::new(path_db)
                .join(&conf.ensembl.regions.path)
                .as_path(),
        )?,
    };

    Ok(result)
}
