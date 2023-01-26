//! Background database overlapping.

use std::{fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use tracing::{debug, info};

use crate::{
    common::CHROMS, sv::conf::BackgroundDbsConf,
    world_flatbuffers::var_fish_server_worker::BackgroundDatabase,
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Code for background database overlappers.
pub struct BgDb {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<BgDbRecord>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl Default for BgDb {
    fn default() -> Self {
        Self {
            records: Default::default(),
            trees: Default::default(),
        }
    }
}

/// Information to store for background database.
pub struct BgDbRecord {
    /// Count associated with the record.
    pub count: u32,
}

impl Default for BgDbRecord {
    fn default() -> Self {
        Self {
            count: Default::default(),
        }
    }
}

/// Load background database from a `.bin` file as created by `sv convert-bgdb`.
#[tracing::instrument]
pub fn load_bg_db_records(path: &Path) -> Result<BgDb, anyhow::Error> {
    debug!("loading binary bg db records from {:?}", path);

    let before_loading = Instant::now();
    let mut result = BgDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let file = File::open(&path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let bg_db = flatbuffers::root::<BackgroundDatabase>(&mmap)?;
    let records = bg_db.records().expect("no records in bg db");

    for record in &records {
        let chrom_no = record.chrom_no() as usize;
        let key = record.begin()..record.end();
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(BgDbRecord {
            count: record.count(),
        });
    }
    debug!(
        "done loading background dbs from {:?} in {:?}",
        path,
        before_loading.elapsed()
    );

    Ok(result)
}

/// Bundle of all background databases (including in-house).
pub struct BgDbBundle {
    pub gnomad: BgDb,
    pub dbvar: BgDb,
    pub dgv: BgDb,
    pub dgv_gs: BgDb,
    pub exac: BgDb,
    pub inhouse: BgDb,
}

// Load all background databases from database given the configuration.
#[tracing::instrument]
pub fn load_bg_dbs(path_db: &str, conf: &BackgroundDbsConf) -> Result<BgDbBundle, anyhow::Error> {
    info!("Loading background dbs");
    let result = BgDbBundle {
        gnomad: load_bg_db_records(Path::new(path_db).join(&conf.gnomad_sv.path).as_path())?,
        dbvar: load_bg_db_records(Path::new(path_db).join(&conf.dbvar.path).as_path())?,
        dgv: load_bg_db_records(Path::new(path_db).join(&conf.dgv.path).as_path())?,
        dgv_gs: load_bg_db_records(Path::new(path_db).join(&conf.dgv_gs.path).as_path())?,
        exac: load_bg_db_records(Path::new(path_db).join(&conf.exac.path).as_path())?,
        inhouse: load_bg_db_records(Path::new(path_db).join(&conf.inhouse.path).as_path())?,
    };

    Ok(result)
}
