//! Background database overlapping.

use std::{collections::HashMap, fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use tracing::{debug, info};

use crate::{
    common::CHROMS,
    sv::conf::BackgroundDbsConf,
    sv_query::dbrecords::{reciprocal_overlap, BeginEnd},
    world_flatbuffers::var_fish_server_worker::BackgroundDatabase,
    world_flatbuffers::var_fish_server_worker::SvType as FlatSvType,
};

use super::schema::{CaseQuery, StructuralVariant, SvType};

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

impl BgDb {
    pub fn count_overlaps(
        &self,
        chrom_map: &HashMap<String, usize>,
        enabled: bool,
        min_overlap: Option<f32>,
        slack_ins: u32,
        slack_bnd: u32,
        sv: &StructuralVariant,
    ) -> u32 {
        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let range = if sv.sv_type == SvType::Ins {
            sv.pos.saturating_sub(slack_ins)..sv.pos.saturating_add(slack_ins)
        } else if sv.sv_type == SvType::Bnd {
            sv.pos.saturating_sub(slack_bnd)..sv.pos.saturating_add(slack_bnd)
        } else {
            sv.pos..sv.end
        };

        self.trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count)
            .sum::<u32>()
    }
}

/// Information to store for background database.
pub struct BgDbRecord {
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
    /// Type of the background database record.
    pub sv_type: SvType,
    /// Count associated with the record.
    pub count: u32,
}

impl BeginEnd for BgDbRecord {
    fn begin(&self) -> u32 {
        self.begin
    }

    fn end(&self) -> u32 {
        self.end
    }
}

impl Default for BgDbRecord {
    fn default() -> Self {
        Self {
            begin: Default::default(),
            end: Default::default(),
            sv_type: Default::default(),
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
            begin: record.begin(),
            end: record.end(),
            sv_type: match record.sv_type() {
                FlatSvType::Del => SvType::Del,
                FlatSvType::Dup => SvType::Dup,
                FlatSvType::Inv => SvType::Inv,
                FlatSvType::Ins => SvType::Ins,
                FlatSvType::Bnd => SvType::Bnd,
                FlatSvType::Cnv => SvType::Cnv,
                _ => panic!("Invalid SV type from flatbuffer"),
            },
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
    pub g1k: BgDb,
    pub inhouse: BgDb,
}

/// Store background database counts for a structural variant.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct BgDbOverlaps {
    pub gnomad: u32,
    pub dbvar: u32,
    pub dgv: u32,
    pub dgv_gs: u32,
    pub exac: u32,
    pub g1k: u32,
    pub inhouse: u32,
}

impl BgDbBundle {
    pub fn count_overlaps(
        &self,
        sv: &StructuralVariant,
        query: &CaseQuery,
        chrom_map: &HashMap<String, usize>,
        slack_ins: u32,
        slack_bnd: u32,
    ) -> BgDbOverlaps {
        BgDbOverlaps {
            gnomad: self.gnomad.count_overlaps(
                chrom_map,
                query.svdb_gnomad_enabled,
                query.svdb_gnomad_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            dbvar: self.dbvar.count_overlaps(
                chrom_map,
                query.svdb_dbvar_enabled,
                query.svdb_dbvar_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            dgv: self.dgv.count_overlaps(
                chrom_map,
                query.svdb_dgv_enabled,
                query.svdb_dgv_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            dgv_gs: self.dgv_gs.count_overlaps(
                chrom_map,
                query.svdb_dgv_gs_enabled,
                query.svdb_dgv_gs_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            exac: self.exac.count_overlaps(
                chrom_map,
                query.svdb_exac_enabled,
                query.svdb_exac_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            g1k: self.g1k.count_overlaps(
                chrom_map,
                query.svdb_g1k_enabled,
                query.svdb_g1k_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
            inhouse: self.inhouse.count_overlaps(
                chrom_map,
                query.svdb_inhouse_enabled,
                query.svdb_inhouse_min_overlap,
                slack_ins,
                slack_bnd,
                sv,
            ),
        }
    }
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
        g1k: load_bg_db_records(Path::new(path_db).join(&conf.g1k.path).as_path())?,
        inhouse: load_bg_db_records(Path::new(path_db).join(&conf.inhouse.path).as_path())?,
    };

    Ok(result)
}
