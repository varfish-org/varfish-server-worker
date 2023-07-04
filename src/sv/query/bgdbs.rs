//! Background database overlapping.

use std::{ops::Range, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use indexmap::IndexMap;
use prost::Message;
use serde::{Deserialize, Serialize};
use strum_macros::{Display, EnumString};
use tracing::info;

use crate::{
    common::{trace_rss_now, CHROMS},
    db::{conf::GenomeRelease, pbs},
};

use super::{
    records::ChromRange,
    schema::{CaseQuery, StructuralVariant, SvType},
};

pub trait BeginEnd {
    /// 0-base begin position
    fn begin(&self) -> i32;
    /// 0-based end position
    fn end(&self) -> i32;
}

pub fn reciprocal_overlap(lhs: &impl BeginEnd, rhs: &Range<i32>) -> f32 {
    let lhs_b = lhs.begin();
    let lhs_e = lhs.end();
    let rhs_b = rhs.start;
    let rhs_e = rhs.end;
    let ovl_b = std::cmp::max(lhs_b, rhs_b);
    let ovl_e = std::cmp::min(lhs_e, rhs_e);
    if ovl_b >= ovl_e {
        0f32
    } else {
        let ovl_len = (ovl_e - ovl_b) as f32;
        let x1 = ovl_len / (lhs_e - lhs_b) as f32;
        let x2 = ovl_len / (rhs_e - rhs_b) as f32;
        x1.min(x2)
    }
}

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

/// Code for background database overlappers.
#[derive(Default, Debug)]
pub struct BgDb {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<BgDbRecord>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl BgDb {
    pub fn fetch_records(
        &self,
        genomic_region: &ChromRange,
        chrom_map: &IndexMap<String, usize>,
    ) -> Vec<BgDbRecord> {
        let chrom_idx = *chrom_map
            .get(&genomic_region.chromosome)
            .expect("invalid chromosome");
        let range = genomic_region.begin..genomic_region.end;

        self.trees[chrom_idx]
            .find(range)
            .iter()
            .map(|e| &self.records[chrom_idx][*e.data() as usize])
            .cloned()
            .collect()
    }

    pub fn count_overlaps(
        &self,
        chrom_map: &IndexMap<String, usize>,
        enabled: bool,
        min_overlap: Option<f32>,
        slack_ins: i32,
        slack_bnd: i32,
        sv: &StructuralVariant,
    ) -> u32 {
        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let range = if sv.sv_type == SvType::Ins {
            (sv.pos - slack_ins)..(sv.pos + slack_ins)
        } else if sv.sv_type == SvType::Bnd {
            (sv.pos - slack_bnd)..(sv.pos + slack_bnd)
        } else {
            (sv.pos - 1)..sv.end
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
#[derive(Serialize, Default, Debug, Clone)]
pub struct BgDbRecord {
    /// 0-based begin position.
    pub begin: i32,
    /// End position.
    pub end: i32,
    /// Type of the background database record.
    pub sv_type: SvType,
    /// Count associated with the record.
    pub count: u32,
}

impl BeginEnd for BgDbRecord {
    fn begin(&self) -> i32 {
        self.begin
    }

    fn end(&self) -> i32 {
        self.end
    }
}

/// Load background database from a `.bin` file as created by `sv convert-bgdb`.
#[tracing::instrument]
pub fn load_bg_db_records(path: &Path) -> Result<BgDb, anyhow::Error> {
    tracing::debug!("loading binary bg db records from {:?}", path);

    let before_loading = Instant::now();
    let mut result = BgDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let fcontents =
        std::fs::read(path).map_err(|e| anyhow::anyhow!("error reading {:?}: {}", &path, e))?;
    let bg_db = pbs::BackgroundDatabase::decode(std::io::Cursor::new(fcontents))
        .map_err(|e| anyhow::anyhow!("error decoding {:?}: {}", &path, e))?;
    let record_count = bg_db.records.len();

    for record in bg_db.records.into_iter() {
        let chrom_no = record.chrom_no as usize;
        let begin = match pbs::SvType::from_i32(record.sv_type).expect("invalid sv_type") {
            pbs::SvType::Bnd | pbs::SvType::Ins => record.start - 2,
            _ => record.start - 1,
        };
        let end = match pbs::SvType::from_i32(record.sv_type).expect("invalid sv_type") {
            pbs::SvType::Bnd | pbs::SvType::Ins => record.start - 1,
            _ => record.stop,
        };
        let key = begin..end;

        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(BgDbRecord {
            begin: record.start - 1,
            end: record.stop,
            sv_type: match pbs::SvType::from_i32(record.sv_type).expect("invalid sv_type") {
                pbs::SvType::Del => SvType::Del,
                pbs::SvType::Dup => SvType::Dup,
                pbs::SvType::Inv => SvType::Inv,
                pbs::SvType::Ins => SvType::Ins,
                pbs::SvType::Bnd => SvType::Bnd,
                pbs::SvType::Cnv => SvType::Cnv,
            },
            count: record.count,
        });
    }
    tracing::debug!(
        "done loading background db with {} records from {:?} in {:?}",
        record_count,
        path,
        before_loading.elapsed()
    );

    let before_building = Instant::now();
    result.trees.iter_mut().for_each(|tree| tree.index());
    tracing::debug!("done building itrees in {:?}", before_building.elapsed());

    trace_rss_now();

    Ok(result)
}

/// Enumeration of background database types.
#[derive(Serialize, Deserialize, Debug, PartialEq, EnumString, Display)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum BgDbType {
    GnomadSv,
    Dgv,
    DgvGs,
    Exac,
    G1k,
    Inhouse,
}

/// Bundle of all background databases (including in-house).
#[derive(Default, Debug)]
pub struct BgDbBundle {
    pub gnomad: BgDb,
    pub dbvar: BgDb,
    pub dgv: BgDb,
    pub dgv_gs: BgDb,
    pub exac: Option<BgDb>,
    pub g1k: Option<BgDb>,
    pub inhouse: Option<BgDb>,
}

/// Store background database counts for a structural variant.
#[derive(Serialize, Clone, Debug, PartialEq, Default)]
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
    pub fn fetch_records(
        &self,
        genome_range: &ChromRange,
        chrom_map: &IndexMap<String, usize>,
        db_type: BgDbType,
    ) -> Vec<BgDbRecord> {
        match db_type {
            BgDbType::GnomadSv => self.gnomad.fetch_records(genome_range, chrom_map),
            BgDbType::Dgv => self.dgv.fetch_records(genome_range, chrom_map),
            BgDbType::DgvGs => self.dgv_gs.fetch_records(genome_range, chrom_map),
            BgDbType::Exac => self
                .exac
                .as_ref()
                .map(|exac| exac.fetch_records(genome_range, chrom_map))
                .unwrap_or_default(),
            BgDbType::G1k => self
                .g1k
                .as_ref()
                .map(|g1k| g1k.fetch_records(genome_range, chrom_map))
                .unwrap_or_default(),
            BgDbType::Inhouse => self
                .inhouse
                .as_ref()
                .map(|inhouse| inhouse.fetch_records(genome_range, chrom_map))
                .unwrap_or_default(),
        }
    }

    pub fn count_overlaps(
        &self,
        sv: &StructuralVariant,
        query: &CaseQuery,
        chrom_map: &IndexMap<String, usize>,
        slack_ins: i32,
        slack_bnd: i32,
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
            exac: self.exac.as_ref().map_or(0, |exac| {
                exac.count_overlaps(
                    chrom_map,
                    query.svdb_exac_enabled,
                    query.svdb_exac_min_overlap,
                    slack_ins,
                    slack_bnd,
                    sv,
                )
            }),
            g1k: self.g1k.as_ref().map_or(0, |g1k| {
                g1k.count_overlaps(
                    chrom_map,
                    query.svdb_g1k_enabled,
                    query.svdb_g1k_min_overlap,
                    slack_ins,
                    slack_bnd,
                    sv,
                )
            }),
            inhouse: self.inhouse.as_ref().map_or(0, |inhouse| {
                inhouse.count_overlaps(
                    chrom_map,
                    query.svdb_inhouse_enabled,
                    query.svdb_inhouse_min_overlap,
                    slack_ins,
                    slack_bnd,
                    sv,
                )
            }),
        }
    }
}

// Load all background databases from database given the configuration.
#[tracing::instrument]
pub fn load_bg_dbs(
    path_db: &str,
    genome_release: GenomeRelease,
) -> Result<BgDbBundle, anyhow::Error> {
    info!("Loading background dbs");

    let path_exac = Path::new(path_db).join(format!("{}/strucvars/bgdbs/exac.bin", genome_release));
    let path_g1k = Path::new(path_db).join(format!("{}/strucvars/bgdbs/g1k.bin", genome_release));
    let path_inhouse = Path::new(path_db).join(format!("{}/strucvars/inhouse.bin", genome_release));

    let result = BgDbBundle {
        gnomad: load_bg_db_records(
            Path::new(path_db)
                .join(format!("{}/strucvars/bgdbs/gnomad.bin", genome_release))
                .as_path(),
        )?,
        dbvar: load_bg_db_records(
            Path::new(path_db)
                .join(format!("{}/strucvars/bgdbs/dbvar.bin", genome_release))
                .as_path(),
        )?,
        dgv: load_bg_db_records(
            Path::new(path_db)
                .join(format!("{}/strucvars/bgdbs/dgv.bin", genome_release))
                .as_path(),
        )?,
        dgv_gs: load_bg_db_records(
            Path::new(path_db)
                .join(format!("{}/strucvars/bgdbs/dgv-gs.bin", genome_release))
                .as_path(),
        )?,
        exac: path_exac
            .exists()
            .then(|| load_bg_db_records(path_exac.as_path()))
            .transpose()?,
        g1k: path_g1k
            .exists()
            .then(|| load_bg_db_records(path_g1k.as_path()))
            .transpose()?,
        inhouse: path_inhouse
            .exists()
            .then(|| load_bg_db_records(path_inhouse.as_path()))
            .transpose()?,
    };

    Ok(result)
}
