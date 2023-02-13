//! Overlap with masked regions (e.g., repeats, segmental duplications)

use std::{collections::HashMap, fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use serde::Serialize;
use tracing::{debug, info};

use crate::{
    common::{trace_rss_now, CHROMS},
    db::conf::MaskedRegionDbs,
    world_flatbuffers::var_fish_server_worker::MaskedDatabase,
};

use super::{
    bgdbs::BeginEnd,
    records::ChromRange,
    schema::{StructuralVariant, SvType},
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Code for masked regions overlappers.
#[derive(Default, Debug)]
pub struct MaskedDb {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<MaskedDbRecord>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

/// Result for `MaskedDb::fetch_records`.
pub struct FetchRecordsResult {
    /// Overlaps with left end.
    pub left: Vec<MaskedDbRecord>,
    /// Overlaps with right end.
    pub right: Vec<MaskedDbRecord>,
}

impl MaskedDb {
    pub fn fetch_records(
        &self,
        genomic_region: &ChromRange,
        chrom_map: &HashMap<String, usize>,
    ) -> FetchRecordsResult {
        let chrom_idx = *chrom_map
            .get(&genomic_region.chromosome)
            .expect("invalid chromosome");
        let range_left = genomic_region.begin..(genomic_region.begin + 1);
        let range_right = genomic_region.end.saturating_sub(1)..genomic_region.end;

        FetchRecordsResult {
            left: self.trees[chrom_idx]
                .find(range_left)
                .iter()
                .map(|e| &self.records[chrom_idx][*e.data() as usize])
                .cloned()
                .collect(),
            right: self.trees[chrom_idx]
                .find(range_right)
                .iter()
                .map(|e| &self.records[chrom_idx][*e.data() as usize])
                .cloned()
                .collect(),
        }
    }

    /// Counts how many of `sv`'s breakpoints (0, 1, 2) overlap with a masked region.
    ///
    /// For insertions and breake-ends, the one primary breakpoint counts as 2.
    pub fn masked_breakpoint_count(
        &self,
        chrom_map: &HashMap<String, usize>,
        sv: &StructuralVariant,
    ) -> u32 {
        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let (range_left, range_right) = if sv.sv_type == SvType::Ins || sv.sv_type == SvType::Bnd {
            (sv.pos..(sv.pos + 1), sv.pos..(sv.pos + 1))
        } else {
            (sv.pos..(sv.pos + 1), sv.end.saturating_sub(1)..sv.end)
        };

        let any_left = !self.trees[chrom_idx].find(range_left).is_empty();
        let any_right = !self.trees[chrom_idx].find(range_right).is_empty();

        (any_left as u32) + (any_right as u32)
    }
}

/// Information to store for background database.
#[derive(Serialize, Default, Debug, Clone)]
pub struct MaskedDbRecord {
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
}

impl BeginEnd for MaskedDbRecord {
    fn begin(&self) -> u32 {
        self.begin
    }

    fn end(&self) -> u32 {
        self.end
    }
}

/// Load background database from a `.bin` file as created by `sv convert-bgdb`.
#[tracing::instrument]
pub fn load_masked_db_records(path: &Path) -> Result<MaskedDb, anyhow::Error> {
    debug!("loading binary masked db records from {:?}", path);

    let before_loading = Instant::now();
    let mut result = MaskedDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let masked_db = flatbuffers::root::<MaskedDatabase>(&mmap)?;
    let records = masked_db.records().expect("no records in masked db");

    for record in &records {
        let chrom_no = record.chrom_no() as usize;
        let key = record.begin()..record.end();
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(MaskedDbRecord {
            begin: record.begin(),
            end: record.end(),
        });
    }
    debug!(
        "done loading masked dbs from {:?} in {:?}",
        path,
        before_loading.elapsed()
    );

    let before_building = Instant::now();
    result.trees.iter_mut().for_each(|tree| tree.index());
    debug!("done building itrees in {:?}", before_building.elapsed());

    trace_rss_now();

    Ok(result)
}

/// Enumeration of all masked region databases.
pub enum MaskedRegionType {
    Repeat,
    SegDup,
}

/// Bundle of all masked region databases (including in-house).
#[derive(Default, Debug)]
pub struct MaskedDbBundle {
    pub repeat: MaskedDb,
    pub segdup: MaskedDb,
}

/// Store masked region database counts for a structural variant.
#[derive(Serialize, Clone, Debug, PartialEq, Default)]
pub struct MaskedBreakpointCount {
    pub repeat: u32,
    pub segdup: u32,
}

impl MaskedDbBundle {
    pub fn fetch_records(
        &self,
        genome_range: &ChromRange,
        chrom_map: &HashMap<String, usize>,
        db_type: MaskedRegionType,
    ) -> FetchRecordsResult {
        match db_type {
            MaskedRegionType::Repeat => self.repeat.fetch_records(genome_range, chrom_map),
            MaskedRegionType::SegDup => self.segdup.fetch_records(genome_range, chrom_map),
        }
    }

    pub fn masked_breakpoint_count(
        &self,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
    ) -> MaskedBreakpointCount {
        MaskedBreakpointCount {
            repeat: self.repeat.masked_breakpoint_count(chrom_map, sv),
            segdup: self.segdup.masked_breakpoint_count(chrom_map, sv),
        }
    }
}

// Load all masked region databases from database given the configuration.
#[tracing::instrument]
pub fn load_bg_dbs(path_db: &str, conf: &MaskedRegionDbs) -> Result<MaskedDbBundle, anyhow::Error> {
    info!("Loading masked region dbs");
    let result = MaskedDbBundle {
        repeat: load_masked_db_records(
            Path::new(path_db)
                .join(
                    conf.repeat
                        .bin_path
                        .as_ref()
                        .expect("no binary path for repeat masked"),
                )
                .as_path(),
        )?,
        segdup: load_masked_db_records(
            Path::new(path_db)
                .join(
                    conf.segdup
                        .bin_path
                        .as_ref()
                        .expect("no binary path for segmental duplications"),
                )
                .as_path(),
        )?,
    };

    Ok(result)
}
