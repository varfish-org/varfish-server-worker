//! Code for annotating with TADs and genes in the same TAD.

use std::{collections::HashMap, path::Path};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use tracing::{debug, info};

use crate::{
    common::{build_chrom_map, open_maybe_gz, CHROMS},
    sv::conf::TadsConf,
};

use super::{
    interpreter::{BND_SLACK, INS_SLACK},
    schema::{StructuralVariant, SvType},
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Enum for selecting TAD set.
#[derive(Debug, PartialEq, Eq)]
pub enum TadSetChoice {
    // HESC data set from Dixon et al. (2012)
    Hesc,
    // IMR90 data set from Dixon et al. (2012)
    Imr90,
}

/// Information to store for a TAD set entry.
#[derive(Default, Clone, Debug)]
pub struct Record {
    /// Chromosome number.
    pub chrom_no: u32,
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
}

/// TAD set overlapping information.
#[derive(Default, Debug)]
pub struct TadSet {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<Record>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl TadSet {
    pub fn overlapping_tads(
        &self,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
    ) -> Vec<Record> {
        let mut result = Vec::new();

        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let queries = match sv.sv_type {
            SvType::Bnd => vec![
                sv.pos.saturating_sub(BND_SLACK)..sv.pos.saturating_add(BND_SLACK),
                sv.end.saturating_sub(BND_SLACK)..sv.end.saturating_add(BND_SLACK),
            ],
            SvType::Ins => vec![sv.pos.saturating_sub(INS_SLACK)..sv.pos.saturating_sub(INS_SLACK)],
            _ => vec![sv.pos.saturating_sub(1)..sv.end],
        };

        for query in queries {
            self.trees[chrom_idx]
                .find(query.clone())
                .iter()
                .for_each(|cursor| {
                    result.push(self.records[chrom_idx][*cursor.data() as usize].clone())
                });
        }

        result
    }
}

/// Bundle of TAD sets packaged with VarFish.
pub struct TadSetBundle {
    pub hesc: TadSet,
    pub imr90: TadSet,
}

impl TadSetBundle {
    pub fn overlapping_tads(
        &self,
        tad_set: TadSetChoice,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
    ) -> Vec<Record> {
        match tad_set {
            TadSetChoice::Hesc => self.hesc.overlapping_tads(sv, chrom_map),
            TadSetChoice::Imr90 => self.imr90.overlapping_tads(sv, chrom_map),
        }
    }
}
/// Module with code for loading data from input.
mod input {
    use serde::Deserialize;

    /// Type for record structs from input.
    #[derive(Deserialize, Debug)]
    pub struct Record {
        /// Chromosome name
        pub chrom: String,
        /// 0-based begin position from BEd.
        pub begin: u32,
        /// 0-based end position from BED.
        pub end: u32,
    }
}

#[tracing::instrument]
fn load_tad_sets(path: &Path) -> Result<TadSet, anyhow::Error> {
    debug!("loading TAD set records from {:?}...", path);
    let chrom_map = build_chrom_map();

    let mut result = TadSet::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false) // BED has no header
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(path.to_str().unwrap())?);
    let mut total_count = 0;
    for record in reader.deserialize() {
        let record: input::Record = record?;
        let chrom_idx = *chrom_map.get(&record.chrom).expect("invalid chromosome");

        let key = record.begin..record.end;
        result.trees[chrom_idx].insert(key, result.records[chrom_idx].len() as u32);
        result.records[chrom_idx].push(Record {
            chrom_no: chrom_idx as u32,
            begin: record.begin,
            end: record.end,
        });

        total_count += 1;
    }
    result.trees.iter_mut().for_each(|tree| tree.index());
    debug!(
        "... done loading {} records and building trees",
        total_count
    );

    Ok(result)
}

// Load all pathogenic SV databases from database given the configuration.
#[tracing::instrument]
pub fn load_tads(path_db: &str, conf: &TadsConf) -> Result<TadSetBundle, anyhow::Error> {
    info!("Loading TAD sets dbs");
    let result = TadSetBundle {
        hesc: load_tad_sets(Path::new(path_db).join(&conf.hesc.path).as_path())?,
        imr90: load_tad_sets(Path::new(path_db).join(&conf.imr90.path).as_path())?,
    };

    Ok(result)
}
