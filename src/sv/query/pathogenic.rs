//! Code for working with lists of known pathogenic variants.

use std::{collections::HashMap, path::Path};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use serde::Serialize;
use tracing::{debug, info, warn};

use crate::{
    common::{build_chrom_map, open_read_maybe_gz, CHROMS},
    sv::conf::KnownPathogenicSvsConf,
};

use super::schema::{StructuralVariant, SvType};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Information to store for known pathogenic SV database.
#[derive(Default, Debug, Serialize, Clone)]
pub struct Record {
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
    /// Type of the known pathogenic SV database.
    pub sv_type: SvType,
    /// Identifier associated with the record.
    pub id: String,
}

/// Code for known pathogenic SV database overlappers.
#[derive(Default, Debug)]
pub struct PathoDb {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<Record>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl PathoDb {
    pub fn overlapping_records(
        &self,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
    ) -> Vec<Record> {
        if sv.sv_type == SvType::Ins || sv.sv_type == SvType::Bnd {
            return Vec::new();
        }

        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let range = sv.pos.saturating_sub(1)..sv.end;

        self.trees[chrom_idx]
            .find(range)
            .iter()
            .map(|e| &self.records[chrom_idx][*e.data() as usize])
            .cloned()
            .collect()
    }
}

/// Bundle of databases of known pathogenic variants.
pub struct PathoDbBundle {
    pub mms: PathoDb,
}

impl PathoDbBundle {
    pub fn overlapping_records(
        &self,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
    ) -> Vec<Record> {
        self.mms.overlapping_records(sv, chrom_map)
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
        /// Identifier of the record.
        pub id: String,
    }
}

#[tracing::instrument]
fn load_patho_db_records(path: &Path) -> Result<PathoDb, anyhow::Error> {
    debug!("loading pathogenic SV records from {:?}...", path);
    let chrom_map = build_chrom_map();

    let mut result = PathoDb::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false) // BED has no header
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path.to_str().unwrap())?);
    let mut total_count = 0;
    for record in reader.deserialize() {
        let record: input::Record = record?;
        let chrom_idx = *chrom_map.get(&record.chrom).expect("invalid chromosome");

        let key = record.begin..record.end;
        result.trees[chrom_idx].insert(key, result.records[chrom_idx].len() as u32);
        result.records[chrom_idx].push(Record {
            begin: record.begin,
            end: record.end,
            sv_type: SvType::Cnv,
            id: record.id,
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
pub fn load_patho_dbs(
    path_db: &str,
    conf: &KnownPathogenicSvsConf,
) -> Result<PathoDbBundle, anyhow::Error> {
    info!("Loading pathogenic SV dbs");
    let result = PathoDbBundle {
        mms: load_patho_db_records(
            Path::new(path_db)
                .join(&conf.mms_wetzel_danbro.path)
                .as_path(),
        )?,
    };

    Ok(result)
}
