//! Code for annotating with TADs and genes in the same TAD.

use std::path::Path;

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use indexmap::IndexMap;
use tracing::info;

use crate::{
    common::{build_chrom_map, open_read_maybe_gz, CHROMS},
    db::conf::{GenomeRelease, TadSet as TadSetChoice},
};

use super::{
    interpreter::{BND_SLACK, INS_SLACK},
    records::ChromRange,
    schema::{StructuralVariant, SvType},
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

/// Information to store for a TAD set entry.
#[derive(Default, Clone, Debug)]
pub struct Record {
    /// Chromosome number.
    pub chrom_no: u32,
    /// 0-based begin position.
    pub begin: i32,
    /// End position.
    pub end: i32,
}

/// TAD set overlapping information.
#[derive(Default, Debug)]
pub struct TadSet {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<Record>>,
    /// Interval trees, stored by chromosome.
    pub records_trees: Vec<IntervalTree>,
    /// Maximal distance to boundary to track.
    pub boundary_max_dist: i32,
    /// Boundaries, stored by chromosome.
    pub boundaries: Vec<Vec<i32>>,
    /// Interval triees to boundaries, stored by chromosome.
    pub boundaries_trees: Vec<IntervalTree>,
}

impl TadSet {
    /// Fetch TADs overlapping with a `ChromRange`.
    pub fn fetch_tads(
        &self,
        chrom_range: &ChromRange,
        chrom_map: &IndexMap<String, usize>,
    ) -> Vec<Record> {
        let mut result = Vec::new();

        let queries = {
            let chrom_idx = *chrom_map
                .get(&chrom_range.chromosome)
                .expect("invalid chromosome");
            vec![(chrom_idx, chrom_range.begin..chrom_range.end)]
        };

        for (chrom_idx, query) in queries {
            self.records_trees[chrom_idx]
                .find(query.clone())
                .iter()
                .for_each(|cursor| {
                    result.push(self.records[chrom_idx][*cursor.data() as usize].clone())
                });
        }

        result
    }

    pub fn overlapping_tads(
        &self,
        sv: &StructuralVariant,
        chrom_map: &IndexMap<String, usize>,
    ) -> Vec<Record> {
        let mut result = Vec::new();

        let queries = {
            let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
            match sv.sv_type {
                SvType::Bnd => {
                    let chrom_idx2 = *chrom_map
                        .get(sv.chrom2.as_ref().expect("no chrom2?"))
                        .expect("invalid chromosome");
                    vec![
                        (
                            chrom_idx,
                            sv.pos.saturating_sub(BND_SLACK)..sv.pos.saturating_add(BND_SLACK),
                        ),
                        (
                            chrom_idx2,
                            sv.end.saturating_sub(BND_SLACK)..sv.end.saturating_add(BND_SLACK),
                        ),
                    ]
                }
                SvType::Ins => vec![(
                    chrom_idx,
                    sv.pos.saturating_sub(INS_SLACK)..sv.pos.saturating_sub(INS_SLACK),
                )],
                _ => vec![(chrom_idx, sv.pos.saturating_sub(1)..sv.end)],
            }
        };

        for (chrom_idx, query) in queries {
            self.records_trees[chrom_idx]
                .find(query.clone())
                .iter()
                .for_each(|cursor| {
                    result.push(self.records[chrom_idx][*cursor.data() as usize].clone())
                });
        }

        result
    }

    pub fn boundary_dist(
        &self,
        sv: &StructuralVariant,
        chrom_map: &IndexMap<String, usize>,
    ) -> Option<u32> {
        let delta = self.boundary_max_dist;

        let queries = {
            let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
            match sv.sv_type {
                SvType::Bnd => {
                    let chrom_idx2 = *chrom_map
                        .get(sv.chrom2.as_ref().expect("no chrom2?"))
                        .expect("invalid chromosome");
                    vec![
                        (
                            chrom_idx,
                            sv.pos.saturating_sub(delta)..sv.pos.saturating_add(delta),
                            sv.pos,
                        ),
                        (
                            chrom_idx2,
                            sv.end.saturating_sub(delta)..sv.end.saturating_add(delta),
                            sv.end,
                        ),
                    ]
                }
                SvType::Ins => vec![(
                    chrom_idx,
                    sv.pos.saturating_sub(delta)..sv.pos.saturating_add(delta),
                    sv.pos,
                )],
                _ => vec![
                    (
                        chrom_idx,
                        sv.pos.saturating_sub(delta)..sv.pos.saturating_add(delta),
                        sv.pos,
                    ),
                    (
                        chrom_idx,
                        sv.end.saturating_sub(delta)..sv.end.saturating_add(delta),
                        sv.end,
                    ),
                ],
            }
        };

        let mut dists = Vec::new();
        for (chrom_idx, r, pos) in queries {
            self.boundaries_trees[chrom_idx]
                .find(r.clone())
                .iter()
                .for_each(|cursor| {
                    let boundary = self.boundaries[chrom_idx][*cursor.data() as usize];
                    dists.push(pos.abs_diff(boundary));
                });
        }
        dists.into_iter().min()
    }
}

/// Bundle of TAD sets packaged with VarFish.
#[derive(Default, Debug)]
pub struct TadSetBundle {
    pub hesc: TadSet,
}

impl TadSetBundle {
    pub fn fetch_tads(
        &self,
        tad_set: TadSetChoice,
        chrom_range: &ChromRange,
        chrom_map: &IndexMap<String, usize>,
    ) -> Vec<Record> {
        match tad_set {
            TadSetChoice::Hesc => self.hesc.fetch_tads(chrom_range, chrom_map),
        }
    }

    pub fn overlapping_tads(
        &self,
        tad_set: TadSetChoice,
        sv: &StructuralVariant,
        chrom_map: &IndexMap<String, usize>,
    ) -> Vec<Record> {
        match tad_set {
            TadSetChoice::Hesc => self.hesc.overlapping_tads(sv, chrom_map),
        }
    }

    pub fn boundary_dist(
        &self,
        tad_set: TadSetChoice,
        sv: &StructuralVariant,
        chrom_map: &IndexMap<String, usize>,
    ) -> Option<u32> {
        match tad_set {
            TadSetChoice::Hesc => self.hesc.boundary_dist(sv, chrom_map),
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
        pub begin: i32,
        /// 0-based end position from BED.
        pub end: i32,
    }
}

#[tracing::instrument]
fn load_tad_sets(path: &Path, boundary_max_dist: i32) -> Result<TadSet, anyhow::Error> {
    tracing::debug!("loading TAD set records from {:?}...", path);
    let chrom_map = build_chrom_map();

    let mut result = TadSet {
        boundary_max_dist,
        ..Default::default()
    };
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.records_trees.push(IntervalTree::new());
        result.boundaries.push(Vec::new());
        result.boundaries_trees.push(IntervalTree::new());
    }

    // Setup CSV reader for BED file - header is written as comment and must be
    // ignored.
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false) // BED has no header
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path.to_str().unwrap())?);
    let mut total_count = 0;
    for (i, record) in reader.deserialize().enumerate() {
        let record: input::Record = record?;
        let chrom_idx = *chrom_map.get(&record.chrom).expect("invalid chromosome");

        // TAD interval
        {
            let key = record.begin..record.end;
            result.records_trees[chrom_idx].insert(key, result.records[chrom_idx].len() as u32);
            result.records[chrom_idx].push(Record {
                chrom_no: chrom_idx as u32,
                begin: record.begin,
                end: record.end,
            });
        }

        // TAD boundary
        {
            if i == 0 {
                let key = record.begin.saturating_sub(1)..(record.begin + 1);
                result.boundaries_trees[chrom_idx]
                    .insert(key, result.boundaries[chrom_idx].len() as u32);
                result.boundaries[chrom_idx].push(record.begin);
            }
            let key = record.end.saturating_sub(1)..(record.end + 1);
            result.boundaries_trees[chrom_idx]
                .insert(key, result.boundaries[chrom_idx].len() as u32);
            result.boundaries[chrom_idx].push(record.end);
        }

        total_count += 1;
    }
    result
        .records_trees
        .iter_mut()
        .for_each(|tree| tree.index());
    result
        .boundaries_trees
        .iter_mut()
        .for_each(|tree| tree.index());
    tracing::debug!(
        "... done loading {} records and building trees",
        total_count
    );

    Ok(result)
}

// Load all pathogenic SV databases from database given the configuration.
#[tracing::instrument]
pub fn load_tads(
    path_db: &str,
    genome_release: GenomeRelease,
    max_tad_distance: i32,
) -> Result<TadSetBundle, anyhow::Error> {
    info!("Loading TAD sets dbs");
    let result = TadSetBundle {
        hesc: load_tad_sets(
            Path::new(path_db)
                .join(format!("{}/tads/hesc.bed", genome_release))
                .as_path(),
            max_tad_distance,
        )?,
    };

    Ok(result)
}
