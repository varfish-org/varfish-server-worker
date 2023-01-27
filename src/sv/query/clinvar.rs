//! Code for working with Clinvar SV.

use std::{collections::HashMap, fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use memmap2::Mmap;
use thousands::Separable;
use tracing::{debug, info, warn};

use crate::{
    common::{reciprocal_overlap, CHROMS},
    sv::conf::ClinvarSvConf,
    world_flatbuffers::var_fish_server_worker::ClinvarSvDatabase,
};

use super::schema::{Pathogenicity, StructuralVariant, SvType, VariationType};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<u32, u32>;

/// Information to store for ClinVar SV database.
#[derive(Default, Debug)]
pub struct Record {
    /// 0-based begin position.
    pub begin: u32,
    /// End position.
    pub end: u32,
    /// Type of the ClinVar SV variation.
    pub variation_type: VariationType,
    /// Pathogenicity annotation.
    pub pathogenicity: Pathogenicity,
}

/// Code for known Clinvar SV database.
#[derive(Default, Debug)]
pub struct ClinvarSv {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<Record>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl ClinvarSv {
    pub fn count_overlaps(
        &self,
        sv: &StructuralVariant,
        chrom_map: &HashMap<String, usize>,
        min_patho: Option<Pathogenicity>,
        min_overlap: Option<f32>,
    ) -> u32 {
        if sv.sv_type == SvType::Ins || sv.sv_type == SvType::Bnd {
            return 0;
        }

        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let range = sv.pos.saturating_sub(1)..sv.end;

        self.trees[chrom_idx]
            .find(range)
            .iter()
            .map(|e| &self.records[chrom_idx][*e.data() as usize])
            .filter(|record| {
                min_overlap.map_or(true, |min_overlap| {
                    reciprocal_overlap(record.begin..record.end, sv.pos.saturating_sub(1)..sv.end)
                        >= min_overlap
                })
            })
            .map(|record| record.pathogenicity)
            .filter(|patho| *patho >= min_patho.unwrap_or(Pathogenicity::Benign))
            .map(|_record| 1)
            .sum::<u32>()
    }
}

// Load the Clinvar SV databases from database given the configuration.
#[tracing::instrument]
pub fn load_clinvar_sv(path_db: &str, conf: &ClinvarSvConf) -> Result<ClinvarSv, anyhow::Error> {
    info!("loading binaryclinvar SV dbs");

    let before_loading = Instant::now();
    let mut result = ClinvarSv::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let path = Path::new(path_db).join(&conf.svs.path);
    let file = File::open(&path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let bg_db = flatbuffers::root::<ClinvarSvDatabase>(&mmap)?;
    let records = bg_db.records().expect("no records in clinvar SV db");
    let mut total_count = 0;

    for record in &records {
        let chrom_no = record.chrom_no() as usize;
        let key = record.begin()..record.end();
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(Record {
            begin: record.begin(),
            end: record.end(),
            variation_type: record.variation_type().try_into()?,
            pathogenicity: record.pathogenicity().try_into()?,
        });
        total_count += 1;
    }
    debug!(
        "done loading {} clinvar SV records from {:?} in {:?}",
        &total_count.separate_with_commas(),
        &path,
        before_loading.elapsed()
    );

    let before_building = Instant::now();
    result.trees.iter_mut().for_each(|tree| tree.index());
    debug!("done building itrees in {:?}", before_building.elapsed());

    Ok(result)
}
