//! Code for working with Clinvar SV.

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use indexmap::IndexMap;
use prost::Message;
use thousands::Separable;
use tracing::{info, warn};

use crate::common::{reciprocal_overlap, GenomeRelease, CHROMS};

use super::{
    schema::ChromRange,
    schema::{Pathogenicity, StructuralVariant, SvType},
};

/// Alias for the interval tree that we use.
type IntervalTree = ArrayBackedIntervalTree<i32, u32>;

/// Code for known Clinvar SV database.
///
/// NB that the records have 1-based coordinates whereas the trees use 0-based coordinates.
#[derive(Default, Debug)]
pub struct ClinvarSv {
    /// Records, stored by chromosome.
    pub records: Vec<Vec<crate::pbs::clinvar::SvRecord>>,
    /// Interval trees, stored by chromosome.
    pub trees: Vec<IntervalTree>,
}

impl ClinvarSv {
    /// Fetches records of overlapping entries.
    pub fn fetch_records(
        &self,
        chrom_range: &ChromRange,
        chrom_map: &IndexMap<String, usize>,
        min_patho: Option<Pathogenicity>,
    ) -> Vec<crate::pbs::clinvar::SvRecord> {
        let chrom_idx = *chrom_map
            .get(&chrom_range.chromosome)
            .expect("invalid chromosome");
        let range = chrom_range.begin..chrom_range.end;

        self.trees[chrom_idx]
            .find(range)
            .iter()
            .map(|e| self.records[chrom_idx][*e.data() as usize].clone())
            .filter(|record| {
                record.pathogenicity >= min_patho.unwrap_or(Pathogenicity::Benign) as i32
            })
            .collect()
    }

    /// Returns the overlapping RCVs
    pub fn overlapping_rcvs(
        &self,
        sv: &StructuralVariant,
        chrom_map: &IndexMap<String, usize>,
        min_patho: Option<Pathogenicity>,
        min_overlap: Option<f32>,
    ) -> Vec<u32> {
        if sv.sv_type == SvType::Ins || sv.sv_type == SvType::Bnd {
            return Vec::new();
        }

        let chrom_idx = *chrom_map.get(&sv.chrom).expect("invalid chromosome");
        let range = sv.pos.saturating_sub(1)..sv.end;

        self.trees[chrom_idx]
            .find(range)
            .into_iter()
            .map(|e| &self.records[chrom_idx][*e.data() as usize])
            .filter(|record| {
                min_overlap.map_or(true, |min_overlap| {
                    reciprocal_overlap((record.start - 1)..record.stop, (sv.pos - 1)..sv.end)
                        >= min_overlap
                })
            })
            .filter(|record| {
                record.pathogenicity >= min_patho.unwrap_or(Pathogenicity::Benign) as i32
            })
            .map(|record| record.rcv)
            .collect()
    }
}

// Load the Clinvar SV databases from database given the configuration.
#[tracing::instrument]
pub fn load_clinvar_sv(
    path_db: &str,
    genome_release: GenomeRelease,
) -> Result<ClinvarSv, anyhow::Error> {
    info!("loading binary ClinVar SV dbs");

    let before_loading = std::time::Instant::now();
    let mut result = ClinvarSv::default();
    for _ in CHROMS {
        result.records.push(Vec::new());
        result.trees.push(IntervalTree::new());
    }

    let path =
        std::path::Path::new(path_db).join(format!("{}/strucvars/clinvar.bin", genome_release));
    let fcontents =
        std::fs::read(&path).map_err(|e| anyhow::anyhow!("error reading {:?}: {}", &path, e))?;
    let bg_db = crate::pbs::clinvar::SvDatabase::decode(std::io::Cursor::new(fcontents))
        .map_err(|e| anyhow::anyhow!("error decoding {:?}: {}", &path, e))?;

    let mut total_count = 0;
    for record in bg_db.records.into_iter() {
        let chrom_no = record.chrom_no as usize;
        let key = (record.start - 1)..record.stop;
        result.trees[chrom_no].insert(key, result.records[chrom_no].len() as u32);
        result.records[chrom_no].push(record);
        total_count += 1;
    }
    tracing::debug!(
        "done loading {} clinvar SV records from {:?} in {:?}",
        &total_count.separate_with_commas(),
        &path,
        before_loading.elapsed()
    );

    let before_building = std::time::Instant::now();
    result.trees.iter_mut().for_each(|tree| tree.index());
    tracing::debug!("done building itrees in {:?}", before_building.elapsed());

    Ok(result)
}
