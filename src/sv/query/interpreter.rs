//! Apply settings from a `CaseQuery` to `StructuralVariant` records.

use std::collections::HashSet;

use anyhow::anyhow;
use tracing::trace;

use super::{
    bgdbs::BgDbOverlaps,
    masked::MaskedBreakpointCount,
    schema::{CaseQuery, GenotypeChoice, Range, StructuralVariant, SvSubType, SvType},
};

/// Slack around break-end positions
pub static BND_SLACK: u32 = 50;

/// Slack around insertion position
pub static INS_SLACK: u32 = 50;

/// Returns whether the intervals `[s1, e1)` and `[s2, e2)` overlap.
fn overlaps(s1: u32, e1: u32, s2: u32, e2: u32) -> bool {
    s1 < e2 && e1 > s2
}

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug)]
pub struct QueryInterpreter {
    pub query: CaseQuery,
}

impl QueryInterpreter {
    /// Construct new `QueryInterpreter` with the given query settings.
    pub fn new(query: CaseQuery) -> Self {
        QueryInterpreter { query }
    }

    /// Determine whether this record passes the genotype criteria.
    pub fn passes_genotype(
        &self,
        sv: &StructuralVariant,
        masked_count: &MaskedBreakpointCount,
    ) -> Result<bool, anyhow::Error> {
        // Ensure that the sample sets in query and sv are the same
        let query_samples = self.query.genotype.keys().collect::<HashSet<&String>>();
        let sv_samples = sv.call_info.keys().collect::<HashSet<&String>>();
        if !query_samples.eq(&sv_samples) {
            return Err(anyhow!(
                "Samples in query and SV are not equal: {:?} vs {:?}",
                &query_samples,
                &sv_samples
            ));
        }

        // Now check whether for each sample, the selected genotype in `self.query.genotype`
        // matches what we have in terms of `CallInfo` for the sample in `sv`.  For this, we
        // go through all `GenotypeCriteria` in `self.query.genotype_criteria` and look
        // for all matching such records (by genotype, sv sub type, size)...
        for sample in query_samples {
            let query_genotype = *self
                .query
                .genotype
                .get(sample)
                .expect("cannot happen: checked for sample equality earlier");
            let call_info = sv
                .call_info
                .get(sample)
                .expect("cannot happen: checked for sample quality earlier");

            if query_genotype == GenotypeChoice::Any {
                // If any genotype is fine for this sample, we still check the breakpoint/masked info.
                // We check all criteria in this case and any passing is OK.
                let mut is_any_criteria_pass = false;
                for criteria in &self.query.genotype_criteria {
                    if criteria.is_masked_pass(masked_count) {
                        is_any_criteria_pass = true;
                        break;
                    }
                }
                if !is_any_criteria_pass {
                    return Ok(false);
                }
            } else {
                let mut is_any_criteria_pass = false;
                for criteria in &self.query.genotype_criteria {
                    if criteria.is_applicable_to(query_genotype, sv.sv_sub_type, sv.size())
                        && criteria.is_call_info_pass(call_info)
                        && criteria.is_masked_pass(masked_count)
                    {
                        is_any_criteria_pass = true;
                        break; // one matching criteria is enough
                    }
                }
                // If none of the criteria is pass for this sample then the sample's
                // CallInfo does not comply with the selected genotype.  This variant
                // does thus not pass the genotype filter.
                if !is_any_criteria_pass {
                    return Ok(false);
                }
            }
        }

        Ok(true)
    }

    /// Determine whether this record passes the selection criteria regarding
    /// size and SV type.
    pub fn passes_selection(&self, sv: &StructuralVariant) -> bool {
        let pass_sv_type =
            self.query.sv_types.is_empty() || self.query.sv_types.contains(&sv.sv_type);

        let pass_sv_sub_type =
            self.query.sv_sub_types.is_empty() || self.query.sv_sub_types.contains(&sv.sv_sub_type);

        let sv_size = sv.size();
        let (pass_sv_size_min, pass_sv_size_max) = if let Some(sv_size) = sv_size {
            let pass_sv_size_min = self
                .query
                .sv_size_min
                .map_or(true, |sv_size_min| sv_size >= sv_size_min);
            let pass_sv_size_max = self
                .query
                .sv_size_max
                .map_or(true, |sv_size_max| sv_size <= sv_size_max);
            (pass_sv_size_min, pass_sv_size_max)
        } else {
            (true, true)
        };

        trace!("does SV pass selection? pass_sv_type={} pass_sv_sub_type={} pass_sv_size_min={} pass_sv_size_max={}", pass_sv_type, pass_sv_sub_type, pass_sv_size_min, pass_sv_size_max);
        pass_sv_type && pass_sv_sub_type && pass_sv_size_min && pass_sv_size_max
    }

    /// Determine whether an SV record passes the genomic region criteria.
    pub fn passes_genomic_region(&self, sv: &StructuralVariant) -> bool {
        if let Some(regions) = &self.query.genomic_region {
            // interpret the allow list, any match is sufficient
            let mut any_match = false;

            if regions.is_empty() {
                trace!("no genomic region allow list given, pass");
                any_match = true;
            }

            if sv.sv_type == SvType::Ins || sv.sv_sub_type.is_ins() {
                // handle case of insertions: overlap position with `INS_SLACK` and region
                for region in regions {
                    // as for all others, the range matches if `None` (whole chrom) or has overlap
                    let range_matches = match region.range {
                        None => true,
                        Some(Range { start, end }) => overlaps(
                            start.saturating_sub(1),
                            end,
                            sv.pos.saturating_sub(INS_SLACK),
                            sv.pos + INS_SLACK,
                        ),
                    };
                    any_match = any_match || (region.chrom.eq(&sv.chrom) && range_matches);
                }
            } else if sv.sv_type == SvType::Bnd || sv.sv_sub_type == SvSubType::Bnd {
                // for break-ends, test both ends and use `BND_SLACK`
                for region in regions {
                    // as for all others, the range matches if `None` (whole chrom) or has overlap
                    let range_matches_chrom = match region.range {
                        None => true,
                        Some(Range { start, end }) => overlaps(
                            start.saturating_sub(1),
                            end,
                            sv.pos.saturating_sub(BND_SLACK),
                            sv.pos + BND_SLACK,
                        ),
                    };
                    let range_matches_chrom2 = match region.range {
                        None => true,
                        Some(Range { start, end }) => overlaps(
                            start.saturating_sub(1),
                            end,
                            sv.end.saturating_sub(BND_SLACK),
                            sv.end + BND_SLACK,
                        ),
                    };
                    any_match = any_match
                        || (region.chrom.eq(&sv.chrom) && range_matches_chrom)
                        || (sv
                            .chrom2
                            .as_ref()
                            .map_or(false, |chrom2| chrom2.eq(&region.chrom))
                            && range_matches_chrom2);
                }
            } else {
                // handle the case of linear structural variants
                for region in regions {
                    // as for all others, the range matches if `None` (whole chrom) or has overlap
                    let range_matches = match region.range {
                        None => true,
                        Some(Range { start, end }) => overlaps(
                            start.saturating_sub(1),
                            end,
                            sv.pos.saturating_sub(1),
                            sv.end,
                        ),
                    };
                    any_match = any_match || (region.chrom.eq(&sv.chrom) && range_matches);
                }
            }

            trace!("does SV pass genomic region? any_match={}", any_match);
            any_match
        } else {
            trace!("no genomic region allow list given, pass");
            true // no allow list given; always pass
        }
    }

    /// Determine whether an SV record with the given overlap counts passes
    /// the criteria.
    pub fn passes_counts(&self, counts: &BgDbOverlaps) -> bool {
        // We simply check for each database separately and pass if the check has not
        // been enabled or no minimal carrier / allele count is given
        let passes_dgv = !self.query.svdb_dgv_enabled
            || counts.dgv <= self.query.svdb_dgv_max_count.unwrap_or(counts.dgv);
        let passes_dgv_gs = !self.query.svdb_dgv_gs_enabled
            || counts.dgv_gs <= self.query.svdb_dgv_gs_max_count.unwrap_or(counts.dgv_gs);
        let passes_gnomad = !self.query.svdb_gnomad_enabled
            || counts.gnomad <= self.query.svdb_gnomad_max_count.unwrap_or(counts.gnomad);
        let passes_exac = !self.query.svdb_exac_enabled
            || counts.exac <= self.query.svdb_exac_max_count.unwrap_or(counts.exac);
        let passes_dbvar = !self.query.svdb_dbvar_enabled
            || counts.dbvar <= self.query.svdb_dbvar_max_count.unwrap_or(counts.dbvar);
        let passes_g1k = !self.query.svdb_g1k_enabled
            || counts.g1k <= self.query.svdb_g1k_max_count.unwrap_or(counts.g1k);
        let passes_inhouse = !self.query.svdb_inhouse_enabled
            || counts.inhouse <= self.query.svdb_inhouse_max_count.unwrap_or(counts.inhouse);

        trace!("does SV pass counts? passes_dgv={}, passes_dgv_gs={}, passes_gnomad={}, passes_exac={}, passes_dbvar={}, passes_g1k={}, passes_inhouse={}",
        passes_dgv
        ,passes_dgv_gs
        ,passes_gnomad
        ,passes_exac
        ,passes_dbvar
        ,passes_g1k
        ,passes_inhouse);

        passes_dgv
            && passes_dgv_gs
            && passes_gnomad
            && passes_exac
            && passes_dbvar
            && passes_g1k
            && passes_inhouse
    }

    /// Determine whether the annotated `StructuralVariant` passes all criteria.
    pub fn passes<CountBg, CountMasked>(
        &self,
        sv: &StructuralVariant,
        count_bg: &mut CountBg,
        count_masked: &mut CountMasked,
    ) -> Result<bool, anyhow::Error>
    where
        CountBg: FnMut(&StructuralVariant) -> BgDbOverlaps,
        CountMasked: FnMut(&StructuralVariant) -> MaskedBreakpointCount,
    {
        // We first check for matching genotype.  If this succeeds then we execute the
        // overlapper for known pathogenic and then for frequency in background.
        trace!("checking whether SV passes filter");
        if !self.passes_selection(sv)
            || !self.passes_genomic_region(sv)
            || !self.passes_genotype(sv, &count_masked(sv))?
        {
            trace!("... SV does not pass filter");
            Ok(false)
        } else {
            trace!("... SV passes filter");
            Ok(self.passes_counts(&count_bg(sv)))
        }
    }

    // TODO: gene allow list
    // TODO: regulatory ensembl/vista
    // TODO: regulatory custom
    // TODO: recessive mdoe
}

#[cfg(test)]
mod tests {
    use indexmap::IndexMap;

    use crate::{
        db::conf::Database,
        sv::query::schema::{
            CallInfo, GenomicRegion, GenotypeChoice, GenotypeCriteria, StrandOrientation,
        },
    };

    use super::*;

    #[test]
    fn test_overlaps() {
        assert!(overlaps(1, 10, 1, 10));
        assert!(overlaps(1, 10, 9, 20));
        assert!(!overlaps(1, 10, 10, 20));
        assert!(overlaps(1, 10, 1, 10));
        assert!(overlaps(9, 20, 1, 10));
        assert!(!overlaps(10, 20, 1, 10));
    }

    #[test]
    fn test_query_interpreter_smoke() {
        let query = CaseQuery::new(Database::RefSeq);
        let _interpreter = QueryInterpreter::new(query);
    }

    #[test]
    fn test_query_interpreter_passes_simple_pass_sv_type() {
        let query = CaseQuery {
            sv_types: vec![SvType::Del],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_selection(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_sv_type() {
        let query = CaseQuery {
            sv_types: vec![SvType::Ins],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_selection(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_sv_sub_type() {
        let query = CaseQuery {
            sv_sub_types: vec![SvSubType::Ins],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_selection(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_pass_size() {
        let query = CaseQuery {
            sv_size_min: Some(50),
            sv_size_max: Some(500),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_selection(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_size_min() {
        let query = CaseQuery {
            sv_size_min: Some(5000),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_selection(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_size_max() {
        let query = CaseQuery {
            sv_size_max: Some(1),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_selection(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_linear_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 150, 160)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_linear_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_linear_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 201, 250)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_linear_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_ins_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 150, 160)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_ins_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_ins_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 201, 250)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_ins_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_overlap_pos() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 150, 160)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 150,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr1".to_owned()),
            end: 1000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_whole_chromosome_pos() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr1".to_owned()),
            end: 1000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_overlap_end() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr2", 150, 160)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr2".to_owned()),
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_whole_chromosome_end() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr2".to_owned()),
            end: 1000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_bnd_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new("chr1", 201, 250)]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr1".to_owned()),
            end: 1000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_bnd_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr1".to_owned()),
            end: 1000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_counts_pass() {
        let query = CaseQuery {
            svdb_dgv_enabled: true,
            svdb_dgv_max_count: Some(10),
            svdb_dgv_gs_enabled: true,
            svdb_dgv_gs_max_count: Some(10),
            svdb_gnomad_enabled: true,
            svdb_gnomad_max_count: Some(10),
            svdb_exac_enabled: true,
            svdb_exac_max_count: Some(10),
            svdb_dbvar_enabled: true,
            svdb_dbvar_max_count: Some(10),
            svdb_g1k_enabled: true,
            svdb_g1k_max_count: Some(10),
            svdb_inhouse_enabled: true,
            svdb_inhouse_max_count: Some(10),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let counts_pass = BgDbOverlaps {
            dgv: 5,
            dgv_gs: 5,
            gnomad: 5,
            exac: 5,
            g1k: 5,
            inhouse: 5,
            dbvar: 5,
        };

        assert!(interpreter.passes_counts(&counts_pass));
    }

    #[test]
    fn test_query_interpreter_passes_counts_fail() {
        let query = CaseQuery {
            svdb_dgv_enabled: true,
            svdb_dgv_max_count: Some(10),
            svdb_dgv_gs_enabled: true,
            svdb_dgv_gs_max_count: Some(10),
            svdb_gnomad_enabled: true,
            svdb_gnomad_max_count: Some(10),
            svdb_exac_enabled: true,
            svdb_exac_max_count: Some(10),
            svdb_dbvar_enabled: true,
            svdb_dbvar_max_count: Some(10),
            svdb_g1k_enabled: true,
            svdb_g1k_max_count: Some(10),
            svdb_inhouse_enabled: true,
            svdb_inhouse_max_count: Some(10),
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let counts_fail = BgDbOverlaps {
            dgv: 11,
            dgv_gs: 11,
            gnomad: 11,
            exac: 11,
            g1k: 11,
            inhouse: 11,
            dbvar: 11,
        };

        assert!(!interpreter.passes_counts(&counts_fail));
    }

    #[test]
    fn test_query_interpreter_pass_genotype_fail_no_match() -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            genotype: IndexMap::from([("sample".to_owned(), GenotypeChoice::Het)]),
            genotype_criteria: vec![GenotypeCriteria {
                select_sv_sub_type: vec![SvSubType::Del],
                select_sv_min_size: Some(1000),
                select_sv_max_size: Some(5000),
                gt_one_of: Some(vec![
                    "0/1".to_owned(),
                    "0|1".to_owned(),
                    "1/0".to_owned(),
                    "1|0".to_owned(),
                ]),
                min_gq: Some(5.0),
                min_pr_cov: Some(10),
                min_pr_var: Some(5),
                min_sr_cov: Some(10),
                min_sr_var: Some(5),
                ..GenotypeCriteria::new(GenotypeChoice::Het)
            }],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let call_info = CallInfo {
            genotype: Some("0/1".to_owned()),
            quality: Some(10.0),
            paired_end_cov: Some(10),
            paired_end_var: Some(5),
            split_read_cov: Some(10),
            split_read_var: Some(5),
            copy_number: Some(1),
            average_normalized_cov: Some(0.5),
            point_count: Some(10),
            average_mapping_quality: Some(60.0),
        };

        let sv_fail = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 1000,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 2000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::from([("sample".to_owned(), call_info.clone())]),
        };

        // The following tests fail because the SV does not match the match criteria for
        // dfiferent reasons.

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                end: 1100,
                ..sv_fail.clone()
            },
            &Default::default()
        )?); // too small to match

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                end: 10000,
                ..sv_fail.clone()
            },
            &Default::default()
        )?); // too large to match

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                sv_type: SvType::Dup,
                sv_sub_type: SvSubType::Dup,
                ..sv_fail.clone()
            },
            &Default::default()
        )?); // wrong sv sub type

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                call_info: IndexMap::from([(
                    "sample".to_owned(),
                    CallInfo {
                        quality: Some(1.0),
                        ..call_info.clone()
                    },
                )]),
                ..sv_fail.clone()
            },
            &Default::default()
        )?); // quality too low

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                call_info: IndexMap::from([(
                    "sample".to_owned(),
                    CallInfo {
                        genotype: Some("0/0".to_owned()),
                        paired_end_cov: Some(1),
                        ..call_info
                    },
                )]),
                ..sv_fail.clone()
            },
            &Default::default()
        )?); // pr coverage too low

        assert!(!interpreter.passes_genotype(
            &StructuralVariant {
                call_info: IndexMap::from([(
                    "sample".to_owned(),
                    CallInfo {
                        genotype: Some("0/0".to_owned()),
                        split_read_cov: Some(1),
                        ..call_info
                    },
                )]),
                ..sv_fail
            },
            &Default::default()
        )?); // sr coverage too low

        Ok(())
    }

    #[test]
    fn test_query_interpreter_ins_min_pr() -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            genotype: IndexMap::from([("sample".to_owned(), GenotypeChoice::Variant)]),
            genotype_criteria: vec![GenotypeCriteria {
                select_sv_sub_type: vec![SvSubType::Ins],
                select_sv_min_size: Some(1000),
                select_sv_max_size: Some(5000),
                gt_one_of: Some(vec![
                    "0/1".to_owned(),
                    "0|1".to_owned(),
                    "1/0".to_owned(),
                    "1|0".to_owned(),
                ]),
                min_gq: Some(5.0),
                // min_pr_cov: Some(10),
                // min_pr_var: Some(1),
                // min_sr_cov: Some(10),
                min_sr_var: Some(1),
                ..GenotypeCriteria::new(GenotypeChoice::Variant)
            }],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let sv = StructuralVariant {
            chrom: "1".to_owned(),
            pos: 12345,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 12345,
            strand_orientation: Some(StrandOrientation::NotApplicable),
            call_info: IndexMap::from([(
                "sample".to_owned(),
                CallInfo {
                    genotype: Some("0/1".to_owned()),
                    quality: Some(16.0),
                    paired_end_cov: Some(0),
                    paired_end_var: Some(0),
                    split_read_cov: Some(5),
                    split_read_var: Some(1),
                    copy_number: None,
                    average_normalized_cov: None,
                    point_count: None,
                    average_mapping_quality: None,
                },
            )]),
        };

        assert!(interpreter.passes_genotype(&sv, &Default::default())?);

        Ok(())
    }

    #[test]
    fn test_query_interpreter_pass_genotype_pass_smoke() -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            genotype: IndexMap::from([("sample".to_owned(), GenotypeChoice::Het)]),
            genotype_criteria: vec![GenotypeCriteria {
                select_sv_sub_type: vec![SvSubType::Del],
                select_sv_min_size: Some(1000),
                select_sv_max_size: Some(5000),
                gt_one_of: Some(vec![
                    "0/1".to_owned(),
                    "0|1".to_owned(),
                    "1/0".to_owned(),
                    "1|0".to_owned(),
                ]),
                min_gq: Some(5.0),
                min_pr_cov: Some(10),
                min_pr_var: Some(5),
                min_sr_cov: Some(10),
                min_sr_var: Some(5),
                ..GenotypeCriteria::new(GenotypeChoice::Het)
            }],
            ..CaseQuery::new(Database::RefSeq)
        };
        let interpreter = QueryInterpreter::new(query);

        let call_info = CallInfo {
            genotype: Some("0/1".to_owned()),
            quality: Some(10.0),
            paired_end_cov: Some(10),
            paired_end_var: Some(5),
            split_read_cov: Some(10),
            split_read_var: Some(5),
            copy_number: Some(1),
            average_normalized_cov: Some(0.5),
            point_count: Some(10),
            average_mapping_quality: Some(60.0),
        };

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 1000,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 2000,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::from([("sample".to_owned(), call_info)]),
        };

        assert!(interpreter.passes_genotype(&sv_pass, &Default::default())?);
        Ok(())
    }

    #[test]
    fn test_query_interpreter_passes_smoke() -> Result<(), anyhow::Error> {
        let query = CaseQuery::new(Database::RefSeq);
        let interpreter = QueryInterpreter::new(query);

        let sv_pass = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::Del,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: IndexMap::new(),
        };
        let counts_pass = BgDbOverlaps {
            dgv: 5,
            dgv_gs: 5,
            gnomad: 5,
            exac: 5,
            g1k: 5,
            inhouse: 5,
            dbvar: 5,
        };

        assert!(
            interpreter.passes(&sv_pass, &mut |_sv| counts_pass.clone(), &mut |_sv| {
                Default::default()
            })?
        );

        Ok(())
    }
}
