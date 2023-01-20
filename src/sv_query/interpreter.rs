// Apply settings from a `CaseQuery` to `StructuralVariant` records.

use super::{
    dbrecords::SvOverlapCounts,
    schema::{CaseQuery, Range, StructuralVariant, SvSubType, SvType},
};

/// Slack around break-end positions
static BND_SLACK: u32 = 50;

/// Slack around insertion position
static INS_SLACK: u32 = 50;

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug)]
struct QueryInterpreter {
    query: CaseQuery,
}

/// Returns whether the intervals `[s1, e1)` and `[s2, e2)` overlap.
fn overlaps(s1: u32, e1: u32, s2: u32, e2: u32) -> bool {
    s1 < e2 && e1 > s2
}

impl QueryInterpreter {
    /// Construct new `QueryInterpreter` with the given query settings.
    pub fn new(query: CaseQuery) -> Self {
        QueryInterpreter { query }
    }

    /// Determine whether this record passes the genotype criteria.
    pub fn passes_genotype(&self, sv: &StructuralVariant) -> bool {
        let _ = sv;
        true
    }

    /// Determine whether this record passes the simple criteria regarding
    /// size and SV type.
    pub fn passes_simple(&self, sv: &StructuralVariant) -> bool {
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

        pass_sv_type && pass_sv_sub_type && pass_sv_size_min && pass_sv_size_max
    }

    /// Determine whether an SV record passes the genomic region criteria.
    pub fn passes_genomic_region(&self, sv: &StructuralVariant) -> bool {
        if let Some(regions) = &self.query.genomic_region {
            // interpret the allow list, any match is sufficient
            let mut any_match = false;
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

            any_match
        } else {
            true // no allow list given; always pass
        }
    }

    /// Determine whether an SV record with the given overlap counts passes
    /// the criteria.
    pub fn passes_counts(&self, counts: &SvOverlapCounts) -> bool {
        // We simply check for each database separately and pass if the check has not
        // been enabled or no minimal carrier / allele count is given
        let passes_dgv = !self.query.svdb_dgv_enabled
            || counts.dgv_carriers
                <= self
                    .query
                    .svdb_dgv_max_carriers
                    .unwrap_or(counts.dgv_carriers);
        let passes_dgv_gs = !self.query.svdb_dgv_gs_enabled
            || counts.dgv_gs_carriers
                <= self
                    .query
                    .svdb_dgv_gs_max_carriers
                    .unwrap_or(counts.dgv_gs_carriers);
        let passes_gnomad = !self.query.svdb_gnomad_enabled
            || counts.gnomad_carriers
                <= self
                    .query
                    .svdb_gnomad_max_carriers
                    .unwrap_or(counts.gnomad_carriers);
        let passes_exac = !self.query.svdb_exac_enabled
            || counts.exac_carriers
                <= self
                    .query
                    .svdb_exac_max_carriers
                    .unwrap_or(counts.exac_carriers);
        let passes_dbvar = !self.query.svdb_dbvar_enabled
            || counts.dbvar_carriers
                <= self
                    .query
                    .svdb_dbvar_max_carriers
                    .unwrap_or(counts.dbvar_carriers);
        let passes_g1k = !self.query.svdb_g1k_enabled
            || counts.g1k_alleles
                <= self
                    .query
                    .svdb_g1k_max_alleles
                    .unwrap_or(counts.g1k_alleles);
        let passes_inhouse = !self.query.svdb_inhouse_enabled
            || counts.inhouse_carriers
                <= self
                    .query
                    .svdb_inhouse_max_carriers
                    .unwrap_or(counts.inhouse_carriers);

        passes_dgv
            && passes_dgv_gs
            && passes_gnomad
            && passes_exac
            && passes_dbvar
            && passes_g1k
            && passes_inhouse
    }

    /// Determine whether the annotated `StructuralVariant` passes all criteria.
    pub fn passes(&self, sv: &StructuralVariant, counts: &SvOverlapCounts) -> bool {
        // simply AND-concatenate all `passes_*` functions
        self.passes_simple(sv)
            && self.passes_counts(counts)
            && self.passes_genomic_region(sv)
            && self.passes_genotype(sv)
    }

    // TODO: gene allow list
    // TODO: regulatory ensembl/vista
    // TODO: regulatory custom
    // TODO: tad boundary
    // TODO: recessive mdoe
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::sv_query::schema::{Database, GenomicRegion, StrandOrientation};

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
        let query = CaseQuery::new(Database::Refseq);
        let _interpreter = QueryInterpreter::new(query);
    }

    #[test]
    fn test_query_interpreter_passes_simple_pass_sv_type() {
        let query = CaseQuery {
            sv_types: vec![SvType::Del],
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_simple(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_sv_type() {
        let query = CaseQuery {
            sv_types: vec![SvType::Ins],
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_simple(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_sv_sub_type() {
        let query = CaseQuery {
            sv_sub_types: vec![SvSubType::Ins],
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_simple(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_pass_size() {
        let query = CaseQuery {
            sv_size_min: Some(50),
            sv_size_max: Some(500),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_simple(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_size_min() {
        let query = CaseQuery {
            sv_size_min: Some(5000),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_simple(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_simple_fail_size_max() {
        let query = CaseQuery {
            sv_size_max: Some(1),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_simple(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_linear_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 150, 160)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_linear_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_linear_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 201, 250)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_linear_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_ins_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 150, 160)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_ins_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_ins_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 201, 250)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_ins_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_overlap_pos() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 150, 160)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_whole_chromosome_pos() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr1")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_overlap_end() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr2", 150, 160)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_pass_bnd_whole_chromosome_end() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(interpreter.passes_genomic_region(&sv_pass));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_bnd_overlap() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::new(&"chr1", 201, 250)]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }

    #[test]
    fn test_query_interpreter_passes_genomic_region_fail_bnd_whole_chromosome() {
        let query = CaseQuery {
            genomic_region: Some(vec![GenomicRegion::whole_chrom("chr2")]),
            ..CaseQuery::new(Database::Refseq)
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
            call_info: HashMap::new(),
        };

        assert!(!interpreter.passes_genomic_region(&sv_fail));
    }
}
