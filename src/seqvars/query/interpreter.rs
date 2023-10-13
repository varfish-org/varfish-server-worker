//! Apply settings from a `strucvar::query::schema::CaseQuery` to `SequenceVariant` records.

use std::collections::HashSet;

use super::schema::{CaseQuery, SequenceVariant};

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug, Default)]
pub struct QueryInterpreter {
    /// The case query settings.
    pub query: CaseQuery,
    /// Gene allowlist with HGNC IDs.
    pub hgvs_allowlist: Option<HashSet<String>>,
}

/// Result type for `QueryInterpreter::passes_genotype()`.
#[derive(Debug, Default)]
pub struct PassesResult {
    /// Whether genotype passes for all samples.
    pub pass_all: bool,
}

impl QueryInterpreter {
    /// Construct new `QueryInterpreter` with the given query settings.
    pub fn new(query: CaseQuery, hgvs_allowlist: Option<HashSet<String>>) -> Self {
        tracing::error!(
            "note well that we will need a second pass for compound heterozygous variants"
        );
        QueryInterpreter {
            query,
            hgvs_allowlist,
        }
    }

    /// Determine whether the annotated `SequenceVariant` passes all criteria.
    pub fn passes(&self, seqvar: &SequenceVariant) -> Result<PassesResult, anyhow::Error> {
        let pass_frequency = self.passes_frequency(seqvar)?;
        let pass_effects = self.passes_effects(seqvar)?;
        let pass_quality = self.passes_quality(seqvar)?;
        let pass_genotype = self.passes_genotype(seqvar)?;
        let pass_genes_regions = self.passes_genes_regions(seqvar)?;
        let pass_clinvar = self.passes_clinvar(seqvar)?;
        let pass_all = pass_frequency
            && pass_effects
            && pass_quality
            && pass_genotype
            && pass_genes_regions
            && pass_clinvar;
        Ok(PassesResult { pass_all })
    }

    /// Determine whether the `SequenceVariant` passes the frequency filter.
    fn passes_frequency(&self, s: &SequenceVariant) -> Result<bool, anyhow::Error> {
        let q = &self.query;
        let is_mtdna = annonars::common::cli::canonicalize(&s.chrom) == "MT";

        if is_mtdna {
            if q.helixmtdb_enabled {
                if q.helixmtdb_frequency.is_some()
                    && s.helixmtdb_af() > q.helixmtdb_frequency.expect("tested before")
                    || q.helixmtdb_heteroplasmic.is_some()
                        && s.helix_het > q.helixmtdb_heteroplasmic.expect("tested before")
                    || q.helixmtdb_homoplasmic.is_some()
                        && s.helix_hom > q.helixmtdb_homoplasmic.expect("tested before")
                {
                    tracing::trace!("variant {:?} fails HelixMtDb frequency filter {:?}", s, &q);
                    return Ok(false);
                }
            }
        } else {
            if q.gnomad_exomes_enabled {
                if q.gnomad_exomes_frequency.is_some()
                    && s.gnomad_exomes_af() > q.gnomad_exomes_frequency.expect("tested before")
                    || q.gnomad_exomes_heterozygous.is_some()
                        && s.gnomad_exomes_het
                            > q.gnomad_exomes_heterozygous.expect("tested before")
                    || q.gnomad_exomes_homozygous.is_some()
                        && s.gnomad_exomes_hom > q.gnomad_exomes_homozygous.expect("tested before")
                    || q.gnomad_exomes_hemizygous.is_some()
                        && s.gnomad_exomes_hemi > q.gnomad_exomes_hemizygous.expect("tested before")
                {
                    tracing::trace!(
                        "variant {:?} fails gnomAD exomes frequency filter {:?}",
                        s,
                        &q.gnomad_exomes_frequency
                    );
                    return Ok(false);
                }
            }
        }

        if q.gnomad_genomes_enabled {
            if q.gnomad_genomes_frequency.is_some()
                && s.gnomad_genomes_af() > q.gnomad_genomes_frequency.expect("tested before")
                || q.gnomad_genomes_heterozygous.is_some()
                    && s.gnomad_genomes_het > q.gnomad_genomes_heterozygous.expect("tested before")
                || q.gnomad_genomes_homozygous.is_some()
                    && s.gnomad_genomes_hom > q.gnomad_genomes_homozygous.expect("tested before")
                || !is_mtdna
                    && q.gnomad_genomes_hemizygous.is_some()
                    && s.gnomad_genomes_hemi > q.gnomad_genomes_hemizygous.expect("tested before")
            {
                tracing::trace!(
                    "variant {:?} fails gnomAD genomes frequency filter {:?}",
                    s,
                    &q.gnomad_genomes_frequency
                );
                return Ok(false);
            }
        }

        Ok(true)
    }

    fn passes_effects(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }

    fn passes_quality(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }

    fn passes_genotype(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }

    fn passes_genes_regions(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }

    fn passes_clinvar(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};

    use super::QueryInterpreter;

    #[rstest]
    // -- frequency ---------------------------------------------------------
    // frequency: pass [het count] (no filter value)
    #[case(1000, 1, 0, 0, true, None, None, None, None, true)]
    // frequency: pass [het count]
    #[case(1000, 1, 0, 0, true, Some(0.001), None, None, None, true)]
    // frequency: fail [het count]
    #[case(1000, 2, 0, 0, true, Some(0.001), None, None, None, false)]
    // frequency: pass [het count] (fail but filter is disabled)
    #[case(1000, 2, 0, 0, false, Some(0.001), None, None, None, true)]
    // frequency: pass [hom count] (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // frequency: pass [hom count]
    #[case(1000, 0, 1, 0, true, Some(0.002), None, None, None, true)]
    // frequency: fail [hom count]
    #[case(1000, 0, 2, 0, true, Some(0.002), None, None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 2, 0, false, Some(0.002), None, None, None, true)]
    // frequency: pass [hemi count] (no filter value)
    #[case(1000, 0, 0, 1, true, None, None, None, None, true)]
    // frequency: pass [hemi count]
    #[case(1000, 0, 0, 1, true, Some(0.001), None, None, None, true)]
    // frequency: fail [hemi count]
    #[case(1000, 0, 0, 2, true, Some(0.001), None, None, None, false)]
    // frequency: pass [hemi count] (fail but filter is disabled)
    #[case(1000, 0, 0, 2, false, Some(0.001), None, None, None, true)]
    // -- heterezygous count ------------------------------------------------
    // het. count: pass (no filter value)
    #[case(1000, 1, 0, 0, true, None, None, None, None, true)]
    // het. count: pass
    #[case(1000, 1, 0, 0, true, None, Some(1), None, None, true)]
    // het. count: fail
    #[case(1000, 2, 0, 0, true, None, Some(1), None, None, false)]
    // het. count: pass (fail but filter is disabled)
    #[case(1000, 2, 0, 0, false, None, Some(1), None, None, true)]
    // -- homozygous count --------------------------------------------------
    // hom. count: pass (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // hom. count: pass
    #[case(1000, 0, 1, 0, true, None, None, Some(1), None, true)]
    // hom. count: fail
    #[case(1000, 0, 2, 0, true, None, None, Some(1), None, false)]
    // hom. count: pass (fail but filter is disabled)
    #[case(1000, 0, 2, 0, false, None, None, Some(1), None, true)]
    // -- hemizygous count --------------------------------------------------
    // hemi. count: pass (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // hemi. count: pass
    #[case(1000, 0, 0, 1, true, None, None, None, Some(1), true)]
    // hemi. count: fail
    #[case(1000, 0, 0, 2, true, None, None, None, Some(1), false)]
    // hemi. count: pass (fail but filter is disabled)
    #[case(1000, 0, 0, 2, false, None, None, None, Some(1), true)]
    fn passes_frequency_gnomad_exomes_nuclear_dna(
        #[case] seqvar_gnomad_exomes_an: i32,
        #[case] seqvar_gnomad_exomes_het: i32,
        #[case] seqvar_gnomad_exomes_hom: i32,
        #[case] seqvar_gnomad_exomes_hemi: i32,
        #[case] query_gnomad_exomes_enabled: bool,
        #[case] query_gnomad_exomes_frequency: Option<f32>,
        #[case] query_gnomad_exomes_heterozygous: Option<i32>,
        #[case] query_gnomad_exomes_homozygous: Option<i32>,
        #[case] query_gnomad_exomes_hemizygous: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let interpreter = QueryInterpreter {
            query: CaseQuery {
                gnomad_exomes_enabled: query_gnomad_exomes_enabled,
                gnomad_exomes_frequency: query_gnomad_exomes_frequency,
                gnomad_exomes_heterozygous: query_gnomad_exomes_heterozygous,
                gnomad_exomes_homozygous: query_gnomad_exomes_homozygous,
                gnomad_exomes_hemizygous: query_gnomad_exomes_hemizygous,
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            gnomad_exomes_an: seqvar_gnomad_exomes_an,
            gnomad_exomes_het: seqvar_gnomad_exomes_het,
            gnomad_exomes_hom: seqvar_gnomad_exomes_hom,
            gnomad_exomes_hemi: seqvar_gnomad_exomes_hemi,
            chrom: "X".to_string(),
            ..Default::default()
        };

        assert_eq!(interpreter.passes(&seq_var)?.pass_all, expected_pass_all);

        Ok(())
    }

    #[rstest]
    // -- frequency ---------------------------------------------------------
    // frequency: pass [het count] (no filter value)
    #[case(1000, 1, 0, 0, true, None, None, None, None, true)]
    // frequency: pass [het count]
    #[case(1000, 1, 0, 0, true, Some(0.001), None, None, None, true)]
    // frequency: fail [het count]
    #[case(1000, 2, 0, 0, true, Some(0.001), None, None, None, false)]
    // frequency: pass [het count] (fail but filter is disabled)
    #[case(1000, 2, 0, 0, false, Some(0.001), None, None, None, true)]
    // frequency: pass [hom count] (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // frequency: pass [hom count]
    #[case(1000, 0, 1, 0, true, Some(0.002), None, None, None, true)]
    // frequency: fail [hom count]
    #[case(1000, 0, 2, 0, true, Some(0.002), None, None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 2, 0, false, Some(0.002), None, None, None, true)]
    // frequency: pass [hemi count] (no filter value)
    #[case(1000, 0, 0, 1, true, None, None, None, None, true)]
    // frequency: pass [hemi count]
    #[case(1000, 0, 0, 1, true, Some(0.001), None, None, None, true)]
    // frequency: fail [hemi count]
    #[case(1000, 0, 0, 2, true, Some(0.001), None, None, None, false)]
    // frequency: pass [hemi count] (fail but filter is disabled)
    #[case(1000, 0, 0, 2, false, Some(0.001), None, None, None, true)]
    // -- heterezygous count ------------------------------------------------
    // het. count: pass (no filter value)
    #[case(1000, 1, 0, 0, true, None, None, None, None, true)]
    // het. count: pass
    #[case(1000, 1, 0, 0, true, None, Some(1), None, None, true)]
    // het. count: fail
    #[case(1000, 2, 0, 0, true, None, Some(1), None, None, false)]
    // het. count: pass (fail but filter is disabled)
    #[case(1000, 2, 0, 0, false, None, Some(1), None, None, true)]
    // -- homozygous count --------------------------------------------------
    // hom. count: pass (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // hom. count: pass
    #[case(1000, 0, 1, 0, true, None, None, Some(1), None, true)]
    // hom. count: fail
    #[case(1000, 0, 2, 0, true, None, None, Some(1), None, false)]
    // hom. count: pass (fail but filter is disabled)
    #[case(1000, 0, 2, 0, false, None, None, Some(1), None, true)]
    // -- hemizygous count --------------------------------------------------
    // hemi. count: pass (no filter value)
    #[case(1000, 0, 1, 0, true, None, None, None, None, true)]
    // hemi. count: pass
    #[case(1000, 0, 0, 1, true, None, None, None, Some(1), true)]
    // hemi. count: fail
    #[case(1000, 0, 0, 2, true, None, None, None, Some(1), false)]
    // hemi. count: pass (fail but filter is disabled)
    #[case(1000, 0, 0, 2, false, None, None, None, Some(1), true)]
    fn passes_frequency_gnomad_genomes_nuclear_dna(
        #[case] seqvar_gnomad_genomes_an: i32,
        #[case] seqvar_gnomad_genomes_het: i32,
        #[case] seqvar_gnomad_genomes_hom: i32,
        #[case] seqvar_gnomad_genomes_hemi: i32,
        #[case] query_gnomad_genomes_enabled: bool,
        #[case] query_gnomad_genomes_frequency: Option<f32>,
        #[case] query_gnomad_genomes_heterozygous: Option<i32>,
        #[case] query_gnomad_genomes_homozygous: Option<i32>,
        #[case] query_gnomad_genomes_hemizygous: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let interpreter = QueryInterpreter {
            query: CaseQuery {
                gnomad_genomes_enabled: query_gnomad_genomes_enabled,
                gnomad_genomes_frequency: query_gnomad_genomes_frequency,
                gnomad_genomes_heterozygous: query_gnomad_genomes_heterozygous,
                gnomad_genomes_homozygous: query_gnomad_genomes_homozygous,
                gnomad_genomes_hemizygous: query_gnomad_genomes_hemizygous,
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            gnomad_genomes_an: seqvar_gnomad_genomes_an,
            gnomad_genomes_het: seqvar_gnomad_genomes_het,
            gnomad_genomes_hom: seqvar_gnomad_genomes_hom,
            gnomad_genomes_hemi: seqvar_gnomad_genomes_hemi,
            chrom: "X".to_string(),
            ..Default::default()
        };

        assert_eq!(interpreter.passes(&seq_var)?.pass_all, expected_pass_all);

        Ok(())
    }

    #[rstest]
    // -- frequency ---------------------------------------------------------
    // frequency: pass [het count] (no filter value)
    #[case(1000, 1, 0, true, None, None, None, true)]
    // frequency: pass [het count]
    #[case(1000, 1, 0, true, Some(0.001), None, None, true)]
    // frequency: fail [het count]
    #[case(1000, 2, 0, true, Some(0.001), None, None, false)]
    // frequency: pass [het count] (fail but filter is disabled)
    #[case(1000, 2, 0, false, Some(0.001), None, None, true)]
    // frequency: pass [hom count] (no filter value)
    #[case(1000, 0, 1, true, None, None, None, true)]
    // frequency: pass [hom count]
    #[case(1000, 0, 1, true, Some(0.002), None, None, true)]
    // frequency: fail [hom count]
    #[case(1000, 0, 2, true, Some(0.002), None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 2, false, Some(0.002), None, None, true)]
    // -- heteroplasmy count ------------------------------------------------
    // het. count: pass (no filter value)
    #[case(1000, 1, 0, true, None, None, None, true)]
    // het. count: pass
    #[case(1000, 1, 0, true, None, Some(1), None, true)]
    // het. count: fail
    #[case(1000, 2, 0, true, None, Some(1), None, false)]
    // het. count: pass (fail but filter is disabled)
    #[case(1000, 2, 0, false, None, Some(1), None, true)]
    // -- homoplasmy count --------------------------------------------------
    // hom. count: pass (no filter value)
    #[case(1000, 0, 1, true, None, None, None, true)]
    // hom. count: pass
    #[case(1000, 0, 1, true, None, None, Some(1), true)]
    // hom. count: fail
    #[case(1000, 0, 2, true, None, None, Some(1), false)]
    // hom. count: pass (fail but filter is disabled)
    #[case(1000, 0, 2, false, None, None, Some(1), true)]
    fn passes_frequency_helix_chrmt(
        #[case] seqvar_helix_an: i32,
        #[case] seqvar_helix_het: i32,
        #[case] seqvar_helix_hom: i32,
        #[case] query_helix_enabled: bool,
        #[case] query_helix_frequency: Option<f32>,
        #[case] query_helix_heteroplasmic: Option<i32>,
        #[case] query_helix_homoplasmic: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let interpreter = QueryInterpreter {
            query: CaseQuery {
                helixmtdb_enabled: query_helix_enabled,
                helixmtdb_frequency: query_helix_frequency,
                helixmtdb_heteroplasmic: query_helix_heteroplasmic,
                helixmtdb_homoplasmic: query_helix_homoplasmic,
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            helix_an: seqvar_helix_an,
            helix_het: seqvar_helix_het,
            helix_hom: seqvar_helix_hom,
            chrom: "MT".to_string(),
            ..Default::default()
        };

        assert_eq!(interpreter.passes(&seq_var)?.pass_all, expected_pass_all);

        Ok(())
    }

    #[rstest]
    // -- frequency ---------------------------------------------------------
    // frequency: pass [het count] (no filter value)
    #[case(1000, 1, 0, true, None, None, None, true)]
    // frequency: pass [het count]
    #[case(1000, 1, 0, true, Some(0.001), None, None, true)]
    // frequency: fail [het count]
    #[case(1000, 2, 0, true, Some(0.001), None, None, false)]
    // frequency: pass [het count] (fail but filter is disabled)
    #[case(1000, 2, 0, false, Some(0.001), None, None, true)]
    // frequency: pass [hom count] (no filter value)
    #[case(1000, 0, 1, true, None, None, None, true)]
    // frequency: pass [hom count]
    #[case(1000, 0, 1, true, Some(0.002), None, None, true)]
    // frequency: fail [hom count]
    #[case(1000, 0, 2, true, Some(0.002), None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 2, false, Some(0.002), None, None, true)]
    // -- heteroplasmy count ------------------------------------------------
    // het. count: pass (no filter value)
    #[case(1000, 1, 0, true, None, None, None, true)]
    // het. count: pass
    #[case(1000, 1, 0, true, None, Some(1), None, true)]
    // het. count: fail
    #[case(1000, 2, 0, true, None, Some(1), None, false)]
    // het. count: pass (fail but filter is disabled)
    #[case(1000, 2, 0, false, None, Some(1), None, true)]
    // -- homoplasmy count --------------------------------------------------
    // hom. count: pass (no filter value)
    #[case(1000, 0, 1, true, None, None, None, true)]
    // hom. count: pass
    #[case(1000, 0, 1, true, None, None, Some(1), true)]
    // hom. count: fail
    #[case(1000, 0, 2, true, None, None, Some(1), false)]
    // hom. count: pass (fail but filter is disabled)
    #[case(1000, 0, 2, false, None, None, Some(1), true)]
    fn passes_frequency_gnomad_genomes_chrmt(
        #[case] seqvar_gnomad_genomes_an: i32,
        #[case] seqvar_gnomad_genomes_het: i32,
        #[case] seqvar_gnomad_genomes_hom: i32,
        #[case] query_gnomad_genomes_enabled: bool,
        #[case] query_gnomad_genomes_frequency: Option<f32>,
        #[case] query_gnomad_genomes_heteroplasmic: Option<i32>,
        #[case] query_gnomad_genomes_homoplasmic: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let interpreter = QueryInterpreter {
            query: CaseQuery {
                gnomad_genomes_enabled: query_gnomad_genomes_enabled,
                gnomad_genomes_frequency: query_gnomad_genomes_frequency,
                gnomad_genomes_heterozygous: query_gnomad_genomes_heteroplasmic,
                gnomad_genomes_homozygous: query_gnomad_genomes_homoplasmic,
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            gnomad_genomes_an: seqvar_gnomad_genomes_an,
            gnomad_genomes_het: seqvar_gnomad_genomes_het,
            gnomad_genomes_hom: seqvar_gnomad_genomes_hom,
            chrom: "MT".to_string(),
            ..Default::default()
        };

        assert_eq!(interpreter.passes(&seq_var)?.pass_all, expected_pass_all);

        Ok(())
    }
}
