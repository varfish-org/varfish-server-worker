use crate::seqvars::query::schema::{
    data::{Af, VariantRecord},
    query::CaseQuery,
};

/// Determine whether the `VariantRecord` passes the frequency filter.
pub fn passes(query: &CaseQuery, s: &VariantRecord) -> Result<bool, anyhow::Error> {
    let frequency = &query.frequency;
    let is_mtdna = annonars::common::cli::canonicalize(&s.vcf_variant.chrom) == "MT";

    if is_mtdna {
        if frequency.helixmtdb.enabled
            && (frequency.helixmtdb.max_af.is_some()
                && s.population_frequencies.helixmtdb.af()
                    > frequency.helixmtdb.max_af.expect("tested before")
                || frequency.helixmtdb.max_het.is_some()
                    && s.population_frequencies.helixmtdb.het
                        > frequency.helixmtdb.max_het.expect("tested before")
                || frequency.helixmtdb.max_hom.is_some()
                    && s.population_frequencies.helixmtdb.hom
                        > frequency.helixmtdb.max_hom.expect("tested before"))
        {
            tracing::trace!(
                "variant {:?} fails HelixMtDb frequency filter {:?}",
                s,
                &frequency.helixmtdb
            );
            return Ok(false);
        }
        if frequency.gnomad_mtdna.enabled
            && (frequency.gnomad_mtdna.max_af.is_some()
                && s.population_frequencies.gnomad_mtdna.af()
                    > frequency.gnomad_mtdna.max_af.expect("tested before")
                || frequency.gnomad_mtdna.max_het.is_some()
                    && s.population_frequencies.gnomad_mtdna.het
                        > frequency.gnomad_mtdna.max_het.expect("tested before")
                || frequency.gnomad_mtdna.max_hom.is_some()
                    && s.population_frequencies.gnomad_mtdna.hom
                        > frequency.gnomad_mtdna.max_hom.expect("tested before"))
        {
            tracing::trace!(
                "variant {:?} fails gnomAD-MT frequency filter {:?}",
                s,
                &frequency.gnomad_mtdna
            );
            return Ok(false);
        }
    } else {
        if frequency.gnomad_exomes.enabled
            && (frequency.gnomad_exomes.max_af.is_some()
                && s.population_frequencies.gnomad_exomes.af()
                    > frequency.gnomad_exomes.max_af.expect("tested before")
                || frequency.gnomad_exomes.max_het.is_some()
                    && s.population_frequencies.gnomad_exomes.het
                        > frequency.gnomad_exomes.max_het.expect("tested before")
                || frequency.gnomad_exomes.max_hom.is_some()
                    && s.population_frequencies.gnomad_exomes.hom
                        > frequency.gnomad_exomes.max_hom.expect("tested before")
                || frequency.gnomad_exomes.max_hemi.is_some()
                    && s.population_frequencies.gnomad_exomes.hemi
                        > frequency.gnomad_exomes.max_hemi.expect("tested before"))
        {
            tracing::trace!(
                "variant {:?} fails gnomAD-exomes frequency filter {:?}",
                s,
                &frequency.gnomad_exomes
            );
            return Ok(false);
        }
        if frequency.gnomad_genomes.enabled
            && (frequency.gnomad_genomes.max_af.is_some()
                && s.population_frequencies.gnomad_genomes.af()
                    > frequency.gnomad_genomes.max_af.expect("tested before")
                || frequency.gnomad_genomes.max_het.is_some()
                    && s.population_frequencies.gnomad_genomes.het
                        > frequency.gnomad_genomes.max_het.expect("tested before")
                || frequency.gnomad_genomes.max_hom.is_some()
                    && s.population_frequencies.gnomad_genomes.hom
                        > frequency.gnomad_genomes.max_hom.expect("tested before")
                || frequency.gnomad_genomes.max_hemi.is_some()
                    && s.population_frequencies.gnomad_genomes.hemi
                        > frequency.gnomad_genomes.max_hemi.expect("tested before"))
        {
            tracing::trace!(
                "variant {:?} fails gnomAD-genomes frequency filter {:?}",
                s,
                &frequency.gnomad_genomes
            );
            return Ok(false);
        }
    }

    Ok(true)
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
mod test {
    use mehari::annotate::seqvars::ann::{AnnField, Consequence};
    use rstest::rstest;

    use crate::seqvars::query::schema::{
        data::{
            MitochondrialFrequencies, NuclearFrequencies, PopulationFrequencies, VariantRecord,
            VcfVariant,
        },
        query::{
            CaseQuery, MitochondrialFrequencySettings, NuclearFrequencySettings,
            QuerySettingsFrequency,
        },
    };

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
        let query = CaseQuery {
            frequency: QuerySettingsFrequency {
                gnomad_exomes: NuclearFrequencySettings {
                    enabled: query_gnomad_exomes_enabled,
                    max_af: query_gnomad_exomes_frequency,
                    max_het: query_gnomad_exomes_heterozygous,
                    max_hom: query_gnomad_exomes_homozygous,
                    max_hemi: query_gnomad_exomes_hemizygous,
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = VariantRecord {
            population_frequencies: PopulationFrequencies {
                gnomad_exomes: NuclearFrequencies {
                    an: seqvar_gnomad_exomes_an,
                    het: seqvar_gnomad_exomes_het,
                    hom: seqvar_gnomad_exomes_hom,
                    hemi: seqvar_gnomad_exomes_hemi,
                },
                ..Default::default()
            },
            vcf_variant: VcfVariant {
                chrom: "X".to_string(),
                pos: 1,
                ref_allele: "G".into(),
                alt_allele: "A".into(),
            },
            ann_fields: vec![AnnField {
                allele: mehari::annotate::seqvars::ann::Allele::Alt {
                    alternative: "A".into(),
                },
                consequences: vec![Consequence::MissenseVariant],
                putative_impact: Consequence::MissenseVariant.impact(),
                gene_symbol: Default::default(),
                gene_id: Default::default(),
                feature_type: mehari::annotate::seqvars::ann::FeatureType::SoTerm {
                    term: mehari::annotate::seqvars::ann::SoFeature::Transcript,
                },
                strand: Default::default(),
                feature_id: Default::default(),
                feature_biotype: vec![mehari::annotate::seqvars::ann::FeatureBiotype::Coding],
                rank: Default::default(),
                hgvs_t: Default::default(),
                hgvs_p: Default::default(),
                tx_pos: Default::default(),
                cds_pos: Default::default(),
                protein_pos: Default::default(),
                distance: Default::default(),
                messages: Default::default(),
            }],
            ..Default::default()
        };

        assert_eq!(
            super::passes(&query, &seq_var)?,
            expected_pass_all,
            "query: {:#?}, seq_var: {:#?}",
            query,
            seq_var
        );

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
        let query = CaseQuery {
            frequency: QuerySettingsFrequency {
                gnomad_genomes: NuclearFrequencySettings {
                    enabled: query_gnomad_genomes_enabled,
                    max_af: query_gnomad_genomes_frequency,
                    max_het: query_gnomad_genomes_heterozygous,
                    max_hom: query_gnomad_genomes_homozygous,
                    max_hemi: query_gnomad_genomes_hemizygous,
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = VariantRecord {
            population_frequencies: PopulationFrequencies {
                gnomad_genomes: NuclearFrequencies {
                    an: seqvar_gnomad_genomes_an,
                    het: seqvar_gnomad_genomes_het,
                    hom: seqvar_gnomad_genomes_hom,
                    hemi: seqvar_gnomad_genomes_hemi,
                },
                ..Default::default()
            },
            vcf_variant: VcfVariant {
                chrom: "X".to_string(),
                pos: 1,
                ref_allele: "G".into(),
                alt_allele: "A".into(),
            },
            ann_fields: vec![AnnField {
                allele: mehari::annotate::seqvars::ann::Allele::Alt {
                    alternative: "A".into(),
                },
                consequences: vec![Consequence::MissenseVariant],
                putative_impact: Consequence::MissenseVariant.impact(),
                gene_symbol: Default::default(),
                gene_id: Default::default(),
                feature_type: mehari::annotate::seqvars::ann::FeatureType::SoTerm {
                    term: mehari::annotate::seqvars::ann::SoFeature::Transcript,
                },
                strand: Default::default(),
                feature_id: Default::default(),
                feature_biotype: vec![mehari::annotate::seqvars::ann::FeatureBiotype::Coding],
                rank: Default::default(),
                hgvs_t: Default::default(),
                hgvs_p: Default::default(),
                tx_pos: Default::default(),
                cds_pos: Default::default(),
                protein_pos: Default::default(),
                distance: Default::default(),
                messages: Default::default(),
            }],
            ..Default::default()
        };

        assert_eq!(
            super::passes(&query, &seq_var)?,
            expected_pass_all,
            "query: {:#?}, seq_var: {:#?}",
            query,
            seq_var
        );

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
    #[case(1000, 0, 3, true, Some(0.002), None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 3, false, Some(0.002), None, None, true)]
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
        #[case] query_helixmtdb_enabled: bool,
        #[case] query_helixmtdb_af: Option<f32>,
        #[case] query_helixmtdb_het: Option<i32>,
        #[case] query_helixmtdb_hom: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            frequency: QuerySettingsFrequency {
                helixmtdb: MitochondrialFrequencySettings {
                    enabled: query_helixmtdb_enabled,
                    max_af: query_helixmtdb_af,
                    max_het: query_helixmtdb_het,
                    max_hom: query_helixmtdb_hom,
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = VariantRecord {
            population_frequencies: PopulationFrequencies {
                helixmtdb: MitochondrialFrequencies {
                    an: seqvar_helix_an,
                    het: seqvar_helix_het,
                    hom: seqvar_helix_hom,
                },
                ..Default::default()
            },
            vcf_variant: VcfVariant {
                chrom: "MT".to_string(),
                pos: 1,
                ref_allele: "G".into(),
                alt_allele: "A".into(),
            },
            ann_fields: vec![AnnField {
                allele: mehari::annotate::seqvars::ann::Allele::Alt {
                    alternative: "A".into(),
                },
                consequences: vec![Consequence::MissenseVariant],
                putative_impact: Consequence::MissenseVariant.impact(),
                gene_symbol: Default::default(),
                gene_id: Default::default(),
                feature_type: mehari::annotate::seqvars::ann::FeatureType::SoTerm {
                    term: mehari::annotate::seqvars::ann::SoFeature::Transcript,
                },
                strand: Default::default(),
                feature_id: Default::default(),
                feature_biotype: vec![mehari::annotate::seqvars::ann::FeatureBiotype::Coding],
                rank: Default::default(),
                hgvs_t: Default::default(),
                hgvs_p: Default::default(),
                tx_pos: Default::default(),
                cds_pos: Default::default(),
                protein_pos: Default::default(),
                distance: Default::default(),
                messages: Default::default(),
            }],
            ..Default::default()
        };

        assert_eq!(
            super::passes(&query, &seq_var)?,
            expected_pass_all,
            "query: {:#?}, seq_var: {:#?}",
            query,
            seq_var
        );

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
    #[case(1000, 0, 3, true, Some(0.002), None, None, false)]
    // frequency: pass [hom count] (fail but filter is disabled)
    #[case(1000, 0, 3, false, Some(0.002), None, None, true)]
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
    #[allow(clippy::too_many_arguments)]
    fn passes_frequency_gnomad_mtdna(
        #[case] seqvar_gnomad_mtdna_an: i32,
        #[case] seqvar_gnomad_mtdna_het: i32,
        #[case] seqvar_gnomad_mtdna_hom: i32,
        #[case] query_gnomad_mtdna_enabled: bool,
        #[case] query_gnomad_mtdna_frequency: Option<f32>,
        #[case] query_gnomad_mtdna_heteroplasmic: Option<i32>,
        #[case] query_gnomad_mtdna_homoplasmic: Option<i32>,
        #[case] expected_pass_all: bool,
    ) -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            frequency: QuerySettingsFrequency {
                gnomad_mtdna: MitochondrialFrequencySettings {
                    enabled: query_gnomad_mtdna_enabled,
                    max_af: query_gnomad_mtdna_frequency,
                    max_het: query_gnomad_mtdna_heteroplasmic,
                    max_hom: query_gnomad_mtdna_homoplasmic,
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = VariantRecord {
            population_frequencies: PopulationFrequencies {
                gnomad_mtdna: MitochondrialFrequencies {
                    an: seqvar_gnomad_mtdna_an,
                    het: seqvar_gnomad_mtdna_het,
                    hom: seqvar_gnomad_mtdna_hom,
                },
                ..Default::default()
            },
            vcf_variant: VcfVariant {
                chrom: "MT".to_string(),
                pos: 1,
                ref_allele: "G".into(),
                alt_allele: "A".into(),
            },
            ann_fields: vec![AnnField {
                allele: mehari::annotate::seqvars::ann::Allele::Alt {
                    alternative: "A".into(),
                },
                consequences: vec![Consequence::MissenseVariant],
                putative_impact: Consequence::MissenseVariant.impact(),
                gene_symbol: Default::default(),
                gene_id: Default::default(),
                feature_type: mehari::annotate::seqvars::ann::FeatureType::SoTerm {
                    term: mehari::annotate::seqvars::ann::SoFeature::Transcript,
                },
                strand: Default::default(),
                feature_id: Default::default(),
                feature_biotype: vec![mehari::annotate::seqvars::ann::FeatureBiotype::Coding],
                rank: Default::default(),
                hgvs_t: Default::default(),
                hgvs_p: Default::default(),
                tx_pos: Default::default(),
                cds_pos: Default::default(),
                protein_pos: Default::default(),
                distance: Default::default(),
                messages: Default::default(),
            }],
            ..Default::default()
        };

        assert_eq!(
            super::passes(&query, &seq_var)?,
            expected_pass_all,
            "query: {:#?}, seq_var: {:#?}",
            query,
            seq_var
        );

        Ok(())
    }
}
