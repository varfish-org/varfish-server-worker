use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the frequency filter.
pub fn passes(query: &CaseQuery, s: &SequenceVariant) -> Result<bool, anyhow::Error> {
    let pop = &query.population_frequency;
    let is_mtdna = annonars::common::cli::canonicalize(&s.chrom) == "MT";

    if is_mtdna {
        if pop.helixmtdb.enabled
            && (pop.helixmtdb.frequency.is_some()
                && s.helixmtdb_af() > pop.helixmtdb.frequency.expect("tested before")
                || pop.helixmtdb.heteroplasmic.is_some()
                    && s.population_frequencies.helixmtdb.het
                        > pop.helixmtdb.heteroplasmic.expect("tested before")
                || pop.helixmtdb.homoplasmic.is_some()
                    && s.population_frequencies.helixmtdb.hom
                        > pop.helixmtdb.homoplasmic.expect("tested before"))
        {
            tracing::trace!(
                "variant {:?} fails HelixMtDb frequency filter {:?}",
                s,
                &pop
            );
            return Ok(false);
        }
    } else if pop.gnomad_exomes.enabled
        && (pop.gnomad_exomes.allele_frequency.is_some()
            && s.gnomad_exomes_af() > pop.gnomad_exomes.allele_frequency.expect("tested before")
            || pop.gnomad_exomes.heterozygous.is_some()
                && s.population_frequencies.gnomad.exomes_het
                    > pop.gnomad_exomes.heterozygous.expect("tested before")
            || pop.gnomad_exomes.homozygous.is_some()
                && s.population_frequencies.gnomad.exomes_hom
                    > pop.gnomad_exomes.homozygous.expect("tested before")
            || pop.gnomad_exomes.hemizygous.is_some()
                && s.population_frequencies.gnomad.exomes_hemi
                    > pop.gnomad_exomes.hemizygous.expect("tested before"))
    {
        tracing::trace!(
            "variant {:?} fails gnomAD exomes frequency filter {:?}",
            s,
            &pop.gnomad_exomes.allele_frequency
        );
        return Ok(false);
    }

    if pop.gnomad_genomes.enabled
        && (pop.gnomad_genomes.allele_frequency.is_some()
            && s.gnomad_genomes_af() > pop.gnomad_genomes.allele_frequency.expect("tested before")
            || pop.gnomad_genomes.heterozygous.is_some()
                && s.population_frequencies.gnomad.genomes_het
                    > pop.gnomad_genomes.heterozygous.expect("tested before")
            || pop.gnomad_genomes.homozygous.is_some()
                && s.population_frequencies.gnomad.genomes_hom
                    > pop.gnomad_genomes.homozygous.expect("tested before")
            || !is_mtdna
                && pop.gnomad_genomes.hemizygous.is_some()
                && s.population_frequencies.gnomad.genomes_hemi
                    > pop.gnomad_genomes.hemizygous.expect("tested before"))
    {
        tracing::trace!(
            "variant {:?} fails gnomAD genomes allele_frequency filter {:?}",
            s,
            &pop.gnomad_genomes.allele_frequency
        );
        return Ok(false);
    }

    Ok(true)
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
mod test {
    use mehari::annotate::seqvars::ann::{AnnField, Consequence};
    use rstest::rstest;

    use crate::seqvars::query::schema::{
        CaseQuery, HelixMtDBs, PopulationFrequencies, SequenceVariant,
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
        use crate::seqvars::query::schema::{
            GnomadNuclearOptions, Gnomads, PopulationFrequencies, PopulationFrequencyOptions,
        };

        let query = CaseQuery {
            population_frequency: PopulationFrequencyOptions {
                gnomad_exomes: GnomadNuclearOptions {
                    enabled: query_gnomad_exomes_enabled,
                    allele_frequency: query_gnomad_exomes_frequency,
                    heterozygous: query_gnomad_exomes_heterozygous,
                    homozygous: query_gnomad_exomes_homozygous,
                    hemizygous: query_gnomad_exomes_hemizygous,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            population_frequencies: PopulationFrequencies {
                gnomad: Gnomads {
                    exomes_an: seqvar_gnomad_exomes_an,
                    exomes_het: seqvar_gnomad_exomes_het,
                    exomes_hom: seqvar_gnomad_exomes_hom,
                    exomes_hemi: seqvar_gnomad_exomes_hemi,
                    ..Default::default()
                },
                ..Default::default()
            },
            chrom: "X".to_string(),
            reference: "G".into(),
            alternative: "A".into(),
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

        assert_eq!(super::passes(&query, &seq_var)?, expected_pass_all);

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
        use crate::seqvars::query::schema::{
            GnomadNuclearOptions, Gnomads, PopulationFrequencies, PopulationFrequencyOptions,
        };

        let query = CaseQuery {
            population_frequency: PopulationFrequencyOptions {
                gnomad_genomes: GnomadNuclearOptions {
                    enabled: query_gnomad_genomes_enabled,
                    allele_frequency: query_gnomad_genomes_frequency,
                    heterozygous: query_gnomad_genomes_heterozygous,
                    homozygous: query_gnomad_genomes_homozygous,
                    hemizygous: query_gnomad_genomes_hemizygous,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            population_frequencies: PopulationFrequencies {
                gnomad: Gnomads {
                    genomes_an: seqvar_gnomad_genomes_an,
                    genomes_het: seqvar_gnomad_genomes_het,
                    genomes_hom: seqvar_gnomad_genomes_hom,
                    genomes_hemi: seqvar_gnomad_genomes_hemi,
                    ..Default::default()
                },
                ..Default::default()
            },
            chrom: "X".to_string(),
            reference: "G".into(),
            alternative: "A".into(),
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

        assert_eq!(super::passes(&query, &seq_var)?, expected_pass_all);

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
        use crate::seqvars::query::schema::{HelixMtDbOptions, PopulationFrequencyOptions};

        let query = CaseQuery {
            population_frequency: PopulationFrequencyOptions {
                helixmtdb: HelixMtDbOptions {
                    enabled: query_helix_enabled,
                    frequency: query_helix_frequency,
                    heteroplasmic: query_helix_heteroplasmic,
                    homoplasmic: query_helix_homoplasmic,
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            population_frequencies: PopulationFrequencies {
                helixmtdb: HelixMtDBs {
                    an: seqvar_helix_an,
                    het: seqvar_helix_het,
                    hom: seqvar_helix_hom,
                },
                ..Default::default()
            },
            chrom: "MT".to_string(),
            reference: "G".into(),
            alternative: "A".into(),
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

        assert_eq!(super::passes(&query, &seq_var)?, expected_pass_all);

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
    #[allow(clippy::too_many_arguments)]
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
        use crate::seqvars::query::schema::{
            GnomadNuclearOptions, Gnomads, PopulationFrequencies, PopulationFrequencyOptions,
        };

        let query = CaseQuery {
            population_frequency: PopulationFrequencyOptions {
                gnomad_genomes: GnomadNuclearOptions {
                    enabled: query_gnomad_genomes_enabled,
                    allele_frequency: query_gnomad_genomes_frequency,
                    heterozygous: query_gnomad_genomes_heteroplasmic,
                    homozygous: query_gnomad_genomes_homoplasmic,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let seq_var = SequenceVariant {
            population_frequencies: PopulationFrequencies {
                gnomad: Gnomads {
                    genomes_an: seqvar_gnomad_genomes_an,
                    genomes_het: seqvar_gnomad_genomes_het,
                    genomes_hom: seqvar_gnomad_genomes_hom,
                    ..Default::default()
                },
                ..Default::default()
            },
            chrom: "MT".to_string(),
            reference: "G".into(),
            alternative: "A".into(),
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

        assert_eq!(super::passes(&query, &seq_var)?, expected_pass_all);

        Ok(())
    }
}
