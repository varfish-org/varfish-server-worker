use crate::seqvars::query::schema::{
    data::VariantRecord,
    query::{
        considered_no_call, CaseQuery, GenotypeChoice, MatchesGenotypeStr as _,
        QuerySettingsGenotype, RecessiveMode, RecessiveParents,
    },
};

/// Determine whether the `VariantRecord` passes the genotype filter.
pub fn passes(query: &CaseQuery, seqvar: &VariantRecord) -> Result<bool, anyhow::Error> {
    let result = if query.genotype.recessive_mode != RecessiveMode::Disabled {
        passes_recessive_modes(&query.genotype, seqvar)?
    } else {
        passes_non_recessive_mode(&query.genotype, seqvar)?
    };

    if !result {
        tracing::trace!(
            "variant {:?} fails for genotype filter {:?}",
            seqvar,
            &query.genotype
        );
    }
    Ok(result)
}

/// Handle case of the mode being one of the recessive modes.
///
/// This means
///
/// - fail on chrMT/chrY
/// - in case of chrX, require het./hom./hemi. in the index, het. in the mother and
///   hom. ref. in the father
/// - in case of autosomal chromosomes, require het. in index and exactly one parent
///   and hom. ref. in other parent OR require hom. in index and het. in both parents
///
/// In the future, we could also provide the sex of the index here and include cases
/// of X inactivation where mother is het., father is hom. ref. and index is het.
fn passes_recessive_modes(
    query_genotype: &QuerySettingsGenotype,
    seqvar: &VariantRecord,
) -> Result<bool, anyhow::Error> {
    // Is/must never be called with disabled recessive mode.
    assert_ne!(query_genotype.recessive_mode, RecessiveMode::Disabled);
    // Get normalized chromosome, short-circuit in case of chrMT/chrY
    // (recessive inheritance does not make sense here).
    let normalized_chrom = annonars::common::cli::canonicalize(seqvar.vcf_variant.chrom.as_str());
    if normalized_chrom == "MT" || normalized_chrom == "Y" {
        tracing::trace!(
            "variant {:?} fails for genotype filter {:?} (chrMT/chrY)",
            seqvar,
            query_genotype
        );
        return Ok(false);
    }

    // Extract genotypes of index and potentially mother/father.
    let (index_gt, father_gt, mother_gt) = extract_trio_genotypes(query_genotype, seqvar)?;

    // Branch into X-linked and autosomal recessive mode.
    Ok(if normalized_chrom == "X" {
        passes_recessive_mode_x_linked(index_gt, father_gt, mother_gt)
    } else {
        passes_recessive_mode_autosomal(
            index_gt,
            father_gt,
            mother_gt,
            query_genotype.recessive_mode,
        )
    })
}

/// Extract genotypes of index and potentially mother/father.
///
/// This function is used to extract the genotypes of the index and the parents
///
/// # Arguments
///
/// * `query_genotype` - The genotype filter to apply
/// * `seqvar` - The variant record to check
///
/// # Returns
///
/// A tuple containing the genotypes of the index, father and mother
///
/// # Errors
///
/// This function returns an error if the parents names could not be extracted
/// from the genotype query settings (more than one mother/father), if the
/// index sample name could not be determined, or the genotypes of the family
/// members could not be extracted from the variant record.
fn extract_trio_genotypes<'a>(
    query_genotype: &QuerySettingsGenotype,
    seqvar: &'a VariantRecord,
) -> Result<(&'a str, Option<&'a str>, Option<&'a str>), anyhow::Error> {
    let index = query_genotype.recessive_index().map_err(|e| {
        anyhow::anyhow!(
            "invalid recessive index in genotype filter {:?}: {}",
            &query_genotype,
            e
        )
    })?;
    let RecessiveParents { father, mother } = query_genotype.recessive_parents().map_err(|e| {
        anyhow::anyhow!(
            "invalid recessive parents in genotype filter {:?}: {}",
            &query_genotype,
            e
        )
    })?;
    let index_gt = seqvar
        .call_infos
        .get(&index)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} not found in call info for {:?}",
                &index,
                &seqvar
            )
        })?
        .genotype
        .as_ref()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} has no genotype in call info for {:?}",
                &index,
                &seqvar
            )
        })?
        .as_str();
    let father_gt = father
        .map(|father| {
            seqvar
                .call_infos
                .get(&father)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "father sample {} not found in call info for {:?}",
                        &father,
                        &seqvar
                    )
                })?
                .genotype
                .as_ref()
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "father sample {} has no genotype in call info for {:?}",
                        &father,
                        &seqvar
                    )
                })
        })
        .transpose()?
        .map(|s| s.as_str());
    let mother_gt = mother
        .map(|mother| {
            seqvar
                .call_infos
                .get(&mother)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "mother sample {} not found in call info for {:?}",
                        &mother,
                        &seqvar
                    )
                })?
                .genotype
                .as_ref()
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "mother sample {} has no genotype in call info for {:?}",
                        &mother,
                        &seqvar
                    )
                })
        })
        .transpose()?
        .map(|s| s.as_str());
    Ok((index_gt, father_gt, mother_gt))
}

/// Handle case of the mode being "recessive" on chromosome X.
fn passes_recessive_mode_x_linked(
    index_gt: &str,
    father_gt: Option<&str>,
    mother_gt: Option<&str>,
) -> bool {
    if GenotypeChoice::Hom
        .matches(index_gt)
        .expect("matches() cannot fail for Hom")
    {
        // Index is hemi. alt.
        //
        // If father/mother is missing, we assume them to have compatible genotypes.
        let father_ref = father_gt
            .map(|father_gt| {
                GenotypeChoice::Ref
                    .matches(father_gt)
                    .expect("matches() cannot fail for Ref")
            })
            .unwrap_or(true);
        let mother_het = mother_gt
            .map(|mother_gt| {
                GenotypeChoice::Het
                    .matches(mother_gt)
                    .expect("matches() cannot fail for Het")
            })
            .unwrap_or(true);
        father_ref && mother_het
    } else {
        false
    }
}

/// Handle case of the mode being "recessive" on autosomal chromosomes.
fn passes_recessive_mode_autosomal(
    index_gt: &str,
    father_gt: Option<&str>,
    mother_gt: Option<&str>,
    mode: RecessiveMode,
) -> bool {
    // Is/must never be called with disabled recessive mode.
    assert_ne!(mode, RecessiveMode::Disabled);

    if GenotypeChoice::Hom
        .matches(index_gt)
        .expect("matches() cannot fail for Hom")
    {
        // Index is hom. alt.
        //
        // Fail if in comp. het. mode.
        if mode == RecessiveMode::CompoundHeterozygous {
            return false;
        }
        assert!(matches!(
            mode,
            RecessiveMode::Any | RecessiveMode::Homozygous
        ));
        // Both parents must be het. if given.
        let father_matches_het = father_gt
            .map(|father_gt| {
                GenotypeChoice::Het
                    .matches(father_gt)
                    .expect("matches() cannot fail for Het")
            })
            .unwrap_or(true);
        let mother_matches_het = mother_gt
            .map(|mother_gt| {
                GenotypeChoice::Het
                    .matches(mother_gt)
                    .expect("matches() cannot fail for Het")
            })
            .unwrap_or(true);
        father_matches_het && mother_matches_het
    } else if GenotypeChoice::Het
        .matches(index_gt)
        .expect("matches() cannot fail for Het")
    {
        // Index is het.
        //
        // Fail if in hom. mode.
        if mode == RecessiveMode::Homozygous {
            return false;
        }
        assert!(matches!(
            mode,
            RecessiveMode::Any | RecessiveMode::CompoundHeterozygous
        ));
        // Exactly one parent must be het., the other hom. ref.  This means
        // when counting het./ref. alleles, each sum must be <= 1.
        let (hom_father, het_father, ref_father) = father_gt
            .map(|father_gt| {
                (
                    GenotypeChoice::Hom
                        .matches(father_gt)
                        .expect("matches() cannot fail for Hom") as i32,
                    GenotypeChoice::Het
                        .matches(father_gt)
                        .expect("matches() cannot fail for Het") as i32,
                    GenotypeChoice::Ref
                        .matches(father_gt)
                        .expect("matches() cannot fail for Ref") as i32,
                )
            })
            .unwrap_or((0, 0, 0));
        let (hom_mother, het_mother, ref_mother) = mother_gt
            .map(|mother_gt| {
                (
                    GenotypeChoice::Hom
                        .matches(mother_gt)
                        .expect("matches() cannot fail for Hom") as i32,
                    GenotypeChoice::Het
                        .matches(mother_gt)
                        .expect("matches() cannot fail for Het") as i32,
                    GenotypeChoice::Ref
                        .matches(mother_gt)
                        .expect("matches() cannot fail for Ref") as i32,
                )
            })
            .unwrap_or((0, 0, 0));
        (hom_father + hom_mother == 0)
            && (het_father + het_mother <= 1)
            && (ref_father + ref_mother <= 1)
    } else {
        // None of the inclusion criteria met.
        false
    }
}

/// Handle case if the mode is not "recessive".
///
/// Note that this includes the homozygous recessive mode.
///
/// # Arguments
///
/// * `query_genotype` - The genotype filter to apply
/// * `seqvar` - The variant record to check
///
/// # Returns
///
/// A boolean indicating whether the variant passes the genotype filter
///
/// # Errors
///
/// This function returns an error if the genotype choice from the settings
/// could not be matched with the genotype (e.g., is recessive indicator).
fn passes_non_recessive_mode(
    query_genotype: &QuerySettingsGenotype,
    seqvar: &VariantRecord,
) -> Result<bool, anyhow::Error> {
    for (sample_name, genotype_choice) in query_genotype.sample_genotypes.iter() {
        // Extract genotype from call info, skip if not present.
        let genotype = if let Some(call_info) = seqvar.call_infos.get(sample_name) {
            if let Some(genotype) = call_info.genotype.as_ref() {
                genotype
            } else {
                tracing::trace!("no GT for sample {} (skip&fail)", sample_name);
                return Ok(false);
            }
        } else {
            tracing::trace!("no call info for sample {} (skip&fail)", sample_name);
            return Ok(false);
        };

        if considered_no_call(genotype) {
            // Handle case of nocall genotype.
            if genotype_choice.include_no_call {
                continue;
            } else {
                tracing::trace!(
                    "variant {:?} fails genotype filter {:?} on sample {} (no call)",
                    seqvar,
                    &query_genotype,
                    sample_name
                );
                return Ok(false);
            }
        } else if !genotype_choice.genotype.matches(genotype)? {
            // Handle case of non-nocall genotype.
            tracing::trace!(
                "variant {:?} fails genotype filter {:?} on sample {}",
                seqvar,
                &query_genotype,
                sample_name
            );
            return Ok(false);
        } else {
            // All good, check next.
        }
    }

    Ok(true) // All good up to the end.
}

#[cfg(test)]
mod test {
    use crate::seqvars::query::schema::data::{CallInfo, VariantRecord, VcfVariant};
    use crate::seqvars::query::schema::query::{
        GenotypeChoice::{self, *},
        QuerySettingsGenotype, RecessiveMode, SampleGenotypeChoice,
    };

    static INDEX_NAME: &str = "sample";
    static FATHER_NAME: &str = "father";
    static MOTHER_NAME: &str = "mother";

    #[rstest::rstest]
    // any: passes
    #[case::any_pass_01("0/0", Any, false, true)]
    #[case::any_pass_02("0|0", Any, false, true)]
    #[case::any_pass_03("0/1", Any, false, true)]
    #[case::any_pass_04("0|1", Any, false, true)]
    #[case::any_pass_05("0", Any, false, true)]
    #[case::any_pass_06("1/0", Any, false, true)]
    #[case::any_pass_07("1|0", Any, false, true)]
    #[case::any_pass_08("1/1", Any, false, true)]
    #[case::any_pass_09("1|1", Any, false, true)]
    #[case::any_pass_10("1", Any, false, true)]
    #[case::any_pass_11(".", Any, true, true)]
    #[case::any_pass_12("./.", Any, true, true)]
    #[case::any_pass_13(".|.", Any, true, true)]
    // any: passes NOT
    #[case::any_nopass_01(".", Any, false, false)]
    #[case::any_nopass_02("./.", Any, false, false)]
    #[case::any_nopass_03(".|.", Any, false, false)]
    // ref: passes
    #[case::ref_pass_01("0/0", Ref, false, true)]
    #[case::ref_pass_02("0|0", Ref, false, true)]
    #[case::ref_pass_03("0", Ref, false, true)]
    #[case::ref_pass_04(".", Ref, true, true)]
    #[case::ref_pass_05("./.", Ref, true, true)]
    #[case::ref_pass_06(".|.", Ref, true, true)]
    // ref: passes NOT
    #[case::ref_nopass_01("0/1", Ref, false, false)]
    #[case::ref_nopass_02("0|1", Ref, false, false)]
    #[case::ref_nopass_03("1/0", Ref, false, false)]
    #[case::ref_nopass_04("1|0", Ref, false, false)]
    #[case::ref_nopass_05("1/1", Ref, false, false)]
    #[case::ref_nopass_06("1|1", Ref, false, false)]
    #[case::ref_nopass_07("1", Ref, false, false)]
    #[case::ref_nopass_08(".", Ref, false, false)]
    #[case::ref_nopass_09("./.", Ref, false, false)]
    #[case::ref_nopass_10(".|.", Ref, false, false)]
    // het: passes
    #[case::het_pass_01("0/1", Het, false, true)]
    #[case::het_pass_02("0|1", Het, false, true)]
    #[case::het_pass_03("1/0", Het, false, true)]
    #[case::het_pass_04("1|0", Het, false, true)]
    #[case::het_pass_05(".", Het, true, true)]
    #[case::het_pass_06("./.", Het, true, true)]
    #[case::het_pass_07(".|.", Het, true, true)]
    // het: passes NOT
    #[case::het_nopass_01("0/0", Het, false, false)]
    #[case::het_nopass_02("0|0", Het, false, false)]
    #[case::het_nopass_03("0", Het, false, false)]
    #[case::het_nopass_04("1/1", Het, false, false)]
    #[case::het_nopass_05("1|1", Het, false, false)]
    #[case::het_nopass_06("1", Het, false, false)]
    #[case::het_nopass_07(".", Het, false, false)]
    #[case::het_nopass_08("./.", Het, false, false)]
    #[case::het_nopass_09(".|.", Het, false, false)]
    // hom: passes
    #[case::hom_pass_01("1/1", Hom, false, true)]
    #[case::hom_pass_02("1|1", Hom, false, true)]
    #[case::hom_pass_03("1", Hom, false, true)]
    #[case::hom_pass_04(".", Hom, true, true)]
    #[case::hom_pass_05("./.", Hom, true, true)]
    #[case::hom_pass_06(".|.", Hom, true, true)]
    // hom: passes NOT
    #[case::hom_nopass_01("0/0", Hom, false, false)]
    #[case::hom_nopass_02("0|0", Hom, false, false)]
    #[case::hom_nopass_03("0/1", Hom, false, false)]
    #[case::hom_nopass_04("0|1", Hom, false, false)]
    #[case::hom_nopass_05("0", Hom, false, false)]
    #[case::hom_nopass_06("1/0", Hom, false, false)]
    #[case::hom_nopass_07("1|0", Hom, false, false)]
    #[case::hom_nopass_08(".", Hom, false, false)]
    #[case::hom_nopass_09("./.", Hom, false, false)]
    #[case::hom_nopass_10(".|.", Hom, false, false)]
    // non-hom: passes
    #[case::nonhom_pass_01("0/0", NonHom, false, true)]
    #[case::nonhom_pass_02("0|0", NonHom, false, true)]
    #[case::nonhom_pass_03("0/1", NonHom, false, true)]
    #[case::nonhom_pass_04("0|1", NonHom, false, true)]
    #[case::nonhom_pass_05("0", NonHom, false, true)]
    #[case::nonhom_pass_06("1/0", NonHom, false, true)]
    #[case::nonhom_pass_07("1|0", NonHom, false, true)]
    #[case::nonhom_pass_08(".", NonHom, true, true)]
    #[case::nonhom_pass_09("./.", NonHom, true, true)]
    #[case::nonhom_pass_10(".|.", NonHom, true, true)]
    // non-hom: passes NOT
    #[case::nonhom_nopass_01("1/1", NonHom, false, false)]
    #[case::nonhom_nopass_02("1|1", NonHom, false, false)]
    #[case::nonhom_nopass_03("1", NonHom, false, false)]
    #[case::nonhom_nopass_04(".", NonHom, false, false)]
    #[case::nonhom_nopass_05("./.", NonHom, false, false)]
    #[case::nonhom_nopass_06(".|.", NonHom, false, false)]
    // variant: passes
    #[case::variant_pass_01("0/1", Variant, false, true)]
    #[case::variant_pass_02("0|1", Variant, false, true)]
    #[case::variant_pass_03("1/0", Variant, false, true)]
    #[case::variant_pass_04("1|0", Variant, false, true)]
    #[case::variant_pass_05("1/1", Variant, false, true)]
    #[case::variant_pass_06("1|1", Variant, false, true)]
    #[case::variant_pass_07("1", Variant, false, true)]
    #[case::variant_pass_08(".", Variant, true, true)]
    #[case::variant_pass_09("./.", Variant, true, true)]
    #[case::variant_pass_10(".|.", Variant, true, true)]
    // variant: passes NOT
    #[case::variant_nopass_01("0/0", Variant, false, false)]
    #[case::variant_nopass_02("0|0", Variant, false, false)]
    #[case::variant_nopass_03("0", Variant, false, false)]
    #[case::variant_nopass_04(".", Variant, false, false)]
    #[case::variant_nopass_05("./.", Variant, false, false)]
    #[case::variant_nopass_06(".|.", Variant, false, false)]

    fn passes_non_recessive_mode_singleton(
        #[case] sample_gt: &str,
        #[case] query_gt: GenotypeChoice,
        #[case] include_no_call: bool,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Any,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: query_gt,
                    include_no_call,
                    ..Default::default()
                }
            },
        };

        let seq_var = VariantRecord {
            call_infos: indexmap::indexmap! {
                INDEX_NAME.into() =>
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes_non_recessive_mode(&query_genotype, &seq_var)?,
            expected,
            "sample_gt = {}, query_gt = {:?}, expected = {}",
            sample_gt,
            query_gt,
            expected
        );

        Ok(())
    }

    #[rstest::rstest]
    // any: passes
    #[case::any_pass_01("0/0,0/0,0/0", Any, false, Any, false, Any, false, true)]
    #[case::any_pass_01(".,.,.", Any, true, Any, true, Any, true, true)]
    // some combinations: passes
    #[case::some_pass_01("0/1,0/1,0/0", Het, false, Het, false, Ref, false, true)]
    #[case::some_pass_02("0/1,0|1,0/0", Het, false, Het, false, Ref, false, true)]
    #[case::some_pass_03("0|1,0/1,0/0", Het, false, Het, false, Ref, false, true)]
    #[case::some_pass_04("0/1,0/0,0/1", Any, false, Ref, false, Het, false, true)]
    #[case::some_pass_05("0|1,0/0,0|1", Any, false, Ref, false, Het, false, true)]
    #[case::some_pass_06("0/1,0/0,0|1", Any, false, Ref, false, Het, false, true)]
    #[case::some_pass_07("0|1,0/0,0/1", Any, false, Ref, false, Het, false, true)]
    #[case::some_pass_08("0/1,0/0,0/0", Variant, false, Ref, false, Ref, false, true)]
    #[case::some_pass_09("0|1,0/0,0/0", Variant, false, Ref, false, Ref, false, true)]
    #[case::some_pass_10("1/0,0/0,0/0", Variant, false, Ref, false, Ref, false, true)]
    #[case::some_pass_11("1|0,0/0,0/0", Variant, false, Ref, false, Ref, false, true)]
    #[case::some_pass_12("1,0,0/0", Variant, false, Ref, false, Ref, false, true)]
    // some combinations: passes NOT
    #[case::some_nopass_01("0/0,0/0,0/0", Het, false, Any, false, Any, false, false)]
    #[case::some_nopass_02("0/0,0/0,0/0", Any, false, Het, false, Any, false, false)]
    #[case::some_nopass_03("0/0,0/0,0/0", Any, false, Any, false, Het, false, false)]
    #[case::some_nopass_04("1,0,0/1", Variant, false, Ref, false, Ref, false, false)]
    #[case::some_nopass_05("1,0,.", Variant, false, Ref, false, Ref, false, false)]
    fn passes_non_recessive_mode_trio(
        #[case] sample_gts: &str,
        #[case] query_gt_index: GenotypeChoice,
        #[case] nocall_index: bool,
        #[case] query_gt_father: GenotypeChoice,
        #[case] nocall_father: bool,
        #[case] query_gt_mother: GenotypeChoice,
        #[case] nocall_mother: bool,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Any,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: query_gt_index,
                    include_no_call: nocall_index,
                    ..Default::default()
                },
                String::from(FATHER_NAME) => SampleGenotypeChoice {
                    sample: String::from(FATHER_NAME),
                    genotype: query_gt_father,
                    include_no_call: nocall_father,
                    ..Default::default()
                },
                String::from(MOTHER_NAME) => SampleGenotypeChoice {
                    sample: String::from(MOTHER_NAME),
                    genotype: query_gt_mother,
                    include_no_call: nocall_mother,
                    ..Default::default()
                },
            },
        };
        let sample_gts = sample_gts
            .split(',')
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let seq_var = VariantRecord {
            call_infos: indexmap::indexmap! {
                String::from(INDEX_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[0].clone()),
                    ..Default::default()
                },
                String::from(FATHER_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[1].clone()),
                    ..Default::default()
                },
                String::from(MOTHER_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[2].clone()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes_non_recessive_mode(&query_genotype, &seq_var)?,
            expected,
            "sample_gt = {:?}, query_gt_index = {:?}, query_gt_father = {:?}, \
            query_gt_mother = {:?}, expected = {}",
            sample_gts,
            query_gt_index,
            query_gt_father,
            query_gt_mother,
            expected
        );

        Ok(())
    }

    #[rstest::rstest]
    #[case::recessive_index("0/0", RecessiveIndex)]
    #[case::recessive_father("0/0", RecessiveFather)]
    #[case::recessive_mother("0/0", RecessiveMother)]
    fn passes_non_recessive_mode_fails_on_recessive_markers(
        #[case] sample_gt: &str,
        #[case] query_gt: GenotypeChoice,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Any,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: query_gt,
                    ..Default::default()
                }
            },
        };
        let seq_var = VariantRecord {
            call_infos: indexmap::indexmap! {
                INDEX_NAME.into() =>
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert!(super::passes_non_recessive_mode(&query_genotype, &seq_var).is_err(),);

        Ok(())
    }

    #[rstest::rstest]
    // any: passes
    #[case::any_pass_01("0/1", RecessiveMode::Any, true)]
    #[case::any_pass_02("0|1", RecessiveMode::Any, true)]
    #[case::any_pass_03("1/0", RecessiveMode::Any, true)]
    #[case::any_pass_04("1|0", RecessiveMode::Any, true)]
    #[case::any_pass_05("1/1", RecessiveMode::Any, true)]
    #[case::any_pass_06("1|1", RecessiveMode::Any, true)]
    #[case::any_pass_07("1", RecessiveMode::Any, true)]
    // any: passes NOT
    #[case::any_nopass_01("0/0", RecessiveMode::Any, false)]
    #[case::any_nopass_02("0|0", RecessiveMode::Any, false)]
    #[case::any_nopass_03("0", RecessiveMode::Any, false)]
    #[case::any_nopass_04(".", RecessiveMode::Any, false)]
    #[case::any_nopass_05("./.", RecessiveMode::Any, false)]
    #[case::any_nopass_06(".|.", RecessiveMode::Any, false)]
    // comp. het.: passes
    #[case::comphet_pass_01("0/1", RecessiveMode::CompoundHeterozygous, true)]
    #[case::comphet_pass_02("0|1", RecessiveMode::CompoundHeterozygous, true)]
    #[case::comphet_pass_03("1/0", RecessiveMode::CompoundHeterozygous, true)]
    #[case::comphet_pass_04("1|0", RecessiveMode::CompoundHeterozygous, true)]
    // comp. het.: passes NOT
    #[case::comphet_nopass_01("0/0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_02("0|0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_03("0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_04("1/1", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_05("1|1", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_06("1", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_07(".", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_08("./.", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_09(".|.", RecessiveMode::CompoundHeterozygous, false)]
    // hom.: passes
    #[case::hom_pass_04("1/1", RecessiveMode::Homozygous, true)]
    #[case::hom_pass_05("1|1", RecessiveMode::Homozygous, true)]
    #[case::hom_pass_06("1", RecessiveMode::Homozygous, true)]
    // hom.: passes NOT
    #[case::hom_nopass_01("0/0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_02("0|0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_03("0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_04("0/1", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_05("0|1", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_06("1/0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_07("1|0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_08(".", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_09("./.", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_10(".|.", RecessiveMode::Homozygous, false)]
    fn passes_recessive_modes_autosomes_singleton(
        #[case] sample_gt: &str,
        #[case] recessive_mode: RecessiveMode,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: GenotypeChoice::RecessiveIndex,
                    ..Default::default()
                }
            },
        };
        let seq_var = VariantRecord {
            call_infos: indexmap::indexmap! {
                INDEX_NAME.into() =>
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, &seq_var)?,
            expected,
            "sample_gt = {}, recessive_mode = {:?}, expected = {}",
            sample_gt,
            recessive_mode,
            expected
        );

        Ok(())
    }

    #[rstest::rstest]
    // any: passes
    #[case::any_pass_01("1/1", RecessiveMode::Any, true)]
    #[case::any_pass_02("1|1", RecessiveMode::Any, true)]
    #[case::any_pass_03("1", RecessiveMode::Any, true)]
    // any: passes NOT
    #[case::any_nopass_01("0/1", RecessiveMode::Any, false)]
    #[case::any_nopass_02("0|1", RecessiveMode::Any, false)]
    #[case::any_nopass_03("1/0", RecessiveMode::Any, false)]
    #[case::any_nopass_04("1|0", RecessiveMode::Any, false)]
    #[case::any_nopass_05("0/0", RecessiveMode::Any, false)]
    #[case::any_nopass_06("0|0", RecessiveMode::Any, false)]
    #[case::any_nopass_07("0", RecessiveMode::Any, false)]
    #[case::any_nopass_08(".", RecessiveMode::Any, false)]
    #[case::any_nopass_09("./.", RecessiveMode::Any, false)]
    #[case::any_nopass_10(".|.", RecessiveMode::Any, false)]
    // comp. het.: passes
    // recessive mode is ignored -- sic!
    #[case::comphet_pass_01("1/1", RecessiveMode::CompoundHeterozygous, true)]
    #[case::comphet_pass_02("1|1", RecessiveMode::CompoundHeterozygous, true)]
    #[case::comphet_pass_03("1", RecessiveMode::CompoundHeterozygous, true)]
    // comp. het.: passes NOT
    #[case::comphet_nopass_01("0/1", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_02("0|1", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_03("1/0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_04("1|0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_05("0/0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_06("0|0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_07("0", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_08(".", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_09("./.", RecessiveMode::CompoundHeterozygous, false)]
    #[case::comphet_nopass_10(".|.", RecessiveMode::CompoundHeterozygous, false)]
    // hom.: passes
    #[case::hom_pass_04("1/1", RecessiveMode::Homozygous, true)]
    #[case::hom_pass_05("1|1", RecessiveMode::Homozygous, true)]
    #[case::hom_pass_06("1", RecessiveMode::Homozygous, true)]
    // hom.: passes NOT
    #[case::hom_nopass_01("0/0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_02("0|0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_03("0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_04("0/1", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_05("0|1", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_06("1/0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_07("1|0", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_08(".", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_09("./.", RecessiveMode::Homozygous, false)]
    #[case::hom_nopass_10(".|.", RecessiveMode::Homozygous, false)]
    fn passes_recessive_modes_x_linked_singleton(
        #[case] sample_gt: &str,
        #[case] recessive_mode: RecessiveMode,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: GenotypeChoice::RecessiveIndex,
                    ..Default::default()
                }
            },
        };
        let seq_var = VariantRecord {
            vcf_variant: VcfVariant {
                chrom: "X".to_string(),
                ..Default::default()
            },
            call_infos: indexmap::indexmap! {
                INDEX_NAME.into() =>
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, &seq_var)?,
            expected,
            "sample_gt = {}, recessive_mode = {:?}, expected = {}",
            sample_gt,
            recessive_mode,
            expected
        );

        Ok(())
    }

    #[rstest::rstest]
    #[case::chry_any_fail("Y", "0/1", RecessiveMode::Any)]
    #[case::chry_any_fail("Y", "1/1", RecessiveMode::Homozygous)]
    #[case::chry_any_fail("Y", "0/1", RecessiveMode::CompoundHeterozygous)]
    #[case::chrmt_any_fail("MT", "0/1", RecessiveMode::Any)]
    #[case::chrmt_any_fail("MT", "1/1", RecessiveMode::Homozygous)]
    #[case::chrmt_any_fail("MT", "0/1", RecessiveMode::CompoundHeterozygous)]
    fn passes_recessive_modes_chry_chrmt_singleton(
        #[case] chrom: &str,
        #[case] sample_gt: &str,
        #[case] recessive_mode: RecessiveMode,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: GenotypeChoice::RecessiveIndex,
                    ..Default::default()
                }
            },
        };
        let seq_var = VariantRecord {
            vcf_variant: VcfVariant {
                chrom: chrom.to_string(),
                ..Default::default()
            },
            call_infos: indexmap::indexmap! {
                INDEX_NAME.into() =>
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert!(
            !(super::passes_recessive_modes(&query_genotype, &seq_var)?),
            "sample_gt = {}, recessive_mode = {:?}",
            sample_gt,
            recessive_mode,
        );

        Ok(())
    }

    #[rstest::rstest]
    // any recessive mode: passes
    #[case::any_pass_01("0/1,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::any_pass_02("0|1,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::any_pass_03("1/0,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::any_pass_04("1|0,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::any_pass_05(
        "1/0,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::any_pass_06(
        "1/0,0/0,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::any_pass_07(
        "1/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::any_pass_08(
        "1/1,0/0,0/1",
        RecessiveIndex,
        Any,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::any_pass_09(
        "1/1,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Any,
        true
    )]
    #[case::any_pass_10(
        "1|1,0|1,1|0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    // any recessive mode: passes NOT
    #[case::any_nopass_01(
        "1/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Any,
        false
    )]
    #[case::any_nopass_02(
        "0/1,0/0,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::any_nopass_03(
        "0/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::any_nopass_04(
        "0/1,0/0,1/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::any_nopass_05(
        "0/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    // homozygous recessive mode: passes
    #[case::homozygous_pass_01(
        "1/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        true
    )]
    #[case::homozygous_pass_02(
        "1/1,0/0,0/1",
        RecessiveIndex,
        Any,
        RecessiveMother,
        RecessiveMode::Homozygous,
        true
    )]
    #[case::homozygous_pass_03(
        "1/1,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Homozygous,
        true
    )]
    #[case::homozygous_pass_04(
        "1|1,0|1,1|0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        true
    )]
    // homozygous recessive mode: passes NOT
    #[case::homozygous_pass_01(
        "0/1,0/0,0/0",
        RecessiveIndex,
        Any,
        Any,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_02(
        "0|1,0/0,0/0",
        RecessiveIndex,
        Any,
        Any,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_03(
        "1/0,0/0,0/0",
        RecessiveIndex,
        Any,
        Any,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_04(
        "1|0,0/0,0/0",
        RecessiveIndex,
        Any,
        Any,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_05(
        "1/0,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_06(
        "1/0,0/0,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_07(
        "1/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_08(
        "0/1,0/0,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_09(
        "0/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_10(
        "0/1,0/0,1/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    #[case::homozygous_pass_11(
        "0/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Homozygous,
        false
    )]
    // comp. het. recessive mode: passes
    #[case::comphet_pass_01("0/1,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::comphet_pass_02("0|1,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::comphet_pass_03("1/0,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::comphet_pass_04("1|0,0/0,0/0", RecessiveIndex, Any, Any, RecessiveMode::Any, true)]
    #[case::comphet_pass_05(
        "1/0,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::comphet_pass_06(
        "1/0,0/0,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::comphet_pass_07(
        "1/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::comphet_pass_08(
        "1/1,0/0,0/1",
        RecessiveIndex,
        Any,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    #[case::comphet_pass_09(
        "1/1,0/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Any,
        true
    )]
    #[case::comphet_pass_10(
        "1|1,0|1,1|0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        true
    )]
    // comp. het. recessive mode: passes NOT
    #[case::comphet_nopass_01(
        "1/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        Any,
        RecessiveMode::Any,
        false
    )]
    #[case::comphet_nopass_02(
        "0/1,0/0,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::comphet_nopass_03(
        "0/1,1/1,0/0",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::comphet_nopass_04(
        "0/1,0/0,1/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    #[case::comphet_nopass_05(
        "0/1,0/1,0/1",
        RecessiveIndex,
        RecessiveFather,
        RecessiveMother,
        RecessiveMode::Any,
        false
    )]
    fn passes_recessive_modes_autosomes_trio(
        #[case] sample_gts: &str,
        #[case] query_gt_index: GenotypeChoice,
        #[case] query_gt_father: GenotypeChoice,
        #[case] query_gt_mother: GenotypeChoice,
        #[case] recessive_mode: RecessiveMode,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode,
            sample_genotypes: indexmap::indexmap! {
                String::from(INDEX_NAME) => SampleGenotypeChoice {
                    sample: String::from(INDEX_NAME),
                    genotype: query_gt_index,
                    ..Default::default()
                },
                String::from(FATHER_NAME) => SampleGenotypeChoice {
                    sample: String::from(FATHER_NAME),
                    genotype: query_gt_father,
                    ..Default::default()
                },
                String::from(MOTHER_NAME) => SampleGenotypeChoice {
                    sample: String::from(MOTHER_NAME),
                    genotype: query_gt_mother,
                    ..Default::default()
                },
            },
        };
        let sample_gts = sample_gts
            .split(',')
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let seq_var = VariantRecord {
            call_infos: indexmap::indexmap! {
                String::from(INDEX_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[0].clone()),
                    ..Default::default()
                },
                String::from(FATHER_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[1].clone()),
                    ..Default::default()
                },
                String::from(MOTHER_NAME) =>
                CallInfo {
                    genotype: Some(sample_gts[2].clone()),
                    ..Default::default()
                },
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, &seq_var)?,
            expected,
            "sample_gt = {:?}, query_gt_index = {:?}, query_gt_father = {:?}, \
            query_gt_mother = {:?}, recessive_mode = {:?}, expected = {}",
            sample_gts,
            query_gt_index,
            query_gt_father,
            query_gt_mother,
            recessive_mode,
            expected
        );

        Ok(())
    }
}
