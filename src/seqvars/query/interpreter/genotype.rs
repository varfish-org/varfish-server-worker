use crate::seqvars::query::schema::{
    data::VariantRecord,
    query::{
        CaseQuery, GenotypeChoice, MatchesGenotypeStr as _, QuerySettingsGenotype, RecessiveMode,
    },
};

/// Determine whether the `VariantRecord` passes the genotype filter.
pub fn passes(
    query: &CaseQuery,
    seqvar: &VariantRecord,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    let result = if query.genotype.recessive_mode != RecessiveMode::Disabled {
        passes_recessive_modes(&query.genotype, seqvar, no_call_samples)?
    } else {
        passes_non_recessive_mode(&query.genotype, seqvar, no_call_samples)?
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
fn passes_recessive_modes(
    query_genotype: &QuerySettingsGenotype,
    seqvar: &VariantRecord,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    // Pick recessive index and parent names from query.
    let index_name = query_genotype.recessive_index().map_err(|e| {
        anyhow::anyhow!(
            "invalid recessive index in genotype filter {:?}: {}",
            &query_genotype,
            e
        )
    })?;
    let parent_names = query_genotype.recessive_parents().map_err(|e| {
        anyhow::anyhow!(
            "invalid recessive parents in genotype filter {:?}: {}",
            &query_genotype,
            e
        )
    })?;

    // Get genotype choice of index.
    let index_gt_choice = query_genotype
        .sample_genotypes
        .get(&index_name)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} not found in genotype filter {:?}",
                &index_name,
                &query_genotype
            )
        })?;
    // For recessive mode, we have to know the samples selected for index and parents.
    let index_gt_string = if no_call_samples.contains(&index_name.as_str()) {
        String::from(".")
    } else {
        let call_info = seqvar.call_infos.get(&index_name).ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} not found in call info for {:?}",
                &index_name,
                &seqvar
            )
        })?;
        call_info.genotype.as_ref().cloned().ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} has no genotype in call info for {:?}",
                &index_name,
                &seqvar
            )
        })?
    };
    let parent_names = query_genotype
        .sample_genotypes
        .iter()
        .flat_map(|(sample, choice)| {
            if choice.genotype == GenotypeChoice::RecessiveParent {
                Some(sample.clone())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let parent_gt_strings = parent_names
        .iter()
        .map(|parent_name| {
            if no_call_samples.contains(&parent_name.as_str()) {
                Ok(String::from(".")) // no-call
            } else {
                seqvar
                    .call_infos
                    .get(parent_name)
                    .cloned()
                    .ok_or_else(|| {
                        anyhow::anyhow!(
                            "parent sample {} not found in call info for {:?}",
                            &parent_name,
                            &seqvar
                        )
                    })?
                    .genotype
                    .ok_or_else(|| {
                        anyhow::anyhow!(
                            "parent sample {} has no genotype in call info for {:?}",
                            &parent_name,
                            &seqvar
                        )
                    })
            }
        })
        .collect::<Result<Vec<_>, _>>()?;

    let compound_recessive_ok_index = GenotypeChoice::Het.matches(&index_gt_string);
    let compound_recessive_parents_ref = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| GenotypeChoice::Ref.matches(parent_gt_string))
        .count();
    let compound_recessive_parents_het = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| GenotypeChoice::Het.matches(parent_gt_string))
        .count();
    let compound_recessive_parents_hom = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| GenotypeChoice::Hom.matches(parent_gt_string))
        .count();

    let compound_recessive_ok = match parent_names.len() {
        0 => compound_recessive_ok_index,
        1 => {
            compound_recessive_ok_index
                && compound_recessive_parents_ref + compound_recessive_parents_het == 1
                && compound_recessive_parents_hom == 0
        }
        2 => {
            compound_recessive_ok_index
                && compound_recessive_parents_ref == 1
                && compound_recessive_parents_het == 1
                && compound_recessive_parents_hom == 0
        }
        _ => anyhow::bail!("more than two recessive parents selected"),
    };

    let homozygous_recessive_ok_index = GenotypeChoice::Hom.matches(&index_gt_string);
    let homozygous_recessive_ok_parents = parent_gt_strings
        .iter()
        .map(|parent_gt_string| GenotypeChoice::Het.matches(parent_gt_string))
        .collect::<Vec<_>>();
    let homozygous_recessive_ok =
        homozygous_recessive_ok_index && homozygous_recessive_ok_parents.iter().all(|&ok| ok);

    match query_genotype.recessive_mode {
        RecessiveMode::CompoundHeterozygous => Ok(compound_recessive_ok),
        RecessiveMode::Homozygous => Ok(homozygous_recessive_ok),
        RecessiveMode::Any => Ok(compound_recessive_ok || homozygous_recessive_ok),
        _ => anyhow::bail!(
            "invalid recessive mode choice for recessive mode: {:?}",
            index_gt_choice
        ),
    }
}

/// Handle case if the mode is not "recessive".  Note that this actually includes the
/// homozygous recessive mode.
fn passes_non_recessive_mode(
    query_genotype: &QuerySettingsGenotype,
    seqvar: &VariantRecord,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    for (sample_name, genotype_choice) in query_genotype.sample_genotypes.iter() {
        let genotype = if no_call_samples.contains(&sample_name.as_str()) {
            "." // no-call
        } else if let Some(call_info) = seqvar.call_infos.get(sample_name) {
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

        if !genotype_choice.genotype.matches(genotype) {
            tracing::trace!(
                "variant {:?} fails genotype filter {:?} on sample {}",
                seqvar,
                &query_genotype,
                sample_name
            );
            return Ok(false);
        }
    }

    Ok(true) // all good up to here
}

#[cfg(test)]
mod test {
    use crate::seqvars::query::schema::data::{CallInfo, VariantRecord};
    use crate::seqvars::query::schema::query::{
        GenotypeChoice::{self, *},
        QuerySettingsGenotype, RecessiveMode, SampleGenotypeChoice,
    };

    static INDEX_NAME: &str = "sample";
    static FATHER_NAME: &str = "father";
    static MOTHER_NAME: &str = "mother";

    #[rstest::rstest]
    // any: passes
    #[case("0/0", Any, true)]
    #[case("0|0", Any, true)]
    #[case("0/1", Any, true)]
    #[case("0|1", Any, true)]
    #[case("0", Any, true)]
    #[case("1/0", Any, true)]
    #[case("1|0", Any, true)]
    #[case("1/1", Any, true)]
    #[case("1|1", Any, true)]
    #[case("1", Any, true)]
    #[case(".", Any, true)]
    #[case("./.", Any, true)]
    #[case(".|.", Any, true)]
    // ref: passes
    #[case("0/0", Ref, true)]
    #[case("0|0", Ref, true)]
    #[case("0", Ref, true)]
    // ref: passes NOT
    #[case("0/1", Ref, false)]
    #[case("0|1", Ref, false)]
    #[case("1/0", Ref, false)]
    #[case("1|0", Ref, false)]
    #[case("1/1", Ref, false)]
    #[case("1|1", Ref, false)]
    #[case("1", Ref, false)]
    #[case(".", Ref, false)]
    #[case("./.", Ref, false)]
    #[case(".|.", Ref, false)]
    // het: passes
    #[case("0/1", Het, true)]
    #[case("0|1", Het, true)]
    #[case("1/0", Het, true)]
    #[case("1|0", Het, true)]
    // het: passes NOT
    #[case("0/0", Het, false)]
    #[case("0|0", Het, false)]
    #[case("0", Het, false)]
    #[case("1/1", Het, false)]
    #[case("1|1", Het, false)]
    #[case("1", Het, false)]
    #[case(".", Het, false)]
    #[case("./.", Het, false)]
    #[case(".|.", Het, false)]
    // hom: passes
    #[case("1/1", Hom, true)]
    #[case("1|1", Hom, true)]
    #[case("1", Hom, true)]
    // hom: passes NOT
    #[case("0/0", Hom, false)]
    #[case("0|0", Hom, false)]
    #[case("0/1", Hom, false)]
    #[case("0|1", Hom, false)]
    #[case("0", Hom, false)]
    #[case("1/0", Hom, false)]
    #[case("1|0", Hom, false)]
    #[case(".", Hom, false)]
    #[case("./.", Hom, false)]
    #[case(".|.", Hom, false)]
    // non-hom: passes
    #[case("0/0", NonHom, true)]
    #[case("0|0", NonHom, true)]
    #[case("0/1", NonHom, true)]
    #[case("0|1", NonHom, true)]
    #[case("0", NonHom, true)]
    #[case("1/0", NonHom, true)]
    #[case("1|0", NonHom, true)]
    #[case(".", NonHom, true)]
    #[case("./.", NonHom, true)]
    #[case(".|.", NonHom, true)]
    // non-hom: passes NOT
    #[case("1/1", NonHom, false)]
    #[case("1|1", NonHom, false)]
    #[case("1", NonHom, false)]
    // variant: passes
    #[case("0/1", Variant, true)]
    #[case("0|1", Variant, true)]
    #[case("1/0", Variant, true)]
    #[case("1|0", Variant, true)]
    #[case("1/1", Variant, true)]
    #[case("1|1", Variant, true)]
    #[case("1", Variant, true)]
    // variant: passes NOT
    #[case("0/0", Variant, false)]
    #[case("0|0", Variant, false)]
    #[case("0", Variant, false)]
    #[case(".", Variant, false)]
    #[case("./.", Variant, false)]
    #[case(".|.", Variant, false)]

    fn passes_non_recessive_mode_singleton(
        #[case] sample_gt: &str,
        #[case] query_gt: GenotypeChoice,
        #[case] expected: bool,
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

        assert_eq!(
            super::passes_non_recessive_mode(&query_genotype, &seq_var, &[])?,
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
    #[case("0/0,0/0,0/0", Any, Any, Any, true)]
    #[case(".,.,.", Any, Any, Any, true)]
    // some combinations: passes
    #[case("0/1,0/1,0/0", Het, Het, Ref, true)]
    #[case("0/1,0|1,0/0", Het, Het, Ref, true)]
    #[case("0|1,0/1,0/0", Het, Het, Ref, true)]
    #[case("0/1,0/0,0/1", Any, Ref, Het, true)]
    #[case("0|1,0/0,0|1", Any, Ref, Het, true)]
    #[case("0/1,0/0,0|1", Any, Ref, Het, true)]
    #[case("0|1,0/0,0/1", Any, Ref, Het, true)]
    #[case("0/1,0/0,0/0", Variant, Ref, Ref, true)]
    #[case("0|1,0/0,0/0", Variant, Ref, Ref, true)]
    #[case("1/0,0/0,0/0", Variant, Ref, Ref, true)]
    #[case("1|0,0/0,0/0", Variant, Ref, Ref, true)]
    #[case("1,0,0/0", Variant, Ref, Ref, true)]
    // some combinations: passes NOT
    #[case("0/0,0/0,0/0", Het, Any, Any, false)]
    #[case("0/0,0/0,0/0", Any, Het, Any, false)]
    #[case("0/0,0/0,0/0", Any, Any, Het, false)]
    #[case("1,0,0/1", Variant, Ref, Ref, false)]
    #[case("1,0,.", Variant, Ref, Ref, false)]
    fn passes_non_recessive_mode_trio(
        #[case] sample_gts: &str,
        #[case] query_gt_index: GenotypeChoice,
        #[case] query_gt_father: GenotypeChoice,
        #[case] query_gt_mother: GenotypeChoice,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Any,
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
            super::passes_non_recessive_mode(&query_genotype, &seq_var, &[])?,
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
    // any: passes
    #[case("0/0", RecessiveIndex)]
    #[case("0/0", RecessiveParent)]
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

        assert!(super::passes_non_recessive_mode(&query_genotype, &seq_var, &[]).is_err(),);

        Ok(())
    }

    #[rstest::rstest]
    // // comphet-index: passes
    // #[case("0/1", ComphetIndex, true)]
    // #[case("0|1", ComphetIndex, true)]
    // #[case("1/0", ComphetIndex, true)]
    // #[case("1|0", ComphetIndex, true)]
    // // comphet-index: passes NOT
    // #[case("0/0", ComphetIndex, false)]
    // #[case("0|0", ComphetIndex, false)]
    // #[case("0", ComphetIndex, false)]
    // #[case("1/1", ComphetIndex, false)]
    // #[case("1|1", ComphetIndex, false)]
    // #[case("1", ComphetIndex, false)]
    // #[case(".", ComphetIndex, false)]
    // #[case("./.", ComphetIndex, false)]
    // #[case(".|.", ComphetIndex, false)]
    // recessive-index: passes
    #[case("0/1", RecessiveIndex, true)]
    #[case("0|1", RecessiveIndex, true)]
    #[case("1/0", RecessiveIndex, true)]
    #[case("1|0", RecessiveIndex, true)]
    #[case("1/1", RecessiveIndex, true)]
    #[case("1|1", RecessiveIndex, true)]
    #[case("1", RecessiveIndex, true)]
    // recessive-index: passes NOT
    #[case("0/0", RecessiveIndex, false)]
    #[case("0|0", RecessiveIndex, false)]
    #[case("0", RecessiveIndex, false)]
    #[case(".", RecessiveIndex, false)]
    #[case("./.", RecessiveIndex, false)]
    #[case(".|.", RecessiveIndex, false)]
    fn passes_recessive_modes_singleton(
        #[case] sample_gt: &str,
        #[case] query_gt: GenotypeChoice,
        #[case] expected: bool,
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

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, &seq_var, &[])?,
            expected,
            "sample_gt = {}, query_gt = {:?}, expected = {}",
            sample_gt,
            query_gt,
            expected
        );

        Ok(())
    }

    #[rstest::rstest]
    // recessive mode: passes
    #[case("0/1,0/0,0/0", RecessiveIndex, Any, Any, true)]
    #[case("0|1,0/0,0/0", RecessiveIndex, Any, Any, true)]
    #[case("1/0,0/0,0/0", RecessiveIndex, Any, Any, true)]
    #[case("1|0,0/0,0/0", RecessiveIndex, Any, Any, true)]
    #[case("1/0,0/1,0/0", RecessiveIndex, RecessiveParent, RecessiveParent, true)]
    #[case("1/0,0/0,0/1", RecessiveIndex, RecessiveParent, RecessiveParent, true)]
    #[case("1/1,0/1,0/1", RecessiveIndex, RecessiveParent, RecessiveParent, true)]
    #[case("1/1,0/0,0/1", RecessiveIndex, Any, RecessiveParent, true)]
    #[case("1/1,0/1,0/0", RecessiveIndex, RecessiveParent, Any, true)]
    #[case("1|1,0|1,1|0", RecessiveIndex, RecessiveParent, RecessiveParent, true)]
    // recessive mode: passes NOT
    #[case("1/1,1/1,0/0", RecessiveIndex, RecessiveParent, Any, false)]
    #[case("0/1,0/0,0/0", RecessiveIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,1/1,0/0", RecessiveIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,0/0,1/1", RecessiveIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,0/1,0/1", RecessiveIndex, RecessiveParent, RecessiveParent, false)]
    // // compound recessive mode: passes
    // #[case("0/1,0/0,0/0", ComphetIndex, Any, Any, true)]
    // #[case("0|1,0/0,0/0", ComphetIndex, Any, Any, true)]
    // #[case("1/0,0/0,0/0", ComphetIndex, Any, Any, true)]
    // #[case("1|0,0/0,0/0", ComphetIndex, Any, Any, true)]
    // #[case("0/1,1/1,0/0", ComphetIndex, Any, RecessiveParent, true)]
    // #[case("0/1,0/0,1/1", ComphetIndex, RecessiveParent, Any, true)]
    // #[case("1/0,0/1,0/0", ComphetIndex, RecessiveParent, RecessiveParent, true)]
    // #[case("1/0,0/0,0/1", ComphetIndex, RecessiveParent, RecessiveParent, true)]
    // // compound recessive mode: passes NOT
    // #[case("1/1,0/1,0/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    // #[case("1|1,0|1,1|0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    // #[case("0/1,0/0,0/0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    // #[case("0/1,1/1,0/0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    // #[case("0/1,0/0,1/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    // #[case("0/1,0/1,0/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    fn passes_recessive_modes_trio(
        #[case] sample_gts: &str,
        #[case] query_gt_index: GenotypeChoice,
        #[case] query_gt_father: GenotypeChoice,
        #[case] query_gt_mother: GenotypeChoice,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Any,
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
            super::passes_recessive_modes(&query_genotype, &seq_var, &[])?,
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
}
