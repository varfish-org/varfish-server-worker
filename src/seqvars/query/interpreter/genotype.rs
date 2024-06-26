use crate::seqvars::query::schema::{CaseQuery, GenotypeChoice, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the genotype filter.
pub fn passes(
    query: &CaseQuery,
    seqvar: &SequenceVariant,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    let result = if query.recessive_mode() {
        let index_sample = query.index_sample().ok_or_else(|| {
            anyhow::anyhow!("recessive mode requires an index sample, but none was found")
        })?;
        passes_recessive_modes(&query.genotype, &index_sample, seqvar, no_call_samples)?
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
    query_genotype: &indexmap::IndexMap<String, Option<GenotypeChoice>>,
    index_name: &str,
    seqvar: &SequenceVariant,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    // Get genotype choice of index.
    let index_gt_choice = query_genotype
        .get(index_name)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} not found in genotype filter {:?}",
                &index_name,
                &query_genotype
            )
        })?
        .ok_or_else(|| {
            anyhow::anyhow!(
                "index sample {} has no genotype choice in genotype filter {:?}",
                &index_name,
                &query_genotype
            )
        })?;
    // For recessive mode, we have to know the samples selected for index and parents.
    let index_gt_string = if no_call_samples.contains(&index_name) {
        String::from(".")
    } else {
        let call_info = seqvar.call_info.get(index_name).ok_or_else(|| {
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
        .iter()
        .flat_map(|(sample, choice)| {
            if matches!(
                choice,
                Some(crate::seqvars::query::schema::GenotypeChoice::RecessiveParent)
            ) {
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
                    .call_info
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

    let compound_recessive_ok_index = crate::seqvars::query::schema::GenotypeChoice::Het
        .matches(&index_gt_string)
        .map_err(|e| anyhow::anyhow!("invalid index genotype: {}", e))?;
    let compound_recessive_parents_ref = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| {
            crate::seqvars::query::schema::GenotypeChoice::Ref
                .matches(parent_gt_string)
                .unwrap_or(false)
        })
        .count();
    let compound_recessive_parents_het = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| {
            crate::seqvars::query::schema::GenotypeChoice::Het
                .matches(parent_gt_string)
                .unwrap_or(false)
        })
        .count();
    let compound_recessive_parents_hom = parent_gt_strings
        .iter()
        .filter(|parent_gt_string| {
            crate::seqvars::query::schema::GenotypeChoice::Hom
                .matches(parent_gt_string)
                .unwrap_or(false)
        })
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

    let homozygous_recessive_ok_index = crate::seqvars::query::schema::GenotypeChoice::Hom
        .matches(&index_gt_string)
        .map_err(|e| anyhow::anyhow!("invalid index genotype: {}", e))?;
    let homozygous_recessive_ok_parents = parent_gt_strings
        .iter()
        .map(|parent_gt_string| {
            crate::seqvars::query::schema::GenotypeChoice::Het
                .matches(parent_gt_string)
                .map_err(|e| anyhow::anyhow!("invalid parent genotype: {}", e))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let homozygous_recessive_ok =
        homozygous_recessive_ok_index && homozygous_recessive_ok_parents.iter().all(|&ok| ok);

    match index_gt_choice {
        GenotypeChoice::ComphetIndex => Ok(compound_recessive_ok),
        GenotypeChoice::RecessiveIndex => Ok(compound_recessive_ok || homozygous_recessive_ok),
        _ => anyhow::bail!(
            "invalid genotype choice for recessive mode: {:?}",
            index_gt_choice
        ),
    }
}

/// Handle case if the mode is not "recessive".  Note that this actually includes the
/// homozygous recessive mode.
fn passes_non_recessive_mode(
    query_genotype: &indexmap::IndexMap<String, Option<GenotypeChoice>>,
    seqvar: &SequenceVariant,
    no_call_samples: &[&str],
) -> Result<bool, anyhow::Error> {
    for (sample_name, genotype) in query_genotype.iter() {
        let genotype_choice = if let Some(genotype_choice) = genotype {
            genotype_choice
        } else {
            tracing::trace!("no genotype choice for sample {} (skip&pass)", sample_name);
            continue;
        };
        let genotype = if no_call_samples.contains(&sample_name.as_str()) {
            "." // no-call
        } else if let Some(call_info) = seqvar.call_info.get(sample_name) {
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

        if !genotype_choice
            .matches(genotype)
            .map_err(|e| anyhow::anyhow!("invalid genotype choice in {:?}: {}", &seqvar, e))?
        {
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
    use rstest::rstest;

    use crate::seqvars::query::schema::{
        CallInfo,
        GenotypeChoice::{self, *},
        SequenceVariant,
    };

    static INDEX_NAME: &str = "sample";
    static FATHER_NAME: &str = "father";
    static MOTHER_NAME: &str = "mother";

    #[rstest]
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
        let query_genotype: indexmap::IndexMap<_, _> =
            vec![(String::from(INDEX_NAME), Some(query_gt))]
                .into_iter()
                .collect();
        let seq_var = SequenceVariant {
            call_info: vec![(
                INDEX_NAME.into(),
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            )]
            .into_iter()
            .collect(),
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

    #[rstest]
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
        let query_genotype: indexmap::IndexMap<_, _> = vec![
            (String::from(INDEX_NAME), Some(query_gt_index)),
            (String::from(FATHER_NAME), Some(query_gt_father)),
            (String::from(MOTHER_NAME), Some(query_gt_mother)),
        ]
        .into_iter()
        .collect();
        let sample_gts = sample_gts
            .split(',')
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let seq_var = SequenceVariant {
            call_info: vec![
                (
                    String::from(INDEX_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[0].clone()),
                        ..Default::default()
                    },
                ),
                (
                    String::from(FATHER_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[1].clone()),
                        ..Default::default()
                    },
                ),
                (
                    String::from(MOTHER_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[2].clone()),
                        ..Default::default()
                    },
                ),
            ]
            .into_iter()
            .collect(),
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

    #[rstest]
    // any: passes
    #[case("0/0", ComphetIndex)]
    #[case("0/0", RecessiveIndex)]
    #[case("0/0", RecessiveParent)]
    fn passes_non_recessive_mode_fails_on_recessive_markers(
        #[case] sample_gt: &str,
        #[case] query_gt: GenotypeChoice,
    ) -> Result<(), anyhow::Error> {
        let query_genotype: indexmap::IndexMap<_, _> =
            vec![(String::from(INDEX_NAME), Some(query_gt))]
                .into_iter()
                .collect();
        let seq_var = SequenceVariant {
            call_info: vec![(
                INDEX_NAME.into(),
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            )]
            .into_iter()
            .collect(),
            ..Default::default()
        };

        assert!(super::passes_non_recessive_mode(&query_genotype, &seq_var, &[]).is_err(),);

        Ok(())
    }

    #[rstest]
    // comphet-index: passes
    #[case("0/1", ComphetIndex, true)]
    #[case("0|1", ComphetIndex, true)]
    #[case("1/0", ComphetIndex, true)]
    #[case("1|0", ComphetIndex, true)]
    // comphet-index: passes NOT
    #[case("0/0", ComphetIndex, false)]
    #[case("0|0", ComphetIndex, false)]
    #[case("0", ComphetIndex, false)]
    #[case("1/1", ComphetIndex, false)]
    #[case("1|1", ComphetIndex, false)]
    #[case("1", ComphetIndex, false)]
    #[case(".", ComphetIndex, false)]
    #[case("./.", ComphetIndex, false)]
    #[case(".|.", ComphetIndex, false)]
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
        let query_genotype: indexmap::IndexMap<_, _> =
            vec![(String::from(INDEX_NAME), Some(query_gt))]
                .into_iter()
                .collect();
        let seq_var = SequenceVariant {
            call_info: vec![(
                INDEX_NAME.into(),
                CallInfo {
                    genotype: Some(sample_gt.into()),
                    ..Default::default()
                },
            )]
            .into_iter()
            .collect(),
            ..Default::default()
        };

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, INDEX_NAME, &seq_var, &[])?,
            expected,
            "sample_gt = {}, query_gt = {:?}, expected = {}",
            sample_gt,
            query_gt,
            expected
        );

        Ok(())
    }

    #[rstest]
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
    // compound recessive mode: passes
    #[case("0/1,0/0,0/0", ComphetIndex, Any, Any, true)]
    #[case("0|1,0/0,0/0", ComphetIndex, Any, Any, true)]
    #[case("1/0,0/0,0/0", ComphetIndex, Any, Any, true)]
    #[case("1|0,0/0,0/0", ComphetIndex, Any, Any, true)]
    #[case("0/1,1/1,0/0", ComphetIndex, Any, RecessiveParent, true)]
    #[case("0/1,0/0,1/1", ComphetIndex, RecessiveParent, Any, true)]
    #[case("1/0,0/1,0/0", ComphetIndex, RecessiveParent, RecessiveParent, true)]
    #[case("1/0,0/0,0/1", ComphetIndex, RecessiveParent, RecessiveParent, true)]
    // compound recessive mode: passes NOT
    #[case("1/1,0/1,0/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    #[case("1|1,0|1,1|0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,0/0,0/0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,1/1,0/0", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,0/0,1/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    #[case("0/1,0/1,0/1", ComphetIndex, RecessiveParent, RecessiveParent, false)]
    fn passes_recessive_modes_trio(
        #[case] sample_gts: &str,
        #[case] query_gt_index: GenotypeChoice,
        #[case] query_gt_father: GenotypeChoice,
        #[case] query_gt_mother: GenotypeChoice,
        #[case] expected: bool,
    ) -> Result<(), anyhow::Error> {
        let query_genotype: indexmap::IndexMap<_, _> = vec![
            (String::from(INDEX_NAME), Some(query_gt_index)),
            (String::from(FATHER_NAME), Some(query_gt_father)),
            (String::from(MOTHER_NAME), Some(query_gt_mother)),
        ]
        .into_iter()
        .collect();
        let sample_gts = sample_gts
            .split(',')
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let seq_var = SequenceVariant {
            call_info: vec![
                (
                    String::from(INDEX_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[0].clone()),
                        ..Default::default()
                    },
                ),
                (
                    String::from(FATHER_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[1].clone()),
                        ..Default::default()
                    },
                ),
                (
                    String::from(MOTHER_NAME),
                    CallInfo {
                        genotype: Some(sample_gts[2].clone()),
                        ..Default::default()
                    },
                ),
            ]
            .into_iter()
            .collect(),
            ..Default::default()
        };

        assert_eq!(
            super::passes_recessive_modes(&query_genotype, INDEX_NAME, &seq_var, &[])?,
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
