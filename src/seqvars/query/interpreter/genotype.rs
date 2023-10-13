use crate::seqvars::query::schema::{CaseQuery, RecessiveMode, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the genotype filter.
pub fn passes(query: &CaseQuery, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
    let result = if let (Some(index_name), Some(mode)) =
        (query.recessive_index.as_ref(), query.recessive_mode)
    {
        passes_recessive_modes(query, mode, index_name, seqvar)?
    } else {
        passes_non_recessive_mode(query, seqvar)?
    };

    tracing::trace!(
        "variant {:?} has result {} for genotype filter {:?}",
        seqvar,
        result,
        &query.genotype
    );
    Ok(result)
}

/// Handle case of the mode being one of the recessive modes.
fn passes_recessive_modes(
    query: &CaseQuery,
    mode: RecessiveMode,
    index_name: &str,
    seqvar: &SequenceVariant,
) -> Result<bool, anyhow::Error> {
    // For recessive mode, we have to know the samples selected for index and parents.
    let index_gt_string = {
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
    let parent_names = query
        .genotype
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
            seqvar
                .call_info
                .get(parent_name)
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "parent sample {} not found in call info for {:?}",
                        &parent_name,
                        &seqvar
                    )
                })?
                .genotype
                .as_ref()
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "parent sample {} has no genotype in call info for {:?}",
                        &parent_name,
                        &seqvar
                    )
                })
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
    let compound_recessive_ok = compound_recessive_ok_index
        && compound_recessive_parents_ref <= 1
        && compound_recessive_parents_het <= 1
        && compound_recessive_parents_hom == 0;

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

    match mode {
        RecessiveMode::Recessive => Ok(compound_recessive_ok && homozygous_recessive_ok),
        RecessiveMode::CompoundRecessive => Ok(compound_recessive_ok),
    }
}

/// Handle case if the mode is not "recessive".  Note that this actually includes the
/// homozygous recessive mode.
fn passes_non_recessive_mode(
    query: &CaseQuery,
    seqvar: &SequenceVariant,
) -> Result<bool, anyhow::Error> {
    for (sample_name, genotype) in query.genotype.iter() {
        let genotype_choice = if let Some(genotype_choice) = genotype {
            genotype_choice
        } else {
            tracing::trace!("no genotype choice for sample {} (skip&pass)", sample_name);
            continue;
        };
        let genotype = if let Some(call_info) = seqvar.call_info.get(sample_name) {
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
                &query.genotype,
                sample_name
            );
            return Ok(false);
        }
    }

    Ok(true) // all good up to here
}
