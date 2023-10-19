use crate::seqvars::query::{
    annonars::Annotator,
    schema::{CaseQuery, SequenceVariant},
};

use annonars::clinvar_minimal::pbs::{self, ClinicalSignificance::*};

/// Determine whether the `SequenceVariant` passes the clinvar filter.
pub fn passes(
    query: &CaseQuery,
    annotator: &Annotator,
    seqvar: &SequenceVariant,
) -> Result<bool, anyhow::Error> {
    if !query.require_in_clinvar {
        return Ok(true);
    }

    if let Some(record) = annotator
        .query_clinvar_minimal(seqvar)
        .map_err(|e| anyhow::anyhow!("problem querying clinvar-minimal: {}", e))?
    {
        if let Some(assertion) = record.reference_assertions.first() {
            let clinical_significance: pbs::ClinicalSignificance = assertion
                .clinical_significance
                .try_into()
                .map_err(|e| anyhow::anyhow!("could not convert clinical significance: {}", e))?;
            let result = match clinical_significance {
                Benign => query.clinvar_include_benign,
                LikelyBenign => query.clinvar_include_likely_benign,
                UncertainSignificance => query.clinvar_include_uncertain_significance,
                LikelyPathogenic => query.clinvar_include_likely_pathogenic,
                Pathogenic => query.clinvar_include_pathogenic,
            };
            if !result {
                tracing::trace!(
                    "variant {:?} fails clinvar filter from query {:?}",
                    seqvar,
                    query
                );
            }
            Ok(result)
        } else {
            unreachable!("no reference clinvar assertion")
        }
    } else {
        // Because of the annonars API, we currently need to swallow any error,
        // as "not found" currently maps to an error.
        Ok(true)
    }
}
