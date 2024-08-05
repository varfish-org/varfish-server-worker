use crate::seqvars::query::{
    annonars::Annotator,
    schema::{
        data::VariantRecord,
        query::{CaseQuery, ClinvarGermlineAggregateDescription},
    },
};

/// Determine whether the `VariantRecord` passes the clinvar filter.
pub fn passes(
    query: &CaseQuery,
    annotator: &Annotator,
    seqvar: &VariantRecord,
) -> Result<bool, anyhow::Error> {
    if !query.clinvar.presence_required {
        return Ok(true);
    }

    if let Some(record) = annotator
        .query_clinvar_minimal(seqvar)
        .map_err(|e| anyhow::anyhow!("problem querying clinvar-minimal: {}", e))?
    {
        if record.records.is_empty() {
            tracing::error!(
                "variant {:?} found with empty list in ClinVar (should not happen",
                seqvar
            );
            return Ok(false);
        } else if record.records.len() > 1 {
            tracing::warn!(
                "variant {:?} found list with {} entries, using first",
                seqvar,
                record.records.len()
            );
        }
        let vcv_record = &record.records[0];

        let description = vcv_record
            .classifications
            .as_ref()
            .and_then(|c| c.germline_classification.as_ref())
            .and_then(|c| c.description.as_ref())
            .cloned()
            .unwrap_or_default();

        use ClinvarGermlineAggregateDescription::*;
        let result = match description.to_lowercase().as_str() {
            "benign" => query.clinvar.germline_descriptions.contains(&Benign),
            "benign/likely benign" => {
                query.clinvar.germline_descriptions.contains(&Benign)
                    || query.clinvar.germline_descriptions.contains(&LikelyBenign)
            }
            "likely benign" => query.clinvar.germline_descriptions.contains(&LikelyBenign),
            "pathogenic" => query.clinvar.germline_descriptions.contains(&Pathogenic),
            "pathogenic/likely pathogenic" => {
                query.clinvar.germline_descriptions.contains(&LikelyBenign)
                    || query.clinvar.germline_descriptions.contains(&Benign)
            }
            "likely pathogenic" => query
                .clinvar
                .germline_descriptions
                .contains(&LikelyPathogenic),
            "uncertain significance" => query
                .clinvar
                .germline_descriptions
                .contains(&UncertainSignificance),
            "conflicting classifications of pathogenicity" => {
                query.clinvar.allow_conflicting_interpretations
            }
            _ => {
                // We could also downtone this to debug.
                tracing::warn!(
                    "variant {:?} has unknown classification: {}",
                    seqvar,
                    description
                );
                false
            }
        };

        if !result {
            tracing::trace!(
                "variant {:?} present in ClinVar but fails clinvar filter from query {:?}",
                seqvar,
                query
            );
        }

        Ok(result)
    } else {
        tracing::trace!(
            "variant {:?} not present in ClinVar and thus fails filter query {:?}",
            seqvar,
            query
        );
        Ok(false)
    }
}
