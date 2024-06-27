use crate::seqvars::query::{
    annonars::Annotator,
    schema::{CaseQuery, SequenceVariant},
};

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

        let result = match description.to_lowercase().as_str() {
            "benign" => query.clinvar.include_benign,
            "benign/likely benign" => {
                query.clinvar.include_benign || query.clinvar.include_likely_benign
            }
            "likely benign" => query.clinvar.include_likely_benign,
            "pathogenic" => query.clinvar.include_pathogenic,
            "pathogenic/likely pathogenic" => {
                query.clinvar.include_pathogenic || query.clinvar.include_likely_pathogenic
            }
            "likely pathogenic" => query.clinvar.include_likely_pathogenic,
            "uncertain significance" => query.clinvar.include_uncertain_significance,
            "conflicting classifications of pathogenicity" => {
                query.clinvar.include_uncertain_significance
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
