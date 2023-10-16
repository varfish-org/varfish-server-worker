use crate::seqvars::query::{
    annonars::AnnonarsDbs,
    schema::{CaseQuery, SequenceVariant},
};

use annonars::clinvar_minimal::pbs::ClinicalSignificance::{self, *};

/// Determine whether the `SequenceVariant` passes the clinvar filter.
pub fn passes(
    query: &CaseQuery,
    annonars_dbs: &AnnonarsDbs,
    seqvar: &SequenceVariant,
) -> Result<bool, anyhow::Error> {
    if !query.require_in_clinvar {
        return Ok(true);
    }

    let variant = annonars::common::spdi::Var::new(
        annonars::common::cli::canonicalize(&seqvar.chrom),
        seqvar.pos,
        seqvar.reference.clone(),
        seqvar.alternative.clone(),
    );
    let cf_data = annonars_dbs
        .clinvar_db
        .cf_handle("clinvar")
        .ok_or_else(|| anyhow::anyhow!("could not get clinvar column family"))?;

    if let Ok(record) = annonars::clinvar_minimal::cli::query::query_for_variant(
        &variant,
        &annonars_dbs.clinvar_meta,
        &annonars_dbs.clinvar_db,
        &cf_data,
    ) {
        let clinical_significance: ClinicalSignificance =
            record
                .clinical_significance
                .try_into()
                .map_err(|e| anyhow::anyhow!("could not convert clinical significance: {}", e))?;
        Ok(match clinical_significance {
            Benign => query.clinvar_include_benign,
            LikelyBenign => query.clinvar_include_likely_benign,
            UncertainSignificance => query.clinvar_include_uncertain_significance,
            LikelyPathogenic => query.clinvar_include_likely_pathogenic,
            Pathogenic => query.clinvar_include_pathogenic,
        })
    } else {
        // Because of the annonars API, we currently need to swallow any error,
        // as "not found" currently maps to an error.
        Ok(true)
    }
}
