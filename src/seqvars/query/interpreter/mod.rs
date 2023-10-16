//! Apply settings from a `strucvar::query::schema::CaseQuery` to `SequenceVariant` records.

use std::collections::HashSet;

mod clinvar;
mod consequences;
mod frequency;
mod genes_allowlist;
mod genotype;
mod quality;
mod regions_allowlist;

use super::{
    annonars::AnnonarsDbs,
    schema::{CaseQuery, SequenceVariant},
};

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug, Default)]
pub struct QueryInterpreter {
    /// The case query settings.
    pub query: CaseQuery,
    /// Gene allowlist with HGNC IDs.
    pub hgnc_allowlist: Option<HashSet<String>>,
}

/// Result type for `QueryInterpreter::passes_genotype()`.
#[derive(Debug, Default)]
pub struct PassesResult {
    /// Whether genotype passes for all samples.
    pub pass_all: bool,
}

impl QueryInterpreter {
    /// Construct new `QueryInterpreter` with the given query settings.
    pub fn new(query: CaseQuery, hgnc_allowlist: Option<HashSet<String>>) -> Self {
        tracing::error!(
            "note well that we will need a second pass for compound heterozygous variants"
        );
        QueryInterpreter {
            query,
            hgnc_allowlist,
        }
    }

    /// Determine whether the annotated `SequenceVariant` passes all criteria.
    pub fn passes(
        &self,
        seqvar: &SequenceVariant,
        annonars_dbs: &AnnonarsDbs,
    ) -> Result<PassesResult, anyhow::Error> {
        // Check the filters first that are cheap to compute.
        let pass_frequency = frequency::passes(&self.query, seqvar)?;
        let pass_consequences = consequences::passes(&self.query, seqvar)?;
        let res_quality = quality::passes(&self.query, seqvar)?;
        let pass_genes_allowlist = genes_allowlist::passes(&self.hgnc_allowlist, seqvar);
        let pass_regions_allowlist = regions_allowlist::passes(&self.query, seqvar);
        if !pass_frequency
            || !pass_consequences
            || !res_quality.pass
            || !pass_genes_allowlist
            || !pass_regions_allowlist
        {
            return Ok(PassesResult { pass_all: false });
        }
        // Now also check the genotype that needs the quality filter output as input.
        if !genotype::passes(
            &self.query,
            seqvar,
            &res_quality
                .no_call_samples
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        )? {
            return Ok(PassesResult { pass_all: false });
        }
        // If we passed until here, check the presence in ClinVar which needs a database lookup.
        Ok(PassesResult {
            pass_all: clinvar::passes(&self.query, annonars_dbs, seqvar)?,
        })
    }
}
