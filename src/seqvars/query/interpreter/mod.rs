//! Apply settings from a `strucvar::query::schema::CaseQuery` to `SequenceVariant` records.

use std::collections::HashSet;

mod consequences;
mod frequency;
mod genes_allowlist;
mod genotype;
mod quality;
mod regions_allowlist;

use super::schema::{CaseQuery, SequenceVariant};

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug, Default)]
pub struct QueryInterpreter {
    /// The case query settings.
    pub query: CaseQuery,
    /// Gene allowlist with HGNC IDs.
    pub hgvs_allowlist: Option<HashSet<String>>,
}

/// Result type for `QueryInterpreter::passes_genotype()`.
#[derive(Debug, Default)]
pub struct PassesResult {
    /// Whether genotype passes for all samples.
    pub pass_all: bool,
}

impl QueryInterpreter {
    /// Construct new `QueryInterpreter` with the given query settings.
    pub fn new(query: CaseQuery, hgvs_allowlist: Option<HashSet<String>>) -> Self {
        tracing::error!(
            "note well that we will need a second pass for compound heterozygous variants"
        );
        QueryInterpreter {
            query,
            hgvs_allowlist,
        }
    }

    /// Determine whether the annotated `SequenceVariant` passes all criteria.
    pub fn passes(&self, seqvar: &SequenceVariant) -> Result<PassesResult, anyhow::Error> {
        let pass_frequency = frequency::passes(&self.query, seqvar)?;
        let pass_consequences = consequences::passes(&self.query, seqvar)?;
        let res_quality = quality::passes(&self.query, seqvar)?;
        let pass_genotype = genotype::passes(
            &self.query,
            seqvar,
            &res_quality
                .no_call_samples
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        )?;
        let pass_genes_allowlist = genes_allowlist::passes(&self.query, seqvar);
        let pass_regions_allowlist = regions_allowlist::passes(&self.query, seqvar);
        let pass_clinvar = self.passes_clinvar(seqvar)?;
        let pass_all = pass_frequency
            && pass_consequences
            && res_quality.pass
            && pass_genotype
            && pass_genes_allowlist
            && pass_regions_allowlist
            && pass_clinvar;
        Ok(PassesResult { pass_all })
    }

    fn passes_clinvar(&self, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
        Ok(true)
    }
}
