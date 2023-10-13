//! Apply settings from a `strucvar::query::schema::CaseQuery` to `SequenceVariant` records.

use std::collections::HashSet;

use super::schema::{CaseQuery, SequenceVariant};

/// Hold data structures that support the interpretation of one `CaseQuery`
/// to multiple `StructuralVariant` records.
#[derive(Debug)]
pub struct QueryInterpreter {
    pub query: CaseQuery,
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
    pub fn passes(&self, _seqvar: &SequenceVariant) -> Result<PassesResult, anyhow::Error> {
        todo!()
        // Ok(PassesResult { pass_all: true })
    }
}
