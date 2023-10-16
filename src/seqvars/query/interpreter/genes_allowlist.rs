use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the genes allowlist filter.
pub fn passes(query: &CaseQuery, seqvar: &SequenceVariant) -> bool {
    if let Some(gene_allowlist) = &query.gene_allowlist {
        if gene_allowlist.is_empty() {
            true
        } else {
            seqvar
                .ann_fields
                .iter()
                .any(|ann_field| gene_allowlist.contains(&ann_field.gene_id))
        }
    } else {
        true
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};
    use mehari::annotate::seqvars::ann::AnnField;

    #[rstest]
    #[case(None, None, true)]
    #[case(
        Some(vec![]),
        None,
        true,
    )]
    #[case(
        Some(vec![String::from("HGNC:1100")]),
        None,
        false,
    )]
    #[case(
        Some(vec![String::from("HGNC:1")]),
        Some(String::from("HGNC:1100")),
        false,
    )]
    #[case(
        Some(vec![String::from("HGNC:1100")]),
        Some(String::from("HGNC:1100")),
        true,
    )]
    #[case(
        Some(vec![String::from("HGNC:1"), String::from("HGNC:1100")]),
        Some(String::from("HGNC:1100")),
        true,
    )]
    fn passes(
        #[case] query_genes: Option<Vec<String>>,
        #[case] seqvar_gene: Option<String>,
        #[case] expected: bool,
    ) {
        let query = CaseQuery {
            gene_allowlist: query_genes.clone(),
            ..Default::default()
        };
        let seqvar = SequenceVariant {
            ann_fields: if let Some(gene_id) = seqvar_gene.as_ref() {
                vec![AnnField {
                    gene_id: gene_id.clone(),
                    ..Default::default()
                }]
            } else {
                vec![]
            },
            ..Default::default()
        };

        assert_eq!(
            super::passes(&query, &seqvar),
            expected,
            "query_genes: {:?}, seqvar_genes: {:?}",
            query_genes,
            seqvar_gene
        )
    }
}
