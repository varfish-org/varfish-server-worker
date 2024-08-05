use std::collections::HashSet;

use crate::seqvars::query::schema::data::VariantRecord;

/// Determine whether the `VariantRecord` passes the genes allowlist filter.
pub fn passes(hgnc_allowlist: &HashSet<String>, seqvar: &VariantRecord) -> bool {
    if hgnc_allowlist.is_empty() {
        true
    } else {
        let res = seqvar
            .ann_fields
            .iter()
            .any(|ann_field| hgnc_allowlist.contains(&ann_field.gene_id));
        if !res {
            tracing::trace!(
                "variant {:?} fails gene allowlist filter {:?}",
                seqvar,
                &hgnc_allowlist
            );
        }
        res
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::seqvars::query::schema::data::VariantRecord;
    use mehari::annotate::seqvars::ann::AnnField;

    #[rstest]
    #[case(vec![], None, true)]
    #[case(
        vec![],
        None,
        true,
    )]
    #[case(
        vec![String::from("HGNC:1100")],
        None,
        false,
    )]
    #[case(
        vec![String::from("HGNC:1")],
        Some(String::from("HGNC:1100")),
        false,
    )]
    #[case(
        vec![String::from("HGNC:1100")],
        Some(String::from("HGNC:1100")),
        true,
    )]
    #[case(
        vec![String::from("HGNC:1"), String::from("HGNC:1100")],
        Some(String::from("HGNC:1100")),
        true,
    )]
    fn passes(
        #[case] hgnc_allowlist: Vec<String>,
        #[case] seqvar_gene: Option<String>,
        #[case] expected: bool,
    ) {
        let hgnc_allowlist = hgnc_allowlist
            .into_iter()
            .map(|hgnc_id| hgnc_id.to_uppercase())
            .collect::<std::collections::HashSet<_>>();
        let seqvar = VariantRecord {
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
            super::passes(&hgnc_allowlist, &seqvar),
            expected,
            "hgnc_allowlist: {:?}, seqvar_genes: {:?}",
            hgnc_allowlist,
            seqvar_gene
        )
    }
}
