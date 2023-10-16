use std::collections::HashSet;

use crate::seqvars::query::schema::SequenceVariant;

/// Determine whether the `SequenceVariant` passes the genes allowlist filter.
pub fn passes(hgnc_allowlist: &Option<HashSet<String>>, seqvar: &SequenceVariant) -> bool {
    if let Some(hgnc_allowlist) = &hgnc_allowlist {
        if hgnc_allowlist.is_empty() {
            true
        } else {
            seqvar
                .ann_fields
                .iter()
                .any(|ann_field| hgnc_allowlist.contains(&ann_field.gene_id))
        }
    } else {
        true
    }
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::seqvars::query::schema::SequenceVariant;
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
        #[case] hgnc_allowlist: Option<Vec<String>>,
        #[case] seqvar_gene: Option<String>,
        #[case] expected: bool,
    ) {
        let hgnc_allowlist = hgnc_allowlist.map(|hgnc_allowlist| {
            hgnc_allowlist
                .into_iter()
                .map(|hgnc_id| hgnc_id.to_uppercase())
                .collect::<std::collections::HashSet<_>>()
        });
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
            super::passes(&hgnc_allowlist, &seqvar),
            expected,
            "hgnc_allowlist: {:?}, seqvar_genes: {:?}",
            hgnc_allowlist,
            seqvar_gene
        )
    }
}
