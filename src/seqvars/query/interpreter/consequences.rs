use mehari::annotate::seqvars::ann;

use crate::seqvars::query::schema::{data::VariantRecord, query::CaseQuery};

/// Determine whether the `VariantRecord` passes the consequences filter.
pub fn passes(query: &CaseQuery, seqvar: &VariantRecord) -> Result<bool, anyhow::Error> {
    // If no consequences are specified, the variant passes.
    if query.consequence.consequences.is_empty() {
        return Ok(true);
    }
    // Variants on chrMT always pass.
    let chrom = annonars::common::cli::canonicalize(&seqvar.vcf_variant.chrom);
    if chrom == "MT" {
        return Ok(true);
    }

    let query_csq: indexmap::IndexSet<ann::Consequence> = indexmap::IndexSet::from_iter(
        query
            .consequence
            .consequences
            .iter()
            .cloned()
            .map(|c| c.into()),
    );
    for ann_field in &seqvar.ann_fields {
        let seqvar_csq: indexmap::IndexSet<ann::Consequence> =
            indexmap::IndexSet::from_iter(ann_field.consequences.iter().cloned());
        let intersection_csq = query_csq.intersection(&seqvar_csq);
        if intersection_csq.count() > 0 {
            return Ok(true);
        }
    }

    tracing::trace!(
        "variant {:?} fails consequence filter {:?}",
        &seqvar,
        &query.consequence
    );
    Ok(false)
}

#[cfg(test)]
mod test {
    use mehari::annotate::seqvars::ann;
    use rstest::rstest;
    use strum::IntoEnumIterator;

    use crate::seqvars::query::schema::{
        data::{VariantRecord, VcfVariant},
        query::{CaseQuery, QuerySettingsConsequence},
    };

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn passes_consequence(#[case] c_equals_csq: bool) -> Result<(), anyhow::Error> {
        use crate::seqvars::query::schema::query::Consequence;

        for csq in ann::Consequence::iter() {
            let query = CaseQuery {
                consequence: QuerySettingsConsequence {
                    consequences: Consequence::iter()
                        .filter(|c| {
                            (<Consequence as Into<ann::Consequence>>::into(*c) == csq)
                                == c_equals_csq
                        })
                        .collect(),
                    ..Default::default()
                },
                ..Default::default()
            };
            let seq_var = VariantRecord {
                vcf_variant: VcfVariant {
                    chrom: "1".into(),
                    pos: 1,
                    ref_allele: "G".into(),
                    alt_allele: "A".into(),
                    ..Default::default()
                },
                ann_fields: vec![ann::AnnField {
                    allele: mehari::annotate::seqvars::ann::Allele::Alt {
                        alternative: "A".into(),
                    },
                    consequences: vec![csq],
                    ..Default::default()
                }],
                ..Default::default()
            };

            assert_eq!(
                super::passes(&query, &seq_var)?,
                c_equals_csq,
                "csq = {:?}",
                &csq
            );
        }

        Ok(())
    }
}
