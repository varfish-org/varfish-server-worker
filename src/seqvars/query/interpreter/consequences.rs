use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};

/// Determine whether the `SequenceVariant` passes the consequences filter.
pub fn passes(query: &CaseQuery, seqvar: &SequenceVariant) -> Result<bool, anyhow::Error> {
    if query.consequences.is_empty() {
        return Ok(true);
    }

    let query_csq = std::collections::BTreeSet::from_iter(query.consequences.iter().cloned());
    for ann_field in &seqvar.ann_fields {
        let seqvar_csq =
            std::collections::BTreeSet::from_iter(ann_field.consequences.iter().cloned());
        let intersection_csq = query_csq.intersection(&seqvar_csq);
        if intersection_csq.count() > 0 {
            return Ok(true);
        }
    }

    Ok(false)
}

#[cfg(test)]
mod test {
    use mehari::annotate::seqvars::ann::{AnnField, Consequence};
    use rstest::rstest;
    use strum::IntoEnumIterator;

    use crate::seqvars::query::schema::{CaseQuery, SequenceVariant};

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn passes_consequence(#[case] c_equals_csq: bool) -> Result<(), anyhow::Error> {
        for csq in Consequence::iter() {
            let query = CaseQuery {
                consequences: Consequence::iter()
                    .filter(|c| (*c == csq) == c_equals_csq)
                    .collect(),
                ..Default::default()
            };
            let seq_var = SequenceVariant {
                reference: "G".into(),
                alternative: "A".into(),
                ann_fields: vec![AnnField {
                    allele: mehari::annotate::seqvars::ann::Allele::Alt {
                        alternative: "A".into(),
                    },
                    consequences: vec![csq],
                    putative_impact: csq.impact(),
                    gene_symbol: Default::default(),
                    gene_id: Default::default(),
                    feature_type: mehari::annotate::seqvars::ann::FeatureType::SoTerm {
                        term: mehari::annotate::seqvars::ann::SoFeature::Transcript,
                    },
                    feature_id: Default::default(),
                    feature_biotype: mehari::annotate::seqvars::ann::FeatureBiotype::Coding,
                    rank: Default::default(),
                    hgvs_t: Default::default(),
                    hgvs_p: Default::default(),
                    tx_pos: Default::default(),
                    cds_pos: Default::default(),
                    protein_pos: Default::default(),
                    distance: Default::default(),
                    messages: Default::default(),
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
