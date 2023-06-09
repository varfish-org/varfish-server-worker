/// Phenotype-related algorithms.

/// Similarity computation using the Phenomizer method.
pub mod phenomizer {
    use hpo::{
        similarity::{Builtins, Similarity},
        term::{HpoGroup, InformationContentKind},
        Ontology,
    };

    /// Compute symmetric similarity score.
    pub fn score(q: &HpoGroup, d: &HpoGroup, o: &Ontology) -> f32 {
        let s = Builtins::Resnik(InformationContentKind::Gene);
        (score_dir(q, d, o, &s) + score_dir(d, q, o, &s)) / 2.0
    }

    /// "Directed" score part of phenomizer score.
    fn score_dir(qs: &HpoGroup, ds: &HpoGroup, o: &Ontology, s: &impl Similarity) -> f32 {
        // Handle case of empty `qs`.
        if qs.is_empty() {
            return 0f32;
        }

        // For each `q in qs` compute max similarity to any `d in ds`.
        let mut tmp: Vec<f32> = Vec::new();
        for q in qs {
            if let Some(q) = o.hpo(q) {
                tmp.push(
                    ds.iter()
                        .filter_map(|d| o.hpo(d).map(|d| q.similarity_score(&d, s)))
                        .max_by(|a, b| a.partial_cmp(b).expect("try to compare NaN"))
                        .unwrap_or_default(),
                );
            }
        }

        tmp.iter().sum::<f32>() / (qs.len() as f32)
    }
}

#[cfg(test)]
mod test {
    use super::phenomizer;
    use hpo::{annotations::OmimDiseaseId, term::HpoGroup, HpoTermId, Ontology};

    fn load_hpo() -> Result<Ontology, anyhow::Error> {
        Ok(Ontology::from_standard("tests/db/hpo")?)
    }

    fn prepare(terms: &[&str]) -> HpoGroup {
        HpoGroup::from(
            terms
                .iter()
                .map(|s| HpoTermId::from(s.to_string()))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn phenomizer_score_gene() -> Result<(), anyhow::Error> {
        let hpo = load_hpo()?;

        let query = &[
            // slender build
            "HP:0001533",
            // high, narrow palate
            "HP:0002705",
        ];
        let omim_marfan = hpo
            .omim_disease(&OmimDiseaseId::from(154700))
            .expect("marfan symdrome must be in HPO");
        let hpo_marfan = HpoGroup::from_iter(
            omim_marfan
                .to_hpo_set(&hpo)
                .child_nodes()
                .without_modifier()
                .into_iter(),
        );

        let score = phenomizer::score(&prepare(query), &hpo_marfan, &hpo);

        assert!((score - 1.7708597).abs() < 0.00001, "score = {}", score);

        Ok(())
    }
}
