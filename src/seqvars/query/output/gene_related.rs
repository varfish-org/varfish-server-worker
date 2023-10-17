//! Gene-related information for the gene.

use mehari::annotate::seqvars::ann::Consequence;

use crate::seqvars::query::schema::SequenceVariant;

/// Gene-related information for a `ResultPayload`.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Record {
    /// Gene identity related (for display of gene symbol).
    pub identity: Identity,
    /// Gene-related consequences of a variant.
    pub consequences: Consequences,
    /// Phenotype-related information, if any.
    pub phenotype: Option<Phenotype>,
    /// Gene-wise constraints on the gene, if any.
    pub constraints: Option<Constraints>,
}

impl Record {
    /// Construct given a `SequenceVariant` if the information is given in the annotation.
    ///
    /// Note that we will only look at the first annotation record as the ingest creates
    /// one `SequenceVariant` record per gene.
    ///
    /// # Error
    ///
    /// Returns an error if `seqvar` does not contain all necessary information.
    pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Option<Self>, anyhow::Error> {
        if let Some(ann) = seqvar.ann_fields.first() {
            if !ann.gene_id.is_empty() && !ann.gene_symbol.is_empty() {
                return Ok(Some(Self {
                    identity: Identity::new(ann.gene_id.clone(), ann.gene_symbol.clone()),
                    consequences: Consequences::new(
                        ann.hgvs_t
                            .clone()
                            .ok_or_else(|| anyhow::anyhow!("missing hgvs_t annotation"))?,
                        ann.hgvs_p.clone(),
                        ann.consequences.clone(),
                    ),
                    // TODO: phenotype
                    // TODO: constraints
                    ..Default::default()
                }));
            }
        }
        Ok(None)
    }
}

/// Result information for gene identity.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Identity {
    /// HGNC gene ID.
    pub hgnc_id: String,
    /// HGNC gene symbol.
    pub hgnc_symbol: String,
}

/// Result information related to phenotype / disease.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Phenotype {
    /// Whether the gene is a known disease gene.
    pub is_disease_gene: bool,
    // TODO: modes of inheritance
    // TODO: disease/phenotype terms?
}

/// Consequences related to a gene.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Consequences {
    /// HGVS.{c,n} code of variant
    pub hgvs_t: String,
    /// HGVS.p code of variant
    pub hgvs_p: Option<String>,

    /// The predicted variant consequences.
    pub consequences: Vec<Consequence>,
}

/// Result gene constraint information.
#[allow(clippy::too_many_arguments)]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Constraints {
    /// gnomAD loeuf score
    pub gnomad_loeuf: f32,
    /// gnomAD mis_z score
    pub gnomad_mis_z: f32,
    /// gnomAD oe_lof score
    pub gnomad_oe_lof: f32,
    /// gnomAD oe_lof_lower score
    pub gnomad_oe_lof_lower: f32,
    /// gnomAD oe_lof_upper score
    pub gnomad_oe_lof_upper: f32,
    /// gnomAD oe_mis score
    pub gnomad_oe_mis: f32,
    /// gnomAD oe_mis_lower score
    pub gnomad_oe_mis_lower: f32,
    /// gnomAD oe_mis_upper score
    pub gnomad_oe_mis_upper: f32,
    /// gnomAD pLI score
    pub gnomad_pli: f32,
    /// gnomAD syn_z score
    pub gnomad_syn_z: f32,
}
