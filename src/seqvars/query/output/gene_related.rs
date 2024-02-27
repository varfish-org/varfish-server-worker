//! Gene-related information for the gene.

use mehari::annotate::seqvars::ann::Consequence;

use crate::seqvars::query::{annonars::Annotator, schema::SequenceVariant};

/// Gene-related information for a `ResultPayload`.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Record {
    /// Gene identity related (for display of gene symbol).
    pub identity: Identity,
    /// Gene-related consequences of a variant.
    pub consequences: Consequences,
    /// Phenotype-related information, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phenotype: Option<Phenotype>,
    /// Gene-wise constraints on the gene, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
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
    pub fn with_seqvar_and_annotator(
        seqvar: &SequenceVariant,
        annotator: &Annotator,
    ) -> Result<Option<Self>, anyhow::Error> {
        if let Some(ann) = seqvar.ann_fields.first() {
            let hgnc_id = ann.gene_id.clone();

            let gene_record = annotator
                .query_genes(&hgnc_id)
                .map_err(|e| anyhow::anyhow!("problem querying genes database: {}", e))?;

            if !ann.gene_id.is_empty() && !ann.gene_symbol.is_empty() {
                return Ok(Some(Self {
                    identity: Identity::new(hgnc_id, ann.gene_symbol.clone()),
                    consequences: Consequences::new(
                        ann.hgvs_t
                            .clone()
                            .ok_or_else(|| anyhow::anyhow!("missing hgvs_t annotation"))?,
                        ann.hgvs_p.clone(),
                        ann.consequences.clone(),
                    ),
                    phenotype: gene_record.as_ref().map(Phenotype::with_gene_record),
                    constraints: gene_record.as_ref().and_then(|gene_record| {
                        gene_record
                            .gnomad_constraints
                            .as_ref()
                            .map(|gnomad_constraints| {
                                Constraints::with_constraints_record(gnomad_constraints)
                            })
                    }),
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
    /// Whether the gene is on the ACMG supplementary finding list.
    pub is_acmg_sf: bool,
    /// Whether the gene is a known disease gene.
    pub is_disease_gene: bool,
}

impl Phenotype {
    /// Construct given a `genes` database record.
    pub fn with_gene_record(gene_record: &annonars::pbs::genes::base::Record) -> Self {
        Self {
            is_acmg_sf: gene_record.acmg_sf.is_some(),
            is_disease_gene: gene_record.omim.is_some() || gene_record.orpha.is_some(),
        }
    }
}

/// Consequences related to a gene.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Consequences {
    /// HGVS.{c,n} code of variant
    pub hgvs_t: String,
    /// HGVS.p code of variant
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hgvs_p: Option<String>,

    /// The predicted variant consequences.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub consequences: Vec<Consequence>,
}

/// Result gene constraint information.
#[allow(clippy::too_many_arguments)]
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Constraints {
    /// gnomAD mis_z score
    pub gnomad_mis_z: f32,
    /// gnomAD oe_lof score
    pub gnomad_oe_lof: f32,
    /// gnomAD oe_lof_lower score
    pub gnomad_oe_lof_lower: f32,
    /// gnomAD oe_lof_upper score (LOEF)
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

impl Constraints {
    /// Construct given a `genes` database record.
    pub fn with_constraints_record(
        constraints: &annonars::pbs::genes::base::GnomadConstraintsRecord,
    ) -> Self {
        Self {
            gnomad_mis_z: constraints.mis_z.unwrap_or_default() as f32,
            gnomad_oe_lof: constraints.oe_lof.unwrap_or_default() as f32,
            gnomad_oe_lof_lower: constraints.oe_lof_lower.unwrap_or_default() as f32,
            gnomad_oe_lof_upper: constraints.oe_lof_upper.unwrap_or_default() as f32,
            gnomad_oe_mis: constraints.oe_mis.unwrap_or_default() as f32,
            gnomad_oe_mis_lower: constraints.oe_mis_lower.unwrap_or_default() as f32,
            gnomad_oe_mis_upper: constraints.oe_mis_upper.unwrap_or_default() as f32,
            gnomad_pli: constraints.pli.unwrap_or_default() as f32,
            gnomad_syn_z: constraints.syn_z.unwrap_or_default() as f32,
        }
    }
}
