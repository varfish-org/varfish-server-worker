//! Data structureds for writing the output.

/// Gene-related information for the gene.
pub mod gene_related {
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
}

/// Variant-related information.
pub mod variant_related {
    use crate::seqvars::query::{annonars::Annotator, schema::SequenceVariant};

    /// Record for variant-related scores.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Record {
        /// Precomputed scores.
        pub precomputed_scores: indexmap::IndexMap<String, f32>,
        /// Database identifiers.
        pub db_ids: DbIds,
        /// Clinvar information.
        pub clinvar: Option<Clinvar>,
        /// Frequency information.
        pub frequency: Frequency,
    }

    impl Record {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Self, anyhow::Error> {
            Ok(Self {
                precomputed_scores: Self::query_precomputed_scores(seqvar, annotator)?,
                db_ids: DbIds::with_seqvar_and_annotator(seqvar, annotator)?,
                clinvar: Clinvar::with_seqvar_and_annotator(seqvar, annotator)?,
                frequency: Frequency::with_seqvar(seqvar)?,
            })
        }

        /// Query precomputed scores for `seqvar` from annonars `annotator`.
        pub fn query_precomputed_scores(
            _seqvar: &SequenceVariant,
            _annotator: &Annotator,
        ) -> Result<indexmap::IndexMap<String, f32>, anyhow::Error> {
            // TODO: implement me!
            Ok(indexmap::IndexMap::default())
        }
    }

    /// Database identifiers.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct DbIds {
        /// The variant's dbSNP identifier present.
        pub dbsnp_rs: Option<String>,
    }

    impl DbIds {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Self, anyhow::Error> {
            // TODO: need to properly separate error from no result in annonars.
            Ok(Self {
                dbsnp_rs: annotator
                    .query_dbsnp(seqvar)
                    .ok()
                    .map(|record| format!("rs{}", record.rs_id)),
            })
        }
    }

    /// ClinVar-related information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Clinvar {
        // TODO: switch to VCV summary or hierarchical?
        /// The RCV accession.
        pub rcv: String,
        /// The clinical significance.
        pub significance: String,
        /// The review status.
        pub review_status: String,
    }

    impl Clinvar {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Option<Self>, anyhow::Error> {
            // TODO: need to properly separate error from no result in annonars.
            annotator
                .query_clinvar_minimal(seqvar)
                .ok()
                .map(|record| {
                    let annonars::clinvar_minimal::pbs::Record {
                        rcv,
                        clinical_significance,
                        review_status,
                        ..
                    } = record;
                    Ok(Self {
                        rcv,
                        significance: match clinical_significance {
                            0 => "Pathogenic",
                            1 => "Likely pathogenic",
                            2 => "Uncertain significance",
                            3 => "Likely benign",
                            4 => "Benign",
                            _ => anyhow::bail!(
                                "invalid clinical significance enum: {}",
                                clinical_significance
                            ),
                        }
                        .into(),
                        review_status: match review_status {
                            0 => "no assertion provided",
                            1 => "no assertion criteria provided",
                            2 => "criteria provided, conflicting interpretations",
                            3 => "criteria provided, single submitter",
                            4 => "criteria provided, multiple submitters, no conflicts",
                            5 => "reviewed by expert panel",
                            6 => "practice guideline",
                            _ => anyhow::bail!("invalid review status enum: {}", review_status),
                        }
                        .into(),
                    })
                })
                .transpose()
        }
    }

    /// Frequency information.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Frequency {
        /// gnomAD-genomes frequency
        pub gnomad_genomes: Option<NuclearFrequency>,
        /// gnomAD-exomes frequency
        pub gnomad_exomes: Option<NuclearFrequency>,
        /// gnomad-mtDNA frequency
        pub gnomad_mtdna: Option<MtdnaFrequency>,
        /// HelixMtDb frequency
        pub helixmtdb: Option<MtdnaFrequency>,
        /// in-house frequency
        pub inhouse: Option<NuclearFrequency>,
    }

    impl Frequency {
        /// Extract frequency information from `seqvar`
        pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Frequency, anyhow::Error> {
            let chrom = annonars::common::cli::canonicalize(&seqvar.chrom);
            let frequency = if chrom == "MT" {
                FrequencyBuilder::default()
                    .gnomad_genomes(Some(NuclearFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                        seqvar.gnomad_genomes_hemi,
                    )))
                    .gnomad_exomes(Some(NuclearFrequency::new(
                        seqvar.gnomad_exomes_af(),
                        seqvar.gnomad_exomes_an,
                        seqvar.gnomad_exomes_het,
                        seqvar.gnomad_exomes_hom,
                        seqvar.gnomad_exomes_hemi,
                    )))
                    .build()
            } else {
                FrequencyBuilder::default()
                    .gnomad_mtdna(Some(MtdnaFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                    )))
                    .helixmtdb(Some(MtdnaFrequency::new(
                        seqvar.helixmtdb_af(),
                        seqvar.helix_an,
                        seqvar.helix_het,
                        seqvar.helix_hom,
                    )))
                    .build()
            }
            .map_err(|e| anyhow::anyhow!("could not build frequency information: {}", e))?;
            Ok(frequency)
        }
    }

    /// Nuclear frequency information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct NuclearFrequency {
        /// Overall allele frequency.
        pub allele_freq: f32,
        /// Number of alleles.
        pub allele_count: i32,
        /// Number of heterozygous carriers.
        pub het_carriers: i32,
        /// Number of homozygous carriers.
        pub hom_carriers: i32,
        /// Number of hemizygous carriers.
        pub hemi_carriers: i32,
    }

    /// Mitochondrial frequency information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct MtdnaFrequency {
        /// Overall allele frequency.
        pub allele_freq: f32,
        /// Number of alleles.
        pub allele_count: i32,
        /// Number of heterplasmic carriers.
        pub het_carriers: i32,
        /// Number of homoplasmic carriers.
        pub hom_carriers: i32,
    }
}

/// Call-related information.
pub mod call_related {
    use crate::seqvars::query::schema::SequenceVariant;

    /// Call-related record.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Record {
        /// The genotype information for each sample.
        pub call_info: indexmap::IndexMap<String, CallInfo>,
    }

    impl Record {
        /// Construct a new `Record` from a `SequenceVariant`.
        ///
        /// # Error
        ///
        /// Returns an error if the `SequenceVariant` does not contain all necessary information.
        pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Self, anyhow::Error> {
            Ok(Self {
                call_info: seqvar
                    .call_info
                    .iter()
                    .map(|(sample_name, call_info)| {
                        (
                            sample_name.clone(),
                            CallInfo::new(
                                call_info.dp,
                                call_info.ad,
                                call_info.quality.map(|q| q as i32),
                                call_info.genotype.clone(),
                            ),
                        )
                    })
                    .collect(),
            })
        }
    }

    /// Genotype information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct CallInfo {
        /// Depth of coverage.
        pub dp: Option<i32>,
        /// Alternate read depth.
        pub ad: Option<i32>,
        /// Genotype quality.
        pub gq: Option<i32>,
        /// Genotype.
        pub gt: Option<String>,
    }
}

/// Result record information.
pub mod result {
    /// A result record from the query.
    ///
    /// These records are written to TSV for import into the database.   They contain the
    /// bare necessary information for sorting etc. in the database.  The actual main
    /// data is in the payload, serialized as JSON.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Record {
        /// UUID for the record.
        pub sodar_uuid: uuid::Uuid,
        /// Genome release for the coordinate.
        pub release: String,
        /// Chromosome name.
        pub chromosome: String,
        /// Chromosome number.
        pub chromosome_no: i32,
        /// Reference allele sequence.
        pub reference: String,
        /// Alternative allele sequence.
        pub alternative: String,
        /// UCSC bin of the record.
        pub bin: u32,
        /// Start position of the record.
        pub start: i32,
        /// End position of the record.
        pub end: i32,
        /// The result set ID as specified on the command line.
        pub smallvariantqueryresultset_id: String,
        /// The JSON-serialized `ResultPayload`.
        pub payload: String,
    }

    /// The structured result information of the result record.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Payload {
        /// Case UUID as specified on the command line.
        pub case_uuid: uuid::Uuid,
        /// The affected gene and consequence, if any.
        pub gene_related: Option<super::gene_related::Record>,
        /// Variant-related information, always present.
        pub variant_related: super::variant_related::Record,
        /// Genotypes call related, always present.
        pub call_related: super::call_related::Record,
    }
}
