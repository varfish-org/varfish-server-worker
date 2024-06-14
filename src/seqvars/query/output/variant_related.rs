//! Variant-related information.

use annonars::pbs::clinvar_data::clinvar_public::AggregateGermlineReviewStatus;

use crate::seqvars::query::{
    annonars::Annotator,
    output::variant_related::score_collection::{
        Collector, ExtremalValueCollector, SingleValueCollector,
    },
    schema::SequenceVariant,
};

/// Helper modules for score collection.
pub mod score_collection {
    /// Trait for score collection.
    pub trait Collector {
        /// Return column names of interest.
        fn column_names(&self) -> Vec<String>;
        /// Register one column value.
        fn register(&mut self, column_name: &str, value: &serde_json::Value);
        /// Write collected scores to the given `indexmap::IndexMap``.
        fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>);
    }

    /// Simple implementation for collecting a single score.
    #[derive(Debug, Clone)]
    pub struct SingleValueCollector {
        /// The column name to collect.
        pub input_column_name: String,
        /// The output column name.
        pub output_column_name: String,
        /// The collected value, if any.
        pub value: Option<serde_json::Value>,
    }

    impl SingleValueCollector {
        /// Construct given column name.
        pub fn new(input_column_name: &str, output_column_name: &str) -> Self {
            Self {
                input_column_name: input_column_name.into(),
                output_column_name: output_column_name.into(),
                value: None,
            }
        }
    }

    impl Collector for SingleValueCollector {
        fn column_names(&self) -> Vec<String> {
            vec![self.input_column_name.clone()]
        }

        fn register(&mut self, column_name: &str, value: &serde_json::Value) {
            if column_name == self.input_column_name && !value.is_null() {
                self.value = Some(value.clone());
            }
        }

        fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>) {
            if let Some(value) = self.value.as_ref() {
                dict.insert(self.input_column_name.clone(), value.clone());
            }
        }
    }

    /// Simple implementation for collecting aggregated min/max score.
    #[derive(Debug, Clone)]
    pub struct ExtremalValueCollector {
        /// Whether is max (else is min).
        pub is_max: bool,
        /// The output column name.
        pub output_column_name: String,
        /// The column names to collect.
        pub column_names: Vec<String>,
        /// The collected values, if any.
        pub values: indexmap::IndexMap<String, serde_json::Value>,
    }

    impl ExtremalValueCollector {
        /// Construct given column name.
        pub fn new(column_names: &[&str], output_column_name: &str, is_max: bool) -> Self {
            Self {
                is_max,
                output_column_name: output_column_name.to_string(),
                column_names: column_names.iter().map(|s| s.to_string()).collect(),
                values: Default::default(),
            }
        }
    }

    impl Collector for ExtremalValueCollector {
        fn column_names(&self) -> Vec<String> {
            self.column_names.clone()
        }

        fn register(&mut self, column_name: &str, value: &serde_json::Value) {
            if self.column_names.iter().any(|s| s.as_str() == column_name) && !value.is_null() {
                self.values.insert(column_name.to_string(), value.clone());
            }
        }

        fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>) {
            let mut sel_name = None;
            let mut sel_value = serde_json::Value::Null;
            for (name, value) in self.values.iter() {
                if value.is_null() || !value.is_number() {
                    // can only select numeric values that are not null
                    continue;
                }

                if sel_value.is_null()
                    || (self.is_max
                        && value.as_f64().unwrap_or_default()
                            > sel_value.as_f64().unwrap_or_default())
                    || (!self.is_max
                        && value.as_f64().unwrap_or_default()
                            < sel_value.as_f64().unwrap_or_default())
                {
                    sel_name = Some(name.clone());
                    sel_value = value.clone();
                }
            }

            if let Some(sel_name) = sel_name {
                dict.insert(self.output_column_name.clone(), sel_value);
                dict.insert(
                    format!("{}_argmax", self.output_column_name),
                    serde_json::Value::String(sel_name),
                );
            }
        }
    }
}

/// Simple

/// Record for variant-related scores.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Record {
    /// Precomputed scores.
    #[serde(skip_serializing_if = "indexmap::IndexMap::is_empty")]
    pub precomputed_scores: indexmap::IndexMap<String, serde_json::Value>,
    /// Database identifiers.
    #[serde(skip_serializing_if = "DbIds::is_empty")]
    pub db_ids: DbIds,
    /// Clinvar information.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub clinvar: Option<Clinvar>,
    /// Frequency information.
    #[serde(skip_serializing_if = "Frequency::is_empty")]
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
        seqvar: &SequenceVariant,
        annotator: &Annotator,
    ) -> Result<indexmap::IndexMap<String, serde_json::Value>, anyhow::Error> {
        let mut result = indexmap::IndexMap::new();

        // Extract values from CADD.
        if let Some(cadd_values) = annotator
            .query_cadd(seqvar)
            .as_ref()
            .map_err(|e| anyhow::anyhow!("problem querying CADD: {}", e))?
        {
            let mut collectors: Vec<Box<dyn Collector>> = vec![
                Box::new(SingleValueCollector::new("PHRED", "cadd_phred")),
                Box::new(SingleValueCollector::new("SIFTval", "sift")),
                Box::new(SingleValueCollector::new("PolyPhenVal", "polyphen")),
                Box::new(ExtremalValueCollector::new(
                    &[
                        "SpliceAI-acc-gain",
                        "SpliceAI-acc-loss",
                        "SpliceAI-don-gain",
                        "SpliceAI-don-loss",
                    ],
                    "spliceai",
                    true,
                )),
            ];

            for (column, value) in annotator
                .annonars_dbs
                .cadd_ctx
                .schema
                .columns
                .iter()
                .zip(cadd_values.iter())
            {
                for collector in collectors.iter_mut() {
                    collector.register(column.name.as_str(), value);
                }
            }

            collectors.iter_mut().for_each(|collector| {
                collector.write_to(&mut result);
            })
        }

        // Extract values from dbNSFP

        if let Some(dbnsfp_values) = annotator
            .query_dbnsfp(seqvar)
            .as_ref()
            .map_err(|e| anyhow::anyhow!("problem querying dbNSFP: {}", e))?
        {
            // REVEL_score
            // BayesDel_addAF_score
            let mut collectors: Vec<Box<dyn Collector>> = vec![
                Box::new(SingleValueCollector::new("REVEL_score", "revel")),
                Box::new(SingleValueCollector::new(
                    "BayesDel_addAF_score",
                    "bayesel_addaf",
                )),
            ];

            for (column, value) in annotator
                .annonars_dbs
                .cadd_ctx
                .schema
                .columns
                .iter()
                .zip(dbnsfp_values.iter())
            {
                for collector in collectors.iter_mut() {
                    collector.register(column.name.as_str(), value);
                }
            }

            collectors.iter_mut().for_each(|collector| {
                collector.write_to(&mut result);
            })
        }

        Ok(result)
    }

    /// Returns whether all fields are empty and would not be serialized.
    pub fn is_empty(&self) -> bool {
        self.precomputed_scores.is_empty()
            && self.db_ids.is_empty()
            && self.clinvar.is_none()
            && self.frequency.is_empty()
    }
}

/// Database identifiers.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct DbIds {
    /// The variant's dbSNP identifier present.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dbsnp_rs: Option<String>,
}

impl DbIds {
    /// Construct given sequence variant and annonars annotator.
    pub fn with_seqvar_and_annotator(
        seqvar: &SequenceVariant,
        annotator: &Annotator,
    ) -> Result<Self, anyhow::Error> {
        Ok(Self {
            dbsnp_rs: annotator
                .query_dbsnp(seqvar)
                .map_err(|e| anyhow::anyhow!("problem querying dbSNP: {}", e))?
                .map(|record| format!("rs{}", record.rs_id)),
        })
    }

    /// Returns whether all database identifiers are unset.
    pub fn is_empty(&self) -> bool {
        self.dbsnp_rs.is_none()
    }
}

/// ClinVar-related information.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Clinvar {
    /// The VCV accession.
    pub vcv: String,
    /// The clinical significance.
    pub germline_significance_description: String,
    /// The review status.
    pub germline_review_status: AggregateGermlineReviewStatus,
}

impl Clinvar {
    /// Construct given sequence variant and annonars annotator.
    pub fn with_seqvar_and_annotator(
        seqvar: &SequenceVariant,
        annotator: &Annotator,
    ) -> Result<Option<Self>, anyhow::Error> {
        let record = annotator
            .query_clinvar_minimal(seqvar)
            .map_err(|e| anyhow::anyhow!("problem querying clinvar-minimal: {}", e))?;
        if let Some(record) = record.as_ref() {
            if record.records.is_empty() {
                tracing::error!(
                    "variant {:?} found with empty list in ClinVar (should not happen)",
                    seqvar
                );
                return Ok(None);
            } else if record.records.len() > 1 {
                tracing::warn!(
                    "variant {:?} found list with {} entries, using first",
                    seqvar,
                    record.records.len()
                );
            }
            let vcv_record = &record.records[0];
            let accession = vcv_record.accession.as_ref().expect("no accession?");
            let vcv = format!("{}.{}", &accession.accession, accession.version);

            if let Some(agc) = vcv_record
                .classifications
                .as_ref()
                .and_then(|c| c.germline_classification.as_ref())
            {
                let germline_significance_description = if let Some(description) =
                    agc.description.as_ref()
                {
                    description.clone()
                } else {
                    tracing::error!("variant {:?} has germline classification without description (should not happen)", &seqvar);
                    return Ok(None);
                };
                let germline_review_status = if let Ok(review_status) =
                    AggregateGermlineReviewStatus::try_from(agc.review_status)
                {
                    review_status
                } else {
                    tracing::error!("variant {:?} has germline classification with invalid review status (should not happen)", &seqvar);
                    return Ok(None);
                };

                Ok(Some(Self {
                    vcv,
                    germline_significance_description,
                    germline_review_status,
                }))
            } else {
                tracing::trace!(
                    "variant {:?} has no germline classification (likely somatic only)",
                    &seqvar
                );
                Ok(None)
            }
        } else {
            Ok(None)
        }
    }
}

/// Frequency information.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder)]
pub struct Frequency {
    /// gnomAD-genomes frequency
    #[serde(skip_serializing_if = "Option::is_none")]
    #[builder(default)]
    pub gnomad_genomes: Option<NuclearFrequency>,
    /// gnomAD-exomes frequency
    #[serde(skip_serializing_if = "Option::is_none")]
    #[builder(default)]
    pub gnomad_exomes: Option<NuclearFrequency>,
    /// gnomad-mtDNA frequency
    #[serde(skip_serializing_if = "Option::is_none")]
    #[builder(default)]
    pub gnomad_mtdna: Option<MtdnaFrequency>,
    /// HelixMtDb frequency
    #[serde(skip_serializing_if = "Option::is_none")]
    #[builder(default)]
    pub helixmtdb: Option<MtdnaFrequency>,
    /// inhouse frequency
    #[serde(skip_serializing_if = "Option::is_none")]
    #[builder(default)]
    pub inhouse: Option<NuclearFrequency>,
}

impl Frequency {
    /// Extract frequency information from `seqvar`
    pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Frequency, anyhow::Error> {
        let chrom = annonars::common::cli::canonicalize(&seqvar.chrom);
        let frequency = if chrom == "MT" {
            FrequencyBuilder::default()
                .gnomad_genomes(
                    NuclearFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                        seqvar.gnomad_genomes_hemi,
                    )
                    .some_unless_empty(),
                )
                .gnomad_exomes(
                    NuclearFrequency::new(
                        seqvar.gnomad_exomes_af(),
                        seqvar.gnomad_exomes_an,
                        seqvar.gnomad_exomes_het,
                        seqvar.gnomad_exomes_hom,
                        seqvar.gnomad_exomes_hemi,
                    )
                    .some_unless_empty(),
                )
                .build()
        } else {
            FrequencyBuilder::default()
                .gnomad_mtdna(
                    MtdnaFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                    )
                    .some_unless_empty(),
                )
                .helixmtdb(
                    MtdnaFrequency::new(
                        seqvar.helixmtdb_af(),
                        seqvar.helix_an,
                        seqvar.helix_het,
                        seqvar.helix_hom,
                    )
                    .some_unless_empty(),
                )
                .build()
        }
        .map_err(|e| anyhow::anyhow!("could not build frequency information: {}", e))?;
        Ok(frequency)
    }

    /// Returns whether all frequencies are empty.
    pub fn is_empty(&self) -> bool {
        self.gnomad_genomes.is_none()
            && self.gnomad_exomes.is_none()
            && self.gnomad_mtdna.is_none()
            && self.helixmtdb.is_none()
            && self.inhouse.is_none()
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

impl NuclearFrequency {
    /// Return self or None if empty.
    pub fn some_unless_empty(self) -> Option<Self> {
        if self.allele_count == 0 {
            None
        } else {
            Some(self)
        }
    }
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

impl MtdnaFrequency {
    /// Return self or None if empty.
    pub fn some_unless_empty(self) -> Option<Self> {
        if self.allele_count == 0 {
            None
        } else {
            Some(self)
        }
    }
}
