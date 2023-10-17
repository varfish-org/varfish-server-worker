//! Variant-related information.

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
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder)]
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
