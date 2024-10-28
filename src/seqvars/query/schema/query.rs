//! Code for representing query definitions.
//!
//! We generally use protobuf JSON for representing queries when stored.
//! After deserialization, they are converted into the data structures defined
//! here.

use crate::pbs::varfish::v1::seqvars::query as pb_query;

/// Enumeration for recessvive mode queries.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
pub enum RecessiveMode {
    /// Disabled recessive mode,
    #[default]
    Disabled,
    /// Compound heterozygous recessive mode.
    CompoundHeterozygous,
    /// Homozygous recessive mode.
    Homozygous,
    /// Generic recessive mode.
    Any,
}

/// Supporting code for `RecessiveMode`.
pub(crate) mod recessive_mode {
    /// Error type for `RecessiveMode::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf RecessiveMode: {0:?}")]
        UnknownRecessiveModeValue(super::pb_query::RecessiveMode),
    }
}

impl TryFrom<pb_query::RecessiveMode> for RecessiveMode {
    type Error = recessive_mode::Error;

    fn try_from(value: pb_query::RecessiveMode) -> Result<Self, Self::Error> {
        match value {
            pb_query::RecessiveMode::Disabled => Ok(RecessiveMode::Disabled),
            pb_query::RecessiveMode::CompoundHeterozygous => {
                Ok(RecessiveMode::CompoundHeterozygous)
            }
            pb_query::RecessiveMode::Homozygous => Ok(RecessiveMode::Homozygous),
            pb_query::RecessiveMode::Any => Ok(RecessiveMode::Any),
            _ => Err(recessive_mode::Error::UnknownRecessiveModeValue(value)),
        }
    }
}

/// Enumeration type for genotype choice.
#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    serde::Serialize,
    serde::Deserialize,
)]
pub enum GenotypeChoice {
    /// Any genotype.
    #[default]
    Any,
    /// Reference genotype.
    Ref,
    /// Heterozygous genotype.
    Het,
    /// Homozygous genotype.
    Hom,
    /// Non-heterozygous genotype.
    NonHet,
    /// Non-homozygous genotype.
    NonHom,
    /// Variant genotype.
    Variant,
    /// Recessive index.
    RecessiveIndex,
    /// Recessive father.
    RecessiveFather,
    /// Recessive mother.
    RecessiveMother,
}

/// Supporting code for `GenotypeChoice`.
pub(crate) mod genotype_choice {
    /// Error type for `GenotypeChoice::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum TryFromError {
        #[error("Cannot convert protobuf GenotypeChoice: {0:?}")]
        UnknownGenotypeChoiceValue(super::pb_query::GenotypeChoice),
    }

    /// Error type for `GenotypeChoice::matches()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum MatchesError {
        #[error("Cannot use genotype matches on recessive indicator: {0:?}")]
        RecessiveIndicator(super::GenotypeChoice),
    }
}

impl TryFrom<pb_query::GenotypeChoice> for GenotypeChoice {
    type Error = genotype_choice::TryFromError;

    fn try_from(value: pb_query::GenotypeChoice) -> Result<Self, Self::Error> {
        match value {
            pb_query::GenotypeChoice::Any => Ok(GenotypeChoice::Any),
            pb_query::GenotypeChoice::Ref => Ok(GenotypeChoice::Ref),
            pb_query::GenotypeChoice::Het => Ok(GenotypeChoice::Het),
            pb_query::GenotypeChoice::Hom => Ok(GenotypeChoice::Hom),
            pb_query::GenotypeChoice::NonHet => Ok(GenotypeChoice::NonHet),
            pb_query::GenotypeChoice::NonHom => Ok(GenotypeChoice::NonHom),
            pb_query::GenotypeChoice::Variant => Ok(GenotypeChoice::Variant),
            pb_query::GenotypeChoice::RecessiveIndex => Ok(GenotypeChoice::RecessiveIndex),
            pb_query::GenotypeChoice::RecessiveFather => Ok(GenotypeChoice::RecessiveFather),
            pb_query::GenotypeChoice::RecessiveMother => Ok(GenotypeChoice::RecessiveMother),
            _ => Err(Self::Error::UnknownGenotypeChoiceValue(value)),
        }
    }
}

/// Returns whether the given genotype script is treated as no-call.
///
/// This is the case if the genotype string contains at least one ".".
pub fn considered_no_call(gt_str: &str) -> bool {
    gt_str.contains('.')
}

/// Trait that describes whether a string matches a value.
///
/// Note that we assume properly ingested VCFs with only one alternate allele.
/// The valid genotype strings have the form "<VAL>/<VAL>", "<VAL>|<VAL>" or
/// "<VAL>" with "<VAL>" being one of "0", "1", and ".".
pub trait MatchesGenotypeStr {
    type Error;

    /// Whether `self` matches `s`.
    ///
    /// Note that no-call genotypes do not match anything.
    ///
    /// # Arguments
    ///
    /// * `gt_str` - The string to match against.
    ///
    /// # Returns
    ///
    /// `true` if `self` matches `gt_str`, `false` otherwise.
    ///
    /// # Errors
    ///
    /// * When the conversion fails.
    fn matches(&self, gt_tr: &str) -> Result<bool, Self::Error>;
}

impl MatchesGenotypeStr for GenotypeChoice {
    type Error = genotype_choice::MatchesError;

    fn matches(&self, gt_str: &str) -> Result<bool, Self::Error> {
        let gt_str = if gt_str.starts_with('/') || gt_str.starts_with('|') {
            &gt_str[1..]
        } else {
            gt_str
        };
        Ok(match self {
            // atoms
            GenotypeChoice::Ref => ["0", "0|0", "0/0"].contains(&gt_str),
            GenotypeChoice::Het => ["0/1", "0|1", "1/0", "1|0"].contains(&gt_str),
            GenotypeChoice::Hom => ["1", "1/1", "1|1"].contains(&gt_str),
            // combinations
            GenotypeChoice::Variant => {
                GenotypeChoice::Het.matches(gt_str)? || GenotypeChoice::Hom.matches(gt_str)?
            }
            GenotypeChoice::Any => {
                GenotypeChoice::Ref.matches(gt_str)? || GenotypeChoice::Variant.matches(gt_str)?
            }
            GenotypeChoice::NonHom => {
                GenotypeChoice::Ref.matches(gt_str)? || GenotypeChoice::Het.matches(gt_str)?
            }
            GenotypeChoice::NonHet => {
                GenotypeChoice::Ref.matches(gt_str)? || GenotypeChoice::Hom.matches(gt_str)?
            }
            // recessive markers
            GenotypeChoice::RecessiveIndex
            | GenotypeChoice::RecessiveFather
            | GenotypeChoice::RecessiveMother => {
                return Err(Self::Error::RecessiveIndicator(*self))
            }
        })
    }
}

/// Query for a single sample.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SampleGenotypeChoice {
    /// Name of the sample filtered for.
    pub sample: String,
    /// Genotype choice.
    pub genotype: GenotypeChoice,
    /// Whether to include no-call (will disable quality filter).
    pub include_no_call: bool,
    /// Whether to enable sample in filtration.
    pub enabled: bool,
}

impl Default for SampleGenotypeChoice {
    fn default() -> Self {
        Self {
            sample: Default::default(),
            genotype: Default::default(),
            include_no_call: false,
            enabled: true,
        }
    }
}

/// Supporting code for `SampleGenotypeChoice`.
pub(crate) mod sample_genotype_choice {
    /// Error type for `SampleGenotypeChoice::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert i32 into protobuf GenotypeChoice: {0}")]
        UnknownGenotypeChoiceInt(i32),
        #[error("Cannot convert protobuf GenotypeChoice: {0:?}")]
        UnknownGenotypeChoiceValue(super::pb_query::GenotypeChoice),
    }
}

impl TryFrom<pb_query::SampleGenotypeChoice> for SampleGenotypeChoice {
    type Error = sample_genotype_choice::Error;

    fn try_from(value: pb_query::SampleGenotypeChoice) -> Result<Self, Self::Error> {
        let pb_genotype = pb_query::GenotypeChoice::try_from(value.genotype)
            .map_err(|_| Self::Error::UnknownGenotypeChoiceInt(value.genotype))?;
        Ok(Self {
            sample: value.sample,
            genotype: GenotypeChoice::try_from(pb_genotype)
                .map_err(|_| Self::Error::UnknownGenotypeChoiceValue(pb_genotype))?,
            include_no_call: value.include_no_call,
            enabled: value.enabled,
        })
    }
}

/// Query settings for genotypes.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsGenotype {
    /// Recessive mode.
    pub recessive_mode: RecessiveMode,
    /// Mapping from sample name to sample genotype choice.
    pub sample_genotypes: indexmap::IndexMap<String, SampleGenotypeChoice>,
}

/// Support code for `QuerySettingsGenotype`.
pub(crate) mod query_settings_genotype {
    /// Error type for `QuerySettingsGenotype::recessive_index()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum RecessiveIndexError {
        #[error("No recessive index sample found")]
        NoRecessiveIndexSample,
        #[error("Multiple recessive index samples found: {0:?}")]
        MultipleRecessiveIndexSamples(Vec<String>),
    }

    /// Error type for `QuerySettingsGenotype::recessive_parents()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum RecessiveParentsError {
        #[error("Found two or more recessive fathers: {0:?}")]
        Fathers(Vec<String>),
        #[error("Found two or more recessive mothers: {0:?}")]
        Mothers(Vec<String>),
    }

    /// Error type for `QuerySettingsGenotype::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert i32 into protobuf RecessiveMode: {0}")]
        UnknownRecessiveModeInt(i32),
        #[error("Cannot convert protobuf RecessiveMode: {0:?}")]
        UnknownRecessiveModeValue(super::pb_query::RecessiveMode),
        #[error("Sample occured twice in sample_genotypes: {0:?}")]
        DuplicateSample(String),
        #[error("Invalid sample genotype choice: {0}")]
        InvalidSampleGenotypeChoice(#[from] super::sample_genotype_choice::Error),
    }
}

/// Struct for storing recessive mother/father.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct RecessiveParents {
    /// Name of the father, if any.
    pub father: Option<String>,
    /// Name of the mother, if any.
    pub mother: Option<String>,
}

impl From<RecessiveParents> for Vec<String> {
    fn from(val: RecessiveParents) -> Self {
        let mut result = Vec::new();
        if let Some(father) = val.father {
            result.push(father);
        }
        if let Some(mother) = val.mother {
            result.push(mother);
        }
        result
    }
}

impl QuerySettingsGenotype {
    /// Return the selected recessive index sample.
    ///
    /// # Returns
    ///
    /// The sample name of the recessive index sample.
    ///
    /// # Errors
    ///
    /// * `RecessiveIndexError::NoRecessiveIndexSample` if no recessive index sample is found.
    /// * `RecessiveIndexError::MultipleRecessiveIndexSamples` if multiple recessive index samples are found.
    pub fn recessive_index(&self) -> Result<String, query_settings_genotype::RecessiveIndexError> {
        let samples = self
            .sample_genotypes
            .values()
            .filter(|sgc| sgc.genotype == GenotypeChoice::RecessiveIndex)
            .map(|sgc| sgc.sample.clone())
            .collect::<Vec<_>>();
        if samples.is_empty() {
            Err(query_settings_genotype::RecessiveIndexError::NoRecessiveIndexSample)
        } else if samples.len() > 1 {
            Err(
                query_settings_genotype::RecessiveIndexError::MultipleRecessiveIndexSamples(
                    samples,
                ),
            )
        } else {
            Ok(samples[0].clone())
        }
    }

    /// Return the selected recessive parent samples.
    ///
    /// # Returns
    ///
    /// Struct holding names of recessive father and mother, if any.
    ///
    /// # Errors
    ///
    /// * `RecessiveParentsError::Fathers` or `RecessiveParentsError::Mothers` if more than two
    ///   recessive fathers/mothers samples are found.
    pub fn recessive_parents(
        &self,
    ) -> Result<RecessiveParents, query_settings_genotype::RecessiveParentsError> {
        let fathers = self
            .sample_genotypes
            .values()
            .filter(|sgc| matches!(sgc.genotype, GenotypeChoice::RecessiveFather))
            .map(|sgc| sgc.sample.clone())
            .collect::<Vec<_>>();
        if fathers.len() > 2 {
            return Err(query_settings_genotype::RecessiveParentsError::Fathers(
                fathers,
            ));
        }

        let mothers = self
            .sample_genotypes
            .values()
            .filter(|sgc| matches!(sgc.genotype, GenotypeChoice::RecessiveMother))
            .map(|sgc| sgc.sample.clone())
            .collect::<Vec<_>>();
        if mothers.len() > 2 {
            return Err(query_settings_genotype::RecessiveParentsError::Mothers(
                mothers,
            ));
        }

        Ok(RecessiveParents {
            father: fathers.into_iter().next(),
            mother: mothers.into_iter().next(),
        })
    }
}

impl TryFrom<pb_query::QuerySettingsGenotype> for QuerySettingsGenotype {
    type Error = query_settings_genotype::Error;

    fn try_from(value: pb_query::QuerySettingsGenotype) -> Result<Self, Self::Error> {
        let pb_recessive_mode = pb_query::RecessiveMode::try_from(value.recessive_mode)
            .map_err(|_| Self::Error::UnknownRecessiveModeInt(value.recessive_mode))?;
        let recessive_mode = RecessiveMode::try_from(pb_recessive_mode)
            .map_err(|_| Self::Error::UnknownRecessiveModeValue(pb_recessive_mode))?;

        let mut sample_genotypes = indexmap::IndexMap::new();
        for pb_sample_genotype in value.sample_genotypes {
            let sample_genotype = SampleGenotypeChoice::try_from(pb_sample_genotype)
                .map_err(Self::Error::InvalidSampleGenotypeChoice)?;
            if sample_genotypes.contains_key(&sample_genotype.sample) {
                return Err(Self::Error::DuplicateSample(sample_genotype.sample));
            }
            sample_genotypes.insert(sample_genotype.sample.clone(), sample_genotype);
        }

        Ok(Self {
            recessive_mode,
            sample_genotypes,
        })
    }
}

/// Quality settings for one sample.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SampleQualitySettings {
    /// Name of the sample filtered for.
    pub sample: String,
    /// Drop whole variant on failure.
    pub filter_active: bool,
    /// Minimal coverage for het. sites.
    pub min_dp_het: Option<i32>,
    /// Minimal coverage for hom. sites.
    pub min_dp_hom: Option<i32>,
    /// Minimal genotype quality.
    pub min_gq: Option<i32>,
    /// Minimal allele balance for het. variants.
    pub min_ab: Option<f32>,
    /// Minimal number of alternate reads.
    pub min_ad: Option<i32>,
    /// Maximal number of alternate reads.
    pub max_ad: Option<i32>,
}

impl Eq for SampleQualitySettings {}

impl From<pb_query::SampleQualitySettings> for SampleQualitySettings {
    fn from(value: pb_query::SampleQualitySettings) -> Self {
        Self {
            sample: value.sample,
            filter_active: value.filter_active,
            min_dp_het: value.min_dp_het,
            min_dp_hom: value.min_dp_hom,
            min_gq: value.min_gq,
            min_ab: value.min_ab,
            min_ad: value.min_ad,
            max_ad: value.max_ad,
        }
    }
}

/// Per-sample quality filter settings.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsQuality {
    /// Mapping from sample name to sample quality settings.
    pub sample_qualities: indexmap::IndexMap<String, SampleQualitySettings>,
}

/// Supporting code for `QuerySettingsQuality`.
pub(crate) mod query_settings_quality {
    /// Error type for `QuerySettingsQuality::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Sample occured twice in sample_qualities: {0:?}")]
        DuplicateSample(String),
    }
}

impl TryFrom<pb_query::QuerySettingsQuality> for QuerySettingsQuality {
    type Error = query_settings_quality::Error;

    fn try_from(value: pb_query::QuerySettingsQuality) -> Result<Self, Self::Error> {
        let mut sample_qualities = indexmap::IndexMap::new();
        for pb_sample_quality in value.sample_qualities {
            let sample_quality = SampleQualitySettings::from(pb_sample_quality);
            if sample_qualities.contains_key(&sample_quality.sample) {
                return Err(Self::Error::DuplicateSample(sample_quality.sample));
            }
            sample_qualities.insert(sample_quality.sample.clone(), sample_quality);
        }
        Ok(Self { sample_qualities })
    }
}

/// gnomAD and In-house nuclear filter options.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct NuclearFrequencySettings {
    /// Whether to enable filtration by 1000 Genomes.
    pub enabled: bool,
    /// Maximal number of in-house heterozygous carriers.
    pub max_het: Option<i32>,
    /// Maximal number of in-house homozygous carriers.
    pub max_hom: Option<i32>,
    /// Maximal number of in-house hemizygous carriers.
    pub max_hemi: Option<i32>,
    /// Maximal allele frequency.
    pub max_af: Option<f32>,
}

impl Eq for NuclearFrequencySettings {}

impl From<pb_query::NuclearFrequencySettings> for NuclearFrequencySettings {
    fn from(value: pb_query::NuclearFrequencySettings) -> Self {
        Self {
            enabled: value.enabled,
            max_het: value.max_het,
            max_hom: value.max_hom,
            max_hemi: value.max_hemi,
            max_af: value.max_af,
        }
    }
}

/// HelixMtDb filter options.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MitochondrialFrequencySettings {
    /// Whether to enable filtration by mtDB.
    pub enabled: bool,
    /// Maximal number of heterozygous carriers in HelixMtDb.
    pub max_het: Option<i32>,
    /// Maximal number of homozygous carriers in HelixMtDb.
    pub max_hom: Option<i32>,
    /// Maximal frequency in HelixMtDb.
    pub max_af: Option<f32>,
}

impl Eq for MitochondrialFrequencySettings {}

impl From<pb_query::MitochondrialFrequencySettings> for MitochondrialFrequencySettings {
    fn from(value: pb_query::MitochondrialFrequencySettings) -> Self {
        Self {
            enabled: value.enabled,
            max_het: value.max_het,
            max_hom: value.max_hom,
            max_af: value.max_af,
        }
    }
}

/// Query settings for population frequencies.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsFrequency {
    /// gnomAD-exomes filter.
    pub gnomad_exomes: NuclearFrequencySettings,
    /// gnomAD-genomes filter.
    pub gnomad_genomes: NuclearFrequencySettings,
    /// gnomAD-mtDNA filter.
    pub gnomad_mtdna: MitochondrialFrequencySettings,
    /// HelixMtDb filter.
    pub helixmtdb: MitochondrialFrequencySettings,
    /// In-house filter.
    pub inhouse: NuclearFrequencySettings,
}

impl From<pb_query::QuerySettingsFrequency> for QuerySettingsFrequency {
    fn from(value: pb_query::QuerySettingsFrequency) -> Self {
        Self {
            gnomad_exomes: NuclearFrequencySettings::from(value.gnomad_exomes.unwrap_or_default()),
            gnomad_genomes: NuclearFrequencySettings::from(
                value.gnomad_genomes.unwrap_or_default(),
            ),
            gnomad_mtdna: MitochondrialFrequencySettings::from(
                value.gnomad_mtdna.unwrap_or_default(),
            ),
            helixmtdb: MitochondrialFrequencySettings::from(value.helixmtdb.unwrap_or_default()),
            inhouse: NuclearFrequencySettings::from(value.inhouse.unwrap_or_default()),
        }
    }
}

/// Variant Types.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum VariantType {
    /// SNV.
    Snv,
    /// Indel.
    Indel,
    /// MNV.
    Mnv,
    /// Complex substitution.
    ComplexSubstitution,
}

/// Supporting code for `VariantType`.
pub(crate) mod variant_type {
    /// Error type for `VariantType::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf VariantType: {0:?}")]
        UnknownVariantTypeValue(super::pb_query::VariantType),
    }
}

impl TryFrom<pb_query::VariantType> for VariantType {
    type Error = variant_type::Error;

    fn try_from(value: pb_query::VariantType) -> Result<Self, Self::Error> {
        match value {
            pb_query::VariantType::Snv => Ok(VariantType::Snv),
            pb_query::VariantType::Indel => Ok(VariantType::Indel),
            pb_query::VariantType::Mnv => Ok(VariantType::Mnv),
            pb_query::VariantType::ComplexSubstitution => Ok(VariantType::ComplexSubstitution),
            _ => Err(variant_type::Error::UnknownVariantTypeValue(value)),
        }
    }
}

/// Transcript types.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum TranscriptType {
    /// Coding transcripts.
    Coding,
    /// Non-coding transcripts.
    NonCoding,
}

/// Supporting code for `TranscriptType`.
pub(crate) mod transcript_type {
    /// Error type for `TranscriptType::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf TranscriptType: {0:?}")]
        UnknownTranscriptTypeValue(super::pb_query::TranscriptType),
    }
}

impl TryFrom<pb_query::TranscriptType> for TranscriptType {
    type Error = transcript_type::Error;

    fn try_from(value: pb_query::TranscriptType) -> Result<Self, Self::Error> {
        match value {
            pb_query::TranscriptType::Coding => Ok(TranscriptType::Coding),
            pb_query::TranscriptType::NonCoding => Ok(TranscriptType::NonCoding),
            _ => Err(transcript_type::Error::UnknownTranscriptTypeValue(value)),
        }
    }
}

/// Consequence types.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    serde::Serialize,
    serde::Deserialize,
    strum::EnumIter,
)]
pub enum Consequence {
    // high impact
    /// "A feature ablation whereby the deleted region includes a transcript feature."
    /// SO:transcript_ablation, VEP:transcript_ablation
    TranscriptAblation,

    /// "A sequence variant whereby an exon is lost from the transcript."
    /// SO:exon_loss_variant, VEP:transcript_ablation
    ExonLossVariant,

    /// "A splice variant that changes the 2 base region at the 3' end of an intron."
    /// SO:splice_acceptor_variant, VEP:splice_acceptor_variant
    SpliceAcceptorVariant,

    /// "A splice variant that changes the 2 base region at the 5' end of an intron."
    /// SO:splice_donor_variant, VEP:splice_donor_variant
    SpliceDonorVariant,

    /// "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript."
    /// SO:stop_gained, VEP:stop_gained
    StopGained,

    /// "A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three."
    /// SO:frameshift_variant, VEP:frameshift_variant
    FrameshiftVariant,

    /// "A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript."
    /// SO:stop_lost, VEP:stop_lost
    StopLost,

    /// "A codon variant that changes at least one base of the canonical start codon."
    /// SO:start_lost, VEP:start_lost
    StartLost,

    /// "A feature amplification of a region containing a transcript."
    /// SO:transcript_amplification, VEP:transcript_amplification
    TranscriptAmplification,

    // Currently never written out (because hgvs::parser::ProteinEdit::Ext not produced)
    /// "A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence."
    /// SO:feature_elongation, VEP:feature_elongation
    FeatureElongation,

    /// "A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence."
    /// SO:feature_truncation, VEP:feature_truncation
    FeatureTruncation,

    // moderate impact
    /// "An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon."
    /// SO:disruptive_inframe_insertion, VEP:inframe_insertion
    DisruptiveInframeInsertion,

    /// "An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon."
    /// SO:disruptive_inframe_deletion, VEP:inframe_deletion
    DisruptiveInframeDeletion,

    /// "An inframe increase in cds length that inserts one or more codons into the coding sequence between existing codons."
    /// SO:conservative_inframe_insertion, VEP:inframe_insertion
    ConservativeInframeInsertion,

    /// "An inframe decrease in cds length that deletes one or more entire codons from the coding sequence but does not change any remaining codons."
    /// SO:conservative_inframe_deletion, VEP:inframe_deletion
    ConservativeInframeDeletion,

    /// "A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved."
    /// SO:missense_variant, VEP:missense_variant
    MissenseVariant,

    // Not used by mehari, but by VEP (we're usually more specific)
    // /// "A sequence_variant which is predicted to change the protein encoded in the coding sequence."
    // /// SO:protein_altering_variant, VEP:missense_variant
    // ProteinAlteringVariant,

    // low impact
    /// "A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript."
    /// SO:splice_donor_5th_base_variant, VEP:splice_donor_5th_base_variant
    SpliceDonorFifthBaseVariant,

    /// "A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron."
    /// SO:splice_region_variant, VEP:splice_region_variant
    SpliceRegionVariant,

    /// "A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron)."
    /// SO:splice_donor_region_variant, VEP:splice_donor_region_variant
    SpliceDonorRegionVariant,

    /// "A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17)."
    /// SO:splice_polypyrimidine_tract_variant, VEP:splice_polypyrimidine_tract_variant
    SplicePolypyrimidineTractVariant,

    // Not used by mehari, but by VEP
    // /// "A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed."
    // /// SO:incomplete_terminal_codon_variant, VEP:incomplete_terminal_codon_variant
    // IncompleteTerminalCodonVariant
    /// "A sequence variant where at least one base in the start codon is changed, but the start remains."
    /// SO:start_retained_variant, VEP:start_retained_variant
    StartRetainedVariant,

    /// "A sequence variant where at least one base in the terminator codon is changed, but the terminator remains."
    /// SO:stop_retained_variant, VEP:stop_retained_variant
    StopRetainedVariant,

    /// "A sequence variant where there is no resulting change to the encoded amino acid."
    /// SO:synonymous_variant, VEP:synonymous_variant
    SynonymousVariant,

    // modifier
    /// "A sequence variant that changes the coding sequence."
    /// SO:coding_sequence_variant, VEP:coding_sequence_variant
    CodingSequenceVariant,

    // Not yet implemented
    /// "A transcript variant located with the sequence of the mature miRNA."
    /// SO:mature_miRNA_variant, VEP:mature_miRNA_variant
    MatureMirnaVariant,

    /// "A UTR variant of exonic sequence of the 5' UTR."
    /// SO:5_prime_UTR_exon_variant, VEP:5_prime_UTR_variant
    FivePrimeUtrExonVariant,

    /// "A UTR variant of intronic sequence of the 5' UTR."
    /// SO:5_prime_UTR_intron_variant, VEP:5_prime_UTR_variant
    FivePrimeUtrIntronVariant,

    /// "A UTR variant of exonic sequence of the 3' UTR."
    /// SO:3_prime_UTR_exon_variant, VEP:3_prime_UTR_variant
    ThreePrimeUtrExonVariant,

    /// "A UTR variant of intronic sequence of the 3' UTR."
    /// SO:3_prime_UTR_intron_variant, VEP:3_prime_UTR_variant
    ThreePrimeUtrIntronVariant,

    /// "A sequence variant that changes non-coding exon sequence in a non-coding transcript."
    /// SO:non_coding_transcript_exon_variant, VEP:non_coding_transcript_variant
    NonCodingTranscriptExonVariant,

    /// "A sequence variant that changes non-coding intron sequence in a non-coding transcript."
    /// SO:non_coding_transcript_intron_variant, VEP:non_coding_transcript_variant
    NonCodingTranscriptIntronVariant,

    // Not used by mehari, but by VEP
    // /// "A transcript variant of a protein coding gene."
    // /// SO:coding_transcript_variant, VEP:coding_transcript_variant
    // CodingTranscriptVariant,
    /// "A sequence variant located 5' of a gene."
    /// SO:upstream_gene_variant, VEP:upstream_gene_variant
    UpstreamGeneVariant,

    /// "A sequence variant located 3' of a gene."
    /// SO:downstream_gene_variant, VEP:downstream_gene_variant
    DownstreamGeneVariant,

    /// "A feature ablation whereby the deleted region includes a transcription factor binding site."
    /// SO:TFBS_ablation, VEP:TFBS_ablation
    TfbsAblation,

    /// "A feature amplification of a region containing a transcription factor binding site."
    /// SO:TFBS_amplification, VEP:TFBS_amplification
    TfbsAmplification,

    /// "A sequence variant located within a transcription factor binding site."
    /// SO:TF_binding_site_variant, VEP:TF_binding_site_variant
    TfBindingSiteVariant,

    /// "A feature ablation whereby the deleted region includes a regulatory region."
    /// SO:regulatory_region_ablation, VEP:regulatory_region_ablation
    RegulatoryRegionAblation,

    /// "A feature amplification of a region containing a regulatory region."
    /// SO:regulatory_region_amplification, VEP:regulatory_region_amplification
    RegulatoryRegionAmplification,

    /// "A sequence variant located within a regulatory region."
    /// SO:regulatory_region_variant, VEP:regulatory_region_variant
    RegulatoryRegionVariant,

    /// "A sequence variant located in the intergenic region, between genes."
    /// SO:intergenic_variant, VEP:intergenic_variant
    IntergenicVariant,

    // Not used by mehari, but by VEP
    // /// "A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration."
    // /// SO:sequence_variant, VEP:sequence_variant
    // SequenceVariant,
    IntronVariant,

    /// "A sequence variant where the structure of the gene is changed."
    /// SO:gene_variant
    GeneVariant,
}

/// Supporting code for `Consequence`.
pub(crate) mod consequence {
    /// Error type for `Consequence::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf Consequence: {0:?}")]
        UnknownConsequenceValue(super::pb_query::Consequence),
    }
}

impl TryFrom<pb_query::Consequence> for Consequence {
    type Error = consequence::Error;

    fn try_from(value: pb_query::Consequence) -> Result<Self, Self::Error> {
        match value {
            pb_query::Consequence::TranscriptAblation => Ok(Consequence::TranscriptAblation),
            pb_query::Consequence::ExonLossVariant => Ok(Consequence::ExonLossVariant),
            pb_query::Consequence::SpliceAcceptorVariant => Ok(Consequence::SpliceAcceptorVariant),
            pb_query::Consequence::SpliceDonorVariant => Ok(Consequence::SpliceDonorVariant),
            pb_query::Consequence::StopGained => Ok(Consequence::StopGained),
            pb_query::Consequence::FrameshiftVariant => Ok(Consequence::FrameshiftVariant),
            pb_query::Consequence::StopLost => Ok(Consequence::StopLost),
            pb_query::Consequence::StartLost => Ok(Consequence::StartLost),
            pb_query::Consequence::TranscriptAmplification => {
                Ok(Consequence::TranscriptAmplification)
            }
            pb_query::Consequence::FeatureElongation => Ok(Consequence::FeatureElongation),
            pb_query::Consequence::FeatureTruncation => Ok(Consequence::FeatureTruncation),
            pb_query::Consequence::DisruptiveInframeInsertion => {
                Ok(Consequence::DisruptiveInframeInsertion)
            }
            pb_query::Consequence::DisruptiveInframeDeletion => {
                Ok(Consequence::DisruptiveInframeDeletion)
            }
            pb_query::Consequence::ConservativeInframeInsertion => {
                Ok(Consequence::ConservativeInframeInsertion)
            }
            pb_query::Consequence::ConservativeInframeDeletion => {
                Ok(Consequence::ConservativeInframeDeletion)
            }
            pb_query::Consequence::MissenseVariant => Ok(Consequence::MissenseVariant),
            pb_query::Consequence::SpliceDonorFifthBaseVariant => {
                Ok(Consequence::SpliceDonorFifthBaseVariant)
            }
            pb_query::Consequence::SpliceRegionVariant => Ok(Consequence::SpliceRegionVariant),
            pb_query::Consequence::SpliceDonorRegionVariant => {
                Ok(Consequence::SpliceDonorRegionVariant)
            }
            pb_query::Consequence::SplicePolypyrimidineTractVariant => {
                Ok(Consequence::SplicePolypyrimidineTractVariant)
            }
            pb_query::Consequence::StartRetainedVariant => Ok(Consequence::StartRetainedVariant),
            pb_query::Consequence::StopRetainedVariant => Ok(Consequence::StopRetainedVariant),
            pb_query::Consequence::SynonymousVariant => Ok(Consequence::SynonymousVariant),
            pb_query::Consequence::CodingSequenceVariant => Ok(Consequence::CodingSequenceVariant),
            pb_query::Consequence::MatureMirnaVariant => Ok(Consequence::MatureMirnaVariant),
            pb_query::Consequence::FivePrimeUtrExonVariant => {
                Ok(Consequence::FivePrimeUtrExonVariant)
            }
            pb_query::Consequence::FivePrimeUtrIntronVariant => {
                Ok(Consequence::FivePrimeUtrIntronVariant)
            }
            pb_query::Consequence::ThreePrimeUtrExonVariant => {
                Ok(Consequence::ThreePrimeUtrExonVariant)
            }
            pb_query::Consequence::ThreePrimeUtrIntronVariant => {
                Ok(Consequence::ThreePrimeUtrIntronVariant)
            }
            pb_query::Consequence::NonCodingTranscriptExonVariant => {
                Ok(Consequence::NonCodingTranscriptExonVariant)
            }
            pb_query::Consequence::NonCodingTranscriptIntronVariant => {
                Ok(Consequence::NonCodingTranscriptIntronVariant)
            }
            pb_query::Consequence::UpstreamGeneVariant => Ok(Consequence::UpstreamGeneVariant),
            pb_query::Consequence::DownstreamGeneVariant => Ok(Consequence::DownstreamGeneVariant),
            pb_query::Consequence::TfbsAblation => Ok(Consequence::TfbsAblation),
            pb_query::Consequence::TfbsAmplification => Ok(Consequence::TfbsAmplification),
            pb_query::Consequence::TfBindingSiteVariant => Ok(Consequence::TfBindingSiteVariant),
            pb_query::Consequence::RegulatoryRegionAblation => {
                Ok(Consequence::RegulatoryRegionAblation)
            }
            pb_query::Consequence::RegulatoryRegionAmplification => {
                Ok(Consequence::RegulatoryRegionAmplification)
            }
            pb_query::Consequence::RegulatoryRegionVariant => {
                Ok(Consequence::RegulatoryRegionVariant)
            }
            pb_query::Consequence::IntergenicVariant => Ok(Consequence::IntergenicVariant),
            pb_query::Consequence::IntronVariant => Ok(Consequence::IntronVariant),
            pb_query::Consequence::GeneVariant => Ok(Consequence::GeneVariant),
            _ => Err(consequence::Error::UnknownConsequenceValue(value)),
        }
    }
}

impl From<Consequence> for mehari::annotate::seqvars::ann::Consequence {
    fn from(val: Consequence) -> Self {
        use mehari::annotate::seqvars::ann;

        match val {
            Consequence::TranscriptAblation => ann::Consequence::TranscriptAblation,
            Consequence::ExonLossVariant => ann::Consequence::ExonLossVariant,
            Consequence::SpliceAcceptorVariant => ann::Consequence::SpliceAcceptorVariant,
            Consequence::SpliceDonorVariant => ann::Consequence::SpliceDonorVariant,
            Consequence::StopGained => ann::Consequence::StopGained,
            Consequence::FrameshiftVariant => ann::Consequence::FrameshiftVariant,
            Consequence::StopLost => ann::Consequence::StopLost,
            Consequence::StartLost => ann::Consequence::StartLost,
            Consequence::TranscriptAmplification => ann::Consequence::TranscriptAmplification,
            Consequence::FeatureElongation => ann::Consequence::FeatureElongation,
            Consequence::FeatureTruncation => ann::Consequence::FeatureTruncation,
            Consequence::DisruptiveInframeInsertion => ann::Consequence::DisruptiveInframeInsertion,
            Consequence::DisruptiveInframeDeletion => ann::Consequence::DisruptiveInframeDeletion,
            Consequence::ConservativeInframeInsertion => {
                ann::Consequence::ConservativeInframeInsertion
            }
            Consequence::ConservativeInframeDeletion => {
                ann::Consequence::ConservativeInframeDeletion
            }
            Consequence::MissenseVariant => ann::Consequence::MissenseVariant,
            Consequence::SpliceDonorFifthBaseVariant => {
                ann::Consequence::SpliceDonorFifthBaseVariant
            }
            Consequence::SpliceRegionVariant => ann::Consequence::SpliceRegionVariant,
            Consequence::SpliceDonorRegionVariant => ann::Consequence::SpliceDonorRegionVariant,
            Consequence::SplicePolypyrimidineTractVariant => {
                ann::Consequence::SplicePolypyrimidineTractVariant
            }
            Consequence::StartRetainedVariant => ann::Consequence::StartRetainedVariant,
            Consequence::StopRetainedVariant => ann::Consequence::StopRetainedVariant,
            Consequence::SynonymousVariant => ann::Consequence::SynonymousVariant,
            Consequence::CodingSequenceVariant => ann::Consequence::CodingSequenceVariant,
            Consequence::MatureMirnaVariant => ann::Consequence::MatureMirnaVariant,
            Consequence::FivePrimeUtrExonVariant => ann::Consequence::FivePrimeUtrExonVariant,
            Consequence::FivePrimeUtrIntronVariant => ann::Consequence::FivePrimeUtrIntronVariant,
            Consequence::ThreePrimeUtrExonVariant => ann::Consequence::ThreePrimeUtrExonVariant,
            Consequence::ThreePrimeUtrIntronVariant => ann::Consequence::ThreePrimeUtrIntronVariant,
            Consequence::NonCodingTranscriptExonVariant => {
                ann::Consequence::NonCodingTranscriptExonVariant
            }
            Consequence::NonCodingTranscriptIntronVariant => {
                ann::Consequence::NonCodingTranscriptIntronVariant
            }
            Consequence::UpstreamGeneVariant => ann::Consequence::UpstreamGeneVariant,
            Consequence::DownstreamGeneVariant => ann::Consequence::DownstreamGeneVariant,
            Consequence::TfbsAblation => ann::Consequence::TfbsAblation,
            Consequence::TfbsAmplification => ann::Consequence::TfbsAmplification,
            Consequence::TfBindingSiteVariant => ann::Consequence::TfBindingSiteVariant,
            Consequence::RegulatoryRegionAblation => ann::Consequence::RegulatoryRegionAblation,
            Consequence::RegulatoryRegionAmplification => {
                ann::Consequence::RegulatoryRegionAmplification
            }
            Consequence::RegulatoryRegionVariant => ann::Consequence::RegulatoryRegionVariant,
            Consequence::IntergenicVariant => ann::Consequence::IntergenicVariant,
            Consequence::IntronVariant => ann::Consequence::IntronVariant,
            Consequence::GeneVariant => ann::Consequence::GeneVariant,
        }
    }
}

/// Query settings for consequence types.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsConsequence {
    /// Variant types.
    pub variant_types: Vec<VariantType>,
    /// Transcript types.
    pub transcript_types: Vec<TranscriptType>,
    /// Consequences to consider.
    pub consequences: Vec<Consequence>,
    /// Maximal distance to next exon.
    pub max_dist_to_exon: Option<i32>,
}

/// Supporting code for `QuerySettingsConsequence`.
pub(crate) mod query_settings_consequence {
    /// Error type for `QuerySettingsConsequence::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert i32 into protobuf VariantType: {0}")]
        VariantTypeInt(i32),
        #[error("Cannot convert protobuf VariantType: {0:?}")]
        VariantTypeValue(super::pb_query::VariantType),
        #[error("Cannot convert i32 into protobuf TranscriptType: {0}")]
        TranscriptTypeInt(i32),
        #[error("Cannot convert protobuf TranscriptType: {0:?}")]
        TranscriptTypeValue(super::pb_query::TranscriptType),
        #[error("Cannot convert i32 into protobuf Consequence: {0}")]
        ConsequenceInt(i32),
        #[error("Cannot convert protobuf Consequence: {0:?}")]
        ConsequenceValue(super::pb_query::Consequence),
    }
}

impl TryFrom<pb_query::QuerySettingsConsequence> for QuerySettingsConsequence {
    type Error = query_settings_consequence::Error;

    fn try_from(value: pb_query::QuerySettingsConsequence) -> Result<Self, Self::Error> {
        let variant_types = value
            .variant_types
            .into_iter()
            .map(|v| {
                let v = pb_query::VariantType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::VariantTypeInt(v))?;
                VariantType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::VariantTypeValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let transcript_types = value
            .transcript_types
            .into_iter()
            .map(|v| {
                let v = pb_query::TranscriptType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::TranscriptTypeInt(v))?;
                TranscriptType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::TranscriptTypeValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let consequences = value
            .consequences
            .into_iter()
            .map(|v| {
                let v = pb_query::Consequence::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::ConsequenceInt(v))?;
                Consequence::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::ConsequenceValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            variant_types,
            transcript_types,
            consequences,
            max_dist_to_exon: value.max_dist_to_exon,
        })
    }
}

/// A 1-based integer range.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Range {
    /// 1-based start position.
    pub start: i32,
    /// 1-based end position.
    pub stop: i32,
}

impl From<pb_query::Range> for Range {
    fn from(value: pb_query::Range) -> Self {
        Self {
            start: value.start,
            stop: value.stop,
        }
    }
}

/// Genomic region.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct GenomicRegion {
    /// Chromosome.
    pub chrom: String,
    /// Range of region.
    pub range: Option<Range>,
}

impl From<pb_query::GenomicRegion> for GenomicRegion {
    fn from(value: pb_query::GenomicRegion) -> Self {
        Self {
            chrom: value.chrom,
            range: value.range.map(Range::from),
        }
    }
}

/// Query settings for locus.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsLocus {
    /// List of HGNC identifiers for filtration to genes.
    pub genes: Vec<String>,
    /// List of genomic regions to limit restrict the resulting variants to.
    pub genome_regions: Vec<GenomicRegion>,
}

impl From<pb_query::QuerySettingsLocus> for QuerySettingsLocus {
    fn from(value: pb_query::QuerySettingsLocus) -> Self {
        Self {
            genes: value.genes,
            genome_regions: value
                .genome_regions
                .into_iter()
                .map(GenomicRegion::from)
                .collect(),
        }
    }
}

// Canonical ClinVar germline aggregate descriptions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ClinvarGermlineAggregateDescription {
    /// Pathogenic.
    Pathogenic,
    /// Likely pathogenic.
    LikelyPathogenic,
    /// Uncertain significance.
    UncertainSignificance,
    /// Likely benign.
    LikelyBenign,
    /// Benign.
    Benign,
}

/// Supporting code for `ClinvarGermlineAggregateDescription`.
pub(crate) mod clinvar_germline_aggregate_description {
    /// Error type for `ClinvarGermlineAggregateDescription::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf ClinvarGermlineAggregateDescription: {0:?}")]
        UnknownClinvarGermlineAggregateDescriptionValue(
            super::pb_query::ClinvarGermlineAggregateDescription,
        ),
    }
}

impl TryFrom<pb_query::ClinvarGermlineAggregateDescription>
    for ClinvarGermlineAggregateDescription
{
    type Error = clinvar_germline_aggregate_description::Error;

    fn try_from(value: pb_query::ClinvarGermlineAggregateDescription) -> Result<Self, Self::Error> {
        match value {
            pb_query::ClinvarGermlineAggregateDescription::Pathogenic => Ok(ClinvarGermlineAggregateDescription::Pathogenic),
            pb_query::ClinvarGermlineAggregateDescription::LikelyPathogenic => Ok(ClinvarGermlineAggregateDescription::LikelyPathogenic),
            pb_query::ClinvarGermlineAggregateDescription::UncertainSignificance => Ok(ClinvarGermlineAggregateDescription::UncertainSignificance),
            pb_query::ClinvarGermlineAggregateDescription::LikelyBenign => Ok(ClinvarGermlineAggregateDescription::LikelyBenign),
            pb_query::ClinvarGermlineAggregateDescription::Benign => Ok(ClinvarGermlineAggregateDescription::Benign),
            _ => Err(clinvar_germline_aggregate_description::Error::UnknownClinvarGermlineAggregateDescriptionValue(value)),
        }
    }
}

/// Query settings for ClinVar.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsClinVar {
    /// Whether to require ClinVar membership.
    pub presence_required: bool,
    /// The ClinVar germline aggregate description to include.
    pub germline_descriptions: Vec<ClinvarGermlineAggregateDescription>,
    /// Whether to include conflicting interpretation ClinVar variants.
    pub allow_conflicting_interpretations: bool,
}

/// Supporting code for `QuerySettingsClinVar`.
pub(crate) mod query_settings_clinvar {
    /// Error type for `QuerySettingsClinVar::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert i32 into protobuf ClinvarGermlineAggregateDescription: {0}")]
        UnknownClinvarGermlineAggregateDescriptionInt(i32),
        #[error("Cannot convert protobuf ClinvarGermlineAggregateDescription: {0:?}")]
        UnknownClinvarGermlineAggregateDescriptionValue(
            super::pb_query::ClinvarGermlineAggregateDescription,
        ),
    }
}

impl TryFrom<pb_query::QuerySettingsClinVar> for QuerySettingsClinVar {
    type Error = query_settings_clinvar::Error;

    fn try_from(value: pb_query::QuerySettingsClinVar) -> Result<Self, Self::Error> {
        let germline_descriptions = value
            .germline_descriptions
            .into_iter()
            .map(|v| {
                pb_query::ClinvarGermlineAggregateDescription::try_from(v)
                    .map_err(|_| Self::Error::UnknownClinvarGermlineAggregateDescriptionInt(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let germline_descriptions = germline_descriptions
            .into_iter()
            .map(|v| {
                ClinvarGermlineAggregateDescription::try_from(v)
                    .map_err(|_| Self::Error::UnknownClinvarGermlineAggregateDescriptionValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            presence_required: value.presence_required,
            germline_descriptions,
            allow_conflicting_interpretations: value.allow_conflicting_interpretations,
        })
    }
}

/// Query settings for one case.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct CaseQuery {
    /// Genotype query settings.
    pub genotype: QuerySettingsGenotype,
    /// Quality query settings.
    pub quality: QuerySettingsQuality,
    /// Frequency query settings.
    pub frequency: QuerySettingsFrequency,
    /// Consequence query settings.
    pub consequence: QuerySettingsConsequence,
    /// Locus query settings.
    pub locus: QuerySettingsLocus,
    /// ClinVar query settings.
    pub clinvar: QuerySettingsClinVar,
}

/// Supporting code for `CaseQuery`.
pub(crate) mod case_query {
    /// Error type for `CaseQuery::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Problem converting protobuf for genotype: {0}")]
        Genotype(#[from] super::query_settings_genotype::Error),
        #[error("Problem converting protobuf for quality: {0}")]
        Quality(#[from] super::query_settings_quality::Error),
        #[error("Problem converting protobuf for consequence: {0}")]
        Consequence(#[from] super::query_settings_consequence::Error),
        #[error("Problem converting protobuf for clinvar: {0}")]
        Clinvar(#[from] super::query_settings_clinvar::Error),
    }
}

impl TryFrom<pb_query::CaseQuery> for CaseQuery {
    type Error = case_query::Error;

    fn try_from(value: pb_query::CaseQuery) -> Result<Self, Self::Error> {
        let pb_query::CaseQuery {
            genotype,
            quality,
            frequency,
            consequence,
            locus,
            clinvar,
        } = value;

        let genotype = QuerySettingsGenotype::try_from(genotype.unwrap_or(Default::default()))
            .map_err(Self::Error::Genotype)?;
        let quality = QuerySettingsQuality::try_from(quality.unwrap_or(Default::default()))
            .map_err(Self::Error::Quality)?;
        let frequency = QuerySettingsFrequency::from(frequency.unwrap_or(Default::default()));
        let consequence =
            QuerySettingsConsequence::try_from(consequence.unwrap_or(Default::default()))
                .map_err(Self::Error::Consequence)?;
        let locus = QuerySettingsLocus::from(locus.unwrap_or(Default::default()));
        let clinvar = QuerySettingsClinVar::try_from(clinvar.unwrap_or(Default::default()))
            .map_err(Self::Error::Clinvar)?;

        Ok(Self {
            genotype,
            quality,
            frequency,
            consequence,
            locus,
            clinvar,
        })
    }
}

#[cfg(test)]
mod tests {
    use query_settings_genotype::RecessiveIndexError;

    use super::*;

    #[rstest::rstest]
    #[case::dot(".", true)]
    #[case::dot_slash_dot("./.", true)]
    #[case::dot_pipe_dot(".|.", true)]
    #[case::zero("0", false)]
    #[case::zero_slash_dot("0/.", true)]
    #[case::zero_pipe_dot("0|.", true)]
    #[case::dot_slash_zero("./0", true)]
    #[case::dot_pipe_zero(".|0", true)]
    #[case::one("1", false)]
    #[case::one_slash_dot("1/.", true)]
    #[case::one_pipe_dot("1|.", true)]
    #[case::dot_slash_one("./1", true)]
    #[case::dot_pipe_one(".|1", true)]
    #[case::zero_slash_one("0/1", false)]
    #[case::zero_pipe_one("0|1", false)]
    #[case::one_slash_zero("1/0", false)]
    #[case::one_pipe_zero("1|0", false)]
    fn test_considered_no_call(#[case] gt: &str, #[case] expected: bool) {
        assert_eq!(super::considered_no_call(gt), expected);
    }

    #[test]
    fn test_recessive_mode_try_from() {
        assert_eq!(
            RecessiveMode::try_from(pb_query::RecessiveMode::Disabled).unwrap(),
            RecessiveMode::Disabled
        );
        assert_eq!(
            RecessiveMode::try_from(pb_query::RecessiveMode::CompoundHeterozygous).unwrap(),
            RecessiveMode::CompoundHeterozygous
        );
        assert_eq!(
            RecessiveMode::try_from(pb_query::RecessiveMode::Homozygous).unwrap(),
            RecessiveMode::Homozygous
        );
        assert_eq!(
            RecessiveMode::try_from(pb_query::RecessiveMode::Any).unwrap(),
            RecessiveMode::Any
        );
        assert!(RecessiveMode::try_from(pb_query::RecessiveMode::Unspecified).is_err());
    }

    #[test]
    fn test_genotype_choice_try_from() {
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::Any).unwrap(),
            GenotypeChoice::Any
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::Ref).unwrap(),
            GenotypeChoice::Ref
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::Het).unwrap(),
            GenotypeChoice::Het
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::Hom).unwrap(),
            GenotypeChoice::Hom
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::NonHet).unwrap(),
            GenotypeChoice::NonHet
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::NonHom).unwrap(),
            GenotypeChoice::NonHom
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::Variant).unwrap(),
            GenotypeChoice::Variant
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::RecessiveIndex).unwrap(),
            GenotypeChoice::RecessiveIndex
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::RecessiveFather).unwrap(),
            GenotypeChoice::RecessiveFather
        );
        assert_eq!(
            GenotypeChoice::try_from(pb_query::GenotypeChoice::RecessiveMother).unwrap(),
            GenotypeChoice::RecessiveMother
        );
        assert!(GenotypeChoice::try_from(pb_query::GenotypeChoice::Unspecified).is_err());
    }

    #[test]
    fn test_sample_genotype_choice_try_from() {
        let pb_sample_genotype_choice = pb_query::SampleGenotypeChoice {
            sample: "sample".to_string(),
            genotype: pb_query::GenotypeChoice::Any as i32,
            include_no_call: true,
            enabled: true,
        };
        let sample_genotype_choice = SampleGenotypeChoice {
            sample: "sample".to_string(),
            genotype: GenotypeChoice::Any,
            include_no_call: true,
            enabled: true,
        };
        assert_eq!(
            SampleGenotypeChoice::try_from(pb_sample_genotype_choice).unwrap(),
            sample_genotype_choice
        );
    }

    pub mod genotype_choice {
        use super::GenotypeChoice::{self, *};

        #[rstest::rstest]
        #[case::any_pass(Any, &["0", "0/0", "0|0", "1", "1/1", "1|1", "0/1", "0|1", "1/0", "1|0"], true)]
        #[case::any_fail(Any, &["0/.", ".", "./.", ".|."], false)]
        #[case::ref_pass(Ref, &["0", "0/0", "0|0"], true)]
        #[case::ref_fail(Ref, &[ "1", "1/1", "1|1", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], false)]
        #[case::het_pass(Het, &["0/1", "0|1", "1/0", "1|0"], true)]
        #[case::het_fail(Het, &["0", "0/0", "0|0", "1", "1/1", "1|1", ".", "./.", ".|."], false)]
        #[case::hom_pass(Hom, &["1", "1/1", "1|1"], true)]
        #[case::hom_fail(Hom, &["0", "0/0", "0|0", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], false)]
        #[case::nonhom_pass(NonHom, &["0", "0/0", "0|0", "0/1", "0|1", "1/0", "1|0"], true)]
        #[case::nonhom_fail(NonHom, &["1", "1/1", "1|1", "0/.", ".", "./.", ".|."], false)]
        #[case::variant_pass(Variant, &["1", "1/1", "1|1", "0/1",  "0|1",  "1/0", "1|0"], true)]
        #[case::variant_fail(Variant, &["0", "0/0", "0|0", ".", "./.", ".|."], false)]
        pub fn matches(
            #[case] genotype_choice: GenotypeChoice,
            #[case] gt_strs: &[&str],
            #[case] expected: bool,
        ) -> Result<(), anyhow::Error> {
            use crate::seqvars::query::schema::query::MatchesGenotypeStr as _;

            for gt_str in gt_strs {
                assert_eq!(
                    genotype_choice.matches(gt_str)?,
                    expected,
                    "{:?} {}",
                    genotype_choice,
                    gt_str
                );
            }

            Ok(())
        }
    }

    #[test]
    fn test_query_settings_genotype_recessive_index_none() {
        let query_settings_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::CompoundHeterozygous,
            sample_genotypes: Default::default(),
        };

        assert_eq!(
            query_settings_genotype.recessive_index(),
            Err(RecessiveIndexError::NoRecessiveIndexSample)
        );
    }

    #[test]
    fn test_query_settings_genotype_recessive_index_single() {
        let query_settings_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::CompoundHeterozygous,
            sample_genotypes: indexmap::indexmap! {
                String::from("sample") => SampleGenotypeChoice {
                    sample: "sample".to_string(),
                    genotype: GenotypeChoice::RecessiveIndex,
                    include_no_call: true,
                    enabled: true,
                }
            },
        };

        assert_eq!(
            query_settings_genotype.recessive_index(),
            Ok(String::from("sample"))
        );
    }

    #[test]
    fn test_query_settings_genotype_recessive_index_two() {
        let query_settings_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::CompoundHeterozygous,
            sample_genotypes: indexmap::indexmap! {
                String::from("sample") => SampleGenotypeChoice {
                    sample: "sample".to_string(),
                    genotype: GenotypeChoice::RecessiveIndex,
                    include_no_call: true,
                    enabled: true,
                },
                String::from("sample2") => SampleGenotypeChoice {
                    sample: "sample2".to_string(),
                    genotype: GenotypeChoice::RecessiveIndex,
                    include_no_call: true,
                    enabled: true,
                }
            },
        };

        assert_eq!(
            query_settings_genotype.recessive_index(),
            Err(RecessiveIndexError::MultipleRecessiveIndexSamples(vec![
                String::from("sample"),
                String::from("sample2")
            ]))
        );
    }

    #[test]
    fn test_query_settings_genotype_try_from() {
        let pb_query_settings_genotype = pb_query::QuerySettingsGenotype {
            recessive_mode: pb_query::RecessiveMode::Disabled as i32,
            sample_genotypes: vec![pb_query::SampleGenotypeChoice {
                sample: "sample".to_string(),
                genotype: pb_query::GenotypeChoice::Any as i32,
                include_no_call: true,
                enabled: true,
            }],
        };
        let query_settings_genotype = QuerySettingsGenotype {
            recessive_mode: RecessiveMode::Disabled,
            sample_genotypes: {
                let mut map = indexmap::IndexMap::new();
                map.insert(
                    "sample".to_string(),
                    SampleGenotypeChoice {
                        sample: "sample".to_string(),
                        genotype: GenotypeChoice::Any,
                        include_no_call: true,
                        enabled: true,
                    },
                );
                map
            },
        };
        assert_eq!(
            QuerySettingsGenotype::try_from(pb_query_settings_genotype).unwrap(),
            query_settings_genotype
        );
    }

    #[test]
    fn test_sample_quality_settings_from() {
        let pb_sample_quality_settings = pb_query::SampleQualitySettings {
            sample: "sample".to_string(),
            filter_active: true,
            min_dp_het: Some(10),
            min_dp_hom: Some(20),
            min_gq: Some(30),
            min_ab: Some(0.1),
            min_ad: Some(40),
            max_ad: Some(50),
        };
        let sample_quality_settings = SampleQualitySettings {
            sample: "sample".to_string(),
            filter_active: true,
            min_dp_het: Some(10),
            min_dp_hom: Some(20),
            min_gq: Some(30),
            min_ab: Some(0.1),
            min_ad: Some(40),
            max_ad: Some(50),
        };
        assert_eq!(
            SampleQualitySettings::from(pb_sample_quality_settings),
            sample_quality_settings
        );
    }

    #[test]
    fn test_query_settings_quality_try_from() {
        let pb_query_settings_quality = pb_query::QuerySettingsQuality {
            sample_qualities: vec![pb_query::SampleQualitySettings {
                sample: "sample".to_string(),
                filter_active: true,
                min_dp_het: Some(10),
                min_dp_hom: Some(20),
                min_gq: Some(30),
                min_ab: Some(0.1),
                min_ad: Some(40),
                max_ad: Some(50),
            }],
        };
        let query_settings_quality = QuerySettingsQuality {
            sample_qualities: {
                let mut map = indexmap::IndexMap::new();
                map.insert(
                    "sample".to_string(),
                    SampleQualitySettings {
                        sample: "sample".to_string(),
                        filter_active: true,
                        min_dp_het: Some(10),
                        min_dp_hom: Some(20),
                        min_gq: Some(30),
                        min_ab: Some(0.1),
                        min_ad: Some(40),
                        max_ad: Some(50),
                    },
                );
                map
            },
        };
        assert_eq!(
            QuerySettingsQuality::try_from(pb_query_settings_quality).unwrap(),
            query_settings_quality
        );
    }

    #[test]
    fn test_gnomad_nuclear_frequency_settings_from() {
        let pb_gnomad_nuclear_frequency_settings = pb_query::NuclearFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_hemi: Some(30),
            max_af: Some(0.1),
        };
        let gnomad_nuclear_frequency_settings = NuclearFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_hemi: Some(30),
            max_af: Some(0.1),
        };
        assert_eq!(
            NuclearFrequencySettings::from(pb_gnomad_nuclear_frequency_settings),
            gnomad_nuclear_frequency_settings
        );
    }

    #[test]
    fn test_gnomad_mitochondrial_frequency_settings_from() {
        let pb_gnomad_mitochondrial_frequency_settings = pb_query::MitochondrialFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_af: Some(0.1),
        };
        let gnomad_mitochondrial_frequency_settings = MitochondrialFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_af: Some(0.1),
        };
        assert_eq!(
            MitochondrialFrequencySettings::from(pb_gnomad_mitochondrial_frequency_settings),
            gnomad_mitochondrial_frequency_settings
        );
    }

    #[test]
    fn test_helix_mtdb_frequency_settings_from() {
        let pb_helix_mtdb_frequency_settings = pb_query::MitochondrialFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_af: Some(0.1),
        };
        let helix_mtdb_frequency_settings = MitochondrialFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_af: Some(0.1),
        };
        assert_eq!(
            MitochondrialFrequencySettings::from(pb_helix_mtdb_frequency_settings),
            helix_mtdb_frequency_settings
        );
    }

    #[test]
    fn test_inhouse_frequency_settings_from() {
        let pb_inhouse_frequency_settings = pb_query::NuclearFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_hemi: Some(30),
            max_af: Some(0.1),
        };
        let inhouse_frequency_settings = NuclearFrequencySettings {
            enabled: true,
            max_het: Some(10),
            max_hom: Some(20),
            max_hemi: Some(30),
            max_af: Some(0.1),
        };
        assert_eq!(
            NuclearFrequencySettings::from(pb_inhouse_frequency_settings),
            inhouse_frequency_settings
        );
    }

    #[test]
    fn test_query_settings_frequency_try_from() {
        let pb_query_settings_frequency = pb_query::QuerySettingsFrequency {
            gnomad_exomes: Some(pb_query::NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            }),
            gnomad_genomes: Some(pb_query::NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            }),
            gnomad_mtdna: Some(pb_query::MitochondrialFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_af: Some(0.1),
            }),
            helixmtdb: Some(pb_query::MitochondrialFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_af: Some(0.1),
            }),
            inhouse: Some(pb_query::NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            }),
        };
        let query_settings_frequency = QuerySettingsFrequency {
            gnomad_exomes: NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            },
            gnomad_genomes: NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            },
            gnomad_mtdna: MitochondrialFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_af: Some(0.1),
            },
            helixmtdb: MitochondrialFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_af: Some(0.1),
            },
            inhouse: NuclearFrequencySettings {
                enabled: true,
                max_het: Some(10),
                max_hom: Some(20),
                max_hemi: Some(30),
                max_af: Some(0.1),
            },
        };
        assert_eq!(
            QuerySettingsFrequency::from(pb_query_settings_frequency),
            query_settings_frequency
        );
    }

    #[test]
    fn test_variant_type_try_from() {
        assert_eq!(
            VariantType::try_from(pb_query::VariantType::Snv).unwrap(),
            VariantType::Snv
        );
        assert_eq!(
            VariantType::try_from(pb_query::VariantType::Indel).unwrap(),
            VariantType::Indel
        );
        assert_eq!(
            VariantType::try_from(pb_query::VariantType::Mnv).unwrap(),
            VariantType::Mnv
        );
        assert_eq!(
            VariantType::try_from(pb_query::VariantType::ComplexSubstitution).unwrap(),
            VariantType::ComplexSubstitution
        );
        assert!(VariantType::try_from(pb_query::VariantType::Unspecified).is_err());
    }

    #[test]
    fn test_transcript_type_try_from() {
        assert_eq!(
            TranscriptType::try_from(pb_query::TranscriptType::Coding).unwrap(),
            TranscriptType::Coding
        );
        assert_eq!(
            TranscriptType::try_from(pb_query::TranscriptType::NonCoding).unwrap(),
            TranscriptType::NonCoding
        );
        assert!(TranscriptType::try_from(pb_query::TranscriptType::Unspecified).is_err());
    }

    #[test]
    fn test_consequence_try_from() {
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::TranscriptAblation).unwrap(),
            Consequence::TranscriptAblation
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::ExonLossVariant).unwrap(),
            Consequence::ExonLossVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SpliceAcceptorVariant).unwrap(),
            Consequence::SpliceAcceptorVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SpliceDonorVariant).unwrap(),
            Consequence::SpliceDonorVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::StopGained).unwrap(),
            Consequence::StopGained
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::FrameshiftVariant).unwrap(),
            Consequence::FrameshiftVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::StopLost).unwrap(),
            Consequence::StopLost
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::StartLost).unwrap(),
            Consequence::StartLost
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::TranscriptAmplification).unwrap(),
            Consequence::TranscriptAmplification
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::DisruptiveInframeInsertion).unwrap(),
            Consequence::DisruptiveInframeInsertion
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::DisruptiveInframeDeletion).unwrap(),
            Consequence::DisruptiveInframeDeletion
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::ConservativeInframeInsertion).unwrap(),
            Consequence::ConservativeInframeInsertion
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::ConservativeInframeDeletion).unwrap(),
            Consequence::ConservativeInframeDeletion
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::MissenseVariant).unwrap(),
            Consequence::MissenseVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SpliceDonorFifthBaseVariant).unwrap(),
            Consequence::SpliceDonorFifthBaseVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SpliceRegionVariant).unwrap(),
            Consequence::SpliceRegionVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SpliceDonorRegionVariant).unwrap(),
            Consequence::SpliceDonorRegionVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SplicePolypyrimidineTractVariant).unwrap(),
            Consequence::SplicePolypyrimidineTractVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::StartRetainedVariant).unwrap(),
            Consequence::StartRetainedVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::StopRetainedVariant).unwrap(),
            Consequence::StopRetainedVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::SynonymousVariant).unwrap(),
            Consequence::SynonymousVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::CodingSequenceVariant).unwrap(),
            Consequence::CodingSequenceVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::FivePrimeUtrExonVariant).unwrap(),
            Consequence::FivePrimeUtrExonVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::FivePrimeUtrIntronVariant).unwrap(),
            Consequence::FivePrimeUtrIntronVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::ThreePrimeUtrExonVariant).unwrap(),
            Consequence::ThreePrimeUtrExonVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::ThreePrimeUtrIntronVariant).unwrap(),
            Consequence::ThreePrimeUtrIntronVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::NonCodingTranscriptExonVariant).unwrap(),
            Consequence::NonCodingTranscriptExonVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::NonCodingTranscriptIntronVariant).unwrap(),
            Consequence::NonCodingTranscriptIntronVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::UpstreamGeneVariant).unwrap(),
            Consequence::UpstreamGeneVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::DownstreamGeneVariant).unwrap(),
            Consequence::DownstreamGeneVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::IntergenicVariant).unwrap(),
            Consequence::IntergenicVariant
        );
        assert_eq!(
            Consequence::try_from(pb_query::Consequence::IntronVariant).unwrap(),
            Consequence::IntronVariant
        );
        assert!(Consequence::try_from(pb_query::Consequence::Unspecified).is_err());
    }

    #[test]
    fn test_query_settings_consequence_try_from() {
        let pb_query_settings_consequence = pb_query::QuerySettingsConsequence {
            variant_types: vec![
                pb_query::VariantType::Snv as i32,
                pb_query::VariantType::Indel as i32,
                pb_query::VariantType::Mnv as i32,
                pb_query::VariantType::ComplexSubstitution as i32,
            ],
            transcript_types: vec![
                pb_query::TranscriptType::Coding as i32,
                pb_query::TranscriptType::NonCoding as i32,
            ],
            consequences: vec![
                pb_query::Consequence::TranscriptAblation as i32,
                pb_query::Consequence::ExonLossVariant as i32,
                pb_query::Consequence::SpliceAcceptorVariant as i32,
                pb_query::Consequence::SpliceDonorVariant as i32,
                pb_query::Consequence::StopGained as i32,
                pb_query::Consequence::FrameshiftVariant as i32,
                pb_query::Consequence::StopLost as i32,
                pb_query::Consequence::StartLost as i32,
            ],
            max_dist_to_exon: Some(10),
        };
        let query_settings_consequence = QuerySettingsConsequence {
            variant_types: vec![
                VariantType::Snv,
                VariantType::Indel,
                VariantType::Mnv,
                VariantType::ComplexSubstitution,
            ],
            transcript_types: vec![TranscriptType::Coding, TranscriptType::NonCoding],
            consequences: vec![
                Consequence::TranscriptAblation,
                Consequence::ExonLossVariant,
                Consequence::SpliceAcceptorVariant,
                Consequence::SpliceDonorVariant,
                Consequence::StopGained,
                Consequence::FrameshiftVariant,
                Consequence::StopLost,
                Consequence::StartLost,
            ],
            max_dist_to_exon: Some(10),
        };
        assert_eq!(
            QuerySettingsConsequence::try_from(pb_query_settings_consequence).unwrap(),
            query_settings_consequence
        );
    }

    #[test]
    fn test_range_from() {
        let pb_range = pb_query::Range { start: 1, stop: 2 };
        let range = Range { start: 1, stop: 2 };
        assert_eq!(Range::from(pb_range), range);
    }

    #[test]
    fn test_genomic_region_from() {
        let pb_genomic_region = pb_query::GenomicRegion {
            chrom: "chrom".to_string(),
            range: Some(pb_query::Range { start: 1, stop: 2 }),
        };
        let genomic_region = GenomicRegion {
            chrom: "chrom".to_string(),
            range: Some(Range { start: 1, stop: 2 }),
        };
        assert_eq!(GenomicRegion::from(pb_genomic_region), genomic_region);
    }

    #[test]
    fn test_query_settings_locus_from() {
        let pb_query_settings_locus = pb_query::QuerySettingsLocus {
            genes: vec!["gene".to_string()],
            genome_regions: vec![pb_query::GenomicRegion {
                chrom: "chrom".to_string(),
                range: Some(pb_query::Range { start: 1, stop: 2 }),
            }],
        };
        let query_settings_locus = QuerySettingsLocus {
            genes: vec!["gene".to_string()],
            genome_regions: vec![GenomicRegion {
                chrom: "chrom".to_string(),
                range: Some(Range { start: 1, stop: 2 }),
            }],
        };
        assert_eq!(
            QuerySettingsLocus::from(pb_query_settings_locus),
            query_settings_locus
        );
    }

    #[test]
    fn test_clinvar_germline_aggregate_description_try_from() {
        assert_eq!(
            ClinvarGermlineAggregateDescription::try_from(
                pb_query::ClinvarGermlineAggregateDescription::Pathogenic
            )
            .unwrap(),
            ClinvarGermlineAggregateDescription::Pathogenic
        );
        assert_eq!(
            ClinvarGermlineAggregateDescription::try_from(
                pb_query::ClinvarGermlineAggregateDescription::LikelyPathogenic
            )
            .unwrap(),
            ClinvarGermlineAggregateDescription::LikelyPathogenic
        );
        assert_eq!(
            ClinvarGermlineAggregateDescription::try_from(
                pb_query::ClinvarGermlineAggregateDescription::UncertainSignificance
            )
            .unwrap(),
            ClinvarGermlineAggregateDescription::UncertainSignificance
        );
        assert_eq!(
            ClinvarGermlineAggregateDescription::try_from(
                pb_query::ClinvarGermlineAggregateDescription::LikelyBenign
            )
            .unwrap(),
            ClinvarGermlineAggregateDescription::LikelyBenign
        );
        assert_eq!(
            ClinvarGermlineAggregateDescription::try_from(
                pb_query::ClinvarGermlineAggregateDescription::Benign
            )
            .unwrap(),
            ClinvarGermlineAggregateDescription::Benign
        );
        assert!(ClinvarGermlineAggregateDescription::try_from(
            pb_query::ClinvarGermlineAggregateDescription::Unspecified
        )
        .is_err());
    }

    #[test]
    fn test_query_settings_clinvar_try_from() {
        let pb_query_settings_clinvar = pb_query::QuerySettingsClinVar {
            presence_required: true,
            germline_descriptions: vec![
                pb_query::ClinvarGermlineAggregateDescription::Pathogenic as i32,
                pb_query::ClinvarGermlineAggregateDescription::LikelyPathogenic as i32,
            ],
            allow_conflicting_interpretations: true,
        };
        let query_settings_clinvar = QuerySettingsClinVar {
            presence_required: true,
            germline_descriptions: vec![
                ClinvarGermlineAggregateDescription::Pathogenic,
                ClinvarGermlineAggregateDescription::LikelyPathogenic,
            ],
            allow_conflicting_interpretations: true,
        };
        assert_eq!(
            QuerySettingsClinVar::try_from(pb_query_settings_clinvar).unwrap(),
            query_settings_clinvar
        );
    }

    #[test]
    fn test_case_query_try_from() {
        let pb_case_query = pb_query::CaseQuery {
            genotype: Some(pb_query::QuerySettingsGenotype {
                recessive_mode: pb_query::RecessiveMode::Disabled as i32,
                sample_genotypes: vec![pb_query::SampleGenotypeChoice {
                    sample: "sample".to_string(),
                    genotype: pb_query::GenotypeChoice::Any as i32,
                    include_no_call: true,
                    enabled: true,
                }],
            }),
            quality: Some(pb_query::QuerySettingsQuality {
                sample_qualities: vec![pb_query::SampleQualitySettings {
                    sample: "sample".to_string(),
                    filter_active: true,
                    min_dp_het: Some(10),
                    min_dp_hom: Some(20),
                    min_gq: Some(30),
                    min_ab: Some(0.1),
                    min_ad: Some(40),
                    max_ad: Some(50),
                }],
            }),
            frequency: Some(pb_query::QuerySettingsFrequency {
                gnomad_exomes: Some(pb_query::NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                }),
                gnomad_genomes: Some(pb_query::NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                }),
                gnomad_mtdna: Some(pb_query::MitochondrialFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_af: Some(0.1),
                }),
                helixmtdb: Some(pb_query::MitochondrialFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_af: Some(0.1),
                }),
                inhouse: Some(pb_query::NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                }),
            }),
            consequence: Some(pb_query::QuerySettingsConsequence {
                variant_types: vec![
                    pb_query::VariantType::Snv as i32,
                    pb_query::VariantType::Indel as i32,
                    pb_query::VariantType::Mnv as i32,
                    pb_query::VariantType::ComplexSubstitution as i32,
                ],
                transcript_types: vec![
                    pb_query::TranscriptType::Coding as i32,
                    pb_query::TranscriptType::NonCoding as i32,
                ],
                consequences: vec![
                    pb_query::Consequence::TranscriptAblation as i32,
                    pb_query::Consequence::ExonLossVariant as i32,
                    pb_query::Consequence::SpliceAcceptorVariant as i32,
                    pb_query::Consequence::SpliceDonorVariant as i32,
                    pb_query::Consequence::StopGained as i32,
                    pb_query::Consequence::FrameshiftVariant as i32,
                    pb_query::Consequence::StopLost as i32,
                    pb_query::Consequence::StartLost as i32,
                ],
                max_dist_to_exon: Some(10),
            }),
            locus: Some(pb_query::QuerySettingsLocus {
                genes: vec!["gene".to_string()],
                genome_regions: vec![pb_query::GenomicRegion {
                    chrom: "chrom".to_string(),
                    range: Some(pb_query::Range { start: 1, stop: 2 }),
                }],
            }),
            clinvar: Some(pb_query::QuerySettingsClinVar {
                presence_required: true,
                germline_descriptions: vec![
                    pb_query::ClinvarGermlineAggregateDescription::Pathogenic as i32,
                    pb_query::ClinvarGermlineAggregateDescription::LikelyPathogenic as i32,
                ],
                allow_conflicting_interpretations: true,
            }),
        };
        let case_query = CaseQuery {
            genotype: QuerySettingsGenotype {
                recessive_mode: RecessiveMode::Disabled,
                sample_genotypes: {
                    let mut map = indexmap::IndexMap::new();
                    map.insert(
                        "sample".to_string(),
                        SampleGenotypeChoice {
                            sample: "sample".to_string(),
                            genotype: GenotypeChoice::Any,
                            include_no_call: true,
                            enabled: true,
                        },
                    );
                    map
                },
            },
            quality: QuerySettingsQuality {
                sample_qualities: {
                    let mut map = indexmap::IndexMap::new();
                    map.insert(
                        "sample".to_string(),
                        SampleQualitySettings {
                            sample: "sample".to_string(),
                            filter_active: true,
                            min_dp_het: Some(10),
                            min_dp_hom: Some(20),
                            min_gq: Some(30),
                            min_ab: Some(0.1),
                            min_ad: Some(40),
                            max_ad: Some(50),
                        },
                    );
                    map
                },
            },
            frequency: QuerySettingsFrequency {
                gnomad_exomes: NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                },
                gnomad_genomes: NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                },
                gnomad_mtdna: MitochondrialFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_af: Some(0.1),
                },
                helixmtdb: MitochondrialFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_af: Some(0.1),
                },
                inhouse: NuclearFrequencySettings {
                    enabled: true,
                    max_het: Some(10),
                    max_hom: Some(20),
                    max_hemi: Some(30),
                    max_af: Some(0.1),
                },
            },
            consequence: QuerySettingsConsequence {
                variant_types: vec![
                    VariantType::Snv,
                    VariantType::Indel,
                    VariantType::Mnv,
                    VariantType::ComplexSubstitution,
                ],
                transcript_types: vec![TranscriptType::Coding, TranscriptType::NonCoding],
                consequences: vec![
                    Consequence::TranscriptAblation,
                    Consequence::ExonLossVariant,
                    Consequence::SpliceAcceptorVariant,
                    Consequence::SpliceDonorVariant,
                    Consequence::StopGained,
                    Consequence::FrameshiftVariant,
                    Consequence::StopLost,
                    Consequence::StartLost,
                ],
                max_dist_to_exon: Some(10),
            },
            locus: QuerySettingsLocus {
                genes: vec!["gene".to_string()],
                genome_regions: vec![GenomicRegion {
                    chrom: "chrom".to_string(),
                    range: Some(Range { start: 1, stop: 2 }),
                }],
            },
            clinvar: QuerySettingsClinVar {
                presence_required: true,
                germline_descriptions: vec![
                    ClinvarGermlineAggregateDescription::Pathogenic,
                    ClinvarGermlineAggregateDescription::LikelyPathogenic,
                ],
                allow_conflicting_interpretations: true,
            },
        };
        assert_eq!(CaseQuery::try_from(pb_case_query).unwrap(), case_query);
    }

    #[rstest::rstest]
    #[case::empty("tests/seqvars/query/empty")]
    // #[case::full("tests/seqvars/query/full")]
    // #[case::with_extra("tests/seqvars/query/with_extra")]
    pub fn smoke_test_load(#[case] path_input: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}.json", path_input.split('/').last().unwrap());
        let query: super::CaseQuery =
            serde_json::from_reader(std::fs::File::open(format!("{}.json", path_input))?)?;

        if !path_input.contains("with_extra") {
            // Check if pbjson loading yields same result
            let query_pb: crate::pbs::varfish::v1::seqvars::query::CaseQuery =
                serde_json::from_reader(std::fs::File::open(format!("{}-pb.json", path_input))?)?;
            let query_pb_conv = super::CaseQuery::try_from(query_pb)?;
            assert_eq!(query, query_pb_conv);
        }
        insta::assert_yaml_snapshot!(&query);

        Ok(())
    }
}
