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
    /// Recessive parent.
    RecessiveParent,
}

/// Supporting code for `GenotypeChoice`.
pub(crate) mod genotype_choice {
    /// Error type for `GenotypeChoice::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("Cannot convert protobuf GenotypeChoice: {0:?}")]
        UnknownGenotypeChoiceValue(super::pb_query::GenotypeChoice),
    }
}

impl TryFrom<pb_query::GenotypeChoice> for GenotypeChoice {
    type Error = genotype_choice::Error;

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
            pb_query::GenotypeChoice::RecessiveParent => Ok(GenotypeChoice::RecessiveParent),
            _ => Err(Self::Error::UnknownGenotypeChoiceValue(value)),
        }
    }
}

/// Trait that describes whether a string matches a value.
///
/// Note that we assume properly ingested VCFs with only one alternate allele.
/// The valid genotype strings have the form "<VAL>/<VAL>", "<VAL>|<VAL>" or
/// "<VAL>" with "<VAL>" being one of "0", "1", and ".".
pub trait MatchesGenotypeStr {
    /// Whether `self` matches `s`.
    ///
    /// # Arguments
    ///
    /// * `gt_str` - The string to match against.
    ///
    /// # Returns
    ///
    /// `true` if `self` matches `gt_str`, `false` otherwise.
    fn matches(&self, gt_tr: &str) -> bool;
}

impl MatchesGenotypeStr for GenotypeChoice {
    fn matches(&self, gt_str: &str) -> bool {
        let gt_str = if gt_str.starts_with('/') || gt_str.starts_with('|') {
            &gt_str[1..]
        } else {
            gt_str
        };
        match self {
            GenotypeChoice::Any => true,
            GenotypeChoice::Ref => ["0", "0|0", "0/0"].contains(&gt_str),
            GenotypeChoice::RecessiveIndex
            | GenotypeChoice::RecessiveParent
            | GenotypeChoice::Het => ["0/1", "0|1", "1/0", "1|0"].contains(&gt_str),
            GenotypeChoice::Hom => ["1", "1/1", "1|1"].contains(&gt_str),
            GenotypeChoice::NonHom => !["1", "1/1", "1|1"].contains(&gt_str),
            GenotypeChoice::NonHet => !["0/1", "0|1", "1/0", "1|0"].contains(&gt_str),
            GenotypeChoice::Variant => {
                ["1", "0/1", "0|1", "1/0", "1|0", "1|1", "1/1"].contains(&gt_str)
            }
        }
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
        #[error("Too many recessive parent samples found: {0:?}")]
        TooManyRecessiveParentsFound(Vec<String>),
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
    /// The sample names of the recessive parent samples (zero to two).
    ///
    /// # Errors
    ///
    /// * `RecessiveParentsError::TooManyRecessiveParentsFound` if more than two recessive parent samples are found.
    pub fn recessive_parents(
        &self,
    ) -> Result<Vec<String>, query_settings_genotype::RecessiveParentsError> {
        let samples = self
            .sample_genotypes
            .values()
            .filter(|sgc| sgc.genotype == GenotypeChoice::RecessiveParent)
            .map(|sgc| sgc.sample.clone())
            .collect::<Vec<_>>();
        if samples.len() > 2 {
            Err(
                query_settings_genotype::RecessiveParentsError::TooManyRecessiveParentsFound(
                    samples,
                ),
            )
        } else {
            Ok(samples)
        }
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
                .map_err(|e| Self::Error::InvalidSampleGenotypeChoice(e))?;
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

/// gnomAD nuclear filter options.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GnomadNuclearFrequencySettings {
    /// Whether to enable filtration by 1000 Genomes.
    pub enabled: bool,
    /// Maximal number of in-house heterozygous carriers.
    pub heterozygous: Option<i32>,
    /// Maximal number of in-house homozygous carriers.
    pub homozygous: Option<i32>,
    /// Maximal number of in-house hemizygous carriers.
    pub hemizygous: Option<i32>,
    /// Maximal allele frequency.
    pub frequency: Option<f32>,
}

impl Eq for GnomadNuclearFrequencySettings {}

impl From<pb_query::GnomadNuclearFreqyencySettings> for GnomadNuclearFrequencySettings {
    fn from(value: pb_query::GnomadNuclearFreqyencySettings) -> Self {
        Self {
            enabled: value.enabled,
            heterozygous: value.heterozygous,
            homozygous: value.homozygous,
            hemizygous: value.hemizygous,
            frequency: value.frequency,
        }
    }
}

/// gnomAD mitochondrial filter options.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GnomadMitochondrialFrequencySettings {
    /// Whether to enable filtration by 1000 Genomes.
    pub enabled: bool,
    /// Maximal number of heteroplasmic carriers.
    pub heteroplasmic: Option<i32>,
    /// Maximal number of homoplasmic carriers.
    pub homoplasmic: Option<i32>,
    /// Maximal allele frequency.
    pub frequency: Option<f32>,
}

impl Eq for GnomadMitochondrialFrequencySettings {}

impl From<pb_query::GnomadMitochondrialFrequencySettings> for GnomadMitochondrialFrequencySettings {
    fn from(value: pb_query::GnomadMitochondrialFrequencySettings) -> Self {
        Self {
            enabled: value.enabled,
            heteroplasmic: value.heteroplasmic,
            homoplasmic: value.homoplasmic,
            frequency: value.frequency,
        }
    }
}

/// HelixMtDb filter options.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct HelixMtDbFrequencySettings {
    /// Whether to enable filtration by mtDB.
    pub enabled: bool,
    /// Maximal number of heterozygous carriers in HelixMtDb.
    pub heteroplasmic: Option<i32>,
    /// Maximal number of homozygous carriers in HelixMtDb.
    pub homoplasmic: Option<i32>,
    /// Maximal frequency in HelixMtDb.
    pub frequency: Option<f32>,
}

impl Eq for HelixMtDbFrequencySettings {}

impl From<pb_query::HelixMtDbFrequencySettings> for HelixMtDbFrequencySettings {
    fn from(value: pb_query::HelixMtDbFrequencySettings) -> Self {
        Self {
            enabled: value.enabled,
            heteroplasmic: value.heteroplasmic,
            homoplasmic: value.homoplasmic,
            frequency: value.frequency,
        }
    }
}

/// In-house frequency filter options.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct InhouseFrequencySettings {
    /// Whether to enable filtration by 1000 Genomes.
    pub enabled: bool,
    /// Maximal number of in-house heterozygous carriers.
    pub heterozygous: Option<i32>,
    /// Maximal number of in-house homozygous carriers.
    pub homozygous: Option<i32>,
    /// Maximal number of in-house hemizygous carriers.
    pub hemizygous: Option<i32>,
    /// Maximal number of in-house carriers.
    pub carriers: Option<i32>,
}

impl From<pb_query::InhouseFrequencySettings> for InhouseFrequencySettings {
    fn from(value: pb_query::InhouseFrequencySettings) -> Self {
        Self {
            enabled: value.enabled,
            heterozygous: value.heterozygous,
            homozygous: value.homozygous,
            hemizygous: value.hemizygous,
            carriers: value.carriers,
        }
    }
}

/// Query settings for population frequencies.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct QuerySettingsFrequency {
    /// gnomAD-exomes filter.
    pub gnomad_exomes: GnomadNuclearFrequencySettings,
    /// gnomAD-genomes filter.
    pub gnomad_genomes: GnomadNuclearFrequencySettings,
    /// gnomAD-MT filter.
    pub gnomad_mt: GnomadMitochondrialFrequencySettings,
    /// HelixMtDb filter.
    pub helixmtdb: HelixMtDbFrequencySettings,
    /// In-house filter.
    pub inhouse: InhouseFrequencySettings,
}

/// Supporting code for `QuerySettingsFrequency`.
pub(crate) mod query_settings_frequency {
    /// Error type for `QuerySettingsFrequency::try_from()`.
    #[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
    pub enum Error {
        #[error("gnomad_exomes missing in protobuf")]
        GnomadExomesMissing,
        #[error("gnomad_genomes missing in protobuf")]
        GnomadGenomesMissing,
        #[error("gnomad_mt missing in protobuf")]
        GnomadMtMissing,
        #[error("helixmtdb missing in protobuf")]
        HelixMtDbMissing,
        #[error("inhouse missing in protobuf")]
        InhouseMissing,
    }
}

impl TryFrom<pb_query::QuerySettingsFrequency> for QuerySettingsFrequency {
    type Error = query_settings_frequency::Error;

    fn try_from(value: pb_query::QuerySettingsFrequency) -> Result<Self, Self::Error> {
        Ok(Self {
            gnomad_exomes: GnomadNuclearFrequencySettings::from(
                value
                    .gnomad_exomes
                    .ok_or(Self::Error::GnomadExomesMissing)?,
            ),
            gnomad_genomes: GnomadNuclearFrequencySettings::from(
                value
                    .gnomad_genomes
                    .ok_or(Self::Error::GnomadGenomesMissing)?,
            ),
            gnomad_mt: GnomadMitochondrialFrequencySettings::from(
                value.gnomad_mt.ok_or(Self::Error::GnomadMtMissing)?,
            ),
            helixmtdb: HelixMtDbFrequencySettings::from(
                value.helixmtdb.ok_or(Self::Error::HelixMtDbMissing)?,
            ),
            inhouse: InhouseFrequencySettings::from(
                value.inhouse.ok_or(Self::Error::InhouseMissing)?,
            ),
        })
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
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize, strum::EnumIter
)]
pub enum Consequence {
    /*
     * High impact.
     */
    /// Transcript ablation.
    TranscriptAblation,
    /// Exon loss variant.
    ExonLossVariant,
    /// Splice acceptor variant.
    SpliceAcceptorVariant,
    /// Splice donor variant.
    SpliceDonorVariant,
    /// Stop gained.
    StopGained,
    /// Frameshift variant.
    FrameshiftVariant,
    /// Stop lost.
    StopLost,
    /// Start lost.
    StartLost,
    /// Transcript amplification.
    TranscriptAmplification,
    /*
     * Moderate impact.
     */
    /// Disruptive inframe insertion.
    DisruptiveInframeInsertion,
    /// Disruptive inframe deletion.
    DisruptiveInframeDeletion,
    /// Conservative inframe insertion.
    ConservativeInframeInsertion,
    /// Conservative inframe deletion.
    ConservativeInframeDeletion,
    /// Missense variant.
    MissenseVariant,
    /*
     * Low impact.
     */
    /// Splice donor 5th base variant.
    SpliceDonorFifthBaseVariant,
    /// Splice region variant.
    SpliceRegionVariant,
    /// Splice donor region variant.
    SpliceDonorRegionVariant,
    /// Splice polypyrimidine tract variant.
    SplicePolypyrimidineTractVariant,
    /// Start retained variant.
    StartRetainedVariant,
    /// Stop retained variant.
    StopRetainedVariant,
    /// Synonymous variant.
    SynonymousVariant,
    /*
     * Modifier.
     */
    /// Coding sequence variant.
    CodingSequenceVariant,
    /// 5' UTR exon variant.
    FivePrimeUtrExonVariant,
    /// 5' UTR intron variant.
    FivePrimeUtrIntronVariant,
    /// 3' UTR exon variant.
    ThreePrimeUtrExonVariant,
    /// 3' UTR intron variant.
    ThreePrimeUtrIntronVariant,
    /// Non-coding transcript exon variant.
    NonCodingTranscriptExonVariant,
    /// Non-coding transcript intron variant.
    NonCodingTranscriptIntronVariant,
    /// Upstream gene variant.
    UpstreamGeneVariant,
    /// Downstream gene variant.
    DownstreamGeneVariant,
    /// Intergenic variant.
    IntergenicVariant,
    /// Intron variant.
    IntronVariant,
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
            pb_query::Consequence::IntergenicVariant => Ok(Consequence::IntergenicVariant),
            pb_query::Consequence::IntronVariant => Ok(Consequence::IntronVariant),
            _ => Err(consequence::Error::UnknownConsequenceValue(value)),
        }
    }
}

impl Into<mehari::annotate::seqvars::ann::Consequence> for Consequence {
    fn into(self) -> mehari::annotate::seqvars::ann::Consequence {
        use mehari::annotate::seqvars::ann;

        match self {
            Consequence::TranscriptAblation => ann::Consequence::TranscriptAblation,
            Consequence::ExonLossVariant => ann::Consequence::ExonLossVariant,
            Consequence::SpliceAcceptorVariant => ann::Consequence::SpliceAcceptorVariant,
            Consequence::SpliceDonorVariant => ann::Consequence::SpliceDonorVariant,
            Consequence::StopGained => ann::Consequence::StopGained,
            Consequence::FrameshiftVariant => ann::Consequence::FrameshiftVariant,
            Consequence::StopLost => ann::Consequence::StopLost,
            Consequence::StartLost => ann::Consequence::StartLost,
            Consequence::TranscriptAmplification => ann::Consequence::TranscriptAmplification,
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
            Consequence::IntergenicVariant => ann::Consequence::IntergenicVariant,
            Consequence::IntronVariant => ann::Consequence::IntronVariant,
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
        UnknownVariantTypeInt(i32),
        #[error("Cannot convert protobuf VariantType: {0:?}")]
        UnknownVariantTypeValue(super::pb_query::VariantType),
        #[error("Cannot convert i32 into protobuf TranscriptType: {0}")]
        UnknownTranscriptTypeInt(i32),
        #[error("Cannot convert protobuf TranscriptType: {0:?}")]
        UnknownTranscriptTypeValue(super::pb_query::TranscriptType),
        #[error("Cannot convert i32 into protobuf Consequence: {0}")]
        UnknownConsequenceInt(i32),
        #[error("Cannot convert protobuf Consequence: {0:?}")]
        UnknownConsequenceValue(super::pb_query::Consequence),
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
                    .map_err(|_| query_settings_consequence::Error::UnknownVariantTypeInt(v))?;
                VariantType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::UnknownVariantTypeValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let transcript_types = value
            .transcript_types
            .into_iter()
            .map(|v| {
                let v = pb_query::TranscriptType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::UnknownTranscriptTypeInt(v))?;
                TranscriptType::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::UnknownTranscriptTypeValue(v))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let consequences = value
            .consequences
            .into_iter()
            .map(|v| {
                let v = pb_query::Consequence::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::UnknownConsequenceInt(v))?;
                Consequence::try_from(v)
                    .map_err(|_| query_settings_consequence::Error::UnknownConsequenceValue(v))
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
    pub end: i32,
}

impl From<pb_query::Range> for Range {
    fn from(value: pb_query::Range) -> Self {
        Self {
            start: value.start,
            end: value.end,
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
        #[error("Missing genotype in protobuf")]
        GenotypeMissing,
        #[error("Problem converting protobuf for genotype: {0}")]
        GenotypeConversion(#[from] super::query_settings_genotype::Error),
        #[error("Missing quality in protobuf")]
        QualityMissing,
        #[error("Problem converting protobuf for quality: {0}")]
        QualityConversion(#[from] super::query_settings_quality::Error),
        #[error("Missing frequency in protobuf")]
        FrequencyMissing,
        #[error("Problem converting protobuf for frequency: {0}")]
        FrequencyConversion(#[from] super::query_settings_frequency::Error),
        #[error("Missing consequence in protobuf")]
        ConsequenceMissing,
        #[error("Problem converting protobuf for consequence: {0}")]
        ConsequenceConversion(#[from] super::query_settings_consequence::Error),
        #[error("Missing locus in protobuf")]
        LocusMissing,
        #[error("Missing clinvar in protobuf")]
        ClinVarMissing,
        #[error("Problem converting protobuf for clinvar: {0}")]
        ClinVarConversion(#[from] super::query_settings_clinvar::Error),
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

        let genotype =
            QuerySettingsGenotype::try_from(genotype.ok_or(Self::Error::GenotypeMissing)?)
                .map_err(|e| Self::Error::GenotypeConversion(e))?;
        let quality = QuerySettingsQuality::try_from(quality.ok_or(Self::Error::QualityMissing)?)
            .map_err(|e| Self::Error::QualityConversion(e))?;
        let frequency =
            QuerySettingsFrequency::try_from(frequency.ok_or(Self::Error::FrequencyMissing)?)
                .map_err(|e| Self::Error::FrequencyConversion(e))?;
        let consequence =
            QuerySettingsConsequence::try_from(consequence.ok_or(Self::Error::ConsequenceMissing)?)
                .map_err(|e| Self::Error::ConsequenceConversion(e))?;
        let locus = QuerySettingsLocus::from(locus.ok_or(Self::Error::LocusMissing)?);
        let clinvar = QuerySettingsClinVar::try_from(clinvar.ok_or(Self::Error::ClinVarMissing)?)
            .map_err(|e| Self::Error::ClinVarConversion(e))?;

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
            GenotypeChoice::try_from(pb_query::GenotypeChoice::RecessiveParent).unwrap(),
            GenotypeChoice::RecessiveParent
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
        #[case(Any, &["0", "0/0", "0|0", "1", "1/1", "1|1", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], true)]
        #[case(Ref, &["0", "0/0", "0|0"], true)]
        #[case(Ref, &[ "1", "1/1", "1|1", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], false)]
        #[case(Het, &["0/1", "0|1", "1/0", "1|0"], true)]
        #[case(Het, &["0", "0/0", "0|0", "1", "1/1", "1|1", ".", "./.", ".|."], false)]
        #[case(Hom, &["1", "1/1", "1|1"], true)]
        #[case(Hom, &["0", "0/0", "0|0", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], false)]
        #[case(NonHom, &["1", "1/1", "1|1"], false)]
        #[case(NonHom, &["0", "0/0", "0|0", "0/1", "0/.", "0|1", "1/0", "1|0", ".", "./.", ".|."], true)]
        #[case(Variant, &["1", "1/1", "1|1", "0/1",  "0|1",  "1/0", "1|0"], true)]
        #[case(Variant, &["0", "0/0", "0|0", ".", "./.", ".|."], false)]
        pub fn matches(
            #[case] genotype_choice: GenotypeChoice,
            #[case] gt_strs: &[&str],
            #[case] expected: bool,
        ) {
            use crate::seqvars::query::schema::query::MatchesGenotypeStr as _;

            for gt_str in gt_strs {
                assert_eq!(
                    genotype_choice.matches(gt_str),
                    expected,
                    "{:?} {}",
                    genotype_choice,
                    gt_str
                );
            }
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
        let pb_gnomad_nuclear_frequency_settings = pb_query::GnomadNuclearFreqyencySettings {
            enabled: true,
            heterozygous: Some(10),
            homozygous: Some(20),
            hemizygous: Some(30),
            frequency: Some(0.1),
        };
        let gnomad_nuclear_frequency_settings = GnomadNuclearFrequencySettings {
            enabled: true,
            heterozygous: Some(10),
            homozygous: Some(20),
            hemizygous: Some(30),
            frequency: Some(0.1),
        };
        assert_eq!(
            GnomadNuclearFrequencySettings::from(pb_gnomad_nuclear_frequency_settings),
            gnomad_nuclear_frequency_settings
        );
    }

    #[test]
    fn test_gnomad_mitochondrial_frequency_settings_from() {
        let pb_gnomad_mitochondrial_frequency_settings =
            pb_query::GnomadMitochondrialFrequencySettings {
                enabled: true,
                heteroplasmic: Some(10),
                homoplasmic: Some(20),
                frequency: Some(0.1),
            };
        let gnomad_mitochondrial_frequency_settings = GnomadMitochondrialFrequencySettings {
            enabled: true,
            heteroplasmic: Some(10),
            homoplasmic: Some(20),
            frequency: Some(0.1),
        };
        assert_eq!(
            GnomadMitochondrialFrequencySettings::from(pb_gnomad_mitochondrial_frequency_settings),
            gnomad_mitochondrial_frequency_settings
        );
    }

    #[test]
    fn test_helix_mtdb_frequency_settings_from() {
        let pb_helix_mtdb_frequency_settings = pb_query::HelixMtDbFrequencySettings {
            enabled: true,
            heteroplasmic: Some(10),
            homoplasmic: Some(20),
            frequency: Some(0.1),
        };
        let helix_mtdb_frequency_settings = HelixMtDbFrequencySettings {
            enabled: true,
            heteroplasmic: Some(10),
            homoplasmic: Some(20),
            frequency: Some(0.1),
        };
        assert_eq!(
            HelixMtDbFrequencySettings::from(pb_helix_mtdb_frequency_settings),
            helix_mtdb_frequency_settings
        );
    }

    #[test]
    fn test_inhouse_frequency_settings_from() {
        let pb_inhouse_frequency_settings = pb_query::InhouseFrequencySettings {
            enabled: true,
            heterozygous: Some(10),
            homozygous: Some(20),
            hemizygous: Some(30),
            carriers: Some(40),
        };
        let inhouse_frequency_settings = InhouseFrequencySettings {
            enabled: true,
            heterozygous: Some(10),
            homozygous: Some(20),
            hemizygous: Some(30),
            carriers: Some(40),
        };
        assert_eq!(
            InhouseFrequencySettings::from(pb_inhouse_frequency_settings),
            inhouse_frequency_settings
        );
    }

    #[test]
    fn test_query_settings_frequency_try_from() {
        let pb_query_settings_frequency = pb_query::QuerySettingsFrequency {
            gnomad_exomes: Some(pb_query::GnomadNuclearFreqyencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                frequency: Some(0.1),
            }),
            gnomad_genomes: Some(pb_query::GnomadNuclearFreqyencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                frequency: Some(0.1),
            }),
            gnomad_mt: Some(pb_query::GnomadMitochondrialFrequencySettings {
                enabled: true,
                heteroplasmic: Some(10),
                homoplasmic: Some(20),
                frequency: Some(0.1),
            }),
            helixmtdb: Some(pb_query::HelixMtDbFrequencySettings {
                enabled: true,
                heteroplasmic: Some(10),
                homoplasmic: Some(20),
                frequency: Some(0.1),
            }),
            inhouse: Some(pb_query::InhouseFrequencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                carriers: Some(40),
            }),
        };
        let query_settings_frequency = QuerySettingsFrequency {
            gnomad_exomes: GnomadNuclearFrequencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                frequency: Some(0.1),
            },
            gnomad_genomes: GnomadNuclearFrequencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                frequency: Some(0.1),
            },
            gnomad_mt: GnomadMitochondrialFrequencySettings {
                enabled: true,
                heteroplasmic: Some(10),
                homoplasmic: Some(20),
                frequency: Some(0.1),
            },
            helixmtdb: HelixMtDbFrequencySettings {
                enabled: true,
                heteroplasmic: Some(10),
                homoplasmic: Some(20),
                frequency: Some(0.1),
            },
            inhouse: InhouseFrequencySettings {
                enabled: true,
                heterozygous: Some(10),
                homozygous: Some(20),
                hemizygous: Some(30),
                carriers: Some(40),
            },
        };
        assert_eq!(
            QuerySettingsFrequency::try_from(pb_query_settings_frequency).unwrap(),
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
                pb_query::VariantType::Unspecified as i32,
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
        let pb_range = pb_query::Range { start: 1, end: 2 };
        let range = Range { start: 1, end: 2 };
        assert_eq!(Range::from(pb_range), range);
    }

    #[test]
    fn test_genomic_region_from() {
        let pb_genomic_region = pb_query::GenomicRegion {
            chrom: "chrom".to_string(),
            range: Some(pb_query::Range { start: 1, end: 2 }),
        };
        let genomic_region = GenomicRegion {
            chrom: "chrom".to_string(),
            range: Some(Range { start: 1, end: 2 }),
        };
        assert_eq!(GenomicRegion::from(pb_genomic_region), genomic_region);
    }

    #[test]
    fn test_query_settings_locus_from() {
        let pb_query_settings_locus = pb_query::QuerySettingsLocus {
            genes: vec!["gene".to_string()],
            genome_regions: vec![pb_query::GenomicRegion {
                chrom: "chrom".to_string(),
                range: Some(pb_query::Range { start: 1, end: 2 }),
            }],
        };
        let query_settings_locus = QuerySettingsLocus {
            genes: vec!["gene".to_string()],
            genome_regions: vec![GenomicRegion {
                chrom: "chrom".to_string(),
                range: Some(Range { start: 1, end: 2 }),
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
                gnomad_exomes: Some(pb_query::GnomadNuclearFreqyencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    frequency: Some(0.1),
                }),
                gnomad_genomes: Some(pb_query::GnomadNuclearFreqyencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    frequency: Some(0.1),
                }),
                gnomad_mt: Some(pb_query::GnomadMitochondrialFrequencySettings {
                    enabled: true,
                    heteroplasmic: Some(10),
                    homoplasmic: Some(20),
                    frequency: Some(0.1),
                }),
                helixmtdb: Some(pb_query::HelixMtDbFrequencySettings {
                    enabled: true,
                    heteroplasmic: Some(10),
                    homoplasmic: Some(20),
                    frequency: Some(0.1),
                }),
                inhouse: Some(pb_query::InhouseFrequencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    carriers: Some(40),
                }),
            }),
            consequence: Some(pb_query::QuerySettingsConsequence {
                variant_types: vec![
                    pb_query::VariantType::Unspecified as i32,
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
                    range: Some(pb_query::Range { start: 1, end: 2 }),
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
                gnomad_exomes: GnomadNuclearFrequencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    frequency: Some(0.1),
                },
                gnomad_genomes: GnomadNuclearFrequencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    frequency: Some(0.1),
                },
                gnomad_mt: GnomadMitochondrialFrequencySettings {
                    enabled: true,
                    heteroplasmic: Some(10),
                    homoplasmic: Some(20),
                    frequency: Some(0.1),
                },
                helixmtdb: HelixMtDbFrequencySettings {
                    enabled: true,
                    heteroplasmic: Some(10),
                    homoplasmic: Some(20),
                    frequency: Some(0.1),
                },
                inhouse: InhouseFrequencySettings {
                    enabled: true,
                    heterozygous: Some(10),
                    homozygous: Some(20),
                    hemizygous: Some(30),
                    carriers: Some(40),
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
                    range: Some(Range { start: 1, end: 2 }),
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
    #[case("tests/seqvars/query/empty")]
    #[case("tests/seqvars/query/full")]
    #[case("tests/seqvars/query/with_extra")]
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
