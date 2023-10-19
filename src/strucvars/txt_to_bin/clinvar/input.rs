//! Module with code supporting the parsing.
//!
//! Note that not the full model is implemented, only the parts that are needed for the
//! conversion of the ClinVar structural variants.

use crate::strucvars::query::clinvar::pbs::{Pathogenicity, VariationType};

/// Accession of a ClinVar record.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct ClinVarAccession {
    /// Accession string of the record.
    pub acc: String,
    /// Version of the accession.
    pub version: u32,
}

/// The values of review status are used to build the 'star ratings' displayed on the
/// ClinVar public site.
///
/// - 0 stars:  a conflict or not classified by submitter
/// - 1 star: classified by single submitter
/// - 2 stars: classified by multiple submitters
/// - 3 stars: reviewed by expert panel
/// - 4 stars: reviewed by professional society
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub enum ReviewStatus {
    /// no assertion provided
    #[serde(alias = "no assertion provided")]
    NoAssertionProvided,
    /// no assertion criteria provided
    #[serde(alias = "no assertion criteria provided")]
    NoAssertionCriteriaProvided,
    /// no assertion provided, single submitter
    #[serde(alias = "criteria provided, single submitter")]
    CriteriaProvidedSingleSubmitter,
    /// no assertion provided, multiple submitters, no conflicts
    #[serde(alias = "criteria provided, multiple submitters, no conflicts")]
    CriteriaProvidedMultipleSubmittersNoConflicts,
    /// no assertion provided, conflicting interpretations
    #[serde(alias = "criteria provided, conflicting interpretations")]
    CriteriaProvidedConflictingInterpretations,
    /// reviewed by expert panel
    #[serde(alias = "reviewed by expert panel")]
    ReviewedByExpertPanel,
    /// practice guideline
    #[serde(alias = "practice guideline")]
    PracticeGuideline,
}

/// Description of the clinical significance of a variant.
#[derive(Debug, Clone, Copy, serde::Deserialize, serde::Serialize)]
pub enum ClinicalSignificanceDescription {
    /// affects
    #[serde(alias = "affects")]
    Affects,
    /// benign
    #[serde(alias = "benign")]
    Benign,
    /// established risk allele
    #[serde(alias = "established risk allele")]
    EstablishedRiskAllele,
    /// likely benign
    #[serde(alias = "likely benign")]
    LikelyBenign,
    /// likely pathogenic
    #[serde(alias = "likely pathogenic")]
    LikelyPathogenic,
    /// likely pathogenic, low penetrance
    #[serde(alias = "likely pathogenic, low penetrance")]
    LikelyPathogenicLowPenetrance,
    /// likely risk allele
    #[serde(alias = "likely risk allele")]
    LikelyRiskAllele,
    /// pathogenic
    #[serde(alias = "pathogenic")]
    Pathogenic,
    /// pathogenic, low penetrance
    #[serde(alias = "pathogenic, low penetrance")]
    PathogenicLowPenetrance,
    /// uncertain risk allele
    #[serde(alias = "uncertain risk allele")]
    UncertainRiskAllele,
    /// uncertain significance
    #[serde(alias = "uncertain significance")]
    UncertainSignificance,
    /// association
    #[serde(alias = "association")]
    Association,
    /// association not found
    #[serde(alias = "association not found")]
    AssociationNotFound,
    /// confers sensitivity
    #[serde(alias = "confers sensitivity")]
    ConfersSensitivity,
    /// conflicting data from submitters
    #[serde(alias = "conflicting data from submitters")]
    ConflictingDataFromSubmitters,
    /// drug response
    #[serde(alias = "drug response")]
    DrugResponse,
    /// not provided
    #[serde(alias = "not provided")]
    NotProvided,
    /// other
    #[serde(alias = "other")]
    Other,
    /// protective
    #[serde(alias = "protective")]
    Protective,
    /// risk factor
    #[serde(alias = "risk factor")]
    RiskFactory,
}

impl TryInto<Pathogenicity> for ClinicalSignificanceDescription {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<Pathogenicity, Self::Error> {
        match self {
            ClinicalSignificanceDescription::LikelyBenign => Ok(Pathogenicity::Benign),
            ClinicalSignificanceDescription::LikelyPathogenic
            | ClinicalSignificanceDescription::LikelyPathogenicLowPenetrance => {
                Ok(Pathogenicity::LikelyPathogenic)
            }
            ClinicalSignificanceDescription::Pathogenic
            | ClinicalSignificanceDescription::PathogenicLowPenetrance => {
                Ok(Pathogenicity::Pathogenic)
            }
            ClinicalSignificanceDescription::UncertainSignificance => Ok(Pathogenicity::Uncertain),
            ClinicalSignificanceDescription::Benign
            | ClinicalSignificanceDescription::Protective => Ok(Pathogenicity::Benign),
            _ => anyhow::bail!("unsupported clinical significance: {:?}", self),
        }
    }
}

/// Description of the clinical significance of a variant.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct ClinicalSignificance {
    /// Review status of the variant.
    pub review_status: ReviewStatus,
    /// Description of the clinical significance of a variant.
    pub description: ClinicalSignificanceDescription,
}

/// Identifier of an assembly.
#[derive(Debug, Clone, PartialEq, Eq, Copy, serde::Deserialize, serde::Serialize)]
pub enum Assembly {
    /// GRCh38
    #[serde(alias = "GRCh38")]
    Grch38,
    /// GRCh37
    #[serde(alias = "GRCh37")]
    Grch37,
    /// NCBI36
    #[serde(alias = "NCBI36")]
    Ncib36,
}

/// Description of a location on a sequence.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct SequenceLocation {
    /// Assembly of the location
    pub assembly: Assembly,
    /// Chromosome
    pub chr: String,
    /// Accession of the chromosome
    pub accession: Option<String>,
    /// Start position
    pub start: Option<i32>,
    /// Stop position,
    pub stop: Option<i32>,
    /// Inner start position
    pub inner_start: Option<i32>,
    /// Inner stop position
    pub inner_stop: Option<i32>,
    /// Outer start position
    pub outer_start: Option<i32>,
    /// Outer stop position
    pub outer_stop: Option<i32>,
    /// Variant length
    pub variant_length: Option<i32>,
    /// Position in VCF format
    pub position_vcf: Option<i32>,
    /// Reference allele in VCF format
    pub reference_allele_vcf: Option<String>,
    /// Alternate allele in VCF format
    pub alternate_allele_vcf: Option<String>,
}

#[derive(Debug, Clone, Copy, serde::Deserialize, serde::Serialize)]
pub enum MeasureType {
    /// Gene
    #[serde(alias = "gene")]
    Gene,
    /// Variation
    #[serde(alias = "variation")]
    Variation,
    /// Insertion
    #[serde(alias = "insertion")]
    Insertion,
    /// Deletion
    #[serde(alias = "deletion")]
    Deletion,
    /// Single nucleotide variant
    #[serde(alias = "single nucleotide variant")]
    SingleNucleotideVariant,
    /// Indel
    #[serde(alias = "indel")]
    Indel,
    /// Duplication
    #[serde(alias = "duplication")]
    Duplication,
    /// Tandem duplication
    #[serde(alias = "tandem duplication")]
    TandemDuplication,
    /// Structural variant
    #[serde(alias = "structural variant")]
    StructuralVariant,
    /// Copy number gain
    #[serde(alias = "copy number gain")]
    CopyNumberGain,
    /// Copy number loss
    #[serde(alias = "copy number loss")]
    CopyNumberLoss,
    /// Protein change only
    #[serde(alias = "protein only")]
    ProteinOnly,
    /// Microsatellite
    #[serde(alias = "microsatellite")]
    Microsatellite,
    // Fusion
    #[serde(alias = "fusion")]
    Fusion,
    /// Inversion
    #[serde(alias = "inversion")]
    Inversion,
    /// Translocation
    #[serde(alias = "translocation")]
    Translocation,
    /// QTL
    #[serde(alias = "qtl")]
    Qtl,
    /// Complex variant
    #[serde(alias = "complex")]
    Complex,
    /// Other variant
    #[serde(alias = "other")]
    Other,
}

impl TryInto<VariationType> for MeasureType {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<VariationType, Self::Error> {
        match self {
            MeasureType::Deletion => Ok(VariationType::Del),
            MeasureType::Duplication => Ok(VariationType::Dup),
            MeasureType::TandemDuplication => Ok(VariationType::Dup),
            MeasureType::CopyNumberGain => Ok(VariationType::Dup),
            MeasureType::CopyNumberLoss => Ok(VariationType::Del),
            MeasureType::Microsatellite => Ok(VariationType::Microsatellite),
            MeasureType::Inversion => Ok(VariationType::Inv),
            MeasureType::Translocation => Ok(VariationType::Bnd),
            MeasureType::Complex => Ok(VariationType::Complex),
            MeasureType::SingleNucleotideVariant
            | MeasureType::Gene
            | MeasureType::Variation
            | MeasureType::Insertion
            | MeasureType::Indel
            | MeasureType::ProteinOnly
            | MeasureType::Fusion
            | MeasureType::Qtl
            | MeasureType::StructuralVariant
            | MeasureType::Other => anyhow::bail!("unsupported measure type: {:?}", self),
        }
    }
}

/// A single measure / variant.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Measure {
    /// Type of the measure.
    pub r#type: MeasureType,
    /// Locations in the sequence, for each assembly.
    pub sequence_locations: Vec<SequenceLocation>,
}

/// A set of measures / variants.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct MeasureSet {
    /// Accession of the record.
    pub acc: String,
    /// Version of the accession.
    pub version: u32,
    /// Measures / variants.
    pub measures: Vec<Measure>,
}

/// A reference ClinVar assertion.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct ReferenceClinVarAssertion {
    /// The ClinVar accession
    pub clinvar_accession: ClinVarAccession,
    /// The ClinVar significance
    pub clinical_significance: ClinicalSignificance,
    /// The ClinVar measure set
    pub measures: MeasureSet,
}

/// A ClinVar set.
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct ClinVarSet {
    /// The reference ClinVar accession
    pub reference_clinvar_assertion: ReferenceClinVarAssertion,
}
