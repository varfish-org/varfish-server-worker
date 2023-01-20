/// Schema for the query definition.
use std::collections::HashMap;

/// Definition of query schemas.
use serde::{Deserialize, Serialize};

/// Range with 1-based positions
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct Range {
    pub start: u32,
    pub end: u32,
}

impl Range {
    pub fn new(start: u32, end: u32) -> Self {
        Range { start, end }
    }
}

/// Genomic region
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct GenomicRegion {
    pub chrom: String,
    pub range: Option<Range>,
}

impl GenomicRegion {
    pub fn new(chrom: &str, start: u32, end: u32) -> Self {
        GenomicRegion {
            chrom: chrom.to_owned(),
            range: Some(Range::new(start, end)),
        }
    }

    pub fn whole_chrom(chrom: &str) -> Self {
        GenomicRegion {
            chrom: chrom.to_owned(),
            range: None,
        }
    }
}

/// Database of transcripts
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum Database {
    Refseq,
    Ensembl,
}

/// Encode the type of an SV
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum SvType {
    /// Deletion
    Del,
    /// Duplication
    Dup,
    /// Inversion
    Inv,
    /// Insertion
    Ins,
    /// Break-end
    Bnd,
    /// Copy number variable region
    Cnv,
}

/// Structural variant sub type
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum SvSubType {
    /// Deletion
    #[serde(rename = "DEL")]
    Del,
    /// Mobile element deletion
    #[serde(rename = "DEL:ME")]
    DelMe,
    /// Mobile element deletion (SVA)
    #[serde(rename = "DEL:ME:SVA")]
    DelMeSva,
    /// Mobile element deletion (L1)
    #[serde(rename = "DEL:ME:L1")]
    DelMeL1,
    /// Mobile element deletion (ALU)
    #[serde(rename = "DEL:ME:ALU")]
    DelMeAlu,
    /// Duplication
    #[serde(rename = "DUP")]
    Dup,
    /// Tandem duplication
    #[serde(rename = "DUP:TANDEM")]
    DupTandem,
    /// Inversion
    #[serde(rename = "INV")]
    Inv,
    /// Insertion
    #[serde(rename = "INS")]
    Ins,
    /// Mobile element insertion
    #[serde(rename = "INS:ME")]
    InsMe,
    /// Mobile element insertion (SVA)
    #[serde(rename = "INS:ME:SVA")]
    InsMeSva,
    /// Mobile element insertion (L1)
    #[serde(rename = "INS:ME:L1")]
    InsMeL1,
    /// Mobile element insertion (ALU)
    #[serde(rename = "INS:ME:ALU")]
    InsMeAlu,
    /// Break-end
    #[serde(rename = "BND")]
    Bnd,
    /// Copy number variable region
    #[serde(rename = "CNV")]
    Cnv,
}

impl SvSubType {
    /// Return whether SV sub type is any insertion
    pub fn is_ins(&self) -> bool {
        matches!(
            self,
            SvSubType::Ins
                | SvSubType::InsMe
                | SvSubType::InsMeSva
                | SvSubType::InsMeL1
                | SvSubType::InsMeAlu
        )
    }
    /// Return whether SV sub type is any deletion
    pub fn is_del(&self) -> bool {
        matches!(
            self,
            SvSubType::Del
                | SvSubType::DelMe
                | SvSubType::DelMeSva
                | SvSubType::DelMeL1
                | SvSubType::DelMeAlu
        )
    }
}

/// Genotype choice for filter dropdowns
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
#[serde(rename_all = "kebab-case")]
pub enum GenotypeChoice {
    /// Any genotype
    Any,
    /// Reference genotype
    Ref,
    /// Heterozygous genotype
    Het,
    /// Homozygous alternative genotype
    Hom,
    /// Not homozygous alternative genotype
    NonHom,
    /// Not wild-type genotype
    Variant,
    /// Not variant genotype
    NonVariant,
    /// Non-reference genotype
    NonReference,
}

/// ENSEMBL regulatory feature
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum EnsemblRegulatoryFeature {
    /// Any fature
    #[serde(rename = "any_feature")]
    AnyFeature,
    /// CTCF binding site
    #[serde(rename = "CTCF_binding_site")]
    CtcfBindingSite,
    /// Enhancer
    #[serde(rename = "enhancer")]
    Enhancer,
    /// Open chromatin region
    #[serde(rename = "open_chromatin_region")]
    OpenChromatinRegion,
    /// Promoter region
    #[serde(rename = "promoter")]
    Promoter,
    /// Promoter flanking region
    #[serde(rename = "promoter_flanking_region")]
    PromoterFlankingRegion,
    /// Transcription factor binding site
    #[serde(rename = "TF_binding_site")]
    TfBindingSite,
}

/// Variant effect description
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum VariantEffect {
    #[serde(rename = "coding_sequence_variant")]
    CodingSequenceVariant,
    #[serde(rename = "coding_transcript_intron_variant")]
    CodingTranscriptIntronVariant,
    #[serde(rename = "coding_transcript_variant")]
    CodingTranscriptVariant,
    #[serde(rename = "copy_number_change")]
    CopyNumberChange,
    #[serde(rename = "direct_tandem_duplication")]
    DirectTandemDuplicatin,
    #[serde(rename = "downstream_gene_variant")]
    DownstreamGeneVariant,
    #[serde(rename = "exon_loss_variant")]
    ExonLossVariant,
    #[serde(rename = "feature_truncation")]
    FeatureTruncation,
    #[serde(rename = "5_prime_UTR_exon_variant")]
    FivePrimeUtrExonVariant,
    #[serde(rename = "5_prime_UTR_intron_variant")]
    FivePrimeUtrIntronVariant,
    #[serde(rename = "5_prime_UTR_truncation")]
    FivePrimeUtrTruncation,
    #[serde(rename = "frameshift_truncation")]
    FrameshiftTruncation,
    #[serde(rename = "insertion")]
    Insertion,
    #[serde(rename = "intron_variant")]
    IntronVariant,
    #[serde(rename = "inversion")]
    Inversion,
    #[serde(rename = "mobile_element_deletion")]
    MobileElementDeletion,
    #[serde(rename = "mobile_element_insertion")]
    MobileElementInsertion,
    #[serde(rename = "non_coding_transcript_exon_variant")]
    NonCodingTranscriptExonVariant,
    #[serde(rename = "non_coding_transcript_intron_variant")]
    NonCodingTranscriptIntronVariant,
    #[serde(rename = "non_coding_transcript_variant")]
    NonCodingTranscriptVariant,
    #[serde(rename = "sequence_variant")]
    SequenceVariant,
    #[serde(rename = "start_lost")]
    StartLost,
    #[serde(rename = "stop_lost")]
    StopLost,
    #[serde(rename = "structural_variant")]
    StructuralVariant,
    #[serde(rename = "3_prime_UTR_exon_variant")]
    ThreePrimeUtrExonVariant,
    #[serde(rename = "3_prime_UTR_intron_variant")]
    ThreePrimeUtrIntronVariant,
    #[serde(rename = "3_prime_UTR_truncation")]
    ThreePrimeUtrTruncation,
    #[serde(rename = "transcript_ablation")]
    TranscriptAblation,
    #[serde(rename = "transcript_amplification")]
    TranscriptAmplification,
    #[serde(rename = "translocation")]
    Translocation,
    #[serde(rename = "upstream_variant")]
    UpstreamVariant,
}

/// VISTA validation
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
#[serde(rename_all = "snake_case")]
pub enum VistaValidation {
    /// Overlap with any VISTA enhancer
    Any,
    /// Overlap with positive VISTSA enhancer
    Positive,
    /// Overlap with negative VISTA enhancer
    Negative,
}

#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct RegulatoryCustomConfig {
    /// The selected cell types
    pub cell_types: Vec<String>,
    /// The selected element type
    pub element_types: Vec<String>,
    /// Alternatively, overlapping with interaction in cell type (includes all elements)
    pub overlaps_interaction: bool,
}

impl RegulatoryCustomConfig {
    pub fn new() -> Self {
        RegulatoryCustomConfig {
            cell_types: vec![],
            element_types: vec![],
            overlaps_interaction: false,
        }
    }
}

/// Enum for recessive mode
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
#[serde(rename_all = "kebab-case")]
pub enum RecessiveMode {
    /// Recessive mode of inheritance
    Recessive,
    /// Compound recessive mode of inheritance
    CompoundRecessive,
}

/// Options for TAD sets to use.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum TadSet {
    /// hESC  cell line
    #[serde(rename = "hesc")]
    Hesc,
    /// IMR90 cell line
    #[serde(rename = "im390")]
    Imr90,
}

/// Define rule to apply to a given sub set of structural variants for matching a genotype.
///
/// See documentation of VarFish Server for full documentation.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct GenotypeCriteria {
    /// The genotype to evaluate for
    pub genotype: GenotypeChoice,

    /// Select SV sub types that this applies to
    pub select_sv_sub_type: Vec<SvSubType>,
    /// Select minimal size of SV to apply to (ignored for BND and INS)
    pub select_sv_min_size: Option<i32>,
    /// Select maximal size of SV to apply to (ignored for BND and INS)
    pub select_sv_max_size: Option<i32>,

    /// The FORMAT/GT field should be one of, unless None
    pub gt_one_of: Option<Vec<String>>,

    /// Minimal genotype quality as returned by caller
    pub min_gq: Option<f32>,
    /// Minimal number of total paired-end reads covering the SV
    pub min_pr_cov: Option<i32>,
    /// Maximal number of total paired-end reads covering the SV
    pub max_pr_cov: Option<i32>,
    /// Minimal number of reference paired-end reads covering the SV
    pub min_pr_ref: Option<i32>,
    /// Maximal number of reference paired-end reads covering the SV
    pub max_pr_ref: Option<i32>,
    /// Minimal number of variant paired-end reads covering the SV
    pub min_pr_var: Option<i32>,
    /// Maximal number of variant paired-end reads covering the SV
    pub max_pr_var: Option<i32>,
    /// Minimal number of total split reads covering the SV
    pub min_sr_cov: Option<i32>,
    /// Maximal number of total split reads covering the SV
    pub max_sr_cov: Option<i32>,
    /// Minimal number of reference split reads covering the SV
    pub min_sr_ref: Option<i32>,
    /// Maximal number of reference split reads covering the SV
    pub max_sr_ref: Option<i32>,
    /// Minimal number of variant split reads covering the SV
    pub min_sr_var: Option<i32>,
    /// Maximal number of variant split reads covering the SV
    pub max_sr_var: Option<i32>,
    /// Minimal sum of total paired-end/split read coverage
    pub min_srpr_cov: Option<i32>,
    /// Maximal sum of total paired-end/split read coverage
    pub max_srpr_cov: Option<i32>,
    /// Minimal sum of reference paired-end/split read coverage
    pub min_srpr_ref: Option<i32>,
    /// Maximal sum of reference paired-end/split read coverage
    pub max_srpr_ref: Option<i32>,
    /// Minimal sum of variant paired-end/split read coverage
    pub min_srpr_var: Option<i32>,
    /// Maximal sum of variant paired-end/split read coverage
    pub max_srpr_var: Option<i32>,
    /// Minimal coverage deviation
    pub min_rd_dev: Option<f32>,
    /// Maximal coverage deviation
    pub max_rd_dev: Option<f32>,
    /// Minimal average mapping quality
    pub min_amq: Option<f32>,
    /// Maximal average mapping quality
    pub max_amq: Option<f32>,

    /// An optional comment
    pub comment: Option<String>,
}

impl GenotypeCriteria {
    /// Return new, empty GenotypeCriteria for the given genotype.
    pub fn new(genotype: GenotypeChoice) -> Self {
        GenotypeCriteria {
            genotype,
            select_sv_sub_type: vec![],
            select_sv_min_size: None,
            select_sv_max_size: None,
            gt_one_of: None,
            min_gq: None,
            min_pr_cov: None,
            max_pr_cov: None,
            min_pr_ref: None,
            max_pr_ref: None,
            min_pr_var: None,
            max_pr_var: None,
            min_sr_cov: None,
            max_sr_cov: None,
            min_sr_ref: None,
            max_sr_ref: None,
            min_sr_var: None,
            max_sr_var: None,
            min_srpr_cov: None,
            max_srpr_cov: None,
            min_srpr_ref: None,
            max_srpr_ref: None,
            min_srpr_var: None,
            max_srpr_var: None,
            min_rd_dev: None,
            max_rd_dev: None,
            min_amq: None,
            max_amq: None,
            comment: None,
        }
    }
}

/// Define a query for structural variants from a case.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct CaseQuery {
    /// The transcript database to use
    pub database: Database,

    /// Whether to enable SVDB overlap queries with DGV.
    pub svdb_dgv_enabled: bool,
    /// The minimal reciprocal overlap for querying DGV.
    pub svdb_dgv_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying DGV.
    pub svdb_dgv_max_carriers: Option<u32>,
    /// Whether to enable SVDB overlap queries with DGV gold standard.
    pub svdb_dgv_gs_enabled: bool,
    /// The minimal reciprocal overlap for querying DGV gold standard.
    pub svdb_dgv_gs_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying DGV gold standard.
    pub svdb_dgv_gs_max_carriers: Option<u32>,
    /// Whether to enable SVDB overlap queries with gnomAD.
    pub svdb_gnomad_enabled: bool,
    /// The minimal reciprocal overlap for querying gnomAD.
    pub svdb_gnomad_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying gnomAD.
    pub svdb_gnomad_max_carriers: Option<u32>,
    /// Whether to enable SVDB overlap queries with ExAC.
    pub svdb_exac_enabled: bool,
    /// The minimal reciprocal overlap for querying ExAC.
    pub svdb_exac_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying ExAC.
    pub svdb_exac_max_carriers: Option<u32>,
    /// Whether to enable SVDB overlap queries with dbVar.
    pub svdb_dbvar_enabled: bool,
    /// The minimal reciprocal overlap for querying dbVar.
    pub svdb_dbvar_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying dbVar.
    pub svdb_dbvar_max_carriers: Option<u32>,
    /// Whether to enable SVDB overlap queries with Thousand Genomes Project.
    pub svdb_g1k_enabled: bool,
    /// The minimal reciprocal overlap for querying Thousand Genomes Project.
    pub svdb_g1k_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying Thousand Genomes Project.
    pub svdb_g1k_max_alleles: Option<u32>,
    /// Whether to enable SVDB overlap queries with in-house DB.
    pub svdb_inhouse_enabled: bool,
    /// The minimal reciprocal overlap for querying in-house DB.
    pub svdb_inhouse_min_overlap: Option<f32>,
    /// The maximal number of alleles for querying in-house DB.
    pub svdb_inhouse_max_carriers: Option<u32>,

    /// The minimal SV size to consider.
    pub sv_size_min: Option<i32>,
    /// The maximal SV size to consider.
    pub sv_size_max: Option<i32>,

    /// The SV types to consider.
    pub sv_types: Vec<SvType>,
    /// The SV subtypes to consider.
    pub sv_sub_types: Vec<SvSubType>,

    /// List of genes to require.
    pub gene_allowlist: Option<Vec<String>>,
    /// Genomic region to limit consideration to.
    pub genomic_region: Option<Vec<GenomicRegion>>,

    /// Regulatory region padding to use.
    pub regulatory_overlap: i32,
    /// Regulatory features to select from ENSEMBL.
    pub regulatory_ensembl_features: Option<Vec<EnsemblRegulatoryFeature>>,
    /// VISTA enhancer validation results.
    pub regulatory_vista_validation: Option<VistaValidation>,

    /// Custom regulatory maps configuration.
    pub regulatory_custom_configs: Vec<RegulatoryCustomConfig>,

    /// Name of the TAD set to use for annotation, if any.
    pub tad_set: Option<TadSet>,

    /// Genotype choices
    pub genotype: HashMap<String, GenotypeChoice>,
    /// Criteria for filtering CNVs.
    pub genotype_criteria: Vec<GenotypeCriteria>,

    /// The mode for recessive inheritance.
    pub recessive_mode: Option<RecessiveMode>,
    /// The index to use for recessive inheritance.
    pub recessive_index: Option<String>,
}

impl CaseQuery {
    pub fn new(database: Database) -> Self {
        CaseQuery {
            database,
            svdb_dgv_enabled: false,
            svdb_dgv_min_overlap: None,
            svdb_dgv_max_carriers: None,
            svdb_dgv_gs_enabled: false,
            svdb_dgv_gs_min_overlap: None,
            svdb_dgv_gs_max_carriers: None,
            svdb_gnomad_enabled: false,
            svdb_gnomad_min_overlap: None,
            svdb_gnomad_max_carriers: None,
            svdb_exac_enabled: false,
            svdb_exac_min_overlap: None,
            svdb_exac_max_carriers: None,
            svdb_dbvar_enabled: false,
            svdb_dbvar_min_overlap: None,
            svdb_dbvar_max_carriers: None,
            svdb_g1k_enabled: false,
            svdb_g1k_min_overlap: None,
            svdb_g1k_max_alleles: None,
            svdb_inhouse_enabled: false,
            svdb_inhouse_min_overlap: None,
            svdb_inhouse_max_carriers: None,
            sv_size_min: None,
            sv_size_max: None,
            sv_types: vec![],
            sv_sub_types: vec![],
            gene_allowlist: None,
            genomic_region: None,
            regulatory_overlap: 100,
            regulatory_ensembl_features: None,
            regulatory_vista_validation: None,
            regulatory_custom_configs: vec![],
            tad_set: None,
            genotype: HashMap::new(),
            genotype_criteria: vec![],
            recessive_mode: None,
            recessive_index: None,
        }
    }
}

/// Strand orientation of
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub enum StrandOrientation {
    #[serde(rename = "3to3")]
    ThreeToThree,
    #[serde(rename = "5to5")]
    FiveToFive,
    #[serde(rename = "3to5")]
    ThreeToFive,
    #[serde(rename = "5to3")]
    FiveToThree,
}

/// Information on the call as combined by the annotator.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct CallInfo {
    /// The genotype, if applicable, e.g., "0/1", "./1", "."
    pub genotype: Option<String>,
    /// Genotype quality score, if applicable
    pub quality: Option<f32>,
    /// Paired-end coverage, if applicable
    pub paired_end_cov: Option<u32>,
    /// Paired-end variant support, if applicable
    pub paired_end_var: Option<u32>,
    /// Split-read coverage, if applicable
    pub split_read_cov: Option<u32>,
    /// Split-read variant support, if applicable
    pub split_read_var: Option<u32>,
    /// Integer copy number estimate, if applicable
    pub copy_number: Option<u32>,
    /// Average normalized coverage, if applicable
    pub average_normalized_cov: Option<f32>,
    /// Number of buckets/targets supporting the CNV call, if applicable
    pub point_count: Option<u32>,
}

/// Definition of a structural variant with per-sample genotype calls.
///
/// This uses a subset/specialization of what is described by the VCF standard
/// for the purpose of running SV queries in `varfish-server-worker`.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct StructuralVariant {
    /// Chromosome name
    pub chrom: String,
    /// 1-based start position of the variant (or position on first chromosome for
    /// break-ends)
    pub pos: u32,
    /// Type of the structural variant
    pub sv_type: SvType,
    /// Sub type of the structural variant
    pub sv_sub_type: SvSubType,
    /// Potentially the second involved chromosome
    pub chrom2: Option<String>,
    /// End position (position on second chromosome for break-ends)
    pub end: u32,
    /// The strand orientation of the structural variant, if applicable.
    pub strand_orientation: Option<StrandOrientation>,

    /// Mapping of sample to genotype information for the SV.
    pub call_info: HashMap<String, CallInfo>,
}

impl StructuralVariant {
    /// Return the size of the structural variant.
    ///
    /// Size is not applicable for insertions and break-ends, so `None` is returned
    /// in this case.
    pub fn size(&self) -> Option<u32> {
        if self.sv_type == SvType::Ins
            || self.sv_type == SvType::Bnd
            || self.sv_sub_type.is_ins()
            || self.sv_sub_type == SvSubType::Bnd
        {
            None
        } else {
            Some(self.end.saturating_sub(self.pos) + 1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use serde_test::{assert_tokens, Token};

    #[test]
    fn test_range_smoke() {
        assert_eq!(Range::new(1, 2), Range { start: 1, end: 2 });
    }

    #[test]
    fn test_range_serde_smoke() {
        assert_tokens(
            &Range { start: 1, end: 2 },
            &[
                Token::Struct {
                    name: "Range",
                    len: 2,
                },
                Token::Str("start"),
                Token::U32(1),
                Token::Str("end"),
                Token::U32(2),
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_genomic_region_smoke() {
        assert_eq!(
            GenomicRegion::new("chr1", 1, 2),
            GenomicRegion {
                chrom: "chr1".to_owned(),
                range: Some(Range::new(1, 2))
            }
        );
        assert_eq!(
            GenomicRegion::whole_chrom("chr1"),
            GenomicRegion {
                chrom: "chr1".to_owned(),
                range: None,
            }
        );
    }

    #[test]
    fn test_genomic_region_serde_smoke() {
        assert_tokens(
            &Database::Refseq,
            &[Token::UnitVariant {
                name: "Database",
                variant: "REFSEQ",
            }],
        );
    }

    #[test]
    fn test_database_smoke() {
        assert_tokens(
            &Database::Refseq,
            &[Token::UnitVariant {
                name: "Database",
                variant: "REFSEQ",
            }],
        );
    }

    #[test]
    fn test_sv_type_serde_smoke() {
        assert_tokens(
            &SvType::Del,
            &[Token::UnitVariant {
                name: "SvType",
                variant: "DEL",
            }],
        );
    }

    #[test]
    fn test_sv_sub_type_serde_smoke() {
        assert_tokens(
            &SvSubType::DelMeL1,
            &[Token::UnitVariant {
                name: "SvSubType",
                variant: "DEL:ME:L1",
            }],
        );
    }

    #[test]
    fn test_genotype_choice_serde_smoke() {
        assert_tokens(
            &GenotypeChoice::Het,
            &[Token::UnitVariant {
                name: "GenotypeChoice",
                variant: "het",
            }],
        );
    }

    #[test]
    fn test_ensembl_regulatory_feature_serde_smoke() {
        assert_tokens(
            &EnsemblRegulatoryFeature::Promoter,
            &[Token::UnitVariant {
                name: "EnsemblRegulatoryFeature",
                variant: "promoter",
            }],
        );
    }

    #[test]
    fn test_variant_effect_serde_smoke() {
        assert_tokens(
            &VariantEffect::FrameshiftTruncation,
            &[Token::UnitVariant {
                name: "VariantEffect",
                variant: "frameshift_truncation",
            }],
        );
    }

    #[test]
    fn test_vista_validation_serde_smoke() {
        assert_tokens(
            &VistaValidation::Positive,
            &[Token::UnitVariant {
                name: "VistaValidation",
                variant: "positive",
            }],
        );
    }

    #[test]
    fn test_regulatory_custom_config_smoke() {
        assert_eq!(
            RegulatoryCustomConfig::new(),
            RegulatoryCustomConfig {
                cell_types: vec![],
                element_types: vec![],
                overlaps_interaction: false,
            }
        )
    }

    #[test]
    fn test_regulatory_custom_config_serde_smoke() {
        assert_tokens(
            &RegulatoryCustomConfig::new(),
            &[
                Token::Struct {
                    name: "RegulatoryCustomConfig",
                    len: 3,
                },
                Token::Str("cell_types"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("element_types"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("overlaps_interaction"),
                Token::Bool(false),
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_recessive_mode_serde_smoke() {
        assert_tokens(
            &RecessiveMode::CompoundRecessive,
            &[Token::UnitVariant {
                name: "RecessiveMode",
                variant: "compound-recessive",
            }],
        );
    }

    #[test]
    fn test_tad_set_serde_smoke() {
        assert_tokens(
            &TadSet::Hesc,
            &[Token::UnitVariant {
                name: "TadSet",
                variant: "hesc",
            }],
        );
    }

    #[test]
    fn test_genotype_criteria_serde_smoke() {
        let crit = GenotypeCriteria::new(GenotypeChoice::Het);
        assert_tokens(
            &crit,
            &[
                Token::Struct {
                    name: "GenotypeCriteria",
                    len: 29,
                },
                Token::Str("genotype"),
                Token::UnitVariant {
                    name: "GenotypeChoice",
                    variant: "het",
                },
                Token::Str("select_sv_sub_type"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("select_sv_min_size"),
                Token::None,
                Token::Str("select_sv_max_size"),
                Token::None,
                Token::Str("gt_one_of"),
                Token::None,
                Token::Str("min_gq"),
                Token::None,
                Token::Str("min_pr_cov"),
                Token::None,
                Token::Str("max_pr_cov"),
                Token::None,
                Token::Str("min_pr_ref"),
                Token::None,
                Token::Str("max_pr_ref"),
                Token::None,
                Token::Str("min_pr_var"),
                Token::None,
                Token::Str("max_pr_var"),
                Token::None,
                Token::Str("min_sr_cov"),
                Token::None,
                Token::Str("max_sr_cov"),
                Token::None,
                Token::Str("min_sr_ref"),
                Token::None,
                Token::Str("max_sr_ref"),
                Token::None,
                Token::Str("min_sr_var"),
                Token::None,
                Token::Str("max_sr_var"),
                Token::None,
                Token::Str("min_srpr_cov"),
                Token::None,
                Token::Str("max_srpr_cov"),
                Token::None,
                Token::Str("min_srpr_ref"),
                Token::None,
                Token::Str("max_srpr_ref"),
                Token::None,
                Token::Str("min_srpr_var"),
                Token::None,
                Token::Str("max_srpr_var"),
                Token::None,
                Token::Str("min_rd_dev"),
                Token::None,
                Token::Str("max_rd_dev"),
                Token::None,
                Token::Str("min_amq"),
                Token::None,
                Token::Str("max_amq"),
                Token::None,
                Token::Str("comment"),
                Token::None,
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_case_query_serde_smoke() {
        let query = CaseQuery::new(Database::Refseq);
        assert_tokens(
            &query,
            &[
                Token::Struct {
                    name: "CaseQuery",
                    len: 37,
                },
                Token::Str("database"),
                Token::UnitVariant {
                    name: "Database",
                    variant: "REFSEQ",
                },
                Token::Str("svdb_dgv_enabled"),
                Token::Bool(false),
                Token::Str("svdb_dgv_min_overlap"),
                Token::None,
                Token::Str("svdb_dgv_max_carriers"),
                Token::None,
                Token::Str("svdb_dgv_gs_enabled"),
                Token::Bool(false),
                Token::Str("svdb_dgv_gs_min_overlap"),
                Token::None,
                Token::Str("svdb_dgv_gs_max_carriers"),
                Token::None,
                Token::Str("svdb_gnomad_enabled"),
                Token::Bool(false),
                Token::Str("svdb_gnomad_min_overlap"),
                Token::None,
                Token::Str("svdb_gnomad_max_carriers"),
                Token::None,
                Token::Str("svdb_exac_enabled"),
                Token::Bool(false),
                Token::Str("svdb_exac_min_overlap"),
                Token::None,
                Token::Str("svdb_exac_max_carriers"),
                Token::None,
                Token::Str("svdb_dbvar_enabled"),
                Token::Bool(false),
                Token::Str("svdb_dbvar_min_overlap"),
                Token::None,
                Token::Str("svdb_dbvar_max_carriers"),
                Token::None,
                Token::Str("svdb_g1k_enabled"),
                Token::Bool(false),
                Token::Str("svdb_g1k_min_overlap"),
                Token::None,
                Token::Str("svdb_g1k_max_alleles"),
                Token::None,
                Token::Str("svdb_inhouse_enabled"),
                Token::Bool(false),
                Token::Str("svdb_inhouse_min_overlap"),
                Token::None,
                Token::Str("svdb_inhouse_max_carriers"),
                Token::None,
                Token::Str("sv_size_min"),
                Token::None,
                Token::Str("sv_size_max"),
                Token::None,
                Token::Str("sv_types"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("sv_sub_types"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("gene_allowlist"),
                Token::None,
                Token::Str("genomic_region"),
                Token::None,
                Token::Str("regulatory_overlap"),
                Token::I32(100),
                Token::Str("regulatory_ensembl_features"),
                Token::None,
                Token::Str("regulatory_vista_validation"),
                Token::None,
                Token::Str("regulatory_custom_configs"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("tad_set"),
                Token::None,
                Token::Str("genotype"),
                Token::Map { len: Some(0) },
                Token::MapEnd,
                Token::Str("genotype_criteria"),
                Token::Seq { len: Some(0) },
                Token::SeqEnd,
                Token::Str("recessive_mode"),
                Token::None,
                Token::Str("recessive_index"),
                Token::None,
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_strand_orientation_serde_smoke() {
        assert_tokens(
            &StrandOrientation::ThreeToFive,
            &[Token::UnitVariant {
                name: "StrandOrientation",
                variant: "3to5",
            }],
        );
    }

    #[test]
    fn test_call_info_serde_smoke() {
        let info = CallInfo {
            genotype: Some("0/1".to_owned()),
            quality: Some(10.0),
            paired_end_cov: Some(10),
            paired_end_var: Some(10),
            split_read_cov: Some(10),
            split_read_var: Some(10),
            copy_number: Some(1),
            average_normalized_cov: Some(0.491),
            point_count: Some(5),
        };
        assert_tokens(
            &info,
            &[
                Token::Struct {
                    name: "CallInfo",
                    len: 9,
                },
                Token::Str("genotype"),
                Token::Some,
                Token::Str("0/1"),
                Token::Str("quality"),
                Token::Some,
                Token::F32(10.0),
                Token::Str("paired_end_cov"),
                Token::Some,
                Token::U32(10),
                Token::Str("paired_end_var"),
                Token::Some,
                Token::U32(10),
                Token::Str("split_read_cov"),
                Token::Some,
                Token::U32(10),
                Token::Str("split_read_var"),
                Token::Some,
                Token::U32(10),
                Token::Str("copy_number"),
                Token::Some,
                Token::U32(1),
                Token::Str("average_normalized_cov"),
                Token::Some,
                Token::F32(0.491),
                Token::Str("point_count"),
                Token::Some,
                Token::U32(5),
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_structural_variant_size_linear() {
        let sv = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::DelMeL1,
            chrom2: None,
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: HashMap::new(),
        };
        assert_eq!(sv.size().unwrap(), 101);
    }

    #[test]
    fn test_structural_variant_size_ins() {
        let sv = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Ins,
            sv_sub_type: SvSubType::Ins,
            chrom2: None,
            end: 100,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: HashMap::new(),
        };
        assert!(sv.size().is_none());
    }

    #[test]
    fn test_structural_variant_size_bnd() {
        let sv = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 100,
            sv_type: SvType::Bnd,
            sv_sub_type: SvSubType::Bnd,
            chrom2: Some("chr2".to_owned()),
            end: 200,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: HashMap::new(),
        };
        assert!(sv.size().is_none());
    }

    #[test]
    fn test_structural_variant_serde_smoke() {
        let sv = StructuralVariant {
            chrom: "chr1".to_owned(),
            pos: 123,
            sv_type: SvType::Del,
            sv_sub_type: SvSubType::DelMeL1,
            chrom2: None,
            end: 245,
            strand_orientation: Some(StrandOrientation::ThreeToFive),
            call_info: HashMap::new(),
        };
        assert_tokens(
            &sv,
            &[
                Token::Struct {
                    name: "StructuralVariant",
                    len: 8,
                },
                Token::Str("chrom"),
                Token::Str("chr1"),
                Token::Str("pos"),
                Token::U32(123),
                Token::Str("sv_type"),
                Token::UnitVariant {
                    name: "SvType",
                    variant: "DEL",
                },
                Token::Str("sv_sub_type"),
                Token::UnitVariant {
                    name: "SvSubType",
                    variant: "DEL:ME:L1",
                },
                Token::Str("chrom2"),
                Token::None,
                Token::Str("end"),
                Token::U32(245),
                Token::Str("strand_orientation"),
                Token::Some,
                Token::UnitVariant {
                    name: "StrandOrientation",
                    variant: "3to5",
                },
                Token::Str("call_info"),
                Token::Map { len: Some(0) },
                Token::MapEnd,
                Token::StructEnd,
            ],
        );
    }

    fn test_sv_sub_type_is_ins() {
        assert_eq!(SvSubType::Del.is_ins(), false);
        assert_eq!(SvSubType::DelMe.is_ins(), false);
        assert_eq!(SvSubType::DelMeSva.is_ins(), false);
        assert_eq!(SvSubType::DelMeL1.is_ins(), false);
        assert_eq!(SvSubType::DelMeAlu.is_ins(), false);
        assert_eq!(SvSubType::Dup.is_ins(), false);
        assert_eq!(SvSubType::DupTandem.is_ins(), false);
        assert_eq!(SvSubType::Inv.is_ins(), false);
        assert_eq!(SvSubType::Ins.is_ins(), true);
        assert_eq!(SvSubType::InsMe.is_ins(), true);
        assert_eq!(SvSubType::InsMeSva.is_ins(), true);
        assert_eq!(SvSubType::InsMeL1.is_ins(), true);
        assert_eq!(SvSubType::InsMeAlu.is_ins(), true);
        assert_eq!(SvSubType::Bnd.is_ins(), false);
        assert_eq!(SvSubType::Cnv.is_ins(), false);
    }

    fn test_sv_sub_type_is_del() {
        assert_eq!(SvSubType::Del.is_ins(), true);
        assert_eq!(SvSubType::DelMe.is_ins(), true);
        assert_eq!(SvSubType::DelMeSva.is_ins(), true);
        assert_eq!(SvSubType::DelMeL1.is_ins(), true);
        assert_eq!(SvSubType::DelMeAlu.is_ins(), true);
        assert_eq!(SvSubType::Dup.is_ins(), false);
        assert_eq!(SvSubType::DupTandem.is_ins(), false);
        assert_eq!(SvSubType::Inv.is_ins(), false);
        assert_eq!(SvSubType::Ins.is_ins(), false);
        assert_eq!(SvSubType::InsMe.is_ins(), false);
        assert_eq!(SvSubType::InsMeSva.is_ins(), false);
        assert_eq!(SvSubType::InsMeL1.is_ins(), false);
        assert_eq!(SvSubType::InsMeAlu.is_ins(), false);
        assert_eq!(SvSubType::Bnd.is_ins(), false);
        assert_eq!(SvSubType::Cnv.is_ins(), false);
    }
}
