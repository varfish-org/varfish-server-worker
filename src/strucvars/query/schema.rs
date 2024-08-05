//! Supporting code for SV query definition.

use crate::common::{genotype_to_string, TadSet};
use indexmap::IndexMap;
use mehari::annotate::strucvars::{
    bnd::Breakend, csq::interface::StrandOrientation, PeOrientation,
};
use noodles::vcf::{self, variant::record::AlternateBases as _};
use regex::Regex;
use serde::{Deserialize, Deserializer, Serialize};
use strum_macros::{Display, EnumIter, EnumString};

use super::masked::MaskedBreakpointCount;

/// Range with 1-based positions
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct Range {
    pub start: i32,
    pub end: i32,
}

impl Range {
    pub fn new(start: i32, end: i32) -> Self {
        Range { start, end }
    }
}

/// Chromosomal interval.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub struct ChromRange {
    /// Chromosome name.
    pub chromosome: String,
    /// 0-based begin position.
    pub begin: i32,
    /// 0-based end position.
    pub end: i32,
}

/// Genomic region
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct GenomicRegion {
    pub chrom: String,
    pub range: Option<Range>,
}

impl GenomicRegion {
    pub fn new(chrom: &str, start: i32, end: i32) -> Self {
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

/// Encode the type of an SV
#[derive(
    Serialize,
    Deserialize,
    EnumIter,
    PartialEq,
    Eq,
    Ord,
    PartialOrd,
    Hash,
    Debug,
    Clone,
    Copy,
    Default,
)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum SvType {
    /// Deletion
    #[default]
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

impl std::str::FromStr for SvType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use SvType::*;
        match s {
            "DEL" => Ok(Del),
            "DUP" => Ok(Dup),
            "INV" => Ok(Inv),
            "INS" => Ok(Ins),
            "BND" => Ok(Bnd),
            "CNV" => Ok(Cnv),
            _ => Err(anyhow::anyhow!("invalid SV type: {}", s)),
        }
    }
}

impl SvType {
    pub fn vec_all() -> Vec<SvType> {
        use SvType::*;
        vec![Del, Dup, Inv, Ins, Bnd, Cnv]
    }

    pub fn is_compatible(&self, other: SvType) -> bool {
        use SvType::*;
        matches!(
            (self, other),
            (Del, Del)
                | (Dup, Dup)
                | (Inv, Inv)
                | (Ins, Ins)
                | (Bnd, Bnd)
                | (Cnv, Cnv)
                | (Del, Cnv)
                | (Cnv, Del)
                | (Dup, Cnv)
                | (Cnv, Dup)
        )
    }
}

/// Structural variant sub type
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy, Default)]
pub enum SvSubType {
    /// Deletion
    #[serde(rename = "DEL")]
    #[default]
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
    /// Return vector with all SV sub types
    pub fn vec_all() -> Vec<SvSubType> {
        use SvSubType::*;
        vec![
            Del, DelMe, DelMeSva, DelMeL1, DelMeAlu, Dup, DupTandem, Inv, Ins, InsMe, InsMeSva,
            InsMeL1, InsMeAlu, Bnd, Cnv,
        ]
    }

    /// Return vector with del/dup/CNV SV sub types
    pub fn vec_cnv() -> Vec<SvSubType> {
        use SvSubType::*;
        vec![Del, DelMe, DelMeSva, DelMeL1, DelMeAlu, Dup, DupTandem, Cnv]
    }

    /// Return vector with deletion SV sub types
    pub fn vec_del() -> Vec<SvSubType> {
        use SvSubType::*;
        vec![Del, DelMe, DelMeSva, DelMeL1, DelMeAlu]
    }

    /// Return vector with insertion SV sub types
    pub fn vec_ins() -> Vec<SvSubType> {
        use SvSubType::*;
        vec![Ins, InsMe, InsMeSva, InsMeL1, InsMeAlu]
    }

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

/// Enumeration for effect on transcript.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum TranscriptEffect {
    /// Affects the full transcript.
    TranscriptVariant,
    /// An exon is affected by the SV.
    ExonVariant,
    /// The splice region is affected by the SV.
    SpliceRegionVariant,
    /// The intron is affected by the SV.
    IntronVariant,
    /// The upstream region of the transcript is affected.
    UpstreamVariant,
    /// The downstream region of the transcript is affected.
    DownstreamVariant,
    /// Only intergenic regions is affected,
    IntergenicVariant,
}

impl TranscriptEffect {
    /// Return vector with all transcript effects.
    pub fn vec_all() -> Vec<TranscriptEffect> {
        vec![
            TranscriptEffect::TranscriptVariant,
            TranscriptEffect::ExonVariant,
            TranscriptEffect::SpliceRegionVariant,
            TranscriptEffect::IntronVariant,
            TranscriptEffect::UpstreamVariant,
            TranscriptEffect::DownstreamVariant,
            TranscriptEffect::IntergenicVariant,
        ]
    }
}

/// Genotype choice for filter dropdowns
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy)]
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

/// Enum for effective genotypes.
///
/// This is very similar to `GenotypeChoice` but has no `Any`.  Also, note that
/// the order of the values corresponds to the priority when the effective
/// genotype is returned.
#[derive(Serialize, Deserialize, PartialEq, PartialOrd, Debug, Clone, Copy)]
#[serde(rename_all = "kebab-case")]
pub enum Genotype {
    /// Homozygous alternative genotype
    Hom,
    /// Heterozygous genotype
    Het,
    /// Not wild-type genotype
    Variant,
    /// Reference genotype
    Ref,
    /// Not variant genotype
    NonVariant,
}

/// ENSEMBL regulatory feature
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy)]
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
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy)]
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
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy)]
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
    /// Alternatively, overlapping with interaction in cell type (includes all
    /// elements)
    pub overlaps_interaction: bool,
}

impl Default for RegulatoryCustomConfig {
    fn default() -> Self {
        Self::new()
    }
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
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Copy)]
#[serde(rename_all = "kebab-case")]
pub enum RecessiveMode {
    /// Recessive mode of inheritance
    Recessive,
    /// Compound recessive mode of inheritance
    CompoundRecessive,
}

fn default_as_true() -> bool {
    true
}

/// Define rule to apply to a given sub set of structural variants for matching
/// a genotype.
///
/// See documentation of VarFish Server for full documentation.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct GenotypeCriteria {
    /// The genotype to evaluate for
    pub genotype: GenotypeChoice,

    /// Select SV sub types that this applies to
    pub select_sv_sub_type: Vec<SvSubType>,
    /// Select minimal size of SV to apply to (ignored for BND and INS)
    pub select_sv_min_size: Option<u32>,
    /// Select maximal size of SV to apply to (ignored for BND and INS)
    pub select_sv_max_size: Option<u32>,

    // Maximal number of ends/breakpoints within segmental duplications
    pub max_brk_segdup: Option<u32>,
    // Maximal number of ends/breakpoints within repeat-masked sequence
    pub max_brk_repeat: Option<u32>,
    // Maximal number of ends/breakpoints within segmental duplications or repeat-masked sequence
    pub max_brk_segduprepeat: Option<u32>,

    /// The FORMAT/GT field should be one of, unless None
    pub gt_one_of: Option<Vec<String>>,

    /// Minimal genotype quality as returned by caller
    pub min_gq: Option<f32>,
    /// Minimal number of total paired-end reads covering the SV
    pub min_pr_cov: Option<u32>,
    /// Maximal number of total paired-end reads covering the SV
    pub max_pr_cov: Option<u32>,
    /// Minimal number of reference paired-end reads covering the SV
    pub min_pr_ref: Option<u32>,
    /// Maximal number of reference paired-end reads covering the SV
    pub max_pr_ref: Option<u32>,
    /// Minimal number of variant paired-end reads covering the SV
    pub min_pr_var: Option<u32>,
    /// Maximal number of variant paired-end reads covering the SV
    pub max_pr_var: Option<u32>,
    /// Minimal allelic balaance of paired-end reads covering the SV
    pub min_pr_ab: Option<f32>,
    /// Maximal allelic balance of paired-end reads covering the SV
    pub max_pr_ab: Option<f32>,
    /// Minimal number of total split reads covering the SV
    pub min_sr_cov: Option<u32>,
    /// Maximal number of total split reads covering the SV
    pub max_sr_cov: Option<u32>,
    /// Minimal number of reference split reads covering the SV
    pub min_sr_ref: Option<u32>,
    /// Maximal number of reference split reads covering the SV
    pub max_sr_ref: Option<u32>,
    /// Minimal number of variant split reads covering the SV
    pub min_sr_var: Option<u32>,
    /// Maximal number of variant split reads covering the SV
    pub max_sr_var: Option<u32>,
    /// Minimal allelic balance of split reads covering the SV
    pub min_sr_ab: Option<f32>,
    /// Maximal allelic balance of split reads covering the SV
    pub max_sr_ab: Option<f32>,
    /// Minimal sum of total paired-end/split read coverage
    pub min_srpr_cov: Option<u32>,
    /// Maximal sum of total paired-end/split read coverage
    pub max_srpr_cov: Option<u32>,
    /// Minimal sum of reference paired-end/split read coverage
    pub min_srpr_ref: Option<u32>,
    /// Maximal sum of reference paired-end/split read coverage
    pub max_srpr_ref: Option<u32>,
    /// Minimal sum of variant paired-end/split read coverage
    pub min_srpr_var: Option<u32>,
    /// Maximal sum of variant paired-end/split read coverage
    pub max_srpr_var: Option<u32>,
    /// Minimal allelic balance of paired-end/split read coverage
    pub min_srpr_ab: Option<f32>,
    /// Maximal allelic balance of paired-end/split read coverage
    pub max_srpr_ab: Option<f32>,
    /// Minimal coverage deviation
    pub min_rd_dev: Option<f32>,
    /// Maximal coverage deviation
    pub max_rd_dev: Option<f32>,
    /// Minimal average mapping quality
    pub min_amq: Option<f32>,
    /// Maximal average mapping quality
    pub max_amq: Option<f32>,

    /// Whether missing genotype call leads to filter out variant
    #[serde(default = "default_as_true")]
    pub missing_gt_ok: bool,
    /// Whether missing genotype quality information leads filter out variant
    #[serde(default = "default_as_true")]
    pub missing_gq_ok: bool,
    /// Whether missing paired-read information leads to filter out variant
    #[serde(default = "default_as_true")]
    pub missing_pr_ok: bool,
    /// Whether missing split read information leads to filter out variant
    #[serde(default = "default_as_true")]
    pub missing_sr_ok: bool,
    /// Whether missing split read or paired read information leads to filter
    /// out variant
    #[serde(default = "default_as_true")]
    pub missing_srpr_ok: bool,
    /// Whether missing read depth information leads to filter out variant
    #[serde(default = "default_as_true")]
    pub missing_rd_dev_ok: bool,
    /// Whether missing mapping quality information leads to filter out variant
    #[serde(default = "default_as_true")]
    pub missing_amq_ok: bool,

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
            max_brk_repeat: None,
            max_brk_segdup: None,
            max_brk_segduprepeat: None,
            gt_one_of: None,
            min_gq: None,
            min_pr_cov: None,
            max_pr_cov: None,
            min_pr_ref: None,
            max_pr_ref: None,
            min_pr_var: None,
            max_pr_var: None,
            min_pr_ab: None,
            max_pr_ab: None,
            min_sr_cov: None,
            max_sr_cov: None,
            min_sr_ref: None,
            max_sr_ref: None,
            min_sr_var: None,
            max_sr_var: None,
            min_sr_ab: None,
            max_sr_ab: None,
            min_srpr_cov: None,
            max_srpr_cov: None,
            min_srpr_ref: None,
            max_srpr_ref: None,
            min_srpr_var: None,
            max_srpr_var: None,
            min_srpr_ab: None,
            max_srpr_ab: None,
            min_rd_dev: None,
            max_rd_dev: None,
            min_amq: None,
            max_amq: None,
            comment: None,
            missing_gt_ok: true,
            missing_gq_ok: true,
            missing_pr_ok: true,
            missing_sr_ok: true,
            missing_srpr_ok: true,
            missing_rd_dev_ok: true,
            missing_amq_ok: true,
        }
    }

    /// Returns whether the `GenotypeCriteria` is applicable to a call in
    /// a structural variant with the given genotype, of the given SV sub
    /// type and with the given size.
    pub fn is_applicable_to(
        &self,
        genotype: GenotypeChoice,
        sv_sub_type: SvSubType,
        sv_size: Option<u32>,
    ) -> bool {
        if (genotype != GenotypeChoice::Any && genotype != self.genotype)
            || !self.select_sv_sub_type.contains(&sv_sub_type)
        {
            false
        } else if sv_sub_type == SvSubType::Bnd || sv_sub_type.is_ins() {
            true // no further size check
        } else if let Some(sv_size) = sv_size {
            // apply the size limits, if defined
            sv_size >= self.select_sv_min_size.unwrap_or(sv_size)
                && sv_size <= self.select_sv_max_size.unwrap_or(sv_size)
        } else {
            false // is not BND or INS, must have size!
        }
    }

    /// Returns whether the `GenotypeCriteria` is pass for the given `CallInfo`.
    ///
    /// Note that this only check the genotype and quality criteria.  Whether
    /// the `GenotypeCriteria` is applicable to the `CallInfo` has to be
    /// checked independently.
    pub fn is_call_info_pass(&self, call_info: &CallInfo) -> bool {
        // The pattern below is always the same: if the constraint in self is
        // None then pass regardlessly of what `call_info` has.  Otherwise
        // fail if the corresponding value of `call_info` has not been set.
        // If both have been set then perform the actual check.

        // gt -- genotype

        let pass_gt_one_of = self.gt_one_of.as_ref().map_or(true, |gt_one_of| {
            call_info
                .genotype
                .as_ref()
                .map_or(self.missing_gt_ok, |gt| gt_one_of.contains(gt))
        });

        // gq -- genotype quality

        let pass_min_gq = self.min_gq.map_or(true, |min_gq| {
            call_info
                .quality
                .map_or(self.missing_gq_ok, |gq| gq >= min_gq)
        });

        // pr -- paired-end reads

        let pass_min_pr_cov = self.min_pr_cov.map_or(true, |min_pr_cov| {
            call_info
                .paired_end_cov
                .map_or(self.missing_pr_ok, |paired_end_cov| {
                    paired_end_cov >= min_pr_cov
                })
        });
        let pass_max_pr_cov = self.max_pr_cov.map_or(true, |max_pr_cov| {
            call_info
                .paired_end_cov
                .map_or(self.missing_pr_ok, |paired_end_cov| {
                    paired_end_cov <= max_pr_cov
                })
        });

        let pass_min_pr_ref = self.min_pr_ref.map_or(true, |min_pr_ref| {
            call_info
                .paired_end_cov
                .map_or(self.missing_pr_ok, |paired_end_cov| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_pr_ok, |paired_end_var| {
                            paired_end_cov.saturating_sub(paired_end_var) >= min_pr_ref
                        })
                })
        });
        let pass_max_pr_ref = self.max_pr_ref.map_or(true, |max_pr_ref| {
            call_info
                .paired_end_cov
                .map_or(self.missing_pr_ok, |paired_end_cov| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_pr_ok, |paired_end_var| {
                            paired_end_cov.saturating_sub(paired_end_var) <= max_pr_ref
                        })
                })
        });

        let pass_min_pr_var = self.min_pr_var.map_or(true, |min_pr_var| {
            call_info
                .paired_end_var
                .map_or(self.missing_pr_ok, |paired_end_var| {
                    paired_end_var >= min_pr_var
                })
        });
        let pass_max_pr_var = self.max_pr_var.map_or(true, |max_pr_var| {
            call_info
                .paired_end_var
                .map_or(self.missing_pr_ok, |paired_end_var| {
                    paired_end_var <= max_pr_var
                })
        });

        let pass_min_pr_ab = self.min_pr_ab.map_or(true, |min_pr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_pr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_pr_ok, |split_read_var| {
                            (split_read_var as f64 / split_read_cov as f64 >= min_pr_ab as f64)
                                || (1.0 - (split_read_var as f64 / split_read_cov as f64)
                                    >= min_pr_ab as f64)
                        })
                })
        });
        let pass_max_pr_ab = self.max_pr_ab.map_or(true, |max_pr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_pr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_pr_ok, |split_read_var| {
                            (split_read_var as f64 / split_read_cov as f64 <= max_pr_ab as f64)
                                || (1.0 - (split_read_var as f64 / split_read_cov as f64)
                                    >= max_pr_ab as f64)
                        })
                })
        });

        // sr -- split reads

        let pass_min_sr_cov = self.min_sr_cov.map_or(true, |min_sr_cov| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    split_read_cov >= min_sr_cov
                })
        });
        let pass_max_sr_cov = self.max_sr_cov.map_or(true, |max_sr_cov| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    split_read_cov <= max_sr_cov
                })
        });

        let pass_min_sr_ref = self.min_sr_ref.map_or(true, |min_sr_ref| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_sr_ok, |split_read_var| {
                            split_read_cov.saturating_sub(split_read_var) >= min_sr_ref
                        })
                })
        });
        let pass_max_sr_ref = self.max_sr_ref.map_or(true, |max_sr_ref| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_sr_ok, |split_read_var| {
                            split_read_cov.saturating_sub(split_read_var) <= max_sr_ref
                        })
                })
        });

        let pass_min_sr_var = self.min_sr_var.map_or(true, |min_sr_var| {
            call_info
                .split_read_var
                .map_or(self.missing_sr_ok, |split_read_var| {
                    split_read_var >= min_sr_var
                })
        });
        let pass_max_sr_var = self.max_sr_var.map_or(true, |max_sr_var| {
            call_info
                .split_read_var
                .map_or(self.missing_sr_ok, |split_read_var| {
                    split_read_var <= max_sr_var
                })
        });

        let pass_min_sr_ab = self.min_sr_ab.map_or(true, |min_sr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_sr_ok, |split_read_var| {
                            (split_read_var as f64 / split_read_cov as f64 >= min_sr_ab as f64)
                                || (1.0 - (split_read_var as f64 / split_read_cov as f64)
                                    >= min_sr_ab as f64)
                        })
                })
        });
        let pass_max_sr_ab = self.max_sr_ab.map_or(true, |max_sr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_sr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_sr_ok, |split_read_var| {
                            (split_read_var as f64 / split_read_cov as f64 <= max_sr_ab as f64)
                                || (1.0 - (split_read_var as f64 / split_read_cov as f64)
                                    >= max_sr_ab as f64)
                        })
                })
        });

        // sr + pr -- split reads + paired-end reads

        let pass_min_srpr_cov = self.min_srpr_cov.map_or(true, |min_srpr_cov| {
            call_info
                .split_read_cov
                .map_or(self.missing_srpr_ok, |split_read_cov| {
                    call_info
                        .paired_end_cov
                        .map_or(self.missing_srpr_ok, |paired_end_cov| {
                            split_read_cov + paired_end_cov >= min_srpr_cov
                        })
                })
        });
        let pass_max_srpr_cov = self.max_srpr_cov.map_or(true, |max_srpr_cov| {
            call_info
                .split_read_cov
                .map_or(self.missing_srpr_ok, |split_read_cov| {
                    call_info
                        .paired_end_cov
                        .map_or(self.missing_srpr_ok, |paired_end_cov| {
                            split_read_cov + paired_end_cov <= max_srpr_cov
                        })
                })
        });

        let pass_min_srpr_ref = self.min_srpr_ref.map_or(true, |min_srpr_ref| {
            call_info
                .paired_end_cov
                .map_or(self.missing_srpr_ok, |paired_end_cov| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_srpr_ok, |paired_end_var| {
                            call_info.split_read_cov.map_or(
                                self.missing_srpr_ok,
                                |split_read_cov| {
                                    call_info.split_read_var.map_or(
                                        self.missing_srpr_ok,
                                        |split_read_var| {
                                            paired_end_cov.saturating_sub(paired_end_var)
                                                + split_read_cov.saturating_sub(split_read_var)
                                                >= min_srpr_ref
                                        },
                                    )
                                },
                            )
                        })
                })
        });
        let pass_max_srpr_ref = self.max_srpr_ref.map_or(true, |max_srpr_ref| {
            call_info
                .paired_end_cov
                .map_or(self.missing_srpr_ok, |paired_end_cov| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_srpr_ok, |paired_end_var| {
                            call_info.split_read_cov.map_or(
                                self.missing_srpr_ok,
                                |split_read_cov| {
                                    call_info.split_read_var.map_or(
                                        self.missing_srpr_ok,
                                        |split_read_var| {
                                            paired_end_cov.saturating_sub(paired_end_var)
                                                + split_read_cov.saturating_sub(split_read_var)
                                                <= max_srpr_ref
                                        },
                                    )
                                },
                            )
                        })
                })
        });

        let pass_min_srpr_var = self.min_srpr_var.map_or(true, |min_srpr_var| {
            call_info
                .split_read_var
                .map_or(self.missing_srpr_ok, |split_read_var| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_srpr_ok, |paired_end_var| {
                            split_read_var + paired_end_var >= min_srpr_var
                        })
                })
        });
        let pass_max_srpr_var = self.max_srpr_var.map_or(true, |max_srpr_var| {
            call_info
                .split_read_var
                .map_or(self.missing_srpr_ok, |split_read_var| {
                    call_info
                        .paired_end_var
                        .map_or(self.missing_srpr_ok, |paired_end_var| {
                            split_read_var + paired_end_var <= max_srpr_var
                        })
                })
        });

        let pass_min_srpr_ab = self.min_srpr_ab.map_or(true, |min_srpr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_srpr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_srpr_ok, |split_read_var| {
                            call_info.paired_end_cov.map_or(
                                self.missing_srpr_ok,
                                |paired_end_cov| {
                                    call_info.paired_end_var.map_or(
                                        self.missing_srpr_ok,
                                        |paired_end_var| {
                                            ((split_read_var as f64 + paired_end_var as f64)
                                                / (split_read_cov as f64 + paired_end_cov as f64)
                                                >= min_srpr_ab as f64)
                                                || (1.0
                                                    - (split_read_var as f64
                                                        + paired_end_var as f64)
                                                        / (split_read_cov as f64
                                                            + paired_end_cov as f64)
                                                    >= min_srpr_ab as f64)
                                        },
                                    )
                                },
                            )
                        })
                })
        });
        let pass_max_srpr_ab = self.max_srpr_ab.map_or(true, |max_srpr_ab| {
            call_info
                .split_read_cov
                .map_or(self.missing_srpr_ok, |split_read_cov| {
                    call_info
                        .split_read_var
                        .map_or(self.missing_srpr_ok, |split_read_var| {
                            call_info.paired_end_cov.map_or(
                                self.missing_srpr_ok,
                                |paired_end_cov| {
                                    call_info.paired_end_var.map_or(
                                        self.missing_srpr_ok,
                                        |paired_end_var| {
                                            ((split_read_var as f64 + paired_end_var as f64)
                                                / (split_read_cov as f64 + paired_end_cov as f64)
                                                <= max_srpr_ab as f64)
                                                || (1.0
                                                    - (split_read_var as f64
                                                        + paired_end_var as f64)
                                                        / (split_read_cov as f64
                                                            + paired_end_cov as f64)
                                                    <= max_srpr_ab as f64)
                                        },
                                    )
                                },
                            )
                        })
                })
        });

        // rd_dev -- read depth deviation

        let pass_min_rd_dev = self.min_rd_dev.map_or(true, |min_rd_dev| {
            call_info
                .average_normalized_cov
                .map_or(self.missing_rd_dev_ok, |average_normalized_cov| {
                    (average_normalized_cov - 1.0).abs() >= min_rd_dev
                })
        });

        let pass_max_rd_dev = self.max_rd_dev.map_or(true, |max_rd_dev| {
            call_info
                .average_normalized_cov
                .map_or(self.missing_rd_dev_ok, |average_normalized_cov| {
                    (average_normalized_cov - 1.0).abs() <= max_rd_dev
                })
        });

        // amq - average mapping quality

        let pass_min_amq = self.min_amq.map_or(true, |min_amq| {
            call_info
                .average_mapping_quality
                .map_or(self.missing_amq_ok, |average_mapping_quality| {
                    average_mapping_quality >= min_amq
                })
        });
        let pass_max_amq = self.max_amq.map_or(true, |max_amq| {
            call_info
                .average_mapping_quality
                .map_or(self.missing_amq_ok, |average_mapping_quality| {
                    average_mapping_quality <= max_amq
                })
        });

        pass_gt_one_of
            && pass_min_gq
            && pass_min_pr_cov
            && pass_max_pr_cov
            && pass_min_pr_ref
            && pass_max_pr_ref
            && pass_min_pr_var
            && pass_max_pr_var
            && pass_min_pr_ab
            && pass_max_pr_ab
            && pass_min_sr_cov
            && pass_max_sr_cov
            && pass_min_sr_ref
            && pass_max_sr_ref
            && pass_min_sr_var
            && pass_max_sr_var
            && pass_min_sr_ab
            && pass_max_sr_ab
            && pass_min_srpr_cov
            && pass_max_srpr_cov
            && pass_min_srpr_ref
            && pass_max_srpr_ref
            && pass_min_srpr_var
            && pass_max_srpr_var
            && pass_min_srpr_ab
            && pass_max_srpr_ab
            && pass_min_rd_dev
            && pass_max_rd_dev
            && pass_min_amq
            && pass_max_amq
    }

    pub fn is_masked_pass(&self, masked_count: &MaskedBreakpointCount) -> bool {
        let pass_max_brk_segdup = self
            .max_brk_segdup
            .map_or(true, |max_brk_segdup| masked_count.segdup <= max_brk_segdup);
        let pass_max_brk_repeat = self
            .max_brk_repeat
            .map_or(true, |max_brk_repeat| masked_count.repeat <= max_brk_repeat);
        let pass_max_brk_segduprepeat =
            self.max_brk_segduprepeat
                .map_or(true, |max_brk_segduprepeat| {
                    masked_count.segdup + masked_count.repeat <= max_brk_segduprepeat
                });

        pass_max_brk_segdup && pass_max_brk_repeat && pass_max_brk_segduprepeat
    }
}

/// Define a query for structural variants from a case.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct CaseQuery {
    /// Whether to enable SVDB overlap queries with DGV.
    pub svdb_dgv_enabled: bool,
    /// The minimal reciprocal overlap for querying DGV.
    pub svdb_dgv_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying DGV.
    pub svdb_dgv_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with DGV gold standard.
    pub svdb_dgv_gs_enabled: bool,
    /// The minimal reciprocal overlap for querying DGV gold standard.
    pub svdb_dgv_gs_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying DGV gold standard.
    pub svdb_dgv_gs_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with gnomAD SV.
    pub svdb_gnomad_genomes_enabled: bool,
    /// The minimal reciprocal overlap for querying gnomAD SV.
    pub svdb_gnomad_genomes_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying gnomAD SV.
    pub svdb_gnomad_genomes_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with gnomAD exomes/ExAC.
    pub svdb_gnomad_exomes_enabled: bool,
    /// The minimal reciprocal overlap for querying gnomAD exomes/ExAC.
    pub svdb_gnomad_exomes_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying gnomAD exomes/ExAC.
    pub svdb_gnomad_exomes_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with dbVar.
    pub svdb_dbvar_enabled: bool,
    /// The minimal reciprocal overlap for querying dbVar.
    pub svdb_dbvar_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying dbVar.
    pub svdb_dbvar_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with Thousand Genomes Project.
    pub svdb_g1k_enabled: bool,
    /// The minimal reciprocal overlap for querying Thousand Genomes Project.
    pub svdb_g1k_min_overlap: Option<f32>,
    /// The maximal number of carriers for querying Thousand Genomes Project.
    pub svdb_g1k_max_count: Option<u32>,
    /// Whether to enable SVDB overlap queries with in-house DB.
    pub svdb_inhouse_enabled: bool,
    /// The minimal reciprocal overlap for querying in-house DB.
    pub svdb_inhouse_min_overlap: Option<f32>,
    /// The maximal number of alleles for querying in-house DB.
    pub svdb_inhouse_max_count: Option<u32>,

    /// Minimal reciprocal overlap when overlapping with ClinVar SVs
    pub clinvar_sv_min_overlap: Option<f32>,
    /// Minimal pathogenicity when overlapping with ClinVar SVs.
    pub clinvar_sv_min_pathogenicity: Option<Pathogenicity>,

    /// The minimal SV size to consider.
    pub sv_size_min: Option<u32>,
    /// The maximal SV size to consider.
    pub sv_size_max: Option<u32>,

    /// The SV types to consider.
    pub sv_types: Vec<SvType>,
    /// The SV subtypes to consider.
    pub sv_sub_types: Vec<SvSubType>,
    /// The transcript effects to consider.
    pub tx_effects: Vec<TranscriptEffect>,

    /// List of genes to require.
    pub gene_allowlist: Option<Vec<String>>,
    /// Genomic region to limit consideration to.
    #[serde(deserialize_with = "deserialize_genomic_region")]
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
    pub genotype: IndexMap<String, GenotypeChoice>,
    /// Criteria for filtering CNVs.
    pub genotype_criteria: Vec<GenotypeCriteria>,

    /// The mode for recessive inheritance.
    pub recessive_mode: Option<RecessiveMode>,
    /// The index to use for recessive inheritance.
    pub recessive_index: Option<String>,
}

fn deserialize_genomic_region<'de, D>(
    deserializer: D,
) -> Result<Option<Vec<GenomicRegion>>, D::Error>
where
    D: Deserializer<'de>,
{
    let re = Regex::new(
        r"^(?P<chrom>(chr)?(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M|MT))(:(?P<start>\d+(,\d+)*)-(?P<stop>\d+(,\d+)*))?$",
    )
    .unwrap();

    let tokens: Option<Vec<String>> = Deserialize::deserialize(deserializer)?;

    Ok(tokens.map(|tokens| {
        tokens
            .iter()
            .filter_map(|token| {
                if let Some(caps) = re.captures(token) {
                    let chrom = caps.name("chrom").unwrap().as_str().to_string();
                    let chrom = if let Some(chrom) = chrom.strip_prefix("chr") {
                        chrom.to_string()
                    } else {
                        chrom
                    };
                    let range = if let (Some(start), Some(stop)) =
                        (caps.name("start"), caps.name("stop"))
                    {
                        let start: i32 =
                            start.as_str().replace(',', "").parse().unwrap_or_else(|_| {
                                panic!("could not parse start position: {}", start.as_str())
                            });
                        let end: i32 =
                            stop.as_str().replace(',', "").parse().unwrap_or_else(|_| {
                                panic!("could not parse stop position: {}", stop.as_str())
                            });
                        Some(Range { start, end })
                    } else {
                        None
                    };
                    Some(GenomicRegion { chrom, range })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    }))
}

impl Default for CaseQuery {
    fn default() -> Self {
        CaseQuery {
            svdb_dgv_enabled: false,
            svdb_dgv_min_overlap: None,
            svdb_dgv_max_count: None,
            svdb_dgv_gs_enabled: false,
            svdb_dgv_gs_min_overlap: None,
            svdb_dgv_gs_max_count: None,
            svdb_gnomad_genomes_enabled: false,
            svdb_gnomad_genomes_min_overlap: None,
            svdb_gnomad_genomes_max_count: None,
            svdb_gnomad_exomes_enabled: false,
            svdb_gnomad_exomes_min_overlap: None,
            svdb_gnomad_exomes_max_count: None,
            svdb_dbvar_enabled: false,
            svdb_dbvar_min_overlap: None,
            svdb_dbvar_max_count: None,
            svdb_g1k_enabled: false,
            svdb_g1k_min_overlap: None,
            svdb_g1k_max_count: None,
            svdb_inhouse_enabled: false,
            svdb_inhouse_min_overlap: None,
            svdb_inhouse_max_count: None,
            sv_size_min: None,
            sv_size_max: None,
            sv_types: SvType::vec_all(),
            sv_sub_types: SvSubType::vec_all(),
            clinvar_sv_min_overlap: None,
            clinvar_sv_min_pathogenicity: None,
            gene_allowlist: None,
            genomic_region: None,
            regulatory_overlap: 100,
            regulatory_ensembl_features: None,
            regulatory_vista_validation: None,
            regulatory_custom_configs: vec![],
            tad_set: None,
            genotype: IndexMap::new(),
            genotype_criteria: vec![],
            recessive_mode: None,
            recessive_index: None,
            tx_effects: TranscriptEffect::vec_all(),
        }
    }
}

/// Information on the call as combined by the annotator.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone, Default)]
pub struct CallInfo {
    /// The genotype, if applicable, e.g., "0/1", "./1", "."
    pub genotype: Option<String>,
    /// The effective genotype, if set.
    pub effective_genotype: Option<Genotype>,
    /// All compatible genotypes by the genotype criteria list.
    pub matched_gt_criteria: Option<Vec<Genotype>>,
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
    /// Average mapping quality, if applicable
    pub average_mapping_quality: Option<f32>,
}

/// Definition of a structural variant with per-sample genotype calls.
///
/// This uses a subset/specialization of what is described by the VCF standard
/// for the purpose of running SV queries in `varfish-server-worker`.
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct StructuralVariant {
    /// Chromosome name
    pub chrom: String,
    /// 1-based start position of the variant (or position on first chromosome
    /// for break-ends)
    pub pos: i32,
    /// Type of the structural variant
    pub sv_type: SvType,
    /// Sub type of the structural variant
    pub sv_sub_type: SvSubType,
    /// Potentially the second involved chromosome
    pub chrom2: Option<String>,
    /// End position (position on second chromosome for break-ends)
    pub end: i32,
    /// The strand orientation of the structural variant.
    pub strand_orientation: StrandOrientation,

    /// The callers of the variant.
    pub callers: Vec<String>,
    /// Mapping of sample to genotype information for the SV.
    pub call_info: IndexMap<String, CallInfo>,
}

impl StructuralVariant {
    /// Return the size of the structural variant.
    ///
    /// Size is not applicable for insertions and break-ends, so `None` is
    /// returned in this case.
    pub fn size(&self) -> Option<u32> {
        if self.sv_type == SvType::Ins
            || self.sv_type == SvType::Bnd
            || self.sv_sub_type.is_ins()
            || self.sv_sub_type == SvSubType::Bnd
        {
            None
        } else {
            Some((self.end - self.pos + 1) as u32)
        }
    }

    /// Convert from VCF record.
    pub fn from_vcf(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<Self, anyhow::Error> {
        use noodles::vcf::variant::record::info::field::key;

        let chrom = record.reference_sequence_name().to_string();
        let pos: usize = record.variant_start().expect("no variant_start?").into();
        let pos = pos as i32;
        let sv_type: SvType =
            if let Some(Some(vcf::variant::record_buf::info::field::Value::String(sv_type))) =
                record.info().get(key::SV_TYPE)
            {
                sv_type
                    .parse()
                    .map_err(|e| anyhow::anyhow!("could not parse SVTYPE {}: {}", &sv_type, e))?
            } else {
                anyhow::bail!("no INFO/SVTYPE in VCF record")
            };
        let sv_sub_type = match sv_type {
            SvType::Del => SvSubType::Del,
            SvType::Dup => SvSubType::Dup,
            SvType::Inv => SvSubType::Inv,
            SvType::Ins => SvSubType::Ins,
            SvType::Bnd => SvSubType::Bnd,
            SvType::Cnv => SvSubType::Cnv,
        };
        let end = if let Some(Some(vcf::variant::record_buf::info::field::Value::Integer(end))) =
            record.info().get(key::END_POSITION)
        {
            *end
        } else {
            anyhow::bail!("no INFO/END in VCF record")
        };
        let chrom2 =
            if let Some(Some(vcf::variant::record_buf::info::field::Value::String(chrom2))) =
                record.info().get("chr2")
            {
                Some(chrom2.clone())
            } else {
                None
            };

        let strand_orientation = match sv_type {
            SvType::Del => StrandOrientation::ThreeToFive,
            SvType::Dup => StrandOrientation::FiveToThree,
            SvType::Inv => StrandOrientation::FiveToFive,
            SvType::Bnd => {
                let pe_orientation = Breakend::from_ref_alt_str(
                    record.reference_bases(),
                    record
                        .alternate_bases()
                        .iter()
                        .next()
                        .ok_or_else(|| anyhow::anyhow!("no alternate allele?"))??,
                )?
                .pe_orientation;
                match pe_orientation {
                    PeOrientation::ThreeToThree => StrandOrientation::ThreeToThree,
                    PeOrientation::FiveToFive => StrandOrientation::FiveToFive,
                    PeOrientation::ThreeToFive => StrandOrientation::ThreeToFive,
                    PeOrientation::FiveToThree => StrandOrientation::FiveToThree,
                    PeOrientation::Other => StrandOrientation::NotApplicable,
                }
            }
            SvType::Ins | SvType::Cnv => StrandOrientation::NotApplicable,
        };

        let callers = if let Some(Some(vcf::variant::record_buf::info::field::Value::Array(
            vcf::variant::record_buf::info::field::value::Array::String(callers),
        ))) = record.info().get("callers")
        {
            callers.iter().flatten().cloned().collect::<Vec<_>>()
        } else {
            anyhow::bail!("no INFO/callers in VCF record")
        };

        let call_info = Self::build_call_info(record, header)?;

        Ok(Self {
            chrom,
            pos,
            sv_type,
            sv_sub_type,
            chrom2,
            end,
            strand_orientation,
            callers,
            call_info,
        })
    }

    /// Build call information.
    fn build_call_info(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<IndexMap<String, CallInfo>, anyhow::Error> {
        let mut result = IndexMap::new();

        for (name, sample) in header.sample_names().iter().zip(record.samples().values()) {
            let mut call_info = CallInfo::default();

            for (key, value) in record.samples().keys().as_ref().iter().zip(sample.values()) {
                let value = if let Some(value) = value.as_ref() {
                    value
                } else {
                    continue; // is empty
                };
                use noodles::vcf::variant::record::samples::keys::key;
                match (key.as_str(), value) {
                    (
                        key::GENOTYPE,
                        vcf::variant::record_buf::samples::sample::value::Value::Genotype(gt),
                    ) => {
                        call_info.genotype =
                            Some((genotype_to_string(&gt)?.as_str()[1..]).to_string());
                    }
                    (
                        key::CONDITIONAL_GENOTYPE_QUALITY,
                        vcf::variant::record_buf::samples::sample::Value::Float(quality),
                    ) => call_info.quality = Some(*quality),
                    (
                        key::CONDITIONAL_GENOTYPE_QUALITY,
                        vcf::variant::record_buf::samples::sample::Value::Integer(quality),
                    ) => call_info.quality = Some(*quality as f32),
                    (
                        "pec",
                        vcf::variant::record_buf::samples::sample::Value::Integer(paired_end_cov),
                    ) => call_info.paired_end_cov = Some(*paired_end_cov as u32),
                    (
                        "pev",
                        vcf::variant::record_buf::samples::sample::Value::Integer(paired_end_var),
                    ) => call_info.paired_end_var = Some(*paired_end_var as u32),
                    (
                        "src",
                        vcf::variant::record_buf::samples::sample::Value::Integer(split_read_cov),
                    ) => call_info.split_read_cov = Some(*split_read_cov as u32),
                    (
                        "srv",
                        vcf::variant::record_buf::samples::sample::Value::Integer(split_read_var),
                    ) => call_info.split_read_var = Some(*split_read_var as u32),
                    (
                        "cn",
                        vcf::variant::record_buf::samples::sample::Value::Integer(copy_number),
                    ) => call_info.copy_number = Some(*copy_number as u32),
                    (
                        "anc",
                        vcf::variant::record_buf::samples::sample::Value::Integer(
                            average_normalized_cov,
                        ),
                    ) => call_info.average_normalized_cov = Some(*average_normalized_cov as f32),
                    (
                        "anc",
                        vcf::variant::record_buf::samples::sample::Value::Float(
                            average_normalized_cov,
                        ),
                    ) => call_info.average_normalized_cov = Some(*average_normalized_cov),
                    (
                        "pc",
                        vcf::variant::record_buf::samples::sample::Value::Integer(point_count),
                    ) => call_info.point_count = Some(*point_count as u32),
                    (
                        "amq",
                        vcf::variant::record_buf::samples::sample::Value::Integer(
                            average_mapping_quality,
                        ),
                    ) => call_info.average_mapping_quality = Some(*average_mapping_quality as f32),
                    _ => {
                        panic!("unknown FORMAT key: {}", key);
                    }
                }
            }

            result.insert(name.clone(), call_info);
        }

        Ok(result)
    }
}

impl From<VariationType> for crate::pbs::varfish::v1::strucvars::bgdb::VariationType {
    fn from(val: VariationType) -> Self {
        match val {
            VariationType::Complex => {
                crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Complex
            }
            VariationType::Microsatellite => {
                crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Microsatellite
            }
            VariationType::Dup => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Dup,
            VariationType::Del => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Del,
            VariationType::Bnd => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Bnd,
            VariationType::Cnv => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Cnv,
            VariationType::Inv => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Inv,
            VariationType::Ins => crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Ins,
        }
    }
}

impl TryInto<VariationType> for crate::pbs::varfish::v1::strucvars::bgdb::VariationType {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<VariationType, anyhow::Error> {
        Ok(match self {
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Complex => {
                VariationType::Complex
            }
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Microsatellite => {
                VariationType::Microsatellite
            }
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Dup => VariationType::Dup,
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Del => VariationType::Del,
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Bnd => VariationType::Bnd,
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Cnv => VariationType::Cnv,
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Inv => VariationType::Inv,
            crate::pbs::varfish::v1::strucvars::bgdb::VariationType::Ins => VariationType::Ins,
        })
    }
}

/// Variation type from Clinvar.
#[derive(
    Serialize, PartialEq, PartialOrd, Eq, Hash, Copy, Clone, Debug, Default, EnumString, Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum VariationType {
    #[default]
    Complex,
    Microsatellite,
    Dup,
    Del,
    Bnd,
    Cnv,
    Inv,
    Ins,
}

/// Clinvar pathogenicity.
#[derive(
    Serialize,
    Deserialize,
    PartialEq,
    PartialOrd,
    Eq,
    Hash,
    Copy,
    Clone,
    Debug,
    Default,
    EnumString,
    Display,
)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum Pathogenicity {
    Benign,
    LikelyBenign,
    #[default]
    Uncertain,
    LikelyPathogenic,
    Pathogenic,
}

impl From<Pathogenicity> for crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity {
    fn from(val: Pathogenicity) -> Self {
        match val {
            Pathogenicity::Benign => {
                crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Benign
            }
            Pathogenicity::LikelyBenign => {
                crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::LikelyBenign
            }
            Pathogenicity::Uncertain => {
                crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Uncertain
            }
            Pathogenicity::LikelyPathogenic => {
                crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::LikelyPathogenic
            }
            Pathogenicity::Pathogenic => {
                crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Pathogenic
            }
        }
    }
}

impl TryInto<Pathogenicity> for crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity {
    type Error = anyhow::Error;

    fn try_into(self) -> Result<Pathogenicity, anyhow::Error> {
        Ok(match self {
            crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Benign => {
                Pathogenicity::Benign
            }
            crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::LikelyBenign => {
                Pathogenicity::LikelyBenign
            }
            crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Uncertain => {
                Pathogenicity::Uncertain
            }
            crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::LikelyPathogenic => {
                Pathogenicity::LikelyPathogenic
            }
            crate::pbs::varfish::v1::strucvars::bgdb::Pathogenicity::Pathogenic => {
                Pathogenicity::Pathogenic
            }
        })
    }
}

/// Clinvar SV as created by clinvar-tsv.
#[derive(Debug)]
pub struct ClinvarSv {
    /// Genome release
    pub release: String,
    /// Chromosome name
    pub chromosome: String,
    /// 1-based start position
    pub start: i32,
    /// 1-based end position
    pub end: i32,
    /// Clinvar variation type
    pub variation_type: VariationType,
    /// Clinvar pathogenicty
    pub pathogenicity: Pathogenicity,
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
                Token::I32(1),
                Token::Str("end"),
                Token::I32(2),
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
        insta::assert_snapshot!(serde_json::to_string_pretty(&crit).unwrap())
    }

    #[test]
    fn test_genotype_criteria_is_applicable_to() {
        let crit = GenotypeCriteria {
            genotype: GenotypeChoice::Het,
            select_sv_sub_type: vec![SvSubType::Del],
            select_sv_min_size: Some(1000),
            select_sv_max_size: Some(5000),
            ..GenotypeCriteria::new(GenotypeChoice::Het)
        };

        assert!(crit.is_applicable_to(GenotypeChoice::Het, SvSubType::Del, Some(1000)));
        assert!(!crit.is_applicable_to(GenotypeChoice::Hom, SvSubType::Ins, Some(1000)));
        assert!(!crit.is_applicable_to(GenotypeChoice::Het, SvSubType::Ins, Some(1000)));
        assert!(!crit.is_applicable_to(GenotypeChoice::Het, SvSubType::Del, Some(100)));
        assert!(!crit.is_applicable_to(GenotypeChoice::Het, SvSubType::Del, Some(10000)));
    }

    #[test]
    fn test_genotype_criteria_is_call_info_pass() {
        let crit = GenotypeCriteria {
            gt_one_of: Some(vec!["0/1".to_owned()]),
            min_gq: Some(10.0),
            min_pr_cov: Some(10),
            max_pr_cov: Some(20),
            min_pr_ref: Some(10),
            max_pr_ref: Some(20),
            min_pr_var: Some(10),
            max_pr_var: Some(20),
            min_sr_cov: Some(10),
            max_sr_cov: Some(20),
            min_sr_var: Some(10),
            max_sr_var: Some(20),
            min_sr_ref: Some(10),
            max_sr_ref: Some(20),
            min_srpr_cov: Some(20),
            max_srpr_cov: Some(50),
            min_srpr_var: Some(20),
            max_srpr_var: Some(50),
            min_srpr_ref: Some(20),
            max_srpr_ref: Some(50),
            min_rd_dev: Some(0.5),
            max_rd_dev: Some(1.5),
            min_amq: Some(60.0),
            max_amq: Some(70.0),
            ..GenotypeCriteria::new(GenotypeChoice::Het)
        };

        let pass_info = CallInfo {
            genotype: Some("0/1".to_owned()),
            quality: Some(10.0),
            paired_end_cov: Some(20),
            paired_end_var: Some(10),
            split_read_cov: Some(20),
            split_read_var: Some(10),
            copy_number: Some(1),
            average_normalized_cov: Some(0.491),
            point_count: Some(5),
            average_mapping_quality: Some(60.0),
            ..Default::default()
        };

        assert!(crit.is_call_info_pass(&pass_info));
    }

    #[test]
    fn test_genotype_criteria_is_call_info_fail() {
        let crit = GenotypeCriteria {
            gt_one_of: Some(vec!["0/1".to_owned()]),
            min_gq: Some(10.0),
            min_pr_cov: Some(10),
            min_pr_var: Some(10),
            min_sr_cov: Some(10),
            min_sr_var: Some(10),
            min_rd_dev: Some(0.5),
            min_amq: Some(60.0),
            ..GenotypeCriteria::new(GenotypeChoice::Het)
        };

        let fail_info = CallInfo {
            genotype: Some("0/1".to_owned()),
            quality: Some(10.0),
            paired_end_cov: Some(10),
            paired_end_var: Some(10),
            split_read_cov: Some(10),
            split_read_var: Some(10),
            copy_number: Some(1),
            average_normalized_cov: Some(0.491),
            point_count: Some(5),
            average_mapping_quality: Some(60.0),
            ..Default::default()
        };

        assert!(!crit.is_call_info_pass(&CallInfo {
            genotype: Some("1/1".to_owned()),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            quality: Some(9.9),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            paired_end_cov: Some(9),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            paired_end_var: Some(9),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            split_read_cov: Some(9),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            split_read_var: Some(9),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            average_normalized_cov: Some(0.6),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            average_normalized_cov: Some(1.4),
            ..fail_info.clone()
        }));
        assert!(!crit.is_call_info_pass(&CallInfo {
            average_mapping_quality: Some(59.0),
            ..fail_info
        }));
    }

    #[test]
    fn test_case_query_serde_smoke() {
        let query: CaseQuery = CaseQuery::default();
        insta::assert_snapshot!(serde_json::to_string_pretty(&query).unwrap());
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
            average_mapping_quality: Some(60.0),
            ..Default::default()
        };
        insta::assert_snapshot!(serde_json::to_string_pretty(&info).unwrap());
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
            strand_orientation: StrandOrientation::ThreeToFive,
            callers: Vec::new(),
            call_info: IndexMap::new(),
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
            strand_orientation: StrandOrientation::ThreeToFive,
            callers: Vec::new(),
            call_info: IndexMap::new(),
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
            strand_orientation: StrandOrientation::ThreeToFive,
            callers: Vec::new(),
            call_info: IndexMap::new(),
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
            strand_orientation: StrandOrientation::ThreeToFive,
            callers: Vec::new(),
            call_info: IndexMap::new(),
        };
        insta::assert_snapshot!(serde_json::to_string_pretty(&sv).unwrap());
    }

    #[test]
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

    #[test]
    fn test_sv_sub_type_is_del() {
        assert_eq!(SvSubType::Del.is_del(), true);
        assert_eq!(SvSubType::DelMe.is_del(), true);
        assert_eq!(SvSubType::DelMeSva.is_del(), true);
        assert_eq!(SvSubType::DelMeL1.is_del(), true);
        assert_eq!(SvSubType::DelMeAlu.is_del(), true);
        assert_eq!(SvSubType::Dup.is_del(), false);
        assert_eq!(SvSubType::DupTandem.is_del(), false);
        assert_eq!(SvSubType::Inv.is_del(), false);
        assert_eq!(SvSubType::Ins.is_del(), false);
        assert_eq!(SvSubType::InsMe.is_del(), false);
        assert_eq!(SvSubType::InsMeSva.is_del(), false);
        assert_eq!(SvSubType::InsMeL1.is_del(), false);
        assert_eq!(SvSubType::InsMeAlu.is_del(), false);
        assert_eq!(SvSubType::Bnd.is_del(), false);
        assert_eq!(SvSubType::Cnv.is_del(), false);
    }
}
