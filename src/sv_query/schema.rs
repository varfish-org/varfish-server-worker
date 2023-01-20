/// Definition of query schemas.
use serde::{Deserialize, Serialize};

/// Range with 1-based positions
#[derive(Serialize, Deserialize, PartialEq, Debug, Clone)]
pub struct Range {
    pub start: i32,
    pub end: i32,
}

impl Range {
    pub fn new(start: i32, end: i32) -> Range {
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
    pub fn new(chrom: &str, start: i32, end: i32) -> GenomicRegion {
        GenomicRegion {
            chrom: chrom.to_owned(),
            range: Some(Range::new(start, end)),
        }
    }

    pub fn whole_chrom(chrom: &str) -> GenomicRegion {
        GenomicRegion {
            chrom: chrom.to_owned(),
            range: None,
        }
    }
}

/// Database of transcripts
#[derive(Serialize, Deserialize, PartialEq, Debug)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum Database {
    Refseq,
    Ensembl,
}

/// Encode the type of an SV
#[derive(Serialize, Deserialize, PartialEq, Debug)]
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

#[derive(Serialize, Deserialize, PartialEq, Debug)]
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

#[derive(Serialize, Deserialize, PartialEq, Debug)]
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

// TODO: VariantEffect
// TODO: EnsemblRegulatoryFeature
// TODO: VistaValidation
// TODO: GenotypeCriteria

/// Define rule to apply to a given sub set of structural variants for matching a genotype.
///
/// See documentation of VarFish Server for full documentation.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
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
    pub fn new(genotype: GenotypeChoice) -> GenotypeCriteria {
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

// /// Define a query for structural variants from a case.
// #[derive(Serialize, Deserialize, PartialEq, Debug)]
// pub struct CaseQuery {
//     /// The transcript database to use
//     pub database: Database,

//     /// Whether to enable SVDB overlap queries with DGV.
//     svdb_dgv_enabled: bool,
//     /// The minimal reciprocal overlap for querying DGV.
//     svdb_dgv_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying DGV.
//     svdb_dgv_max_carriers: Option<i32>,
//     /// Whether to enable SVDB overlap queries with DGV gold standard.
//     svdb_dgv_gs_enabled: bool,
//     /// The minimal reciprocal overlap for querying DGV gold standard.
//     svdb_dgv_gs_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying DGV gold standard.
//     svdb_dgv_gs_max_carriers: Option<i32>,
//     /// Whether to enable SVDB overlap queries with gnomAD.
//     svdb_gnomad_enabled: bool,
//     /// The minimal reciprocal overlap for querying gnomAD.
//     svdb_gnomad_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying gnomAD.
//     svdb_gnomad_max_carriers: Option<i32>,
//     /// Whether to enable SVDB overlap queries with ExAC.
//     svdb_exac_enabled: bool,
//     /// The minimal reciprocal overlap for querying ExAC.
//     svdb_exac_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying ExAC.
//     svdb_exac_max_carriers: Option<i32>,
//     /// Whether to enable SVDB overlap queries with dbVar.
//     svdb_dbvar_enabled: bool,
//     /// The minimal reciprocal overlap for querying dbVar.
//     svdb_dbvar_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying dbVar.
//     svdb_dbvar_max_carriers: Option<i32>,
//     /// Whether to enable SVDB overlap queries with Thousand Genomes Project.
//     svdb_g1k_enabled: bool,
//     /// The minimal reciprocal overlap for querying Thousand Genomes Project.
//     svdb_g1k_min_overlap: Option<f32>,
//     /// The maximal number of carriers for querying Thousand Genomes Project.
//     svdb_g1k_max_alleles: Option<i32>,
//     /// Whether to enable SVDB overlap queries with in-house DB.
//     svdb_inhouse_enabled: bool,
//     /// The minimal reciprocal overlap for querying in-house DB.
//     svdb_inhouse_min_overlap: Option<f32>,
//     /// The maximal number of alleles for querying in-house DB.
//     svdb_inhouse_max_carriers: Option<i32>,

//     /// The minimal SV size to consider.
//     sv_size_min: Option<i32>,
//     /// The maximal SV size to consider.
//     sv_size_max: Option<i32>,

//     /// The SV types to consider.
//     sv_types: Vec<SvType>,
//     /// The SV subtypes to consider.
//     sv_sub_types: Vec<SvSubType>,

//     /// List of genes to require.
//     gene_allowlist: Option<List<String>>,
//     /// Genomic region to limit consideration to.
//     genomic_region: typing.Optional[typing.List[GenomicRegionV1]] = None

//     /// Regulatory region padding to use.
//     regulatory_overlap: int = 100
//     /// Regulatory features to select from ENSEMBL.
//     regulatory_ensembl_features: typing.Optional[typing.List[EnsemblRegulatorFeature]] = None
//     /// VISTA enhancer validation results.
//     regulatory_vista_validation: typing.Optional[VistaValidation] = None

//     /// Custom regulatory maps configuration.
//     regulatory_custom_configs: typing.List[RegulatoryCustomConfig] = attrs.field(factory=list)

//     /// Name of the TAD set to use for annotation, if any.
//     tad_set: typing.Optional[str] = None

//     /// Genotype choices
//     genotype: typing.Dict[str, typing.Optional[GenotypeChoice]] = attrs.field(factory=dict)
//     /// Criteria for filtering CNVs.
//     genotype_criteria: typing.List[GenotypeCriteria] = attrs.field(factory=list)

//     /// The mode for recessive inheritance.
//     recessive_mode: typing.Optional[RecessiveModeV1] = None
//     /// The index to use for recessive inheritance.
//     recessive_index: typing.Optional[str] = None

// }

// impl CaseQuery {
//     pub fn new(database: Database) -> CaseQuery {
//         CaseQuery {
//             database,
//         }
//     }
// }

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
}
