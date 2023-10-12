//! Supporting code for seqvar query definition.

use noodles_vcf as vcf;
use strum::IntoEnumIterator;

/// Variant effects.
#[derive(
    serde::Serialize,
    serde::Deserialize,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Debug,
    Clone,
    Copy,
    strum::EnumIter,
)]
pub enum VariantEffect {
    /// 3' UTR exon variant.
    #[serde(rename = "3_prime_UTR_exon_variant")]
    ThreePrimeUtrExonVariant,
    /// 3' UTR intron variant.
    #[serde(rename = "3_prime_UTR_intron_variant")]
    ThreePrimeUtrIntronVariant,
    /// 5' UTR exon variant.
    #[serde(rename = "5_prime_UTR_exon_variant")]
    FivePrimeUtrExonVariant,
    /// 5' UTR intron variant.
    #[serde(rename = "5_prime_UTR_intron_variant")]
    FivePrimeUtrIntronVariant,
    /// Coding transcript intron variant.
    #[serde(rename = "coding_transcript_intron_variant")]
    CodingTranscriptIntronVariant,
    /// Complex substitution.
    #[serde(rename = "complex_substitution")]
    ComplexSubstitution,
    /// Direct tandem duplication.
    #[serde(rename = "direct_tandem_duplication")]
    DirectTandemDuplication,
    /// Disruptive in-frame deletion.
    #[serde(rename = "disruptive_inframe_deletion")]
    DisruptiveInframeDeletion,
    /// Disruptive in-frame insertion.
    #[serde(rename = "disruptive_inframe_insertion")]
    DisruptiveInframeInsertion,
    /// Downstream gene variant.
    #[serde(rename = "downstream_gene_variant")]
    DownstreamGeneVariant,
    /// Exon loss variant.
    #[serde(rename = "exon_loss_variant")]
    ExonLossVariant,
    /// Feature truncation.
    #[serde(rename = "feature_truncation")]
    FeatureTruncation,
    /// Frameshift elongation.
    #[serde(rename = "frameshift_elongation")]
    FrameshiftElongation,
    /// Frameshift truncation.
    #[serde(rename = "frameshift_truncation")]
    FrameshiftTruncation,
    /// Frameshift variant.
    #[serde(rename = "frameshift_variant")]
    FrameshiftVariant,
    /// In-frame deletion.
    #[serde(rename = "inframe_deletion")]
    InframeDeletion,
    /// In-frame insertion.
    #[serde(rename = "inframe_insertion")]
    InframeInsertion,
    /// Intergenic variant.
    #[serde(rename = "intergenic_variant")]
    IntergenicVariant,
    /// Internal feature elongation.
    #[serde(rename = "internal_feature_elongation")]
    InternalFeatureElongation,
    /// Missense variant.
    #[serde(rename = "missense_variant")]
    MissenseVariant,
    /// MNV.
    #[serde(rename = "mnv")]
    Mnv,
    /// Non-coding transcript exon variant.
    #[serde(rename = "non_coding_transcript_exon_variant")]
    NonCodingTranscriptExonVariant,
    /// Non-coding transcript intron variant.
    #[serde(rename = "non_coding_transcript_intron_variant")]
    NonCodingTranscriptIntronVariant,
    /// Splice acceptor variant.
    #[serde(rename = "splice_acceptor_variant")]
    SpliceAcceptorVariant,
    /// Splice donor variant.
    #[serde(rename = "splice_donor_variant")]
    SpliceDonorVariant,
    /// Splice region variant.
    #[serde(rename = "splice_region_variant")]
    SpliceRegionVariant,
    /// Start lost.
    #[serde(rename = "start_lost")]
    StartLost,
    /// Stop gained.
    #[serde(rename = "stop_gained")]
    StopGained,
    /// Stop lost.
    #[serde(rename = "stop_lost")]
    StopLost,
    /// Stop retained variant.
    #[serde(rename = "stop_retained_variant")]
    StopRetainedVariant,
    /// Structural variant.
    #[serde(rename = "structural_variant")]
    StructuralVariant,
    /// Synonymous variant.
    #[serde(rename = "synonymous_variant")]
    SynonymousVariant,
    /// Transcript ablation.
    #[serde(rename = "transcript_ablation")]
    TranscriptAblation,
    /// Upstream gene variant.
    #[serde(rename = "upstream_gene_variant")]
    UpstreamGeneVariant,
}

impl VariantEffect {
    /// Return vector of all values of `VariantEffect`.
    pub fn all() -> Vec<Self> {
        Self::iter().collect()
    }
}

/// Enumeration for recessive mode queries.
#[derive(
    serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Clone, Copy,
)]
pub enum RecessiveMode {
    /// Recessive.
    #[serde(rename = "recessive")]
    Recessive,
    /// Compound recessive.
    #[serde(rename = "compound-recessive")]
    CompoundRecessive,
}

/// Choices for failing quality thresholds on genotypes.
#[derive(
    serde::Serialize,
    serde::Deserialize,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Debug,
    Clone,
    Copy,
    Default,
)]
pub enum FailChoice {
    /// Ignore failure.
    #[default]
    #[serde(rename = "ignore")]
    Ignore,
    /// Drop whole variant.
    #[serde(rename = "drop-variant")]
    Drop,
    /// Interpret as "no-call".
    #[serde(rename = "no-call")]
    NoCall,
}

/// Choice for genotype.
#[derive(
    serde::Serialize,
    serde::Deserialize,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Debug,
    Clone,
    Copy,
    Default,
)]
pub enum GenotypeChoice {
    /// Any genotype.
    #[default]
    #[serde(rename = "any")]
    Any,
    /// Ref. genotype.
    #[serde(rename = "ref")]
    Ref,
    /// Het. genotype.
    #[serde(rename = "het")]
    Het,
    /// Hom. genotype.
    #[serde(rename = "hom")]
    Hom,
    /// Non-hom. genotype.
    #[serde(rename = "non-hom")]
    NonHom,
    /// Variant genotype.
    #[serde(rename = "variant")]
    Variant,
    /// Non-variant genotype.
    #[serde(rename = "non-variant")]
    NonVariant,
    /// Non-reference genotype.
    #[serde(rename = "non-reference")]
    NonReference,
    /// Index in comp. het. recessive inheritance.
    #[serde(rename = "comphet-index")]
    ComphetIndex,
    /// Index in recessive inheritance.
    #[serde(rename = "recessive-index")]
    RecessiveIndex,
    /// Parent in recessive inheritance.
    #[serde(rename = "recessive-parent")]
    RecessiveParent,
}

/// Quality settings for one sample.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone, Default)]
pub struct QualitySettings {
    /// Minimal coverage for het. sites.
    pub dp_het: Option<i32>,
    /// Minimal coverage for hom. sites.
    pub dp_hom: Option<i32>,
    /// Minimal genotype quality.
    pub gq: Option<i32>,
    /// Minimal allele balance.
    pub ab: Option<f32>,
    /// Minimal number of alternate reads.
    pub ad: Option<i32>,
    /// Maximal number of alternate reads
    pub ad_max: Option<i32>,
    /// Behaviour on failing quality thresholds.
    pub fail: FailChoice,
}

/// Data structure to hold a range.
#[derive(
    serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Clone,
)]
pub struct Range {
    /// Start of range.
    pub start: i32,
    /// End of range.
    pub end: i32,
}

/// Data struture to hold a genomic region.
#[derive(
    serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Clone,
)]
pub struct GenomicRegion {
    /// Chromosome.
    pub chrom: String,
    /// Range of region.
    pub range: Option<Range>,
}

/// Data structure with a single query.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone)]
#[serde(default)]
pub struct CaseQuery {
    /// Effects to consider.
    pub effects: Vec<VariantEffect>,

    /// Whether to enable filtration by gnomAD exomes.
    pub gnomad_exomes_enabled: bool,
    /// Whether to enable filtration by gnomAD genomes
    pub gnomad_genomes_enabled: bool,
    /// Whether to enable filtration by 1000 Genomes.
    pub inhouse_enabled: bool,
    /// Whether to enable filtration by mtDB.
    pub helixmtdb_enabled: bool,

    /// Quality settings for each individual.
    pub quality: indexmap::IndexMap<String, QualitySettings>,
    /// Genotype choice for each individual.
    pub genotype: indexmap::IndexMap<String, Option<GenotypeChoice>>,
    /// List of selected variants.
    pub selected_variants: Option<Vec<String>>, // TODO: remove?

    /// Whether to include coding transcripts.
    pub transcripts_coding: bool,
    /// Whether to include non-coding transcripts.
    pub transcripts_noncoding: bool,

    /// Whether to include SNVs.
    pub var_type_snv: bool,
    /// Whether to include indels.
    pub var_type_indel: bool,
    /// Whether to include MNVs.
    pub var_type_mnv: bool,

    /// Maximal distance to next exon, if any.
    pub max_exon_dist: Option<i32>,

    /// List of HGNC symbols, HGNC:<ID>s, ENSG<ID>s, or NCBI Gene IDs to restrict
    /// the resulting variants to.
    pub gene_allowlist: Option<Vec<String>>,
    /// List of genomic regions to limit restrict the resulting variants to.
    pub genomic_regions: Option<Vec<GenomicRegion>>,

    /// Wether to require ClinVar membership.
    pub require_in_clinvar: bool,
    /// Whether to include benign ClinVar variants.
    pub clinvar_include_benign: bool,
    /// Whether to include pathogenic ClinVar variants.
    pub clinvar_include_pathogenic: bool,
    /// Whether to include likely benign ClinVar variants.
    pub clinvar_include_likely_benign: bool,
    /// Whether to include likely pathogenic ClinVar variants.
    pub clinvar_include_likely_pathogenic: bool,
    /// Whether to include uncertain significance ClinVar variants.
    pub clinvar_include_uncertain_significance: bool,

    /// Maximal frequency in gnomAD exomes.
    pub gnomad_exomes_frequency: Option<f32>,
    /// Maximal number of heterozygous carriers in gnomAD exomes.
    pub gnomad_exomes_heterozygous: Option<i32>,
    /// Maximal number of homozygous carriers in gnomAD exomes.
    pub gnomad_exomes_homozygous: Option<i32>,
    /// Maximal number of hemizygous carriers in gnomAD exomes.
    pub gnomad_exomes_hemizygous: Option<i32>,

    /// Maximal frequency in gnomAD genomes.
    pub gnomad_genomes_frequency: Option<f32>,
    /// Maximal number of heterozygous carriers in gnomAD genomes.
    pub gnomad_genomes_heterozygous: Option<i32>,
    /// Maximal number of homozygous carriers in gnomAD genomes.
    pub gnomad_genomes_homozygous: Option<i32>,
    /// Maximal number of hemizygous carriers in gnomAD genomes.
    pub gnomad_genomes_hemizygous: Option<i32>,

    /// Maximal number of in-house carriers.
    pub inhouse_carriers: Option<i32>,
    /// Maximal number of in-house heterozygous carriers.
    pub inhouse_heterozygous: Option<i32>,
    /// Maximal number of in-house homozygous carriers.
    pub inhouse_homozygous: Option<i32>,
    /// Maximal number of in-house hemizygous carriers.
    pub inhouse_hemizygous: Option<i32>,

    /// Maximal frequency in HelixMtDb
    pub helixmtdb_frequency: Option<f32>,
    /// Maximal number of heterozygous carriers in HelixMtDb.
    pub helixmtdb_heterozygous: Option<i32>,
    /// Maximal number of homozygous carriers in HelixMtDb.
    pub helixmtdb_homozygous: Option<i32>,
}

impl Default for CaseQuery {
    /// Returns default values for a `CaseQuery` which makes all variants pass.
    fn default() -> Self {
        Self {
            effects: VariantEffect::all(),
            gnomad_exomes_enabled: Default::default(),
            gnomad_genomes_enabled: Default::default(),
            inhouse_enabled: Default::default(),
            helixmtdb_enabled: Default::default(),
            quality: Default::default(),
            genotype: Default::default(),
            selected_variants: Default::default(),
            transcripts_coding: true,
            transcripts_noncoding: true,
            var_type_snv: true,
            var_type_indel: true,
            var_type_mnv: true,
            max_exon_dist: Default::default(),
            gene_allowlist: Default::default(),
            genomic_regions: Default::default(),
            require_in_clinvar: Default::default(),
            clinvar_include_benign: true,
            clinvar_include_pathogenic: true,
            clinvar_include_likely_benign: true,
            clinvar_include_likely_pathogenic: true,
            clinvar_include_uncertain_significance: true,
            gnomad_exomes_frequency: Default::default(),
            gnomad_exomes_heterozygous: Default::default(),
            gnomad_exomes_homozygous: Default::default(),
            gnomad_exomes_hemizygous: Default::default(),
            gnomad_genomes_frequency: Default::default(),
            gnomad_genomes_heterozygous: Default::default(),
            gnomad_genomes_homozygous: Default::default(),
            gnomad_genomes_hemizygous: Default::default(),
            inhouse_carriers: Default::default(),
            inhouse_heterozygous: Default::default(),
            inhouse_homozygous: Default::default(),
            inhouse_hemizygous: Default::default(),
            helixmtdb_frequency: Default::default(),
            helixmtdb_heterozygous: Default::default(),
            helixmtdb_homozygous: Default::default(),
        }
    }
}

/// Information on the call as combined by the annotator.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone, Default)]
pub struct CallInfo {
    /// The genotype, if applicable, e.g., "0/1", "./1", "."
    pub genotype: Option<String>,
    /// Genotype quality score, if applicable
    pub quality: Option<f32>,
}

/// Definition of a sequence variant with per-sample genotype calls.
///
/// This uses a subset/specialization of what is described by the VCF standard
/// for the purpose of running SV queries in `varfish-server-worker`.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone)]
pub struct SequenceVariant {
    /// Chromosome name
    pub chrom: String,
    /// 1-based start position of the variant (or position on first chromosome
    /// for break-ends)
    pub pos: i32,
    /// Reference allele.
    pub reference: String,
    /// Alternative allele.
    pub alternative: String,

    /// Mapping of sample to genotype information for the SV.
    pub call_info: indexmap::IndexMap<String, CallInfo>,
}

impl SequenceVariant {
    /// Convert from VCF record.
    pub fn from_vcf(record: &vcf::Record, header: &vcf::Header) -> Result<Self, anyhow::Error> {
        let chrom = record.chromosome().to_string();
        let pos: usize = record.position().into();
        let pos = pos as i32;

        let reference = record.reference_bases().to_string();
        let alternative = record.alternate_bases()[0].to_string();

        let call_info = Self::build_call_info(record, header)?;

        Ok(Self {
            chrom,
            pos,
            call_info,
            reference,
            alternative,
        })
    }

    /// Build call information.
    fn build_call_info(
        record: &vcf::Record,
        header: &vcf::Header,
    ) -> Result<indexmap::IndexMap<String, CallInfo>, anyhow::Error> {
        let mut result = indexmap::IndexMap::new();

        for (name, sample) in header
            .sample_names()
            .iter()
            .zip(record.genotypes().values())
        {
            let genotype = if let Some(Some(vcf::record::genotypes::sample::Value::String(gt))) =
                sample.get(&vcf::record::genotypes::keys::key::GENOTYPE)
            {
                Some(gt.clone())
            } else {
                None
            };
            let quality = if let Some(Some(vcf::record::genotypes::sample::Value::Float(quality))) =
                sample.get(&vcf::record::genotypes::keys::key::CONDITIONAL_GENOTYPE_QUALITY)
            {
                Some(*quality)
            } else {
                None
            };

            result.insert(name.clone(), CallInfo { genotype, quality });
        }

        Ok(result)
    }
}

#[cfg(test)]
pub mod test {
    use rstest::rstest;

    #[rstest]
    #[case("tests/seqvars/query/empty.json")]
    #[case("tests/seqvars/query/full.json")]
    #[case("tests/seqvars/query/with_extra.json")]
    pub fn smoke_test_load(#[case] path: &str) -> Result<(), anyhow::Error> {
        let query: super::CaseQuery = serde_json::from_reader(std::fs::File::open(path)?)?;

        insta::assert_yaml_snapshot!(&query);

        Ok(())
    }
}
