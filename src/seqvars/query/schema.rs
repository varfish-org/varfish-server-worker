//! Supporting code for seqvar query definition.

use noodles::vcf;

use crate::common::genotype_to_string;

/// Enumeration for recessive mode queries.
#[derive(
    serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Clone, Copy,
)]
#[serde(rename_all = "kebab-case")]
pub enum RecessiveMode {
    /// Recessive.
    Recessive,
    /// Compound recessive.
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
#[serde(rename_all = "kebab-case")]
pub enum FailChoice {
    /// Ignore failure.
    #[default]
    Ignore,
    /// Drop whole variant.
    #[serde(rename = "drop-variant")]
    Drop,
    /// Interpret as "no-call".
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
#[serde(rename_all = "kebab-case")]
pub enum GenotypeChoice {
    /// Any genotype.
    #[default]
    Any,
    /// Ref. genotype.
    Ref,
    /// Het. genotype.
    Het,
    /// Hom. genotype.
    Hom,
    /// Non-hom. genotype.
    NonHom,
    /// Variant genotype.
    Variant,
    /// Index in comp. het. recessive inheritance.
    ComphetIndex,
    /// Index in recessive inheritance.
    RecessiveIndex,
    /// Parent in recessive inheritance.
    RecessiveParent,
}

impl GenotypeChoice {
    /// Return wehther the genotype choice matches the genotype string.
    ///
    /// Note that we assume properly ingested VCFs with only one alternate allele.
    /// The valid genotype strings have the form "<VAL>/<VAL>", "<VAL>|<VAL>" or
    /// "<VAL>" with "<VAL>" being one of "0", "1", and ".".
    pub fn matches(&self, gt_str: &str) -> Result<bool, anyhow::Error> {
        let gt_str = if gt_str.starts_with('/') || gt_str.starts_with('|') {
            &gt_str[1..]
        } else {
            gt_str
        };
        Ok(match self {
            GenotypeChoice::Any => true,
            GenotypeChoice::Ref => ["0", "0|0", "0/0"].contains(&gt_str),
            GenotypeChoice::Het => ["0/1", "0|1", "1/0", "1|0"].contains(&gt_str),
            GenotypeChoice::Hom => ["1", "1/1", "1|1"].contains(&gt_str),
            GenotypeChoice::NonHom => !["1", "1/1", "1|1"].contains(&gt_str),
            GenotypeChoice::Variant => {
                ["1", "0/1", "0|1", "1/0", "1|0", "1|1", "1/1"].contains(&gt_str)
            }
            GenotypeChoice::ComphetIndex
            | GenotypeChoice::RecessiveIndex
            | GenotypeChoice::RecessiveParent => {
                anyhow::bail!("recessive marker is not a genotype choice")
            }
        })
    }
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
    /// Minimal allele balance for het. variants.
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

serde_with::with_prefix!(prefix_clinvar "clinvar_");
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone)]
#[serde(default)]
pub struct ClinVarOptions {
    /// Whether to include benign ClinVar variants.
    pub include_benign: bool,
    /// Whether to include pathogenic ClinVar variants.
    pub include_pathogenic: bool,
    /// Whether to include likely benign ClinVar variants.
    pub include_likely_benign: bool,
    /// Whether to include likely pathogenic ClinVar variants.
    pub include_likely_pathogenic: bool,
    /// Whether to include uncertain significance ClinVar variants.
    pub include_uncertain_significance: bool,
    /// Whether to include conflicting interpretation ClinVar variants.
    pub include_conflicting_classifications: bool,
}

serde_with::with_prefix!(prefix_gnomad "gnomad_");
#[derive(serde::Serialize, serde::Deserialize, Default, PartialEq, Debug, Clone)]
#[serde(default)]
pub struct GnomadOptions {
    /// Whether to enable filtration by gnomAD exomes.
    pub exomes_enabled: bool,
    /// Whether to enable filtration by gnomAD genomes
    pub genomes_enabled: bool,

    /// Maximal frequency in gnomAD exomes.
    pub exomes_frequency: Option<f32>,
    /// Maximal number of heterozygous carriers in gnomAD exomes.
    pub exomes_heterozygous: Option<i32>,
    /// Maximal number of homozygous carriers in gnomAD exomes.
    pub exomes_homozygous: Option<i32>,
    /// Maximal number of hemizygous carriers in gnomAD exomes.
    pub exomes_hemizygous: Option<i32>,

    /// Maximal frequency in gnomAD genomes.
    pub genomes_frequency: Option<f32>,
    /// Maximal number of heterozygous carriers in gnomAD genomes.
    pub genomes_heterozygous: Option<i32>,
    /// Maximal number of homozygous carriers in gnomAD genomes.
    pub genomes_homozygous: Option<i32>,
    /// Maximal number of hemizygous carriers in gnomAD genomes.
    pub genomes_hemizygous: Option<i32>,
}

/// Data structure with a single query.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone)]
#[serde(default)]
pub struct CaseQuery {
    /// Molecular consequences to consider.
    pub consequences: Vec<mehari::annotate::seqvars::ann::Consequence>,

    /// Quality settings for each individual.
    pub quality: indexmap::IndexMap<String, QualitySettings>,
    /// Genotype choice for each individual.
    pub genotype: indexmap::IndexMap<String, Option<GenotypeChoice>>,

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

    #[serde(flatten, with = "prefix_clinvar")]
    pub clinvar: ClinVarOptions,

    #[serde(flatten, with = "prefix_gnomad")]
    pub gnomad: GnomadOptions,

    /// Whether to enable filtration by 1000 Genomes.
    pub inhouse_enabled: bool,
    /// Whether to enable filtration by mtDB.
    pub helixmtdb_enabled: bool,

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
    pub helixmtdb_heteroplasmic: Option<i32>,
    /// Maximal number of homozygous carriers in HelixMtDb.
    pub helixmtdb_homoplasmic: Option<i32>,
}

impl Default for ClinVarOptions {
    fn default() -> Self {
        Self {
            include_benign: true,
            include_pathogenic: true,
            include_likely_benign: true,
            include_likely_pathogenic: true,
            include_uncertain_significance: true,
            include_conflicting_classifications: true,
        }
    }
}

impl Default for CaseQuery {
    /// Returns default values for a `CaseQuery` which makes all variants pass.
    fn default() -> Self {
        Self {
            consequences: mehari::annotate::seqvars::ann::Consequence::all(),
            inhouse_enabled: Default::default(),
            helixmtdb_enabled: Default::default(),
            quality: Default::default(),
            genotype: Default::default(),
            transcripts_coding: true,
            transcripts_noncoding: true,
            var_type_snv: true,
            var_type_indel: true,
            var_type_mnv: true,
            max_exon_dist: Default::default(),
            gene_allowlist: Default::default(),
            genomic_regions: Default::default(),
            require_in_clinvar: Default::default(),
            clinvar: Default::default(),
            gnomad: Default::default(),
            inhouse_carriers: Default::default(),
            inhouse_heterozygous: Default::default(),
            inhouse_homozygous: Default::default(),
            inhouse_hemizygous: Default::default(),
            helixmtdb_frequency: Default::default(),
            helixmtdb_heteroplasmic: Default::default(),
            helixmtdb_homoplasmic: Default::default(),
        }
    }
}

impl CaseQuery {
    /// Return whether recessive mode has been enabled.
    pub fn recessive_mode(&self) -> bool {
        self.genotype.values().any(|gt_choice| {
            matches!(
                gt_choice,
                Some(
                    GenotypeChoice::ComphetIndex
                        | GenotypeChoice::RecessiveIndex
                        | GenotypeChoice::RecessiveParent
                )
            )
        })
    }

    /// Returns name of the recessive index sample.
    pub fn index_sample(&self) -> Option<String> {
        for (name, gt_choice) in &self.genotype {
            match gt_choice {
                Some(GenotypeChoice::ComphetIndex | GenotypeChoice::RecessiveIndex) => {
                    return Some(name.clone())
                }
                _ => continue,
            }
        }

        None
    }
}
/// Information on the call as written out by ingest.
///
/// Note that the ingested files have exactly one alternate allele.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone, Default)]
pub struct CallInfo {
    /// The genotype, if applicable, e.g., "0/1"
    pub genotype: Option<String>,
    /// Genotype quality score, if applicable
    pub quality: Option<f32>,
    /// Total read coverage at site in the sample.
    pub dp: Option<i32>,
    /// Alternate allele depth for the single allele in the sample.
    pub ad: Option<i32>,
    /// Physical phasing ID for this sample.
    pub phasing_id: Option<i32>,
}

/// Definition of a sequence variant with per-sample genotype calls.
///
/// This uses a subset/specialization of what is described by the VCF standard
/// for the purpose of running SV queries in `varfish-server-worker`.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Debug, Clone, Default)]
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

    /// Variant effect annotation.
    pub ann_fields: Vec<mehari::annotate::seqvars::ann::AnnField>,

    /// Number of alleles in gnomAD exomes (not for chrMT).
    pub gnomad_exomes_an: i32,
    /// Number of homozygous carriers in gnomAD exomes (not for chrMT).
    pub gnomad_exomes_hom: i32,
    /// Number of heterozygous carriers in gnomAD exomes (not for chrMT).
    pub gnomad_exomes_het: i32,
    /// Number of hemizygous carriers in gnomAD exomes (not for chrMT).
    pub gnomad_exomes_hemi: i32,

    /// Number of alleles in gnomAD genomes (also for chrMT).
    pub gnomad_genomes_an: i32,
    /// Number of homozygous carriers in gnomAD genomes (also for chrMT).
    pub gnomad_genomes_hom: i32,
    /// Number of heterozygous carriers in gnomAD genomes (also for chrMT).
    pub gnomad_genomes_het: i32,
    /// Number of hemizygous carriers in gnomAD genomes (not for chrMT).
    pub gnomad_genomes_hemi: i32,

    /// Number of alleles in HelixMtDb cohort (only chrMT).
    pub helix_an: i32,
    /// Number of homoplasmic carriers in HelixMtDb cohort (only chrMT).
    pub helix_hom: i32,
    /// Number of heteroplasmic carriers in HelixMtDb cohort (only chrMT).
    pub helix_het: i32,

    /// Number of in-house alleles (also for chrMT).
    pub inhouse_an: i32,
    /// Number of homozygous carriers in in-house cohort (also for chrMT).
    pub inhouse_hom: i32,
    /// Number of heterozygous carriers in in-house cohort (also for chrMT).
    pub inhouse_het: i32,
    /// Number of hemizygous carriers in in-house cohort (not for chrMT).
    pub inhouse_hemi: i32,

    /// Mapping of sample to genotype information for the SV.
    pub call_info: indexmap::IndexMap<String, CallInfo>,
}

impl SequenceVariant {
    /// Convert from VCF record.
    pub fn from_vcf(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<Self, anyhow::Error> {
        let chrom = record.reference_sequence_name().to_string();
        let pos: usize = record.variant_start().expect("empty variant_start?").into();
        let pos = pos as i32;

        let reference = record.reference_bases().to_string();
        let alternative = record
            .alternate_bases()
            .as_ref()
            .iter()
            .next()
            .expect("no alternate base?")
            .to_string();

        let call_info = Self::build_call_info(record, header)?;
        let ann_fields = Self::extract_ann_fields(record)?;

        let result = Self {
            chrom,
            pos,
            reference,
            alternative,
            call_info,
            ann_fields,
            ..Default::default()
        };

        Self::with_freqs(result, record)
    }

    /// Build call information.
    fn build_call_info(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<indexmap::IndexMap<String, CallInfo>, anyhow::Error> {
        use vcf::variant::record::samples::keys::key;
        let mut result = indexmap::IndexMap::new();

        for (name, sample) in header.sample_names().iter().zip(record.samples().values()) {
            let genotype = if let Some(Some(
                vcf::variant::record_buf::samples::sample::value::Value::Genotype(gt),
            )) = sample.get(key::GENOTYPE)
            {
                Some(genotype_to_string(&gt)?)
            } else {
                None
            };
            let quality = if let Some(Some(
                vcf::variant::record_buf::samples::sample::value::Value::Integer(quality),
            )) = sample.get(key::CONDITIONAL_GENOTYPE_QUALITY)
            {
                Some(*quality as f32)
            } else {
                None
            };
            let dp = if let Some(Some(
                vcf::variant::record_buf::samples::sample::value::Value::Integer(dp),
            )) = sample.get(key::READ_DEPTH)
            {
                Some(*dp)
            } else {
                None
            };
            let ad =
                if let Some(Some(vcf::variant::record_buf::samples::sample::value::Value::Array(
                    vcf::variant::record_buf::samples::sample::value::Array::Integer(ad),
                ))) = sample.get(key::READ_DEPTHS)
                {
                    Some(ad[1].expect("empty AD?"))
                } else {
                    None
                };
            let phase_set = if let Some(Some(
                vcf::variant::record_buf::samples::sample::value::Value::Integer(id),
            )) = sample.get(key::PHASE_SET)
            {
                Some(*id)
            } else {
                None
            };

            result.insert(
                name.clone(),
                CallInfo {
                    genotype,
                    quality,
                    dp,
                    ad,
                    phasing_id: phase_set,
                },
            );
        }

        Ok(result)
    }

    /// Extract `INFO/ANN` entries
    fn extract_ann_fields(
        record: &vcf::variant::RecordBuf,
    ) -> Result<Vec<mehari::annotate::seqvars::ann::AnnField>, anyhow::Error> {
        if let Some(Some(ann)) = record.info().get("ANN") {
            if let vcf::variant::record_buf::info::field::Value::Array(
                vcf::variant::record_buf::info::field::value::Array::String(ann),
            ) = ann
            {
                ann.iter()
                    .flatten()
                    .map(|s| s.parse::<mehari::annotate::seqvars::ann::AnnField>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| anyhow::anyhow!("problem parsing ANN: {}", e))
            } else {
                anyhow::bail!("invalid type of INFO/ANN")
            }
        } else {
            Ok(Vec::default())
        }
    }

    /// Copy the frequencies from `record` to `result`.
    fn with_freqs(
        result: SequenceVariant,
        record: &vcf::variant::RecordBuf,
    ) -> Result<SequenceVariant, anyhow::Error> {
        use vcf::variant::record_buf::info::field::Value;

        macro_rules! extract_key {
            ($key:ident) => {
                let $key =
                    if let Some(Some(Value::Integer($key))) = record.info().get(stringify!($key)) {
                        *$key
                    } else {
                        0
                    };
            };
        }

        extract_key!(gnomad_exomes_an);
        extract_key!(gnomad_exomes_hom);
        extract_key!(gnomad_exomes_het);
        extract_key!(gnomad_exomes_hemi);

        extract_key!(gnomad_genomes_an);
        extract_key!(gnomad_genomes_hom);
        extract_key!(gnomad_genomes_het);
        extract_key!(gnomad_genomes_hemi);

        extract_key!(helix_an);
        extract_key!(helix_hom);
        extract_key!(helix_het);

        Ok(SequenceVariant {
            gnomad_exomes_an,
            gnomad_exomes_hom,
            gnomad_exomes_het,
            gnomad_exomes_hemi,
            gnomad_genomes_an,
            gnomad_genomes_hom,
            gnomad_genomes_het,
            gnomad_genomes_hemi,
            helix_an,
            helix_hom,
            helix_het,
            ..result
        })
    }

    /// Return allele frequency in gnomAD exomes.
    pub fn gnomad_exomes_af(&self) -> f32 {
        if self.gnomad_exomes_an == 0 {
            return 0f32;
        }
        let an = self.gnomad_exomes_an as f32;
        let hom = self.gnomad_exomes_hom as f32;
        let het = self.gnomad_exomes_het as f32;
        let hemi = self.gnomad_exomes_hemi as f32;
        (2.0 * hom + het + hemi) / an
    }

    /// Return allele frequency in gnomAD genomes.
    pub fn gnomad_genomes_af(&self) -> f32 {
        if self.gnomad_genomes_an == 0 {
            return 0f32;
        }
        let an = self.gnomad_genomes_an as f32;
        let hom = self.gnomad_genomes_hom as f32;
        let het = self.gnomad_genomes_het as f32;
        let hemi = self.gnomad_genomes_hemi as f32;
        (2.0 * hom + het + hemi) / an
    }

    /// Return allele frequency in HelixMtDb.
    pub fn helixmtdb_af(&self) -> f32 {
        if self.helix_an == 0 {
            return 0f32;
        }
        let an = self.helix_an as f32;
        let hom = self.helix_hom as f32;
        let het = self.helix_het as f32;
        (2.0 * hom + het) / an
    }
}

#[cfg(test)]
pub mod test {
    use noodles::vcf;
    use rstest::rstest;

    pub mod genotype_choice {
        use super::super::GenotypeChoice::{self, *};
        use rstest::rstest;

        #[rstest]
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
            for gt_str in gt_strs {
                assert_eq!(
                    genotype_choice.matches(gt_str).unwrap(),
                    expected,
                    "{:?} {}",
                    genotype_choice,
                    gt_str
                );
            }
        }
    }

    #[rstest]
    #[case("tests/seqvars/query/empty.json")]
    #[case("tests/seqvars/query/full.json")]
    #[case("tests/seqvars/query/with_extra.json")]
    pub fn smoke_test_load(#[case] path_input: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path_input.split('/').last().unwrap());

        let query: super::CaseQuery = serde_json::from_reader(std::fs::File::open(path_input)?)?;

        insta::assert_yaml_snapshot!(&query);

        Ok(())
    }

    #[rstest::rstest]
    #[case("tests/seqvars/query/Case_1.ingested.vcf")]
    #[case("tests/seqvars/query/dragen.ingested.vcf")]
    pub fn sequence_variant_from_vcf(#[case] path_input: &str) -> Result<(), anyhow::Error> {
        use noodles::vcf::variant::RecordBuf;

        mehari::common::set_snapshot_suffix!("{}", path_input.split('/').last().unwrap());

        let mut vcf_reader = vcf::io::reader::Builder::default()
            .build_from_path(path_input)
            .unwrap();
        let header = vcf_reader.read_header()?;

        for record in vcf_reader.records() {
            let record = record?;
            let seqvar = super::SequenceVariant::from_vcf(
                &RecordBuf::try_from_variant_record(&header, &record)?,
                &header,
            )?;

            insta::assert_yaml_snapshot!(&seqvar);
        }

        Ok(())
    }
}
