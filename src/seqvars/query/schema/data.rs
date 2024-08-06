//! Code for representing variants internally, corresponds to what is written
//! out by `seqvars ingest`.

use noodles::vcf;

use crate::common::genotype_to_string;

/// Trait for attempting conversion from VCF record.
pub trait TryFromVcf: Sized {
    /// Error type to use.
    type Error;

    /// Convert from VCF record.
    ///
    /// # Arguments
    ///
    /// * `record` - VCF record.
    /// * `header` - VCF header.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be converted.
    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<Self, Self::Error>;
}

/// Information on the call as written out by ingest.
///
/// Corresponds to `FORMAT/*` in VCF.
///
/// Note that the ingested files have exactly one alternate allele.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CallInfo {
    /// The sample name.
    pub sample: String,
    /// The genotype, if applicable, e.g., "0/1"
    pub genotype: Option<String>,
    /// Genotype quality score, if applicable
    pub gq: Option<f32>,
    /// Total read coverage at site in the sample.
    pub dp: Option<i32>,
    /// Alternate allele depth for the single allele in the sample.
    pub ad: Option<i32>,
    /// Physical phasing ID for this sample.
    pub ps: Option<i32>,
}

impl Eq for CallInfo {}

/// pub trait for total allele counts.
pub trait An {
    /// Number of covered alleles.
    fn an(&self) -> i32;
}

/// pub trait for variant alternate allele counts.
pub trait Ac {
    /// Number of homozygous/homoplasmic carriers.
    fn hom(&self) -> i32;
    /// Number of heterozygous/heteroplasmic carriers.
    fn het(&self) -> i32;
    /// Number of total alternate alleles.
    fn ac(&self) -> i32;
}

/// pub trait for allele frequency.
pub trait Af {
    /// Allele frequency.
    fn af(&self) -> f32;
}

/// Blanket implementation of `Af` for any type that implements `Ac` and `An`.
impl<T: Ac + An> Af for T {
    fn af(&self) -> f32 {
        if self.an() == 0 {
            0.0
        } else {
            self.ac() as f32 / self.an() as f32
        }
    }
}

/// pub trait for carrier count.
pub trait Carriers {
    /// Total number of carriers.
    fn carriers(&self) -> i32;
}

/// pub trait for nuclear variant alternate allele including hemizygous.
pub trait Hemi: Ac {
    /// Number of hemizygous carriers.
    fn hemi(&self) -> i32;
}

/// gnomAD Population frequencies for variants on nuclear chromosomes.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct NuclearFrequencies {
    /// Number of alleles.
    pub an: i32,
    /// Number of homozygous carriers.
    pub hom: i32,
    /// Number of heterozygous carriers.
    pub het: i32,
    /// Number of hemizygous carriers.
    pub hemi: i32,
}

impl An for NuclearFrequencies {
    fn an(&self) -> i32 {
        self.an
    }
}

impl Ac for NuclearFrequencies {
    fn hom(&self) -> i32 {
        self.hom
    }
    fn het(&self) -> i32 {
        self.het
    }
    fn ac(&self) -> i32 {
        2 * self.hom + self.het + self.hemi
    }
}

impl Hemi for NuclearFrequencies {
    fn hemi(&self) -> i32 {
        self.hemi
    }
}

impl Carriers for NuclearFrequencies {
    fn carriers(&self) -> i32 {
        self.hom + self.het + self.hemi
    }
}

/// Population frequencies for variants on the mitochondrial chromosome.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct MitochondrialFrequencies {
    /// Number of alleles.
    pub an: i32,
    /// Number of homoplasmic carriers.
    pub hom: i32,
    /// Number of heteroplasmic carriers.
    pub het: i32,
}

impl An for MitochondrialFrequencies {
    fn an(&self) -> i32 {
        self.an
    }
}

impl Ac for MitochondrialFrequencies {
    fn hom(&self) -> i32 {
        self.hom
    }
    fn het(&self) -> i32 {
        self.het
    }
    fn ac(&self) -> i32 {
        self.hom + self.het
    }
}

impl Carriers for MitochondrialFrequencies {
    fn carriers(&self) -> i32 {
        self.hom + self.het
    }
}

/// In-house frequencies.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct InHouseFrequencies {
    // Number of alleles.
    pub an: i32,
    /// Number of homozygous carriers.
    pub hom: i32,
    /// Number of heterozygous carriers.
    pub het: i32,
    /// Number of hemizygous carriers.
    pub hemi: i32,
}

impl Ac for InHouseFrequencies {
    fn hom(&self) -> i32 {
        self.hom
    }
    fn het(&self) -> i32 {
        self.het
    }
    fn ac(&self) -> i32 {
        2 * self.hom + self.het + self.hemi
    }
}

impl An for InHouseFrequencies {
    fn an(&self) -> i32 {
        self.an
    }
}

impl Carriers for InHouseFrequencies {
    fn carriers(&self) -> i32 {
        self.hom + self.het + self.hemi
    }
}

/// Population frequencies for a variant.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct PopulationFrequencies {
    /// Frequencies in gnomAD-exomes.
    pub gnomad_exomes: NuclearFrequencies,
    /// Frequencies in gnomAD-genomes.
    pub gnomad_genomes: NuclearFrequencies,
    /// Frequencies in gnomAD-MT.
    pub gnomad_mtdna: MitochondrialFrequencies,
    /// Frequencies in HelixMtDb.
    pub helixmtdb: MitochondrialFrequencies,
    /// Frequencies in In-house database.
    pub inhouse: InHouseFrequencies,
}

/// Supporting code for `PopulationFrequencies`.
pub(crate) mod population_frequencies {
    /// Error type for `FromVcf` implementation.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        // infallibe
    }
}

impl TryFromVcf for PopulationFrequencies {
    type Error = population_frequencies::Error;

    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        _header: &vcf::Header,
    ) -> Result<Self, Self::Error> {
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

        extract_key!(gnomad_mtdna_an);
        extract_key!(gnomad_mtdna_hom);
        extract_key!(gnomad_mtdna_het);

        extract_key!(helix_an);
        extract_key!(helix_hom);
        extract_key!(helix_het);

        Ok(PopulationFrequencies {
            gnomad_exomes: NuclearFrequencies {
                an: gnomad_exomes_an,
                hom: gnomad_exomes_hom,
                het: gnomad_exomes_het,
                hemi: gnomad_exomes_hemi,
            },
            gnomad_genomes: NuclearFrequencies {
                an: gnomad_genomes_an,
                hom: gnomad_genomes_hom,
                het: gnomad_genomes_het,
                hemi: gnomad_genomes_hemi,
            },
            gnomad_mtdna: MitochondrialFrequencies {
                an: gnomad_mtdna_an,
                hom: gnomad_mtdna_hom,
                het: gnomad_mtdna_het,
            },
            helixmtdb: MitochondrialFrequencies {
                an: helix_an,
                hom: helix_hom,
                het: helix_het,
            },
            inhouse: Default::default(),
        })
    }
}

/// Sequence variant representation VCF-style.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct VcfVariant {
    /// Chromosome.
    pub chrom: String,
    /// Position.
    pub pos: i32,
    /// Reference allele.
    pub ref_allele: String,
    /// Alternate allele.
    pub alt_allele: String,
}

impl From<VcfVariant> for annonars::common::spdi::Var {
    fn from(val: VcfVariant) -> Self {
        annonars::common::spdi::Var::new(
            annonars::common::cli::canonicalize(&val.chrom),
            val.pos,
            val.ref_allele,
            val.alt_allele,
        )
    }
}

/// Supporting code for `VcfVariant`.
pub(crate) mod vcf_variant {
    /// Error type for `FromVcf` implementation.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        #[error("Missing POS value")]
        MissingVariantStart,
        #[error("Missing ALT values")]
        MissingAlternateBases,
    }
}

impl TryFromVcf for VcfVariant {
    type Error = vcf_variant::Error;

    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        _header: &vcf::Header,
    ) -> Result<Self, Self::Error> {
        let chrom = record.reference_sequence_name().to_string();
        let pos = usize::from(
            record
                .variant_start()
                .ok_or(Self::Error::MissingVariantStart)?,
        ) as i32;

        let ref_allele = record.reference_bases().to_string();
        let alt_allele = record
            .alternate_bases()
            .as_ref()
            .iter()
            .next()
            .ok_or(Self::Error::MissingAlternateBases)?
            .to_string();

        Ok(Self {
            chrom,
            pos,
            ref_allele,
            alt_allele,
        })
    }
}

/// Helper type with mapping from sample name to `CallInfo`.
#[derive(Debug, Clone, Default)]
pub struct CallInfos {
    /// Mapping from sample name to `CallInfo`.
    pub call_infos: indexmap::IndexMap<String, CallInfo>,
}

/// Supporting code for `CallInfos`.
pub(crate) mod call_infos {
    use crate::common::GenotypeToStringError;

    /// Error type for `FromVcf` implementation.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        #[error("Empty FORMAT/AD")]
        EmptyFormatAd,
        #[error("Problem parsing genotype: {0}")]
        GenotypeParsing(#[from] GenotypeToStringError),
    }
}

impl TryFromVcf for CallInfos {
    type Error = call_infos::Error;

    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<Self, Self::Error> {
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
                    Some(ad[1].ok_or(Self::Error::EmptyFormatAd)?)
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
                    sample: name.clone(),
                    genotype,
                    gq: quality,
                    dp,
                    ad,
                    ps: phase_set,
                },
            );
        }

        Ok(CallInfos { call_infos: result })
    }
}

/// Helper type with `AnnField` records.
#[derive(Debug, Clone, Default)]
pub struct AnnFields {
    /// Mapping from sample name to `CallInfo`.
    pub ann_fields: Vec<mehari::annotate::seqvars::ann::AnnField>,
}

/// Supporting code for `AnnFields`.
pub(crate) mod ann_fields {
    /// Error type for `FromVcf` implementation.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        #[error("Problem parsing INFO/ANN: {0}")]
        Parsing(String),
        #[error("Invalid type of INFO/ANN")]
        InvalidTypeInfoAnn,
    }
}

impl TryFromVcf for AnnFields {
    type Error = ann_fields::Error;

    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        _header: &vcf::Header,
    ) -> Result<AnnFields, ann_fields::Error> {
        if let Some(Some(ann)) = record.info().get("ANN") {
            if let vcf::variant::record_buf::info::field::Value::Array(
                vcf::variant::record_buf::info::field::value::Array::String(ann),
            ) = ann
            {
                Ok(AnnFields {
                    ann_fields: ann
                        .iter()
                        .flatten()
                        .map(|s| s.parse::<mehari::annotate::seqvars::ann::AnnField>())
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|e| Self::Error::Parsing(format!("{}", e)))?,
                })
            } else {
                Err(ann_fields::Error::InvalidTypeInfoAnn)
            }
        } else {
            Ok(Default::default())
        }
    }
}

/// Sequence variant record.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct VariantRecord {
    /// VCF-style variant.
    pub vcf_variant: VcfVariant,
    /// Call information.
    pub call_infos: indexmap::IndexMap<String, CallInfo>,
    /// Annotation fields.
    pub ann_fields: Vec<mehari::annotate::seqvars::ann::AnnField>,
    /// Population frequencies.
    pub population_frequencies: PopulationFrequencies,
}

/// Supporting code for `VariantRecord`.
pub(crate) mod variant_record {
    /// Error type for `FromVcf` implementation.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        #[error("Problem with variant description: {0:?}")]
        VcfVariant(#[from] super::vcf_variant::Error),
        #[error("Problem with call information: {0:?}")]
        CallInfos(#[from] super::call_infos::Error),
        #[error("Problem with annotation fields: {0:?}")]
        AnnFields(#[from] super::ann_fields::Error),
        #[error("Problem with population frequencies: {0:?}")]
        PopulationFrequencies(#[from] super::population_frequencies::Error),
    }
}

impl TryFromVcf for VariantRecord {
    type Error = variant_record::Error;

    /// Convert from VCF record.
    ///
    /// # Arguments
    ///
    /// * `record` - VCF record.
    /// * `header` - VCF header.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be extracted.
    fn try_from_vcf(
        record: &vcf::variant::RecordBuf,
        header: &vcf::Header,
    ) -> Result<Self, Self::Error> {
        let vcf_variant = VcfVariant::try_from_vcf(record, header)?;
        let CallInfos { call_infos } = CallInfos::try_from_vcf(record, header)?;
        let AnnFields { ann_fields } = AnnFields::try_from_vcf(record, header)?;
        let population_frequencies = PopulationFrequencies::try_from_vcf(record, header)?;

        Ok(Self {
            vcf_variant,
            call_infos,
            ann_fields,
            population_frequencies,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[rstest::fixture]
    fn nuclear_frequencies() -> NuclearFrequencies {
        NuclearFrequencies {
            an: 100,
            hom: 10,
            het: 20,
            hemi: 5,
        }
    }

    #[rstest::fixture]
    fn mitochondrial_frequencies() -> MitochondrialFrequencies {
        MitochondrialFrequencies {
            an: 100,
            hom: 10,
            het: 20,
        }
    }

    #[rstest::fixture]
    fn inhouse_frequencies() -> InHouseFrequencies {
        InHouseFrequencies {
            hom: 10,
            het: 20,
            hemi: 5,
            an: 65,
        }
    }

    #[rstest::fixture]
    fn population_frequencies(
        nuclear_frequencies: NuclearFrequencies,
        mitochondrial_frequencies: MitochondrialFrequencies,
        inhouse_frequencies: InHouseFrequencies,
    ) -> PopulationFrequencies {
        PopulationFrequencies {
            gnomad_exomes: nuclear_frequencies.clone(),
            gnomad_genomes: nuclear_frequencies.clone(),
            gnomad_mtdna: mitochondrial_frequencies.clone(),
            helixmtdb: mitochondrial_frequencies.clone(),
            inhouse: inhouse_frequencies.clone(),
        }
    }

    #[rstest::rstest]
    fn test_nuclear_frequencies_an(nuclear_frequencies: NuclearFrequencies) {
        assert_eq!(nuclear_frequencies.an(), 100);
    }

    #[rstest::rstest]
    fn test_nuclear_frequencies_ac(nuclear_frequencies: NuclearFrequencies) {
        assert_eq!(nuclear_frequencies.hom(), 10);
        assert_eq!(nuclear_frequencies.het(), 20);
        assert_eq!(nuclear_frequencies.ac(), 45);
    }

    #[rstest::rstest]
    fn test_nuclear_frequencies_af(nuclear_frequencies: NuclearFrequencies) {
        assert_eq!(nuclear_frequencies.af(), 0.45);
    }

    #[rstest::rstest]
    fn test_nuclear_frequencies_hemi(nuclear_frequencies: NuclearFrequencies) {
        assert_eq!(nuclear_frequencies.hemi(), 5);
    }

    #[rstest::rstest]
    fn test_nuclear_frequencies_carriers(nuclear_frequencies: NuclearFrequencies) {
        assert_eq!(nuclear_frequencies.carriers(), 35);
    }

    #[rstest::rstest]
    fn test_mitochondrial_frequencies_an(mitochondrial_frequencies: MitochondrialFrequencies) {
        assert_eq!(mitochondrial_frequencies.an(), 100);
    }

    #[rstest::rstest]
    fn test_mitochondrial_frequencies_ac(mitochondrial_frequencies: MitochondrialFrequencies) {
        assert_eq!(mitochondrial_frequencies.hom(), 10);
        assert_eq!(mitochondrial_frequencies.het(), 20);
        assert_eq!(mitochondrial_frequencies.ac(), 30);
    }

    #[rstest::rstest]
    fn test_mitochondrial_frequencies_af(mitochondrial_frequencies: MitochondrialFrequencies) {
        assert_eq!(mitochondrial_frequencies.af(), 0.3);
    }

    #[rstest::rstest]
    fn test_mitochondrial_frequencies_carriers(
        mitochondrial_frequencies: MitochondrialFrequencies,
    ) {
        assert_eq!(mitochondrial_frequencies.carriers(), 30);
    }

    #[rstest::rstest]
    fn test_inhouse_frequencies_ac(inhouse_frequencies: InHouseFrequencies) {
        assert_eq!(inhouse_frequencies.hom(), 10);
        assert_eq!(inhouse_frequencies.het(), 20);
        assert_eq!(inhouse_frequencies.ac(), 45);
    }

    #[rstest::rstest]
    fn test_inhouse_frequencies_carriers(inhouse_frequencies: InHouseFrequencies) {
        assert_eq!(inhouse_frequencies.carriers(), 35);
    }

    #[rstest::rstest]
    fn test_population_frequencies(population_frequencies: PopulationFrequencies) {
        assert_eq!(population_frequencies.gnomad_exomes.an, 100);
        assert_eq!(population_frequencies.gnomad_genomes.an, 100);
        assert_eq!(population_frequencies.gnomad_mtdna.an, 100);
        assert_eq!(population_frequencies.helixmtdb.an, 100);
        assert_eq!(population_frequencies.inhouse.an, 65);
    }

    #[rstest::rstest]
    #[case::case_1_ingested_vcf("tests/seqvars/query/Case_1.ingested.vcf")]
    #[case::dragen_ingested_vcf("tests/seqvars/query/dragen.ingested.vcf")]
    pub fn sequence_variant_from_vcf(#[case] path_input: &str) -> Result<(), anyhow::Error> {
        use noodles::vcf::variant::RecordBuf;

        mehari::common::set_snapshot_suffix!("{}", path_input.split('/').last().unwrap());

        let mut vcf_reader = vcf::io::reader::Builder::default()
            .build_from_path(path_input)
            .unwrap();
        let header = vcf_reader.read_header()?;

        for record in vcf_reader.records() {
            let record = record?;
            let seqvar = super::VariantRecord::try_from_vcf(
                &RecordBuf::try_from_variant_record(&header, &record)?,
                &header,
            )?;

            insta::assert_yaml_snapshot!(&seqvar);
        }

        Ok(())
    }
}
