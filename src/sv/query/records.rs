//! Records for I/O in `sv query`.

use std::collections::HashMap;

use serde::{
    de::{self, IntoDeserializer},
    Deserialize, Deserializer, Serialize,
};

use super::schema::{
    CallInfo as SchemaCallInfo, StrandOrientation, StructuralVariant as SchemaStructuralVariant,
    SvSubType, SvType,
};

/// Structural variant as written out by VarFish Annotator.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub struct StructuralVariant {
    /// Genome release
    pub release: String,
    /// Chromosome name
    pub chromosome: String,
    /// Chromosome number
    pub chromosome_no: u32,
    /// UCSC bin on chromosome
    pub bin: u32,
    /// 1-based start position of the variant (or position on first chromosome for
    /// break-ends)
    pub start: u32,
    /// Type of the structural variant
    #[serde(deserialize_with = "deserialize_sv_type")]
    pub sv_type: SvType,
    /// Sub type of the structural variant
    #[serde(deserialize_with = "deserialize_sv_sub_type")]
    pub sv_sub_type: SvSubType,
    /// Potentially the second involved chromosome
    pub chromosome2: Option<String>,
    /// Chromosome number of second chromosome
    pub chromosome_no2: u32,
    /// UCSC bin on chromosome2
    pub bin2: u32,
    /// End position (position on second chromosome for break-ends)
    pub end: u32,
    /// The strand orientation of the structural variant, if applicable.
    #[serde(deserialize_with = "deserialize_pe_orientation")]
    pub pe_orientation: StrandOrientation,
    /// Mapping of sample to genotype information for the SV.
    #[serde(deserialize_with = "deserialize_genotype")]
    pub genotype: HashMap<String, Genotype>,
}

impl From<StructuralVariant> for SchemaStructuralVariant {
    fn from(val: StructuralVariant) -> Self {
        fn to_u32(val: Option<i32>) -> Option<u32> {
            val.and_then(|val| if val >= 0 { Some(val as u32) } else { None })
        }

        SchemaStructuralVariant {
            chrom: val.chromosome,
            pos: val.start,
            sv_type: val.sv_type,
            sv_sub_type: val.sv_sub_type,
            chrom2: val.chromosome2,
            end: val.end,
            strand_orientation: Some(val.pe_orientation),
            call_info: HashMap::from_iter(val.genotype.into_iter().map(|(k, v)| {
                (
                    k,
                    SchemaCallInfo {
                        genotype: v.gt,
                        quality: v.gq,
                        paired_end_cov: to_u32(v.pec),
                        paired_end_var: to_u32(v.pev),
                        split_read_cov: to_u32(v.src),
                        split_read_var: to_u32(v.srv),
                        copy_number: to_u32(v.cn),
                        average_normalized_cov: v.anc,
                        point_count: to_u32(v.pc),
                        average_mapping_quality: v.amq,
                    },
                )
            })),
        }
    }
}

fn deserialize_sv_type<'de, D>(deserializer: D) -> Result<SvType, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    let end = s.find('_').unwrap_or(s.len());
    SvType::deserialize(s[..end].into_deserializer())
}

fn deserialize_sv_sub_type<'de, D>(deserializer: D) -> Result<SvSubType, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    SvSubType::deserialize(s.replace('_', ":").into_deserializer())
}

fn deserialize_pe_orientation<'de, D>(deserializer: D) -> Result<StrandOrientation, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    if s.eq(".") {
        Ok(StrandOrientation::NotApplicable)
    } else {
        StrandOrientation::deserialize(s.into_deserializer())
    }
}

fn deserialize_genotype<'de, D>(deserializer: D) -> Result<HashMap<String, Genotype>, D::Error>
where
    D: Deserializer<'de>,
{
    let buf = String::deserialize(deserializer)?;
    let json = buf.replace("\"\"\"", "\"");
    serde_json::from_str(&json).map_err(de::Error::custom)
}

/// Information on the call as combined by the annotator.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub struct Genotype {
    /// The genotype, if applicable, e.g., "0/1", "./1", "."
    pub gt: Option<String>,
    /// Genotype quality score, if applicable
    pub gq: Option<f32>,
    /// Paired-end coverage, if applicable
    pub pec: Option<i32>,
    /// Paired-end variant support, if applicable
    pub pev: Option<i32>,
    /// Split-read coverage, if applicable
    pub src: Option<i32>,
    /// Split-read variant support, if applicable
    pub srv: Option<i32>,
    /// Integer copy number estimate, if applicable
    pub cn: Option<i32>,
    /// Average normalized coverage, if applicable
    pub anc: Option<f32>,
    /// Number of buckets/targets supporting the CNV call, if applicable
    pub pc: Option<i32>,
    /// Average mapping quality, if applicable
    pub amq: Option<f32>,
}
