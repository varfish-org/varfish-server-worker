//! Records for I/O in `sv query`.

use std::collections::HashMap;

use serde::{de, Deserialize, Deserializer, Serialize};

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
    /// 1-based start position of the variant (or position on first chromosome for
    /// break-ends)
    pub start: u32,
    /// Type of the structural variant
    pub sv_type: SvType,
    /// Sub type of the structural variant
    pub sv_sub_type: SvSubType,
    /// Potentially the second involved chromosome
    pub chromosome2: Option<String>,
    /// End position (position on second chromosome for break-ends)
    pub end: u32,
    /// The strand orientation of the structural variant, if applicable.
    pub strand_orientation: StrandOrientation,
    /// Mapping of sample to genotype information for the SV.
    #[serde(deserialize_with = "deserialize_genotype")]
    pub genotype: HashMap<String, Genotype>,
}

impl Into<SchemaStructuralVariant> for StructuralVariant {
    fn into(self) -> SchemaStructuralVariant {
        SchemaStructuralVariant {
            chrom: self.chromosome,
            pos: self.start,
            sv_type: self.sv_type,
            sv_sub_type: self.sv_sub_type,
            chrom2: self.chromosome2,
            end: self.end,
            strand_orientation: Some(self.strand_orientation),
            call_info: HashMap::from_iter(self.genotype.into_iter().map(|(k, v)| {
                (
                    k,
                    SchemaCallInfo {
                        genotype: v.gt,
                        quality: v.gq,
                        paired_end_cov: v.pec,
                        paired_end_var: v.pev,
                        split_read_cov: v.src,
                        split_read_var: v.srv,
                        copy_number: v.cn,
                        average_normalized_cov: v.anc,
                        point_count: v.pc,
                        average_mapping_quality: v.amq,
                    },
                )
            })),
        }
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
    pub pec: Option<u32>,
    /// Paired-end variant support, if applicable
    pub pev: Option<u32>,
    /// Split-read coverage, if applicable
    pub src: Option<u32>,
    /// Split-read variant support, if applicable
    pub srv: Option<u32>,
    /// Integer copy number estimate, if applicable
    pub cn: Option<u32>,
    /// Average normalized coverage, if applicable
    pub anc: Option<f32>,
    /// Number of buckets/targets supporting the CNV call, if applicable
    pub pc: Option<u32>,
    /// Average mapping quality, if applicable
    pub amq: Option<f32>,
}
