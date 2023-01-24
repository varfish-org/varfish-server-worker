//! Implementation of I/O for `sv build-bgdb`.

use std::collections::HashMap;

use serde::{de::IntoDeserializer, Deserialize, Deserializer};

use crate::sv_query::schema::{SvSubType, SvType};

/// Representation of the per-sample genotype information from VarFish server files.
#[derive(Debug, Deserialize)]
pub struct FileGenotypeInfo {
    /// The genotype, if applicable, e.g., "0/1", "./1", "."
    #[serde(rename = "gt")]
    pub genotype: Option<String>,
    /// Genotype quality score, if applicable
    #[serde(rename = "gq")]
    pub quality: Option<f32>,
    /// Paired-end coverage, if applicable
    #[serde(rename = "pec")]
    pub paired_end_cov: Option<i32>,
    /// Paired-end variant support, if applicable
    #[serde(rename = "pev")]
    pub paired_end_var: Option<i32>,
    /// Split-read coverage, if applicable
    #[serde(rename = "src")]
    pub split_read_cov: Option<i32>,
    /// Split-read variant support, if applicable
    #[serde(rename = "srv")]
    pub split_read_var: Option<i32>,
    /// Integer copy number estimate, if applicable
    #[serde(rename = "cn")]
    pub copy_number: Option<i32>,
    /// Average normalized coverage, if applicable
    #[serde(rename = "anc")]
    pub average_normalized_cov: Option<f32>,
    /// Number of buckets/targets supporting the CNV call, if applicable
    #[serde(rename = "pc")]
    pub point_count: Option<i32>,
    /// Average mapping quality, if applicable
    #[serde(rename = "amq")]
    pub average_mapping_quality: Option<f32>,
}

/// Representation of the fields from the `StructuralVariant` table from VarFish Server
/// that we need for building the background records.
#[derive(Debug, Deserialize)]
pub struct FileRecord {
    /// genome build
    pub release: String,
    /// chromosome name
    pub chromosome: String,
    /// start position, 1-based
    pub start: i32,
    /// chromosome2 name
    pub chromosome2: String,
    /// end position, 1-based
    pub end: i32,
    /// list of variant callers
    pub caller: String,
    /// SV type of the record
    #[serde(deserialize_with = "from_varfish_sv_type")]
    pub sv_type: SvType,
    /// SV sub type of the record
    #[serde(deserialize_with = "from_varfish_sv_sub_type")]
    pub sv_sub_type: SvSubType,
    /// number of hom. alt. carriers
    pub num_hom_alt: u32,
    /// number of hom. ref. carriers
    pub num_hom_ref: u32,
    /// number of het. carriers
    pub num_het: u32,
    /// number of hemi. alt. carriers
    pub num_hemi_alt: u32,
    /// number of hemi. ref. carriers
    pub num_hemi_ref: u32,
    /// genotypes as postgres triple-quoted JSON
    #[serde(deserialize_with = "from_varfish_genotype")]
    pub genotype: HashMap<String, FileGenotypeInfo>,
}

/// Deserialize "sv_type" from VarFish database.
///
/// This function will strip everything after the first underscore.
fn from_varfish_sv_type<'de, D>(deserializer: D) -> Result<SvType, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    let end = s.find('_').unwrap_or(s.len());
    SvType::deserialize(s[..end].into_deserializer())
}

/// Deserialize "sv_sub_type" from VarFish database.
///
/// This function will replace underscores with colons.
fn from_varfish_sv_sub_type<'de, D>(deserializer: D) -> Result<SvSubType, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    SvSubType::deserialize(s.replace('_', ":").into_deserializer())
}

/// Deserialize "genotype" from VarFish database.
///
/// This function will transmogrify the PostgreSQL JSONB syntax with triple quoting
/// to valid JSON and then use `serde_json`.
fn from_varfish_genotype<'de, D>(
    deserializer: D,
) -> Result<HashMap<String, FileGenotypeInfo>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    let json_data = s.replace("\"\"\"", "\"");
    match serde_json::from_str(&json_data) {
        Ok(value) => Ok(value),
        Err(_) => Err(serde::de::Error::custom("Problem deserializing genotype")),
    }
}
