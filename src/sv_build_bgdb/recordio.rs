//! Implementation of I/O for `sv build-bgdb`.

use serde::{de::IntoDeserializer, Deserialize, Deserializer};

use crate::sv_query::schema::{SvSubType, SvType};

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
    pub genotype: String,
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
