//! Implementation of I/O for `sv build-bgdb`.

use serde::{de::IntoDeserializer, Deserialize, Deserializer, Serialize};

use crate::sv_query::schema::{StrandOrientation, SvType};

/// Representation of the fields from the `StructuralVariant` table from VarFish Server
/// that we need for building the background records.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct FileRecord {
    /// genome build
    pub release: String,
    /// chromosome name
    pub chromosome: String,
    /// UCSC bin
    pub bin: i32,
    /// start position, 1-based
    pub start: i32,
    /// chromosome2 name
    pub chromosome2: String,
    /// end position, 1-based
    pub end: i32,
    /// paired-end orientation
    #[serde(deserialize_with = "from_varfish_pe_orientation")]
    pub pe_orientation: StrandOrientation,
    /// SV type of the record
    #[serde(deserialize_with = "from_varfish_sv_type")]
    pub sv_type: SvType,
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
}

impl FileRecord {
    /// Compute reciprocal overlap between `self` and `other`.
    pub fn overlap(&self, other: &FileRecord) -> f32 {
        let s1 = if self.start > 0 { self.start - 1 } else { 0 };
        let e1 = self.end + 1;
        let s2 = if other.start > 0 { other.start - 1 } else { 0 };
        let e2 = other.end;

        let ovl_s = std::cmp::max(s1, s2);
        let ovl_e = std::cmp::min(e1, e2);
        if ovl_e <= ovl_s {
            0.0
        } else {
            let len1 = (e1 - s1) as f32;
            let len2 = (e2 - s2) as f32;
            let ovl_len = (ovl_e - ovl_s) as f32;
            (ovl_len / len1).min(ovl_len / len2)
        }
    }

    pub fn merge_into(&mut self, other: &FileRecord) {
        self.num_hom_alt += other.num_hom_alt;
        self.num_hom_ref += other.num_hom_alt;
        self.num_het += other.num_het;
        self.num_hemi_alt += other.num_hemi_alt;
        self.num_hemi_ref += other.num_hemi_ref;
    }
}

impl Default for FileRecord {
    fn default() -> Self {
        Self {
            release: "".to_owned(),
            chromosome: "".to_owned(),
            bin: 0,
            start: 0,
            chromosome2: "".to_owned(),
            end: 0,
            pe_orientation: StrandOrientation::NotApplicable,
            sv_type: SvType::Bnd,
            num_hom_alt: 0,
            num_hom_ref: 0,
            num_het: 0,
            num_hemi_alt: 0,
            num_hemi_ref: 0,
        }
    }
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

/// Deserialize "pe_orientation" from VarFish database.
///
/// This function will convert `"."` to `StrandOrientation::NotApplicable`
fn from_varfish_pe_orientation<'de, D>(deserializer: D) -> Result<StrandOrientation, D::Error>
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
