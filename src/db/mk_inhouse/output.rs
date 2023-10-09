//! Code for writing the output file.

use serde::{Deserialize, Serialize};

use crate::strucvars::query::schema::{StrandOrientation, SvType};

use super::input::Record as InputRecord;

/// Representation of the fields for the in-house background database.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Record {
    /// chromosome name
    pub chromosome: String,
    /// start position, 1-based
    pub begin: i32,
    /// chromosome2 name
    pub chromosome2: String,
    /// end position, 1-based
    pub end: i32,
    /// paired-end orientation
    pub pe_orientation: StrandOrientation,
    /// type of the SV
    pub sv_type: SvType,
    /// number of overall carriers
    pub carriers: u32,
    /// number of het. carriers
    pub carriers_het: u32,
    /// number of hom. carriers
    pub carriers_hom: u32,
    /// number of hemi. carriers
    pub carriers_hemi: u32,
}

impl Record {
    /// Compute reciprocal overlap between `self` and `other`.
    pub fn overlap(&self, other: &Record) -> f32 {
        let s1 = if self.begin > 0 { self.begin - 1 } else { 0 };
        let e1 = self.end + 1;
        let s2 = if other.begin > 0 { other.begin - 1 } else { 0 };
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

    pub fn merge_into(&mut self, other: &Record) {
        self.carriers += other.carriers;
        self.carriers_het += other.carriers_het;
        self.carriers_hom += other.carriers_hom;
        self.carriers_hemi += other.carriers_hemi;
    }

    pub fn from_db_record(record: InputRecord) -> Self {
        Record {
            chromosome: record.chromosome,
            begin: record.start - 1,
            chromosome2: record.chromosome2,
            end: record.end,
            pe_orientation: record.pe_orientation,
            sv_type: record.sv_type,
            carriers: record.num_het + record.num_hom_alt + record.num_hemi_alt,
            carriers_het: record.num_het,
            carriers_hom: record.num_hom_alt,
            carriers_hemi: record.num_hemi_alt,
        }
    }
}
