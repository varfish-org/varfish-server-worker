use std::ops::Range;

/// This module provides the code for accessing database records.

/// Provide a chromosome-wise coordinate.
pub trait ChromosomeCoordinate {
    fn chromosome(&self) -> &String;
    /// 0-based begin position
    fn begin(&self) -> u32;
    /// 1-based begin position
    fn start(&self) -> u32;
    /// 0/1-based end position
    fn end(&self) -> u32;
}

pub trait BeginEnd {
    /// 0-base begin position
    fn begin(&self) -> u32;
    /// 0-based end position
    fn end(&self) -> u32;
}

pub trait ToInMemory<InMemory> {
    fn to_in_memory(&self) -> Result<Option<InMemory>, anyhow::Error>;
}

pub fn reciprocal_overlap(lhs: &impl BeginEnd, rhs: &Range<u32>) -> f32 {
    let lhs_b = lhs.begin() as u32;
    let lhs_e = lhs.end() as u32;
    let rhs_b = rhs.start.saturating_sub(1);
    let rhs_e = rhs.end;
    let ovl_b = std::cmp::max(lhs_b, rhs_b);
    let ovl_e = std::cmp::min(lhs_e, rhs_e);
    if ovl_b >= ovl_e {
        0f32
    } else {
        let ovl_len = (ovl_e - ovl_b) as f32;
        let x1 = (lhs_e - lhs_b) as f32 / ovl_len;
        let x2 = (rhs_e - rhs_b) as f32 / ovl_len;
        x1.min(x2)
    }
}

/// Store background database counts for a structural variant.
#[derive(Clone, Debug, PartialEq)]
pub struct SvOverlapCounts {
    /// Number of carriers in DGV
    pub dgv_carriers: u32,
    /// Number of carriers in DGV gold standard
    pub dgv_gs_carriers: u32,
    /// Number of carriers in gnomAD SV
    pub gnomad_carriers: u32,
    /// Number of carriers in ExAC
    pub exac_carriers: u32,
    /// Number of carriers in dbVar
    pub dbvar_carriers: u32,
    /// Number of alleles in Thousand Genomes
    pub g1k_alleles: u32,
    /// Number of carriers in inhouse database
    pub inhouse_carriers: u32,
}

/// Count attached to a background record.
pub trait Count {
    fn count(&self) -> usize;
}

/// Records for in-house SV background database.
pub mod bg_sv {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use serde::Deserialize;

    /// Background SV database record to be kept in memory.
    #[derive(Debug)]
    pub struct Record {
        /// The 0-based begin position.
        pub begin: u32,
        /// The 0-based end position.
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// Total number of carriers.
        pub carriers: u32,
        /// Number of het. carriers.
        pub carriers_het: u32,
        /// Number of hom. carriers.
        pub carriers_hom: u32,
        /// Number of hemi. carriers.
        pub carriers_hemi: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// Background SV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub start: u32,
        /// chromosome2 name
        pub chromosome2: String,
        /// end position, 1-based
        pub end: u32,
        /// paired-end orientation
        pub pe_orientation: String,
        /// type of the SV
        pub sv_type: String,
        /// number of overall carriers
        pub carriers: u32,
        /// number of het. carriers
        pub carriers_het: u32,
        /// number of hom. carriers
        pub carriers_hom: u32,
        /// number of hemi. carriers
        pub carriers_hemi: u32,
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            Ok(Some(Record {
                begin: self.start.saturating_sub(1),
                end: self.end,
                sv_type: serde_json::from_str(&format!(
                    "\"{}\"",
                    &self.sv_type.split(':').next().unwrap()
                ))?,
                carriers: self.carriers,
                carriers_het: self.carriers_het,
                carriers_hom: self.carriers_hom,
                carriers_hemi: self.carriers_hemi,
            }))
        }
    }
}

/// Records for the dbVar
pub mod dbvar {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// dbVar database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall carriers
        pub carriers: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// dbVar database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub start: u32,
        /// end position, 1-based
        pub end: u32,
        /// number of overall carriers
        pub num_carriers: u32,
        /// type of the SV
        pub sv_type: String,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.sv_type.split(';').next().unwrap() {
                "alu_insertion"
                | "herv_insertion"
                | "insertion"
                | "line1_insertion"
                | "mobile_element_insertion"
                | "novel_sequence_insertion"
                | "sva_insertion" => SvType::Ins,
                "copy_number_gain" | "duplication" | "tandem_duplication" => SvType::Dup,
                "alu_deletion" | "copy_number_loss" | "deletion" | "herv_deletion"
                | "line1_deletion" | "sva_deletion" => SvType::Del,
                "copy_number_variation" => SvType::Cnv,
                _ => return Err(anyhow!("Unknown SV type: {}", self.sv_type)),
            };
            Ok(Some(Record {
                begin: self.start - 1,
                end: self.end,
                sv_type,
                carriers: self.num_carriers,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }
}

/// Records for gnomAD SV
pub mod gnomad_sv {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// gnomAD SV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall carriers
        pub carriers: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// gnomAD SV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub start: u32,
        /// end position, 1-based
        pub end: u32,
        /// The structural vairant type
        pub svtype: String,
        /// Number of homozygous alternative carriers
        pub n_homalt: u32,
        /// Number of heterozygous carriers
        pub n_het: u32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.svtype.as_str() {
                "CPX" => return Ok(None), // no correspondence
                "CTX" | "BND" => SvType::Bnd,
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                "INS" => SvType::Ins,
                "INV" => SvType::Inv,
                "MCNV" => SvType::Cnv,
                _ => return Err(anyhow!("Unknown SV type: {}", &self.svtype)),
            };
            Ok(Some(Record {
                begin: self.start - 1,
                end: self.end,
                sv_type,
                carriers: self.n_homalt + self.n_het,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }
}
/// Records for Thousand Genomes SV
pub mod g1k_sv {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// gnomAD SV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall variant alleles
        pub alleles: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.alleles as usize
        }
    }

    /// Thousand Genomes SV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub start: u32,
        /// end position, 1-based
        pub end: u32,
        /// The structural vairant type
        sv_type: String,
        /// Number of variant alleles
        num_var_alleles: u32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.sv_type.as_str() {
                "CNV" => SvType::Cnv,
                "DEL" => SvType::Del,
                "DEL_ALU" | "DEL_HERV" | "DEL_LINE1" | "DEL_SVA" => SvType::Del,
                "DUP" => SvType::Dup,
                "INV" => SvType::Inv,
                "ALU" | "INS" | "LINE1" | "SVA" => SvType::Ins,
                _ => return Err(anyhow!("Unknown SV type {}", &self.sv_type)),
            };
            Ok(Some(Record {
                begin: self.start - 1,
                end: self.end,
                sv_type,
                alleles: self.num_var_alleles,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }
}

/// Records for DGV
pub mod dgv {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// gnomAD SV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall carriers
        pub carriers: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// dbVar database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub start: u32,
        /// end position, 1-based
        pub end: u32,
        /// The structural variant type
        sv_type: String,
        /// Number of observed gains.
        observed_gains: u32,
        /// Number of observed losses
        observed_losses: u32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.sv_type.as_ref() {
                "alu deletion"
                | "deletion"
                | "herv deletion"
                | "line1 deletion"
                | "mobile element deletion"
                | "loss"
                | "sva deletion" => SvType::Del,
                "alu insertion"
                | "herv insertion"
                | "insertion"
                | "line1 insertion"
                | "mobile element insertion"
                | "novel sequence insertion"
                | "sva insertion" => SvType::Ins,
                "duplication" | "gain" | "tandem duplication" => SvType::Dup,
                "sequence alteration" | "complex" => return Ok(None),
                "gain+loss" | "CNV" => SvType::Cnv,
                "inversion" => SvType::Inv,
                "OTHER" => return Ok(None),
                _ => return Err(anyhow!("Unknown sv_type {}", &self.sv_type)),
            };
            Ok(Some(Record {
                begin: self.start - 1,
                end: self.end,
                sv_type,
                carriers: self.observed_gains + self.observed_losses,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }
}

/// Records for DGV Gold Standard
pub mod dgv_gs {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// DGV gold standard database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall carriers
        pub carriers: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// DGV gold standard database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// outer start position, 1-based
        pub start_outer: u32,
        /// outer end position, 1-based
        pub end_outer: u32,
        /// The structural variant type
        pub sv_sub_type: String,
        /// Number of carriers.
        pub num_carriers: u32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.sv_sub_type.as_ref() {
                "Gain" => SvType::Dup,
                "Loss" => SvType::Del,
                _ => return Err(anyhow!("Invalid SV type {}", &self.sv_sub_type)),
            };
            Ok(Some(Record {
                begin: self.start_outer - 1,
                end: self.end_outer,
                sv_type,
                carriers: self.num_carriers,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start_outer - 1
        }

        fn start(&self) -> u32 {
            self.start_outer
        }

        fn end(&self) -> u32 {
            self.end_outer
        }
    }
}

/// Records for ExAC CNV
pub mod exac_cnv {
    use crate::sv_query::schema::SvType;

    use super::{BeginEnd, ChromosomeCoordinate, Count, ToInMemory};
    use anyhow::anyhow;
    use serde::Deserialize;

    /// ExAC CNV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: u32,
        /// end position, 0-based
        pub end: u32,

        /// type of the SV
        pub sv_type: SvType,

        /// number of overall carriers
        pub carriers: u32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> u32 {
            self.begin
        }
        fn end(&self) -> u32 {
            self.end
        }
    }

    impl Count for Record {
        fn count(&self) -> usize {
            self.carriers as usize
        }
    }

    /// ExAC CNV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct FileRecord {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// outer start position, 1-based
        pub start: u32,
        /// outer end position, 1-based
        pub end: u32,
        /// The structural vairant type
        pub sv_type: String,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Result<Option<Record>, anyhow::Error> {
            let sv_type = match self.sv_type.as_ref() {
                "duplication" => SvType::Dup,
                "deletion" => SvType::Del,
                _ => return Err(anyhow!("Invalid SV type {}", &self.sv_type)),
            };
            Ok(Some(Record {
                begin: self.start - 1,
                end: self.end,
                sv_type,
                carriers: 1,
            }))
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> u32 {
            self.start - 1
        }

        fn start(&self) -> u32 {
            self.start
        }

        fn end(&self) -> u32 {
            self.end
        }
    }
}
