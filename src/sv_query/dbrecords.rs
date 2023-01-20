/// This module provides the code for accessing database records.

/// Provide a chromosome-wise coordinate.
pub trait ChromosomeCoordinate {
    fn chromosome(&self) -> &String;
    /// 0-based begin position
    fn begin(&self) -> i32;
    /// 1-based begin position
    fn start(&self) -> i32;
    /// 0/1-based end position
    fn end(&self) -> i32;
}

pub trait BeginEnd {
    fn begin(&self) -> i32;
    fn end(&self) -> i32;
}

pub trait ToInMemory<InMemory> {
    fn to_in_memory(&self) -> InMemory;
}

/// Store background database counts for a structural variant.
#[derive(Debug, PartialEq)]
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

/// Records for in-house SV background database.
pub mod bg_sv {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// Background SV database record to be kept in memory.
    #[derive(Debug)]
    pub struct Record {
        /// The 0-based begin position.
        pub begin: i32,
        /// The 0-based end position.
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// Total number of carriers.
        pub carriers: i32,
        /// Number of het. carriers.
        pub carriers_het: i32,
        /// Number of hom. carriers.
        pub carriers_hom: i32,
        /// Number of hemi. carriers.
        pub carriers_hemi: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start: i32,
        /// chromosome2 name
        pub chromosome2: String,
        /// end position, 1-based
        pub end: i32,
        /// paired-end orientation
        pub pe_orientation: String,
        /// type of the SV
        pub sv_type: String,
        /// number of overall carriers
        pub carriers: i32,
        /// number of het. carriers
        pub carriers_het: i32,
        /// number of hom. carriers
        pub carriers_hom: i32,
        /// number of hemi. carriers
        pub carriers_hemi: i32,
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start - 1
        }

        fn start(&self) -> i32 {
            self.start
        }

        fn end(&self) -> i32 {
            self.end
        }
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start - 1,
                end: self.end,
                sv_type: self.sv_type.clone(),
                carriers: self.carriers,
                carriers_het: self.carriers_het,
                carriers_hom: self.carriers_hom,
                carriers_hemi: self.carriers_hemi,
            }
        }
    }
}

/// Records for the dbVar
pub mod dbvar {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// dbVar database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: i32,
        /// end position, 0-based
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// number of overall carriers
        pub carriers: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start: i32,
        /// end position, 1-based
        pub end: i32,
        /// number of overall carriers
        pub num_carriers: i32,
        /// type of the SV
        pub sv_type: String,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start - 1,
                end: self.end,
                sv_type: self.sv_type.clone(),
                carriers: self.num_carriers,
            }
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start - 1
        }

        fn start(&self) -> i32 {
            self.start
        }

        fn end(&self) -> i32 {
            self.end
        }
    }
}

/// Records for gnomAD SV
pub mod gnomad_sv {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// gnomAD SV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: i32,
        /// end position, 0-based
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// number of overall carriers
        pub carriers: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start: i32,
        /// end position, 1-based
        pub end: i32,
        /// The structural vairant type
        svtype: String,
        /// Number of homozygous alternative carriers
        n_homalt: i32,
        /// Number of heterozygous carriers
        n_het: i32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start - 1,
                end: self.end,
                sv_type: self.svtype.clone(),
                carriers: self.n_homalt + self.n_het,
            }
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start - 1
        }

        fn start(&self) -> i32 {
            self.start
        }

        fn end(&self) -> i32 {
            self.end
        }
    }
}

/// Records for gnomAD SV
pub mod dgv {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// gnomAD SV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: i32,
        /// end position, 0-based
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// number of overall carriers
        pub carriers: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start: i32,
        /// end position, 1-based
        pub end: i32,
        /// The structural vairant type
        sv_type: String,
        /// Number of observed gains.
        observed_gains: i32,
        /// Number of observed losses
        observed_losses: i32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start - 1,
                end: self.end,
                sv_type: self.sv_type.clone(),
                carriers: self.observed_gains + self.observed_losses,
            }
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start - 1
        }

        fn start(&self) -> i32 {
            self.start
        }

        fn end(&self) -> i32 {
            self.end
        }
    }
}

/// Records for DGV Gold Standard
pub mod dgv_gs {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// DGV gold standard database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: i32,
        /// end position, 0-based
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// number of overall carriers
        pub carriers: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start_outer: i32,
        /// outer end position, 1-based
        pub end_outer: i32,
        /// The structural vairant type
        sv_type: String,
        /// Number of carriers.
        num_carriers: i32,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start_outer - 1,
                end: self.end_outer,
                sv_type: self.sv_type.clone(),
                carriers: self.num_carriers,
            }
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start_outer - 1
        }

        fn start(&self) -> i32 {
            self.start_outer
        }

        fn end(&self) -> i32 {
            self.end_outer
        }
    }
}

/// Records for ExAC CNV
pub mod exac_cnv {
    use super::{BeginEnd, ChromosomeCoordinate, ToInMemory};
    use serde::Deserialize;

    /// ExAC CNV database record to be kept in memor
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// start position, 0-based
        pub begin: i32,
        /// end position, 0-based
        pub end: i32,

        /// type of the SV
        pub sv_type: String,

        /// number of overall carriers
        pub carriers: i32,
    }

    impl BeginEnd for Record {
        fn begin(&self) -> i32 {
            self.begin
        }
        fn end(&self) -> i32 {
            self.end
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
        pub start: i32,
        /// outer end position, 1-based
        pub end: i32,
        /// The structural vairant type
        sv_type: String,
    }

    impl ToInMemory<Record> for FileRecord {
        fn to_in_memory(&self) -> Record {
            Record {
                begin: self.start - 1,
                end: self.end,
                sv_type: self.sv_type.clone(),
                carriers: 1,
            }
        }
    }

    impl ChromosomeCoordinate for FileRecord {
        fn chromosome(&self) -> &String {
            &self.chromosome
        }

        fn begin(&self) -> i32 {
            self.start - 1
        }

        fn start(&self) -> i32 {
            self.start
        }

        fn end(&self) -> i32 {
            self.end
        }
    }
}
