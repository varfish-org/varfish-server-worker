//! Code supporting the I/O of public database records and a common
//! `InputRecord` for common representation.

use serde::Deserialize;
use tracing::error;

use crate::strucvars::aggregate::output::Record as InhouseDbRecord;
use crate::strucvars::query::schema::SvType;

/// dbVar database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct DbVarRecord {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// number of overall carriers
    pub num_carriers: u32,
    /// type of the SV
    pub sv_type: String,
}

/// DGV database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct DgvRecord {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// The structural variant type
    sv_type: String,
    /// Number of observed gains.
    observed_gains: u32,
    /// Number of observed losses
    observed_losses: u32,
}

/// DGV gold standard database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct DgvGsRecord {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin_outer: i32,
    /// outer end position, 0-based
    pub end_outer: i32,
    /// The structural variant type
    pub sv_sub_type: String,
    /// Number of carriers.
    pub num_carriers: u32,
}

/// ExAC CNV database record as read from TSV file for deserialization
/// from TSV.
#[derive(Deserialize, Debug)]
pub struct ExacRecord {
    /// chromosome name
    pub chromosome: String,
    /// outer start position, 0-based
    pub begin: i32,
    /// outer end position, 0-based
    pub end: i32,
    /// The structural vairant type
    pub sv_type: String,
}

/// Thousand Genomes SV database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct G1kRecord {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// The structural vairant type
    pub sv_type: String,
    /// Number of hom. alt. alleles.
    pub n_homalt: u32,
    /// Number of het. alt. alleles.
    pub n_het: u32,
}

/// gnomAD SV v2 database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct GnomadSv2Record {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// The structural vairant type
    pub svtype: String,
    /// Number of homozygous alternative carriers
    pub n_homalt: u32,
    /// Number of heterozygous carriers
    pub n_het: u32,
}

/// gnomAD SV v4 database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct GnomadSv4Record {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// The structural vairant type
    pub svtype: String,
    /// Number of male homozygous reference allele carriers.
    pub male_n_homref: u32,
    /// Number of male heterozygous alternate allele carriers.
    pub male_n_het: u32,
    /// Number of male homozygous alternate allele carriers.
    pub male_n_homalt: u32,
    /// Number of male hemizygous alternate allele carriers.
    pub male_n_hemiref: u32,
    /// Number of male hemizygous reference allele carriers.
    pub male_n_hemialt: u32,
    /// Number of female homozygous reference allele carriers.
    pub female_n_homref: u32,
    /// Number of female heterozygous alternate allele carriers.
    pub female_n_het: u32,
    /// Number of female homozygous alternate allele carriers.
    pub female_n_homalt: u32,
    /// Number of samples at this site (CNV only).
    pub cnv_n_total: u32,
    /// Number of samples with a CNV at this site (CNV only).
    pub cnv_n_var: u32,
}

/// gnomAD CNV v$ database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct GnomadCnv4Record {
    /// chromosome name
    pub chromosome: String,
    /// begin position, 0-based
    pub begin: i32,
    /// end position, 0-based
    pub end: i32,
    /// The structural vairant type
    pub svtype: String,
    /// Number of samples at this site (passing QC).
    pub n_total: u32,
    /// Number of samples with a CNV at this site (passing QC).
    pub n_var: u32,
}

/// Common type to convert input data to.
pub struct InputRecord {
    /// Chromosome of start position.
    pub chromosome: String,
    /// Chromosome of end position.
    pub chromosome2: String,
    /// SV type
    pub sv_type: SvType,
    /// 0-based begin position
    pub begin: i32,
    /// 0-based end position
    pub end: i32,
    /// Number of carriers (or alleles), depending on database.
    pub count: u32,
}

impl TryInto<Option<InputRecord>> for InhouseDbRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        Ok(Some(InputRecord {
            chromosome: self.chromosome,
            chromosome2: self.chromosome2,
            sv_type: self.sv_type,
            begin: self.begin,
            end: self.end,
            count: self.carriers,
        }))
    }
}

impl TryInto<Option<InputRecord>> for DbVarRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
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
            _ => {
                error!("sv_type = {}", &self.sv_type);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type,
            count: 1,
        }))
    }
}

impl TryInto<Option<InputRecord>> for DgvRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
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
            "sequence alteration" | "complex" => return Ok(None), // skip
            "gain+loss" | "CNV" => SvType::Cnv,
            "inversion" => SvType::Inv,
            "OTHER" => return Ok(None), // skip
            _ => {
                error!("sv_type = {}", &self.sv_type);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type,
            count: self.observed_gains + self.observed_losses,
        }))
    }
}

impl TryInto<Option<InputRecord>> for DgvGsRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        let sv_type = match self.sv_sub_type.as_ref() {
            "Gain" => SvType::Dup,
            "Loss" => SvType::Del,
            _ => {
                error!("sv_type = {}", &self.sv_sub_type);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin_outer,
            end: self.end_outer,
            sv_type,
            count: self.num_carriers,
        }))
    }
}

impl TryInto<Option<InputRecord>> for ExacRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        let sv_type = match self.sv_type.as_ref() {
            "DUP" => SvType::Dup,
            "DEL" => SvType::Del,
            _ => {
                error!("sv_type = {}", &self.sv_type);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type,
            count: 1,
        }))
    }
}

impl TryInto<Option<InputRecord>> for GnomadSv2Record {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        let sv_type = match self.svtype.as_str() {
            "CPX" => return Ok(None), // no correspondence
            "CTX" | "BND" => SvType::Bnd,
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INS" => SvType::Ins,
            "INV" => SvType::Inv,
            "MCNV" => SvType::Cnv,
            _ => {
                error!("sv_type = {}", &self.svtype);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin - 1,
            end: self.end,
            sv_type,
            count: self.n_homalt + self.n_het,
        }))
    }
}

impl TryInto<Option<InputRecord>> for GnomadCnv4Record {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type: match self.svtype.as_str() {
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                _ => {
                    error!("sv_type = {}", &self.svtype);
                    return Err("unknown SV type");
                }
            },
            count: self.n_var,
        }))
    }
}

impl TryInto<Option<InputRecord>> for GnomadSv4Record {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type: match self.svtype.as_str() {
                "BND" => SvType::Bnd,
                "CNV" => SvType::Cnv,
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                "INS" => SvType::Ins,
                "INV" => SvType::Inv,
                _ => {
                    error!("sv_type = {}", &self.svtype);
                    return Err("unknown SV type");
                }
            },
            count: self.male_n_het
                + self.male_n_homalt
                + self.male_n_hemialt
                + self.female_n_het
                + self.female_n_homalt
                + self.cnv_n_var,
        }))
    }
}

impl TryInto<Option<InputRecord>> for G1kRecord {
    type Error = &'static str;

    fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
        let sv_type = match self.sv_type.as_str() {
            "CN0" | "CNV" => SvType::Cnv,
            "DEL" => SvType::Del,
            "DEL_ALU" | "DEL_HERV" | "DEL_LINE1" | "DEL_SVA" => SvType::Del,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "INS" | "INS:ME:ALU" | "INS:ME:LINE1" | "INS:ME:SVA" => SvType::Ins,
            _ => {
                error!("sv_type = {}", &self.sv_type);
                return Err("unknown SV type");
            }
        };
        Ok(Some(InputRecord {
            chromosome: self.chromosome.clone(),
            chromosome2: self.chromosome,
            begin: self.begin,
            end: self.end,
            sv_type,
            count: self.n_homalt + self.n_het,
        }))
    }
}
