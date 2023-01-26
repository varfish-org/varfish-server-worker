//! Code for converting from background DB TSV file to binary background DB format.

use std::fs::File;
use std::io::Write;
use std::time::Instant;

use anyhow::anyhow;
use clap::{command, Parser, ValueEnum};
use thousands::Separable;
use tracing::{debug, info};

use crate::common::{build_chrom_map, open_maybe_gz, trace_rss_now};
use crate::sv::inhouse_db_build::output::Record as InhouseDbRecord;
use crate::sv::bg_db_to_bin::records::{
    DbVarRecord, DgvGsRecord, DgvRecord, ExacRecord, G1kRecord, GnomadRecord,
};
use crate::sv::query::schema::SvType;
use crate::world_flatbuffers::var_fish_server_worker::{
    BackgroundDatabase, BackgroundDatabaseArgs, BgDbRecord, SvType as FlatSvType,
};

use self::records::InputRecord;

/// Code supporting the I/O of public database records and a common `InputRecord` for
/// common representation.
mod records {
    use serde::Deserialize;

    use crate::sv::inhouse_db_build::output::Record as InhouseDbRecord;
    use crate::sv::query::schema::SvType;

    /// dbVar database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct DbVarRecord {
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

    /// DGV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct DgvRecord {
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

    /// DGV gold standard database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct DgvGsRecord {
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

    /// ExAC CNV database record as read from TSV file for deserialization from TSV.
    #[derive(Deserialize, Debug)]
    pub struct ExacRecord {
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

    /// Thousand Genomes SV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct G1kRecord {
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

    /// gnomAD SV database record as read from TSV file.
    #[derive(Debug, Deserialize)]
    pub struct GnomadRecord {
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

    /// Common type to convert input data to.
    pub struct InputRecord {
        /// Chromosome of start position.
        pub chromosome: String,
        /// Chromosome of end position.
        pub chromosome2: String,
        /// SV type
        pub sv_type: SvType,
        /// 0-based begin position
        pub begin: u32,
        /// 0-based end position
        pub end: u32,
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
                begin: self.start.saturating_sub(1),
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
                _ => return Err("unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start.saturating_sub(1),
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
                _ => return Err("unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start.saturating_sub(1),
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
                _ => return Err("unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start_outer.saturating_sub(1),
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
                "duplication" => SvType::Dup,
                "deletion" => SvType::Del,
                _ => return Err("unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start.saturating_sub(1),
                end: self.end,
                sv_type,
                count: 1,
            }))
        }
    }

    impl TryInto<Option<InputRecord>> for GnomadRecord {
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
                _ => return Err("unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start.saturating_sub(1),
                end: self.end,
                sv_type,
                count: self.n_homalt + self.n_het,
            }))
        }
    }

    impl TryInto<Option<InputRecord>> for G1kRecord {
        type Error = &'static str;

        fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
            let sv_type = match self.sv_type.as_str() {
                "CNV" => SvType::Cnv,
                "DEL" => SvType::Del,
                "DEL_ALU" | "DEL_HERV" | "DEL_LINE1" | "DEL_SVA" => SvType::Del,
                "DUP" => SvType::Dup,
                "INV" => SvType::Inv,
                "ALU" | "INS" | "LINE1" | "SVA" => SvType::Ins,
                _ => return Err("Unknown SV type"),
            };
            Ok(Some(InputRecord {
                chromosome: self.chromosome.clone(),
                chromosome2: self.chromosome,
                begin: self.start.saturating_sub(1),
                end: self.end,
                sv_type,
                count: self.num_var_alleles,
            }))
        }
    }
}

/// Known input file types.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, ValueEnum)]
pub enum InputFileType {
    Dbvar,
    Dgv,
    DgvGs,
    Exac,
    G1k,
    Gnomad,
    InhouseDb,
}

/// Command line arguments for `sv build-inhouse-db` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Convert background TSV to binary file", long_about = None)]
pub struct Args {
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub input_tsv: String,
    /// Type of the input record to expect.
    #[arg(long, required = true)]
    pub input_type: InputFileType,
    /// Path to output binary file.
    #[arg(long, required = true)]
    pub output_bin: String,
}

/// Deserialize from CSV reader to an `Option<records::InputRecord>`
fn deserialize_loop<Rec>(
    reader: &mut csv::Reader<Box<dyn std::io::Read>>,
) -> Result<Vec<BgDbRecord>, anyhow::Error>
where
    Rec: TryInto<Option<InputRecord>> + for<'de> serde::Deserialize<'de>,
    <Rec as TryInto<std::option::Option<InputRecord>>>::Error: core::fmt::Debug,
    <Rec as TryInto<std::option::Option<InputRecord>>>::Error: Send,
    <Rec as TryInto<std::option::Option<InputRecord>>>::Error: std::marker::Sync,
{
    let chrom_map = build_chrom_map();
    let mut result = Vec::new();

    for record in reader.deserialize() {
        let record: Rec = record?;
        let maybe_record: Option<InputRecord> = record
            .try_into()
            .map_err(|err| anyhow!("problem with parsing: {:?}", &err))?;
        if let Some(record) = maybe_record {
            result.push(BgDbRecord::new(
                *chrom_map.get(&record.chromosome).expect("unknown chrom") as u8,
                *chrom_map.get(&record.chromosome2).expect("unknown chrom2") as u8,
                match record.sv_type {
                    SvType::Del => FlatSvType::Del,
                    SvType::Dup => FlatSvType::Dup,
                    SvType::Inv => FlatSvType::Inv,
                    SvType::Ins => FlatSvType::Ins,
                    SvType::Bnd => FlatSvType::Bnd,
                    SvType::Cnv => FlatSvType::Cnv,
                },
                record.begin,
                record.end,
                record.count,
            ));
        }
    }

    Ok(result)
}

/// Perform conversion to flatbuffers `.bin` file.
pub fn convert_to_bin(args: &Args) -> Result<(), anyhow::Error> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.input_tsv)?);
    let before_parsing = Instant::now();

    let output_records = match args.input_type {
        InputFileType::Dbvar => deserialize_loop::<DbVarRecord>(&mut reader)?,
        InputFileType::Dgv => deserialize_loop::<DgvRecord>(&mut reader)?,
        InputFileType::DgvGs => deserialize_loop::<DgvGsRecord>(&mut reader)?,
        InputFileType::Exac => deserialize_loop::<ExacRecord>(&mut reader)?,
        InputFileType::G1k => deserialize_loop::<G1kRecord>(&mut reader)?,
        InputFileType::Gnomad => deserialize_loop::<GnomadRecord>(&mut reader)?,
        InputFileType::InhouseDb => deserialize_loop::<InhouseDbRecord>(&mut reader)?,
    };

    debug!(
        "total time spent reading {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut builder = flatbuffers::FlatBufferBuilder::new();
    let records = builder.create_vector(output_records.as_slice());
    let bg_db = BackgroundDatabase::create(
        &mut builder,
        &BackgroundDatabaseArgs {
            records: Some(records),
        },
    );
    builder.finish_minimal(bg_db);
    let mut output_file = File::create(&args.output_bin)?;
    output_file.write_all(builder.finished_data())?;
    output_file.flush()?;
    debug!(
        "total time spent writing {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}

/// Main entry point for the `sv convert-bgdb` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting sv convert-bgdb");
    info!("common_args = {:?}", &common_args);
    info!("args = {:?}", &args);

    convert_to_bin(args)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use clap_verbosity_flag::Verbosity;
    use temp_testdir::TempDir;

    use super::{run, Args};
    use crate::{common::Args as CommonArgs, sv::bg_db_to_bin::InputFileType};

    #[test]
    fn test_run_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 0),
        };
        let output_bin = tmp_dir.join("output.bin").to_str().unwrap().to_owned();
        let args = Args {
            input_tsv: "tests/sv/convert_bgdb/input.tsv".to_owned(),
            output_bin: output_bin.clone(),
            input_type: InputFileType::InhouseDb,
        };

        run(&common_args, &args)?;

        assert!(file_diff::diff(
            "tests/sv/convert_bgdb/output.bin",
            &output_bin
        ));

        Ok(())
    }
}
