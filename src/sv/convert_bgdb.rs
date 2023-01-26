//! Code for converting from background DB TSV file to binary background DB format.

use std::fs::File;
use std::io::Write;
use std::time::Instant;

use clap::{command, Parser};
use thousands::Separable;
use tracing::{debug, info};

use crate::common::{build_chrom_map, open_maybe_gz, trace_rss_now};
use crate::sv_query::schema::SvType;
use crate::world_flatbuffers::var_fish_server_worker::{
    BgDbRecord, BgDbRecordArgs, SvType as FlatSvType,
};

use super::build_inhouse_db::output::Record as InputRecord;

/// Command line arguments for `sv build-inhouse-db` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Convert background TSV to binary file", long_about = None)]
pub struct Args {
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub input_tsv: String,
    /// Path to output binary file.
    #[arg(long, required = true)]
    pub output_bin: String,
}

/// Perform conversion to flatbuffers `.bin` file.
pub fn convert_to_bin(args: &Args) -> Result<(), anyhow::Error> {
    let chrom_map = build_chrom_map();

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.input_tsv)?);
    let before_parsing = Instant::now();
    let mut count_records = 0;

    let mut builder = flatbuffers::FlatBufferBuilder::new();
    let mut output_records = Vec::new();
    for result in rdr.deserialize() {
        let record: InputRecord = result?;
        output_records.push(BgDbRecord::create(
            &mut builder,
            &BgDbRecordArgs {
                chrom_no: *chrom_map.get(&record.chromosome).expect("unknown chrom") as u16,
                chrom_no2: *chrom_map.get(&record.chromosome2).expect("unknown chrom2") as u16,
                start: record.start as u32,
                end: record.end as u32,
                sv_type: match record.sv_type {
                    SvType::Del => FlatSvType::Del,
                    SvType::Dup => FlatSvType::Dup,
                    SvType::Inv => FlatSvType::Inv,
                    SvType::Ins => FlatSvType::Ins,
                    SvType::Bnd => FlatSvType::Bnd,
                    SvType::Cnv => FlatSvType::Cnv,
                },
                count: record.carriers,
            },
        ));
        count_records += 1;
    }

    debug!(
        "total time spent reading {} records: {:?}",
        count_records.separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let output_vector = builder.create_vector(output_records.as_slice());
    builder.finish_minimal(output_vector);
    let mut output_file = File::create(&args.output_bin)?;
    output_file.write_all(builder.finished_data())?;
    output_file.flush()?;
    debug!(
        "total time spent writing {} records: {:?}",
        count_records.separate_with_commas(),
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
    
    use crate::common::Args as CommonArgs;
    use super::{Args, run};

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
        };

        run(&common_args, &args)?;

        assert!(file_diff::diff("tests/sv/convert_bgdb/output.bin", &output_bin));

        Ok(())
    }
}
