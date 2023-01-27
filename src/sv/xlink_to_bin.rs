//! Conversion from gene ids to binary.

use std::{fs::File, io::Write, time::Instant};

use clap::Parser;
use thousands::Separable;
use tracing::{debug, info};

use crate::{
    common::{open_maybe_gz, trace_rss_now},
    sv::gene_region_to_bin::numeric_gene_id,
    world_flatbuffers::var_fish_server_worker::{
        XlinkDatabase, XlinkDatabaseArgs, XlinkRecord as FlatXlinkRecord,
    },
};

/// Command line arguments for `sv gene-ids-to-bin` sub command.
#[derive(Parser, Debug)]
#[command(about = "Convert gene id mapping to binary", long_about = None)]
pub struct Args {
    /// Path to interlink TSV file.
    #[arg(long, required = true)]
    pub path_input_tsv: String,
    /// Path to output binary file.
    #[arg(long, required = true)]
    pub path_output_bin: String,
}

/// Module with code for parsing the TSVs.
pub mod input {
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct XlinkRecord {
        pub entrez_id: u32,
        pub gene_symbol: String,
        pub ensembl_gene_id: String,
    }
}

/// Perform conversion to flatbuffers `.bin` file.
pub fn convert_to_bin(args: &Args) -> Result<(), anyhow::Error> {
    let mut output_records = Vec::new();
    let mut output_strings = Vec::new();

    let before_parsing = Instant::now();

    let mut builder = flatbuffers::FlatBufferBuilder::new();

    debug!("parsing xlink TSV file from {}", &args.path_input_tsv);
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.path_input_tsv)?);
    for record in reader.deserialize() {
        let record: input::XlinkRecord = record?;
        output_records.push(FlatXlinkRecord::new(
            record.entrez_id,
            numeric_gene_id(&record.ensembl_gene_id)?,
        ));
        let flat_str = builder.create_shared_string(&record.gene_symbol);
        output_strings.push(flat_str);
    }

    debug!(
        "total time spent reading {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let records = builder.create_vector(&output_records);
    let strings = builder.create_vector(&output_strings);
    let xlink_db = XlinkDatabase::create(
        &mut builder,
        &XlinkDatabaseArgs {
            records: Some(records),
            symbols: Some(strings),
        },
    );
    builder.finish_minimal(xlink_db);
    let mut output_file = File::create(&args.path_output_bin)?;
    output_file.write_all(builder.finished_data())?;
    output_file.flush()?;
    debug!(
        "total time spent writing {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_writing.elapsed()
    );
    trace_rss_now();

    Ok(())
}

/// Main entry point for the `sv gene-ids-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting sv gene-ids-to-bin");
    info!("common_args = {:?}", &common_args);
    info!("args = {:?}", &args);

    convert_to_bin(args)?;

    Ok(())
}
