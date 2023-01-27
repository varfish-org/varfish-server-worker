//! Code for converting gene region TSVs to binary format.

use std::{fs::File, io::Write, time::Instant};

use clap::{arg, command, Parser};
use thousands::Separable;
use tracing::{debug, info};

use crate::{
    common::{build_chrom_map, open_maybe_gz, trace_rss_now},
    world_flatbuffers::var_fish_server_worker::{
        GeneRegionDatabase, GeneRegionDatabaseArgs, GeneRegionRecord,
    },
};

/// Command line arguments for `sv gene-region-to-bin` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Convert gene region TSV to binary file", long_about = None)]
pub struct Args {
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub path_input_tsv: String,
    /// Path to output binary file.
    #[arg(long, required = true)]
    pub path_output_bin: String,
}

/// Module with code supporting the parsing.
mod input {
    use serde::Deserialize;

    /// Record as created by VarFish DB Downloader.
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// Genome release
        pub release: String,
        /// Chromosome name
        pub chromosome: String,
        /// 1-based start position
        pub start: u32,
        /// 1-based end position
        pub end: u32,
        /// ENSEMBL or Entrez gene ID
        pub gene_id: String,
    }
}

/// Helper to convert ENSEMBL and RefSeq gene ID to u32.
pub fn numeric_gene_id(raw_id: &str) -> Result<u32, anyhow::Error> {
    let clean_id = if raw_id.starts_with("ENSG") {
        // Strip "ENSG" prefix and as many zeroes as follow
        raw_id
            .chars()
            .skip("ENSG".len())
            .skip_while(|c| *c == '0')
            .collect()
    } else {
        raw_id.to_owned()
    };

    clean_id
        .parse::<u32>()
        .map_err(|e| anyhow::anyhow!("could not parse gene id {:?}: {}", &clean_id, &e))
}

/// Perform conversion to flatbuffers `.bin` file.
pub fn convert_to_bin(args: &Args) -> Result<(), anyhow::Error> {
    let chrom_map = build_chrom_map();

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_maybe_gz(&args.path_input_tsv)?);
    let before_parsing = Instant::now();

    let mut output_records = Vec::new();
    for record in reader.deserialize() {
        let record: input::Record = record?;
        output_records.push(GeneRegionRecord::new(
            *chrom_map.get(&record.chromosome).expect("unknown chrom") as u8,
            record.start.saturating_sub(1),
            record.end,
            numeric_gene_id(&record.gene_id)?,
        ));
    }

    debug!(
        "total time spent reading {:?} records: {:?}",
        output_records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut builder = flatbuffers::FlatBufferBuilder::new();
    let records = builder.create_vector(output_records.as_slice());
    let gene_region_db = GeneRegionDatabase::create(
        &mut builder,
        &GeneRegionDatabaseArgs {
            records: Some(records),
        },
    );
    builder.finish_minimal(gene_region_db);
    let mut output_file = File::create(&args.path_output_bin)?;
    output_file.write_all(builder.finished_data())?;
    output_file.flush()?;
    debug!(
        "total time spent writing {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}

/// Main entry point for the `sv gene-region-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting sv gene-region-to-bin");
    info!("common_args = {:?}", &common_args);
    info!("args = {:?}", &args);

    convert_to_bin(args)?;

    Ok(())
}
