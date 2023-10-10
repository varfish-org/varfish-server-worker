//! Code for converting gene region from text-based to binary format.

use std::{fs::File, io::Write, path::Path, time::Instant};

use mehari::common::open_read_maybe_gz;
use prost::Message;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, numeric_gene_id, trace_rss_now},
    strucvars::pbs::{GeneRegionDatabase, GeneRegionRecord},
};

/// Module with code supporting the parsing.
mod input {
    use serde::Deserialize;

    /// Record as created by VarFish DB Downloader.
    #[derive(Debug, Deserialize)]
    pub struct Record {
        /// Chromosome name
        pub chromosome: String,
        /// 0-based begin position
        pub begin: i32,
        /// 1-based end position
        pub end: i32,
        /// ENSEMBL or Entrez gene ID
        pub gene_id: String,
    }
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(path_input_tsv: P, path_output_bin: Q) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    tracing::debug!(
        "Converting gene regions from BED {:?} to binary {:?}",
        path_input_tsv.as_ref(),
        path_output_bin.as_ref()
    );
    let chrom_map = build_chrom_map();

    // Setup CSV reader for BED file - header is written as comment and must be
    // ignored.
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .from_reader(open_read_maybe_gz(path_input_tsv.as_ref())?);
    let before_parsing = Instant::now();

    let mut records = Vec::new();
    for record in reader.deserialize() {
        let record: input::Record = record?;
        records.push(GeneRegionRecord {
            chrom_no: *chrom_map
                .get(&record.chromosome)
                .unwrap_or_else(|| panic!("unknown chrom {:?}", &record.chromosome))
                as i32,
            start: record.begin + 1,
            stop: record.end,
            gene_id: numeric_gene_id(&record.gene_id)?,
        });
    }
    let gene_region_db = GeneRegionDatabase { records };

    tracing::debug!(
        "total time spent reading {:?} records: {:?}",
        gene_region_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output_bin)?;
    output_file.write_all(&gene_region_db.encode_to_vec())?;
    output_file.flush()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        gene_region_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}
