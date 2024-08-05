//! Code for converting masked regions from text-based to binary format.

use std::{fs::File, io::Write, path::Path, time::Instant};

use mehari::common::io::std::open_read_maybe_gz;
use prost::Message;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, trace_rss_now},
    pbs::varfish::v1::strucvars::bgdb::{MaskedDatabase, MaskedDbRecord},
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
        /// Masked region label
        #[allow(dead_code)]
        pub label: String,
    }
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(path_input_tsv: P, path_output: Q) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    tracing::debug!(
        "Converting masked region from BED {:?} to binary {:?}",
        path_input_tsv.as_ref(),
        path_output.as_ref()
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
        records.push(MaskedDbRecord {
            chrom_no: *chrom_map
                .get(&record.chromosome)
                .unwrap_or_else(|| panic!("unknown chrom {:?}", &record.chromosome))
                as i32,
            start: record.begin + 1,
            stop: record.end,
        });
    }
    let masked_region_db = MaskedDatabase { records };

    tracing::debug!(
        "total time spent reading {:?} records: {:?}",
        masked_region_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output)?;
    output_file.write_all(&masked_region_db.encode_to_vec())?;
    output_file.sync_all()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        masked_region_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}
