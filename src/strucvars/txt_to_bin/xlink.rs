//! Code for converting xlink from text-based to binary format.

use std::{fs::File, io::Write, path::Path, time::Instant};

use mehari::common::io::std::open_read_maybe_gz;
use prost::Message;
use thousands::Separable;

use crate::{
    common::{numeric_gene_id, trace_rss_now},
    strucvars::pbs::{XlinkDatabase, XlinkRecord},
};

/// Module with code for parsing the TSVs.
pub mod input {
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub struct Record {
        pub hgnc_id: Option<String>,
        pub gene_symbol: Option<String>,
        pub ensembl_gene_id: Option<String>,
        pub entrez_id: Option<u32>,
    }
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(path_input_tsv: P, path_output: Q) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut records = Vec::new();

    let before_parsing = Instant::now();

    tracing::debug!("parsing xlink TSV file from {:?}", path_input_tsv.as_ref());
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path_input_tsv)?);
    for record in reader.deserialize() {
        let record: input::Record = record?;
        if let (Some(entrez_id), Some(ensembl_gene_id), Some(gene_symbol), Some(hgnc_id)) = (
            record.entrez_id,
            record.ensembl_gene_id,
            record.gene_symbol,
            record.hgnc_id,
        ) {
            records.push(XlinkRecord {
                entrez_id,
                hgnc_id,
                ensembl_id: numeric_gene_id(&ensembl_gene_id)?,
                symbol: gene_symbol,
            });
        }
    }
    let xlink_db = XlinkDatabase { records };

    tracing::debug!(
        "total time spent reading {} records: {:?}",
        xlink_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output)?;
    output_file.write_all(&xlink_db.encode_to_vec())?;
    output_file.flush()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        xlink_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );
    trace_rss_now();

    Ok(())
}
