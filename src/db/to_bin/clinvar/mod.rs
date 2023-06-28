//! Code for converting ClinVar database to binary.

use std::{fs::File, io::Write, path::Path, time::Instant};

use prost::Message;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, open_read_maybe_gz, trace_rss_now},
    sv::query::clinvar::pbs::{SvDatabase, SvRecord},
};

mod input;

/// Helper to convert VCV IDs to numbers.
fn numeric_vcv_id(raw_id: &str) -> Result<u32, anyhow::Error> {
    let clean_id: String = raw_id
        .chars()
        .skip("VCV".len())
        .skip_while(|c| *c == '0')
        .collect();
    clean_id
        .parse::<u32>()
        .map_err(|e| anyhow::anyhow!("could not parse VCV id {:?}: {}", &clean_id, &e))
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(path_input_tsv: P, path_output_bin: Q) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let chrom_map = build_chrom_map();

    let mut reader = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(open_read_maybe_gz(path_input_tsv)?);
    let before_parsing = Instant::now();

    let mut records = Vec::new();
    for record in reader.deserialize() {
        match record {
            Err(e) => {
                tracing::warn!("error parsing ClinVar record: {}", e);
            }
            Ok(record) => {
                let record: input::Record = record;
                let variation_type: crate::sv::query::clinvar::pbs::VariationType =
                    record.variation_type.into();
                let pathogenicity: crate::sv::query::clinvar::pbs::Pathogenicity =
                    record.pathogenicity.into();
                if let Some(chrom_no) = chrom_map.get(&record.chromosome) {
                    records.push(SvRecord {
                        chrom_no: *chrom_no as i32,
                        start: record.begin + 1,
                        stop: record.end,
                        variation_type: variation_type as i32,
                        pathogenicity: pathogenicity as i32,
                        vcv: numeric_vcv_id(&record.vcv)?,
                    });
                } else {
                    tracing::warn!("unknown chromosome {}", &record.chromosome);
                }
            }
        }
    }
    let clinvar_db = SvDatabase { records };

    tracing::debug!(
        "total time spent reading {} records: {:?}",
        clinvar_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output_bin)?;
    output_file.write_all(&clinvar_db.encode_to_vec())?;
    output_file.flush()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        clinvar_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}
