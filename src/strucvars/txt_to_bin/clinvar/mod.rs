//! Code for converting ClinVar database to binary.

use std::{fs::File, io::BufRead, io::Write, path::Path, time::Instant};

use mehari::common::open_read_maybe_gz;
use prost::Message;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, trace_rss_now},
    strucvars::query::clinvar::pbs::{Pathogenicity, SvDatabase, SvRecord},
};

pub mod input;

/// Helper to convert RCV IDs to numbers.
fn numeric_rcv_id(raw_id: &str) -> Result<u32, anyhow::Error> {
    let clean_id: String = raw_id
        .chars()
        .skip("RCV".len())
        .skip_while(|c| *c == '0')
        .collect();
    clean_id
        .parse::<u32>()
        .map_err(|e| anyhow::anyhow!("could not parse RCV id {:?}: {}", &clean_id, &e))
}

/// Read JSONL file and convert to protobuf records.
fn convert_jsonl_to_protobuf(
    reader: Box<dyn BufRead>,
    assembly: input::Assembly,
) -> Result<Vec<SvRecord>, anyhow::Error> {
    let chrom_map = build_chrom_map();

    let mut records = Vec::new();
    for line in reader.lines() {
        // get next line into a String
        let line = if let Ok(line) = line {
            line
        } else {
            anyhow::bail!("error reading line from input file")
        };

        // deserialize JSONL record from line
        let record = serde_json::from_str(&line);
        let record = match record {
            Err(e) => {
                tracing::warn!("error deserializing JSONL record: \"{}\" in {}", e, &line);
                continue;
            }
            Ok(record) => {
                let record: input::ClinVarSet = record;
                record
            }
        };

        let rcv = numeric_rcv_id(&record.reference_clinvar_assertion.clinvar_accession.acc)?;

        // convert from JSONL to protocolbuffers: pathogenicity
        let pathogenicity: Result<Pathogenicity, anyhow::Error> = record
            .reference_clinvar_assertion
            .clinical_significance
            .description
            .try_into();
        let pathogenicity = if let Ok(pathogenicity) = pathogenicity {
            pathogenicity as i32
        } else {
            continue;
        };

        // there can be multiple measures, we consider them all
        for measure in &record.reference_clinvar_assertion.measures.measures {
            // convert from JSONL to protocolbuffers: variation type
            let variation_type: Result<
                crate::strucvars::query::clinvar::pbs::VariationType,
                anyhow::Error,
            > = measure.r#type.try_into();
            let variation_type = if let Ok(variation_type) = variation_type {
                variation_type as i32
            } else {
                continue;
            };

            // we process sequence locations on the selected assembly
            for sl in &measure.sequence_locations {
                if sl.assembly != assembly {
                    continue;
                }
                let chrom_no = if let Some(chrom_no) = chrom_map.get(&sl.chr) {
                    *chrom_no as i32
                } else {
                    tracing::warn!("unknown chromosome {}", &sl.chr);
                    continue;
                };

                if let (Some(start), Some(stop)) = (sl.start, sl.stop) {
                    records.push(SvRecord {
                        chrom_no,
                        start,
                        stop,
                        variation_type,
                        pathogenicity,
                        rcv,
                    });
                } else if let (Some(inner_start), Some(inner_stop)) =
                    (sl.inner_start, sl.inner_stop)
                {
                    records.push(SvRecord {
                        chrom_no,
                        start: inner_start,
                        stop: inner_stop,
                        variation_type,
                        pathogenicity,
                        rcv,
                    });
                } else if let (Some(outer_start), Some(outer_stop)) =
                    (sl.outer_start, sl.outer_stop)
                {
                    records.push(SvRecord {
                        chrom_no,
                        start: outer_start,
                        stop: outer_stop,
                        variation_type,
                        pathogenicity,
                        rcv,
                    });
                } else if let (Some(position_vcf), Some(reference_allele_vcf), Some(_)) = (
                    sl.position_vcf,
                    sl.reference_allele_vcf.as_ref(),
                    sl.alternate_allele_vcf.as_ref(),
                ) {
                    records.push(SvRecord {
                        chrom_no,
                        start: position_vcf + 1,
                        stop: position_vcf + reference_allele_vcf.len() as i32,
                        variation_type,
                        pathogenicity,
                        rcv,
                    });
                }
            }
        }
    }
    Ok(records)
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(
    path_input_jsonl: P,
    path_output: Q,
    assembly: input::Assembly,
) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reader = open_read_maybe_gz(path_input_jsonl)?;
    let before_parsing = Instant::now();

    let records = convert_jsonl_to_protobuf(reader, assembly)?;

    let clinvar_db = SvDatabase { records };

    tracing::debug!(
        "total time spent reading {} records: {:?}",
        clinvar_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output)?;
    output_file.write_all(&clinvar_db.encode_to_vec())?;
    output_file.flush()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        clinvar_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}

#[cfg(test)]
mod test {
    use crate::strucvars::txt_to_bin::clinvar::input::Assembly;
    use mehari::common::open_read_maybe_gz;

    #[rstest::rstest]
    #[case(Assembly::Grch37)]
    #[case(Assembly::Grch38)]
    fn run_convert_jsonl_to_protobuf(#[case] assembly: Assembly) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{:?}", assembly);
        let reader = open_read_maybe_gz(
            "tests/db/to-bin/varfish-db-downloader/vardbs/clinvar/clinvar-svs.jsonl",
        )?;
        let records = super::convert_jsonl_to_protobuf(reader, assembly)?;

        insta::assert_yaml_snapshot!(records);

        Ok(())
    }
}
