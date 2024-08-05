//! Code for converting ClinVar database to binary.

use std::{fs::File, io::BufRead, io::Write, path::Path, time::Instant};

use mehari::common::io::std::open_read_maybe_gz;
use prost::Message;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, trace_rss_now},
    pbs::{
        self,
        varfish::v1::strucvars::bgdb::{Pathogenicity, SvDatabase, SvRecord},
    },
};

// pub mod input;

/// Helper to convert RCV IDs to numbers.
fn numeric_id(raw_id: &str, prefix: &str) -> Result<u32, anyhow::Error> {
    let clean_id: String = raw_id
        .chars()
        .skip(prefix.len())
        .skip_while(|c| *c == '0')
        .collect();
    clean_id
        .parse::<u32>()
        .map_err(|e| anyhow::anyhow!("could not parse RCV id {:?}: {}", &clean_id, &e))
}

/// Read JSONL file and convert to protobuf records.
fn convert_jsonl_to_protobuf(reader: Box<dyn BufRead>) -> Result<Vec<SvRecord>, anyhow::Error> {
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
        let record: annonars::pbs::clinvar_data::extracted_vars::ExtractedVcvRecord = match record {
            Err(e) => {
                tracing::warn!("error deserializing JSONL record: \"{}\" in {}", e, &line);
                continue;
            }
            Ok(record) => record,
        };

        // Extract VCV and RCV accessions as numeric.
        let vcv = numeric_id(&record.accession.expect("no VCV?").accession, "VCV")?;
        // TODO: drop RCV eventually as it is not so useful but VCV was not readily available earlier.
        let rcv = if record.rcvs.is_empty() {
            0
        } else {
            numeric_id(
                &record.rcvs[0]
                    .accession
                    .as_ref()
                    .expect("no RCV?")
                    .accession,
                "RCV",
            )?
        };

        // Obtain germline classification and skip if not set.
        // TODO: later also support somatic classification
        let agc = if let Some(agc) = record
            .classifications
            .as_ref()
            .expect("no classifications?")
            .germline_classification
            .as_ref()
        {
            agc
        } else {
            continue; // no germline classification, skip
        };

        // Convert pathogenicity from upstream to internal protocolbuffers
        let pathogenicity = match agc
            .description
            .clone()
            .unwrap_or_default()
            .to_lowercase()
            .as_str()
        {
            "benign" => Pathogenicity::Benign,
            "benign/likely benign" => Pathogenicity::Benign,
            "likely benign" => Pathogenicity::LikelyBenign,
            "pathogenic" => Pathogenicity::Pathogenic,
            "pathogenic/likely pathogenic" => Pathogenicity::LikelyPathogenic,
            "likely pathogenic" => Pathogenicity::LikelyPathogenic,
            "uncertain significance" => Pathogenicity::Uncertain,
            // NB: we need to improve the protobuf enum
            "conflicting classifications of pathogenicity" => Pathogenicity::Uncertain,
            _ => {
                continue;
            }
        };

        // Convert variation type from upstream to internal protocolbuffers.
        let variation_type = match annonars::pbs::clinvar_data::extracted_vars::VariationType::try_from(record.variation_type) {
            Ok(variation_type) => match variation_type {
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Unspecified => continue,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Insertion => pbs::varfish::v1::strucvars::bgdb::VariationType::Ins,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Deletion => pbs::varfish::v1::strucvars::bgdb::VariationType::Del,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Snv => continue,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Indel => continue,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Duplication => pbs::varfish::v1::strucvars::bgdb::VariationType::Dup,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::TandemDuplication => pbs::varfish::v1::strucvars::bgdb::VariationType::Dup,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::StructuralVariant => pbs::varfish::v1::strucvars::bgdb::VariationType::Complex,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::CopyNumberGain => pbs::varfish::v1::strucvars::bgdb::VariationType::Dup,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::CopyNumberLoss => pbs::varfish::v1::strucvars::bgdb::VariationType::Del,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::ProteinOnly => continue,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Microsatellite => pbs::varfish::v1::strucvars::bgdb::VariationType::Microsatellite,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Inversion => pbs::varfish::v1::strucvars::bgdb::VariationType::Inv,
                annonars::pbs::clinvar_data::extracted_vars::VariationType::Other => continue,
            },
            Err(_) => {
                tracing::warn!("unknown variation type: {}", record.variation_type);
                continue;
            }
        };

        // Finally, extract sequence location and add record.
        let sl = record
            .sequence_location
            .as_ref()
            .expect("no sequence_location");
        let chrom = annonars::pbs::clinvar_data::clinvar_public::Chromosome::try_from(sl.chr)
            .map_err(|e| anyhow::anyhow!("invalid chromosome: {}: {}", sl.chr, e))?
            .as_str_name()
            .strip_prefix("CHROMOSOME_")
            .map(|s| s.to_string())
            .unwrap_or_default();
        let chrom_no = if let Some(chrom_no) = chrom_map.get(&chrom) {
            *chrom_no as i32
        } else {
            tracing::warn!("unknown chromosome {}", &sl.chr);
            continue;
        };

        if let (Some(start), Some(stop)) = (sl.start, sl.stop) {
            records.push(SvRecord {
                chrom_no,
                start: start as i32,
                stop: stop as i32,
                variation_type: variation_type as i32,
                pathogenicity: pathogenicity as i32,
                rcv,
                vcv,
            });
        } else if let (Some(inner_start), Some(inner_stop)) = (sl.inner_start, sl.inner_stop) {
            records.push(SvRecord {
                chrom_no,
                start: inner_start as i32,
                stop: inner_stop as i32,
                variation_type: variation_type as i32,
                pathogenicity: pathogenicity as i32,
                rcv,
                vcv,
            });
        } else if let (Some(outer_start), Some(outer_stop)) = (sl.outer_start, sl.outer_stop) {
            records.push(SvRecord {
                chrom_no,
                start: outer_start as i32,
                stop: outer_stop as i32,
                variation_type: variation_type as i32,
                pathogenicity: pathogenicity as i32,
                rcv,
                vcv,
            });
        } else if let (Some(position_vcf), Some(reference_allele_vcf), Some(_)) = (
            sl.position_vcf,
            sl.reference_allele_vcf.as_ref(),
            sl.alternate_allele_vcf.as_ref(),
        ) {
            records.push(SvRecord {
                chrom_no,
                start: position_vcf as i32 + 1,
                stop: position_vcf as i32 + reference_allele_vcf.len() as i32,
                variation_type: variation_type as i32,
                pathogenicity: pathogenicity as i32,
                rcv,
                vcv,
            });
        }
    }
    Ok(records)
}

/// Perform conversion to protocolbuffers `.bin` file.
pub fn convert_to_bin<P, Q>(path_input_jsonl: P, path_output: Q) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let reader = open_read_maybe_gz(path_input_jsonl)?;
    let before_parsing = Instant::now();

    let records = convert_jsonl_to_protobuf(reader)?;

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
    output_file.sync_all()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        clinvar_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}

#[cfg(test)]
mod test {
    #[tracing_test::traced_test]
    #[rstest::rstest]
    fn run_convert_jsonl_to_protobuf() -> Result<(), anyhow::Error> {
        let reader = mehari::common::io::std::open_read_maybe_gz(
            "tests/db/to-bin/varfish-db-downloader/vardbs/clinvar/clinvar-svs.jsonl",
        )?;
        let records = super::convert_jsonl_to_protobuf(reader)?;

        insta::assert_yaml_snapshot!(records);

        Ok(())
    }
}
