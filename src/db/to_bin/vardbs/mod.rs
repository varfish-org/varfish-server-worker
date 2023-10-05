//! Code for converting other structural variant database to binary (incl. in-house).

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use anyhow::anyhow;
use prost::Message;
use thousands::Separable;

use crate::common::{build_chrom_map, open_read_maybe_gz, trace_rss_now};
use crate::db;
use crate::db::mk_inhouse::output::Record as InhouseDbRecord;
use crate::db::pbs::{BackgroundDatabase, BgDbRecord};
use crate::sv::query::schema::SvType;

use self::input::InputRecord;

mod input;

/// Known input file types.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum InputFileType {
    Dbvar,
    Dgv,
    DgvGs,
    Exac,
    G1k,
    Gnomad,
    InhouseDb,
}
/// Deserialize from CSV reader to an `Option<records::InputRecord>`
fn deserialize_loop<Rec>(
    reader: &mut csv::Reader<Box<dyn std::io::BufRead>>,
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
            result.push(BgDbRecord {
                chrom_no: *chrom_map
                    .get(&record.chromosome)
                    .unwrap_or_else(|| panic!("unknown chrom: {:?}", &record.chromosome))
                    as i32,
                chrom_no2: *chrom_map
                    .get(&record.chromosome2)
                    .unwrap_or_else(|| panic!("unknown chrom2: {:?}", &record.chromosome2))
                    as i32,
                sv_type: match record.sv_type {
                    SvType::Del => db::pbs::SvType::Del,
                    SvType::Dup => db::pbs::SvType::Dup,
                    SvType::Inv => db::pbs::SvType::Inv,
                    SvType::Ins => db::pbs::SvType::Ins,
                    SvType::Bnd => db::pbs::SvType::Bnd,
                    SvType::Cnv => db::pbs::SvType::Cnv,
                } as i32,
                start: record.begin + 1,
                stop: record.end,
                count: record.count,
            });
        }
    }

    Ok(result)
}

/// Perform conversion to protobuf `.bin` file.
pub fn convert_to_bin<P, Q>(
    path_input_tsv: P,
    path_output_bin: Q,
    input_type: InputFileType,
) -> Result<(), anyhow::Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    // Setup CSV reader for BED file - header is written as comment and must be
    // ignored.
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(path_input_tsv.as_ref())?);
    let before_parsing = Instant::now();

    let records = match input_type {
        InputFileType::Dbvar => deserialize_loop::<input::DbVarRecord>(&mut reader)?,
        InputFileType::Dgv => deserialize_loop::<input::DgvRecord>(&mut reader)?,
        InputFileType::DgvGs => deserialize_loop::<input::DgvGsRecord>(&mut reader)?,
        InputFileType::Exac => deserialize_loop::<input::ExacRecord>(&mut reader)?,
        InputFileType::G1k => deserialize_loop::<input::G1kRecord>(&mut reader)?,
        InputFileType::Gnomad => deserialize_loop::<input::GnomadRecord>(&mut reader)?,
        InputFileType::InhouseDb => deserialize_loop::<InhouseDbRecord>(&mut reader)?,
    };
    let bg_db = BackgroundDatabase { records };

    tracing::debug!(
        "total time spent reading {} records: {:?}",
        bg_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output_bin)?;
    output_file.write_all(&bg_db.encode_to_vec())?;
    output_file.flush()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        bg_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}
