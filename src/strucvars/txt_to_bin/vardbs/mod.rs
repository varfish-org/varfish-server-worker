//! Code for converting other structural variant database to binary (incl. in-house).

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use anyhow::anyhow;
use prost::Message;
use thousands::Separable;

use crate::common::{build_chrom_map, trace_rss_now};
use crate::pbs::svs::{BackgroundDatabase, BgDbRecord};
use crate::strucvars::aggregate::output::Record as InhouseDbRecord;
use crate::strucvars::query::schema::SvType;

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
    GnomadSv2,
    GnomadCnv4,
    GnomadSv4,
    InhouseDb,
}
/// Deserialize from CSV reader to an `Option<records::InputRecord>`
fn deserialize_loop<Rec>(
    reader: &mut csv::Reader<Box<dyn std::io::BufRead>>,
) -> Result<Vec<BgDbRecord>, anyhow::Error>
where
    Rec: core::fmt::Debug + TryInto<Option<InputRecord>> + for<'de> serde::Deserialize<'de>,
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
                    SvType::Del => crate::pbs::svs::SvType::Del,
                    SvType::Dup => crate::pbs::svs::SvType::Dup,
                    SvType::Inv => crate::pbs::svs::SvType::Inv,
                    SvType::Ins => crate::pbs::svs::SvType::Ins,
                    SvType::Bnd => crate::pbs::svs::SvType::Bnd,
                    SvType::Cnv => crate::pbs::svs::SvType::Cnv,
                } as i32,
                start: record.begin + 1,
                stop: record.end,
                count: record.count,
            });
        }
    }

    Ok(result)
}

/// Branch around `deserialize_loop`.
pub fn deserialize_branch(
    input_type: InputFileType,
    reader: &mut csv::Reader<Box<dyn std::io::BufRead>>,
) -> Result<Vec<BgDbRecord>, anyhow::Error> {
    match input_type {
        InputFileType::Dbvar => deserialize_loop::<input::DbVarRecord>(reader),
        InputFileType::Dgv => deserialize_loop::<input::DgvRecord>(reader),
        InputFileType::DgvGs => deserialize_loop::<input::DgvGsRecord>(reader),
        InputFileType::Exac => deserialize_loop::<input::ExacRecord>(reader),
        InputFileType::G1k => deserialize_loop::<input::G1kRecord>(reader),
        InputFileType::InhouseDb => deserialize_loop::<InhouseDbRecord>(reader),
        InputFileType::GnomadSv2 => deserialize_loop::<input::GnomadSv2Record>(reader),
        InputFileType::GnomadCnv4 => deserialize_loop::<input::GnomadCnv4Record>(reader),
        InputFileType::GnomadSv4 => deserialize_loop::<input::GnomadSv4Record>(reader),
    }
}

/// Perform conversion to protobuf `.bin` file.
pub fn convert_to_bin<P, Q>(
    path_input_tsv: P,
    path_output: Q,
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
        .from_reader(mehari::common::io::std::open_read_maybe_gz(
            path_input_tsv.as_ref(),
        )?);
    let before_parsing = Instant::now();

    let records = deserialize_branch(input_type, &mut reader)?;
    let bg_db = BackgroundDatabase { records };

    tracing::debug!(
        "total time spent reading {} records: {:?}",
        bg_db.records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut output_file = File::create(&path_output)?;
    output_file.write_all(&bg_db.encode_to_vec())?;
    output_file.sync_all()?;
    tracing::debug!(
        "total time spent writing {} records: {:?}",
        bg_db.records.len().separate_with_commas(),
        before_writing.elapsed()
    );

    Ok(())
}

#[cfg(test)]
mod test {
    use super::InputFileType;

    #[rstest::rstest]
    #[case::dbvar(
        InputFileType::Dbvar,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dbvar.bed.gz"
    )]
    #[case::dgv(
        InputFileType::Dgv,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dgv.bed.gz"
    )]
    #[case::dgv_gs(
        InputFileType::DgvGs,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dgv_gs.bed.gz"
    )]
    #[case::exac(
        InputFileType::Exac,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/exac.bed.gz"
    )]
    #[case::g1k(
        InputFileType::G1k,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/g1k.bed.gz"
    )]
    #[case::gnomad_sv2(
        InputFileType::GnomadSv2,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/gnomad_sv.bed.gz"
    )]
    #[case::gnomad_cnv4(
        InputFileType::GnomadCnv4,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch38/strucvar/gnomad-cnv.bed.gz"
    )]
    #[case::gnomad_sv4(
        InputFileType::GnomadSv4,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch38/strucvar/gnomad-sv.bed.gz"
    )]
    #[case::inhouse_db(
        InputFileType::InhouseDb,
        "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/inhouse.tsv"
    )]
    fn test_deserialize_branch(
        #[case] input_type: InputFileType,
        #[case] path_input: &str,
    ) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!(
            "{:?}-{}",
            input_type,
            path_input
                .split('/')
                .last()
                .unwrap()
                .split('.')
                .next()
                .unwrap()
        );

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .from_reader(mehari::common::io::std::open_read_maybe_gz(path_input)?);

        let records = super::deserialize_branch(input_type, &mut reader)?;
        insta::assert_yaml_snapshot!(records);

        Ok(())
    }
}
