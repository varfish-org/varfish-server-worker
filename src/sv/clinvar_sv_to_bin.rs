//! Code supporting conversion from ClinVar SV TSV to flatbuffers `.bin`.

use std::{fs::File, io::Write, time::Instant};

use clap::Parser;
use thousands::Separable;
use tracing::{debug, info};

use crate::{
    common::{build_chrom_map, open_maybe_gz, trace_rss_now},
    world_flatbuffers::var_fish_server_worker::{
        ClinvarSvDatabase, ClinvarSvDatabaseArgs, ClinvarSvRecord,
    },
};

/// Command line arguments for `sv clinvar-sv-to-bin` sub command.
#[derive(Parser, Debug)]
#[command(about = "Convert ClinVar SV TSV to binary file", long_about = None)]
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
    use std::collections::HashMap;

    use crate::sv::query::schema::{Pathogenicity, VariationType};
    use serde::{de, Deserialize, Deserializer};
    use tracing::warn;

    lazy_static::lazy_static! {
        static ref VARIATION_TYPE_LABELS: HashMap<&'static str, VariationType> = {
            let mut m = HashMap::new();
            m.insert("Complex", VariationType::Complex);
            m.insert("copy number gain", VariationType::Dup);
            m.insert("copy number loss", VariationType::Del);
            m.insert("Deletion", VariationType::Del);
            m.insert("Duplication", VariationType::Dup);
            m.insert("fusion", VariationType::Bnd);
            m.insert("Indel", VariationType::Cnv);
            m.insert("Insertion", VariationType::Ins);
            m.insert("Inversion", VariationType::Inv);
            m.insert("Microsatellite", VariationType::Microsatellite);
            m.insert("Tandem duplication", VariationType::Dup);
            m.insert("Translocation", VariationType::Bnd);
            m
        };
    }

    impl VariationType {
        pub fn from_label(label: &str) -> Result<VariationType, anyhow::Error> {
            if let Some(result) = VARIATION_TYPE_LABELS.get(label) {
                Ok(*result)
            } else {
                Err(anyhow::anyhow!("Invalid VariationType label: {}", label))
            }
        }
    }

    lazy_static::lazy_static! {
        static ref PATHOGENICITY_LABELS: HashMap<&'static str, Pathogenicity> = {
            let mut m = HashMap::new();
            m.insert("{\"benign\"}", Pathogenicity::Benign);
            m.insert("{\"benign\",\"likely benign\"}", Pathogenicity::LikelyBenign);
            m.insert("{\"likely benign\"}", Pathogenicity::LikelyBenign);
            m.insert("{\"likely pathogenic\"}", Pathogenicity::LikelyPathogenic);
            m.insert("{\"likely pathogenic\",\"pathogenic\"}", Pathogenicity::LikelyPathogenic);
            m.insert("{\"pathogenic\"}", Pathogenicity::Pathogenic);
            m.insert("{\"uncertain significance\"}", Pathogenicity::Uncertain);

            m
        };
    }
    impl Pathogenicity {
        pub fn from_label(label: &str) -> Result<Self, anyhow::Error> {
            if let Some(pathogenicity) = PATHOGENICITY_LABELS.get(label) {
                Ok(*pathogenicity)
            } else {
                warn!("Cannot decode pathogenicity from {}", label);
                Ok(Pathogenicity::Uncertain)
            }
        }
    }

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
        /// ClinVar SV variation type
        #[serde(deserialize_with = "from_variation_type_label")]
        pub variation_type: VariationType,
        /// Pathogenicity
        #[serde(
            alias = "summary_clinvar_pathogenicity",
            deserialize_with = "from_pathogenicity_summary"
        )]
        pub pathogenicity: Pathogenicity,
    }

    /// Deserialize "VariationType" from ClinVar TSV file
    ///
    /// This function will strip everything after the first underscore.
    fn from_variation_type_label<'de, D>(deserializer: D) -> Result<VariationType, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        VariationType::from_label(s).map_err(de::Error::custom)
    }
    /// Deserialize "Pathogenicity" from ClinVar TSV file
    ///
    /// This function will strip everything after the first underscore.
    fn from_pathogenicity_summary<'de, D>(deserializer: D) -> Result<Pathogenicity, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        Pathogenicity::from_label(s).map_err(de::Error::custom)
    }
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
        output_records.push(ClinvarSvRecord::new(
            *chrom_map.get(&record.chromosome).expect("unknown chrom") as u8,
            record.start.saturating_sub(1),
            record.end,
            record.variation_type.into(),
            record.pathogenicity.into(),
        ));
    }

    debug!(
        "total time spent reading {} records: {:?}",
        output_records.len().separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();

    let before_writing = Instant::now();
    let mut builder = flatbuffers::FlatBufferBuilder::new();
    let records = builder.create_vector(output_records.as_slice());
    let clinvar_sv_db = ClinvarSvDatabase::create(
        &mut builder,
        &ClinvarSvDatabaseArgs {
            records: Some(records),
        },
    );
    builder.finish_minimal(clinvar_sv_db);
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

/// Main entry point for the `sv clinvar-sv-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting sv clinvar-sv-to-bin");
    info!("common_args = {:?}", &common_args);
    info!("args = {:?}", &args);

    convert_to_bin(args)?;

    Ok(())
}
