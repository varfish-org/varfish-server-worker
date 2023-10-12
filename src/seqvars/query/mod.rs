//! Code implementing the "seqvars query" sub command.

mod interpreter;
mod schema;

use std::time::Instant;

use clap::{command, Parser};
use noodles_vcf as vcf;

use mehari::{annotate::seqvars::CHROM_TO_CHROM_NO, common::open_read_maybe_gz};
use rand_core::{RngCore, SeedableRng};
use thousands::Separable;
use uuid::Uuid;

use crate::{common::trace_rss_now, common::GenomeRelease};

use self::schema::SequenceVariant;

/// Command line arguments for `seqvars query` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for seqvars", long_about = None)]
pub struct Args {
    /// Genome release to assume.
    #[arg(long, value_enum)]
    pub genome_release: GenomeRelease,
    /// Result set ID.
    #[arg(long)]
    pub result_set_id: Option<String>,
    /// Path to worker database to use for querying.
    #[arg(long)]
    pub path_db: String,
    /// Path to query JSON file.
    #[arg(long)]
    pub path_query_json: String,
    /// Path to input TSV file.
    #[arg(long)]
    pub path_input: String,
    /// Path to the output TSV file.
    #[arg(long)]
    pub path_output: String,

    /// Optional maximal number of total records to write out.
    #[arg(long)]
    pub max_results: Option<usize>,
    /// Optional seed for RNG.
    #[arg(long)]
    pub rng_seed: Option<u64>,
    /// Maximal distance to TAD to consider (unused, but required when loading database).
    #[arg(long, default_value_t = 10_000)]
    pub max_tad_distance: i32,
}

/// The structured result information of the result record.
#[derive(Debug, Default, serde::Serialize, serde::Deserialize)]
struct ResultPayload {}

/// A result record from the query.
#[derive(Debug, Default, serde::Serialize, serde::Deserialize)]
struct ResultRecord {
    sodar_uuid: Uuid,
    release: String,
    chromosome: String,
    chromosome_no: i32,
    reference: String,
    alternative: String,
    bin: u32,
    start: i32,
    end: i32,
    smallvariantqueryresultset_id: String,
    payload: String,
}

/// Utility struct to store statistics about counts.
#[derive(Debug, Default)]
struct QueryStats {
    pub count_passed: usize,
    pub count_total: usize,
    pub by_effect: indexmap::IndexMap<schema::VariantEffect, usize>,
}

/// Run the `args.path_input` VCF file and run through the given `interpreter` writing to
/// `args.path_output`.
fn run_query(
    interpreter: &interpreter::QueryInterpreter,
    args: &Args,
    _dbs: &crate::strucvars::query::Databases,
    rng: &mut rand::rngs::StdRng,
) -> Result<QueryStats, anyhow::Error> {
    let chrom_to_chrom_no = &CHROM_TO_CHROM_NO;
    let mut stats = QueryStats::default();

    // Open VCF file, create reader, and read header.
    let mut input_reader = open_read_maybe_gz(&args.path_input).map_err(|e| {
        anyhow::anyhow!("could not open file {} for reading: {}", args.path_input, e)
    })?;
    let mut input_reader = vcf::Reader::new(&mut input_reader);
    let input_header = input_reader.read_header()?;

    // Create output TSV writer.
    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .quote_style(csv::QuoteStyle::Never)
        .from_path(&args.path_output)?;

    // Read through input records using the query interpreter as a filter
    for input_record in input_reader.records(&input_header) {
        stats.count_total += 1;
        let record_seqvar = SequenceVariant::from_vcf(&input_record?, &input_header)
            .map_err(|e| anyhow::anyhow!("could not parse VCF record: {}", e))?;
        tracing::debug!("processing record {:?}", record_seqvar);

        let passes = interpreter.passes(&record_seqvar)?;

        let start = record_seqvar.pos;
        let end = start + record_seqvar.reference.len() as i32 - 1;
        let bin = mehari::annotate::seqvars::binning::bin_from_range(start - 1, end)? as u32;

        let result_payload = ResultPayload::default();

        let SequenceVariant {
            chrom: chromosome,
            reference,
            alternative,
            ..
        } = record_seqvar;
        let chromosome_no = *chrom_to_chrom_no
            .get(&chromosome)
            .expect("invalid chromosome") as i32;

        if passes.pass_all {
            // Finally, write out the record.
            let mut uuid_buf = [0u8; 16];
            rng.fill_bytes(&mut uuid_buf);
            csv_writer
                .serialize(&ResultRecord {
                    smallvariantqueryresultset_id: args.result_set_id.clone().unwrap_or(".".into()),
                    sodar_uuid: Uuid::from_bytes(uuid_buf),
                    release: match args.genome_release {
                        GenomeRelease::Grch37 => "GRCh37".into(),
                        GenomeRelease::Grch38 => "GRCh38".into(),
                    },
                    chromosome,
                    chromosome_no,
                    start,
                    end,
                    bin,
                    reference,
                    alternative,
                    payload: serde_json::to_string(&result_payload)
                        .map_err(|e| anyhow::anyhow!("could not serialize payload: {}", e))?,
                })
                .map_err(|e| anyhow::anyhow!("could not write record: {}", e))?;
        }
    }

    Ok(stats)
}

/// Main entry point for `seqvars query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = Instant::now();
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    // Initialize the random number generator from command line seed if given or local entropy
    // source.
    let mut rng = if let Some(rng_seed) = args.rng_seed {
        rand::rngs::StdRng::seed_from_u64(rng_seed)
    } else {
        rand::rngs::StdRng::from_entropy()
    };

    tracing::info!("Loading query...");
    let query: schema::CaseQuery =
        serde_json::from_reader(std::fs::File::open(&args.path_query_json)?)?;
    tracing::info!(
        "... done loading query = {}",
        &serde_json::to_string(&query)?
    );

    tracing::info!("Loading worker databases...");
    let before_loading = Instant::now();
    let path_worker_db = format!("{}/worker", &args.path_db);
    let dbs = crate::strucvars::query::load_databases(
        &path_worker_db,
        args.genome_release,
        args.max_tad_distance,
    )?;
    tracing::info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    tracing::info!("Translating gene allow list...");
    let hgvs_allowlist = if let Some(gene_allowlist) = &query.gene_allowlist {
        if gene_allowlist.is_empty() {
            None
        } else {
            Some(crate::strucvars::query::translate_gene_allowlist(
                gene_allowlist,
                &dbs,
            ))
        }
    } else {
        None
    };

    tracing::info!("Running queries...");
    let before_query = Instant::now();
    let query_stats = run_query(
        &interpreter::QueryInterpreter::new(query, hgvs_allowlist),
        args,
        &dbs,
        &mut rng,
    )?;
    tracing::info!("... done running query in {:?}", before_query.elapsed());
    tracing::info!(
        "summary: {} records passed out of {}",
        query_stats.count_passed.separate_with_commas(),
        query_stats.count_total.separate_with_commas()
    );
    tracing::info!("passing records by effect type");
    for (effect, count) in query_stats.by_effect.iter() {
        tracing::info!("{:?} -- {}", effect, count);
    }

    trace_rss_now();

    tracing::info!(
        "All of `seqvars query` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("tests/seqvars/query/Case_1.ingested.vcf")]
    #[case("tests/seqvars/query/dragen.ingested.vcf")]
    fn smoke_test(#[case] path_input: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path_input.split('/').last().unwrap());

        let tmpdir = temp_testdir::TempDir::default();
        let path_output = format!("{}/out.tsv", tmpdir.to_string_lossy());
        let path_input: String = path_input.into();
        let path_query_json = path_input.replace(".ingested.vcf", ".query.json");

        let args_common = Default::default();
        let args = super::Args {
            genome_release: crate::common::GenomeRelease::Grch37,
            path_db: "tests/seqvars/query/db".into(),
            path_query_json,
            path_input,
            path_output: path_output,
            max_results: None,
            rng_seed: Some(42),
            max_tad_distance: 10_000,
            result_set_id: None,
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(args.path_output.as_str())?);

        Ok(())
    }
}
