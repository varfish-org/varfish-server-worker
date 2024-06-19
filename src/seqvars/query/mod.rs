//! Code implementing the "seqvars query" sub command.

pub mod annonars;
pub mod interpreter;
pub mod output;
pub mod schema;
pub mod sorting;

use std::io::{BufRead, Write};
use std::time::Instant;

use clap::{command, Parser};
use ext_sort::LimitedBufferBuilder;
use ext_sort::{ExternalSorter, ExternalSorterBuilder};
use itertools::Itertools;

use mehari::annotate::seqvars::CHROM_TO_CHROM_NO;
use noodles::vcf;
use rand_core::{RngCore, SeedableRng};
use thousands::Separable;
use uuid::Uuid;

use crate::common;
use crate::seqvars::query::schema::GenotypeChoice;
use crate::{common::trace_rss_now, common::GenomeRelease};
use mehari::common::noodles::open_vcf_reader;

use self::annonars::Annotator;
use self::schema::CaseQuery;
use self::schema::SequenceVariant;
use self::sorting::{ByCoordinate, ByHgncId};

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
    /// The case UUID.
    #[arg(long)]
    pub case_uuid_id: Option<uuid::Uuid>,
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

/// Utility struct to store statistics about counts.
#[derive(Debug, Default)]
struct QueryStats {
    pub count_passed: usize,
    pub count_total: usize,
    pub by_consequence: indexmap::IndexMap<mehari::annotate::seqvars::ann::Consequence, usize>,
}

/// Checks whether the variants pass through the query interpreter.
///
/// This function is only relevant if the query uses recessive mode.
fn passes_for_gene(
    query: &CaseQuery,
    seqvars: &Vec<SequenceVariant>,
) -> Result<bool, anyhow::Error> {
    #[derive(Debug, Clone, PartialEq, Eq)]
    enum Mode {
        ComphetRecessive,
        Recessive,
        Other,
    }

    let mut mode = Mode::Other;
    let mut index_name = String::default();
    let mut parents = Vec::new();
    query
        .genotype
        .iter()
        .for_each(|(sample_name, genotype_choice)| match genotype_choice {
            Some(GenotypeChoice::ComphetIndex) => {
                index_name.clone_from(sample_name);
                mode = Mode::ComphetRecessive;
            }
            Some(GenotypeChoice::RecessiveIndex) => {
                index_name.clone_from(sample_name);
                mode = Mode::Recessive;
            }
            Some(GenotypeChoice::RecessiveParent) => {
                parents.push(sample_name.clone());
            }
            _ => (),
        });

    eprintln!(
        "mode = {:?}, index_name = {:?}, parents = {:?}",
        mode, &index_name, &parents
    );

    // No special handling for non-recessive mode.
    if mode == Mode::Other {
        return Ok(true);
    }

    let mut seen_het_parents = Vec::new();
    for seqvar in seqvars {
        // Get parsed index genotype.
        let index_gt: common::Genotype = seqvar
            .call_info
            .get(&index_name)
            .expect("no call info for index")
            .genotype
            .as_ref()
            .expect("no GT for index")
            .parse()
            .map_err(|e| anyhow::anyhow!("could not parse index genotype: {}", e))?;

        eprintln!(
            "seqvar = {:?}, index_gt = {:?}, mode = {:?}",
            &seqvar, &index_gt, mode
        );
        if mode == Mode::Recessive && index_gt == common::Genotype::HomAlt {
            // if hom. recessive is allowed then we are done
            return Ok(true);
        } else if mode != Mode::Recessive && index_gt != common::Genotype::Het {
            // it only makese sense to continue in recessive mode if the index is het.
            return Ok(false);
        }

        // Otherwise, the index must be Het and we have to check which parent is also het.
        // At this point, only one parent can be het.
        let parent_gts = parents
            .iter()
            .map(|parent_name| {
                seqvar
                    .call_info
                    .get(parent_name)
                    .expect("no call info for parent")
                    .genotype
                    .as_ref()
                    .expect("no GT for parent")
                    .parse::<common::Genotype>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let het_parents = parents
            .iter()
            .zip(parent_gts.iter())
            .filter(|(_, gt)| **gt == common::Genotype::Het)
            .map(|(name, _)| name.clone())
            .collect::<Vec<_>>();
        assert!(het_parents.len() <= 1);
        eprintln!("het_parents = {:?}", &het_parents);
        if let Some(parent) = het_parents.first() {
            if !seen_het_parents.contains(parent) {
                seen_het_parents.push(parent.clone());
            }
        }

        eprintln!("seen_het_parents = {:?}", &seen_het_parents);

        // If the number of seen het. parents is equal to the number of parents, we are done.
        if seen_het_parents.len() == parents.len() {
            return Ok(true);
        }
    }

    // one of the recessive modes was active and we did not find all parents
    Ok(false)
}

/// Run the `args.path_input` VCF file and run through the given `interpreter` writing to
/// `args.path_output`.
async fn run_query(
    interpreter: &interpreter::QueryInterpreter,
    args: &Args,
    annotator: &annonars::Annotator,
    rng: &mut rand::rngs::StdRng,
) -> Result<QueryStats, anyhow::Error> {
    let tmp_dir = tempfile::TempDir::new()?;

    let chrom_to_chrom_no = &CHROM_TO_CHROM_NO;
    let mut stats = QueryStats::default();

    // Buffer for generating UUIDs.
    let mut uuid_buf = [0u8; 16];

    // Open VCF file, create reader, and read header.
    let mut input_reader = open_vcf_reader(&args.path_input).await.map_err(|e| {
        anyhow::anyhow!("could not open file {} for reading: {}", args.path_input, e)
    })?;
    let input_header = input_reader.read_header().await?;

    let path_unsorted = tmp_dir.path().join("unsorted.jsonl");
    let path_by_hgnc = tmp_dir.path().join("by_hgnc_filtered.jsonl");
    let path_by_coord = tmp_dir.path().join("by_coord.jsonl");

    // Read through input records using the query interpreter as a filter and write to
    // temporary file for unsorted records.
    {
        // Create temporary output file.
        let mut tmp_unsorted = std::fs::File::create(&path_unsorted)
            .map(std::io::BufWriter::new)
            .map_err(|e| anyhow::anyhow!("could not create temporary unsorted file: {}", e))?;

        let mut record_buf = vcf::variant::RecordBuf::default();
        loop {
            let bytes_read = input_reader
                .read_record_buf(&input_header, &mut record_buf)
                .await
                .map_err(|e| {
                    anyhow::anyhow!("problem reading VCF file {}: {}", &args.path_input, e)
                })?;
            if bytes_read == 0 {
                break; // EOF
            }

            stats.count_total += 1;
            let record_seqvar = SequenceVariant::from_vcf(&record_buf, &input_header)
                .map_err(|e| anyhow::anyhow!("could not parse VCF record: {}", e))?;
            tracing::debug!("processing record {:?}", record_seqvar);

            if interpreter.passes(&record_seqvar, annotator)?.pass_all {
                stats.count_passed += 1;
                if let Some(ann) = record_seqvar.ann_fields.first() {
                    ann.consequences.iter().for_each(|csq| {
                        stats
                            .by_consequence
                            .entry(*csq)
                            .and_modify(|e| *e += 1)
                            .or_insert(1);
                    })
                }
                writeln!(
                    tmp_unsorted,
                    "{}",
                    serde_json::to_string(&sorting::ByHgncId::from(record_seqvar))?
                )
                .map_err(|e| anyhow::anyhow!("could not write record to unsorted: {}", e))?;
            }
        }
        tmp_unsorted.into_inner()?.sync_all().map_err(|e| {
            anyhow::anyhow!("could not flush temporary output file unsorted: {}", e)
        })?;
    }

    let elem_count = 10_000; // at most 10k records in memory

    // Now:
    //
    // - sort the records by HGNC ID using external sorting
    // - group by HGNC id
    // - keep the groups where the recessive criteria are met according to query
    // - write out the records again for later sorting by coordinate
    {
        let tmp_unsorted = std::fs::File::open(&path_unsorted)
            .map(std::io::BufReader::new)
            .map_err(|e| anyhow::anyhow!("could not open temporary unsorted file: {}", e))?;
        let mut tmp_by_hgnc_filtered = std::fs::File::create(&path_by_hgnc)
            .map(std::io::BufWriter::new)
            .map_err(|e| {
                anyhow::anyhow!("could not create temporary by_hgnc_filtered file: {}", e)
            })?;

        let sorter: ExternalSorter<sorting::ByHgncId, std::io::Error, LimitedBufferBuilder> =
            ExternalSorterBuilder::new()
                .with_tmp_dir(tmp_dir.as_ref())
                .with_buffer(LimitedBufferBuilder::new(elem_count, false))
                .build()
                .map_err(|e| anyhow::anyhow!("problem creating external sorter: {}", e))?;
        let sorted_iter = sorter
            .sort(tmp_unsorted.lines().map(|res| {
                Ok(serde_json::from_str(&res.expect("problem reading line"))
                    .expect("problem deserializing"))
            }))
            .map_err(|e| anyhow::anyhow!("problem sorting temporary unsorted file: {}", e))?;

        sorted_iter
            .map(|res| res.expect("problem reading line after sorting by HGNC ID"))
            .group_by(|by_hgnc_id| by_hgnc_id.hgnc_id.clone())
            .into_iter()
            .map(|(_, group)| {
                group
                    .map(|ByHgncId { seqvar, .. }| seqvar)
                    .collect::<Vec<_>>()
            })
            .filter(|seqvars| passes_for_gene(&interpreter.query, seqvars).unwrap())
            .for_each(|seqvars| {
                seqvars.into_iter().for_each(|seqvar| {
                    writeln!(
                        tmp_by_hgnc_filtered,
                        "{}",
                        serde_json::to_string(&sorting::ByCoordinate::from(seqvar)).unwrap()
                    )
                    .expect("could not write record to by_hgnc_filtered");
                })
            });
        tmp_by_hgnc_filtered.flush().map_err(|e| {
            anyhow::anyhow!(
                "could not flush temporary output file by_hgnc_filtered: {}",
                e
            )
        })?;
    }

    // Finally:
    // - sort surviving records by coordinate
    // - generate payload with annotations
    {
        let tmp_by_hgnc_filtered = std::fs::File::open(&path_by_hgnc)
            .map(std::io::BufReader::new)
            .map_err(|e| {
                anyhow::anyhow!("could not open temporary tmp_by_hgnc_filtered file: {}", e)
            })?;
        let mut tmp_by_coord = std::fs::File::create(&path_by_coord)
            .map(std::io::BufWriter::new)
            .map_err(|e| anyhow::anyhow!("could not create temporary by_coord file: {}", e))?;

        let sorter: ExternalSorter<sorting::ByCoordinate, std::io::Error, LimitedBufferBuilder> =
            ExternalSorterBuilder::new()
                .with_tmp_dir(tmp_dir.as_ref())
                .with_buffer(LimitedBufferBuilder::new(elem_count, false))
                .build()
                .map_err(|e| anyhow::anyhow!("problem creating external sorter: {}", e))?;
        let sorted_iter = sorter
            .sort(tmp_by_hgnc_filtered.lines().map(|res| {
                Ok(serde_json::from_str(&res.expect("problem reading line"))
                    .expect("problem deserializing"))
            }))
            .map_err(|e| anyhow::anyhow!("problem sorting temporary unsorted file: {}", e))?;

        sorted_iter
            .map(|res| res.expect("problem reading line after sorting by HGNC ID"))
            .for_each(|ByCoordinate { seqvar, .. }| {
                writeln!(tmp_by_coord, "{}", serde_json::to_string(&seqvar).unwrap())
                    .expect("could not write record to by_coord");
            });

        tmp_by_coord.flush().map_err(|e| {
            anyhow::anyhow!(
                "could not flush temporary output file by_hgnc_filtered: {}",
                e
            )
        })?;
    }

    // Finally, perform annotation of the record using the annonars library and write it
    // in TSV format, ready for import into the database.  However, in recessive mode, we
    // have to do a second pass to properly collect compound heterozygous variants.

    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .quote_style(csv::QuoteStyle::Never)
        .from_path(&args.path_output)?;

    let tmp_by_coord = std::fs::File::open(&path_by_coord)
        .map(std::io::BufReader::new)
        .map_err(|e| anyhow::anyhow!("could not open temporary by_coord file: {}", e))?;

    for line in tmp_by_coord.lines() {
        // get next line into a String
        let line = if let Ok(line) = line {
            line
        } else {
            anyhow::bail!("error reading line from input file")
        };
        let seqvar: SequenceVariant = serde_json::from_str(&line).map_err(|e| {
            anyhow::anyhow!(
                "error parsing line from input file: {:?} (line: {:?})",
                e,
                &line
            )
        })?;

        create_payload_and_write_record(
            seqvar,
            annotator,
            chrom_to_chrom_no,
            &mut csv_writer,
            args,
            rng,
            &mut uuid_buf,
        )?;
    }

    Ok(stats)
}

/// Create output payload and write the record to the output file.
fn create_payload_and_write_record(
    seqvar: SequenceVariant,
    annotator: &Annotator,
    chrom_to_chrom_no: &CHROM_TO_CHROM_NO,
    csv_writer: &mut csv::Writer<std::fs::File>,
    args: &Args,
    rng: &mut rand::rngs::StdRng,
    uuid_buf: &mut [u8; 16],
) -> Result<(), anyhow::Error> {
    let result_payload = output::PayloadBuilder::default()
        .case_uuid(args.case_uuid_id.unwrap_or_default())
        .gene_related(
            output::gene_related::Record::with_seqvar_and_annotator(&seqvar, annotator)
                .map_err(|e| anyhow::anyhow!("problem creating gene-related payload: {}", e))?,
        )
        .variant_related(
            output::variant_related::Record::with_seqvar_and_annotator(&seqvar, annotator)
                .map_err(|e| anyhow::anyhow!("problem creating variant-related payload: {}", e))?,
        )
        .call_related(
            output::call_related::Record::with_seqvar(&seqvar)
                .map_err(|e| anyhow::anyhow!("problem creating call-related payload: {}", e))?,
        )
        .build()
        .map_err(|e| anyhow::anyhow!("could not build payload: {}", e))?;
    eprintln!("result_payload = {:?}", &result_payload);
    let start = seqvar.pos;
    let end = start + seqvar.reference.len() as i32 - 1;
    let bin = mehari::annotate::seqvars::binning::bin_from_range(start - 1, end)? as u32;
    let SequenceVariant {
        chrom: chromosome,
        reference,
        alternative,
        ..
    } = seqvar;
    csv_writer
        .serialize(
            &output::RecordBuilder::default()
                .smallvariantqueryresultset_id(args.result_set_id.clone().unwrap_or(".".into()))
                .sodar_uuid(Uuid::from_bytes({
                    rng.fill_bytes(uuid_buf);
                    *uuid_buf
                }))
                .release(match args.genome_release {
                    GenomeRelease::Grch37 => "GRCh37".into(),
                    GenomeRelease::Grch38 => "GRCh38".into(),
                })
                .chromosome_no(
                    *chrom_to_chrom_no
                        .get(&chromosome)
                        .expect("invalid chromosome") as i32,
                )
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .bin(bin)
                .reference(reference)
                .alternative(alternative)
                .payload(
                    serde_json::to_string(&result_payload)
                        .map_err(|e| anyhow::anyhow!("could not serialize payload: {}", e))?,
                )
                .build()
                .map_err(|e| anyhow::anyhow!("could not build record: {}", e))?,
        )
        .map_err(|e| anyhow::anyhow!("could not write record: {}", e))?;
    Ok(())
}

/// Main entry point for `seqvars query` sub command.
pub async fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
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
    let in_memory_dbs = crate::strucvars::query::load_databases(
        &path_worker_db,
        args.genome_release,
        args.max_tad_distance,
    )
    .map_err(|e| {
        anyhow::anyhow!(
            "could not load worker databases from {}: {}",
            path_worker_db,
            e
        )
    })?;
    let annotator = annonars::Annotator::with_path(&args.path_db, args.genome_release)?;
    tracing::info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    tracing::info!("Translating gene allow list...");
    let hgnc_allowlist = if let Some(gene_allowlist) = &query.gene_allowlist {
        if gene_allowlist.is_empty() {
            None
        } else {
            Some(crate::strucvars::query::translate_gene_allowlist(
                gene_allowlist,
                &in_memory_dbs,
            ))
        }
    } else {
        None
    };

    tracing::info!("Running queries...");
    let before_query = Instant::now();
    let query_stats = run_query(
        &interpreter::QueryInterpreter::new(query, hgnc_allowlist),
        args,
        &annotator,
        &mut rng,
    )
    .await?;
    tracing::info!("... done running query in {:?}", before_query.elapsed());
    tracing::info!(
        "summary: {} records passed out of {}",
        query_stats.count_passed.separate_with_commas(),
        query_stats.count_total.separate_with_commas()
    );
    tracing::info!("passing records by effect type");
    for (effect, count) in query_stats.by_consequence.iter() {
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
    use rstest::rstest;

    use super::schema::{CallInfo, SequenceVariant};
    use crate::seqvars::query::schema::{CaseQuery, GenotypeChoice};

    #[rstest]
    #[case(
        GenotypeChoice::ComphetIndex,
        vec!["0/1,0/1,0/0"],
        false,
    )]
    #[case(
        GenotypeChoice::ComphetIndex,
        vec!["1/1,0/0,0/0"],
        false,
    )]
    #[case(
        GenotypeChoice::ComphetIndex,
        vec!["0/1,0/1,0/0","0/1,0/0,0/1"],
        true,
    )]
    #[case(
        GenotypeChoice::RecessiveIndex,
        vec!["1/1,0/0,0/0","0/1,0/0,0/1"],
        true,
    )]
    #[case(
        GenotypeChoice::RecessiveIndex,
        vec!["1/0,1/0,0/0","0/1,0/0,0/1"],
        true,
    )]
    #[case(
        GenotypeChoice::RecessiveIndex,
        vec!["1/0,1/0,0/0","1/0,0/1,0/0"],
        false,
    )]
    fn passes_for_gene_full_trio(
        #[case] gt_choice_index: GenotypeChoice,
        #[case] trio_gts: Vec<&str>,
        #[case] passes: bool,
    ) -> Result<(), anyhow::Error> {
        let query = CaseQuery {
            genotype: vec![
                ("index".into(), Some(gt_choice_index)),
                ("father".into(), Some(GenotypeChoice::RecessiveParent)),
                ("mother".into(), Some(GenotypeChoice::RecessiveParent)),
            ]
            .into_iter()
            .collect(),
            ..Default::default()
        };
        let seqvars = trio_gts
            .iter()
            .map(|gts| {
                let gts: Vec<&str> = gts.split(',').collect();
                SequenceVariant {
                    call_info: vec![
                        (
                            String::from("index"),
                            CallInfo {
                                genotype: Some(gts[0].into()),
                                ..Default::default()
                            },
                        ),
                        (
                            String::from("father"),
                            CallInfo {
                                genotype: Some(gts[1].into()),
                                ..Default::default()
                            },
                        ),
                        (
                            String::from("mother"),
                            CallInfo {
                                genotype: Some(gts[2].into()),
                                ..Default::default()
                            },
                        ),
                    ]
                    .into_iter()
                    .collect(),
                    ..Default::default()
                }
            })
            .collect::<Vec<_>>();

        assert_eq!(super::passes_for_gene(&query, &seqvars)?, passes);

        Ok(())
    }

    #[tracing_test::traced_test]
    #[rstest::rstest]
    #[case("tests/seqvars/query/Case_1.ingested.vcf")]
    #[case("tests/seqvars/query/dragen.ingested.vcf")]
    #[tokio::test]
    async fn smoke_test(#[case] path_input: &str) -> Result<(), anyhow::Error> {
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
            path_output,
            max_results: None,
            rng_seed: Some(42),
            max_tad_distance: 10_000,
            result_set_id: None,
            case_uuid_id: None,
        };
        super::run(&args_common, &args).await?;

        insta::assert_snapshot!(std::fs::read_to_string(args.path_output.as_str())?);

        Ok(())
    }
}
