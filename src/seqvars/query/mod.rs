//! Code implementing the "seqvars query" sub command.

pub mod annonars;
pub mod hpo;
pub mod interpreter;
pub mod schema;
pub mod sorting;

use std::collections::BTreeSet;
use std::io::{BufRead, Write};
use std::time::Instant;

use clap::{command, Parser};
use ext_sort::LimitedBufferBuilder;
use ext_sort::{ExternalSorter, ExternalSorterBuilder};

use futures::TryStreamExt as _;
use itertools::Itertools as _;
use mehari::annotate::seqvars::CHROM_TO_CHROM_NO;
use mehari::common::noodles::NoodlesVariantReader as _;
use rand_core::{RngCore, SeedableRng};
use schema::data::{TryFromVcf as _, VariantRecord};
use schema::query::{CaseQuery, GenotypeChoice, RecessiveMode, SampleGenotypeChoice};
use thousands::Separable;
use tokio::io::AsyncWriteExt as _;
use uuid::Uuid;

use crate::common;
use crate::pbs::varfish::v1::seqvars::output as pbs_output;
use crate::pbs::varfish::v1::seqvars::query as pbs_query;
use crate::{common::trace_rss_now, common::GenomeRelease};

use self::annonars::Annotator;
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
    pub case_uuid: Option<uuid::Uuid>,
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
    pub passed_by_consequences:
        indexmap::IndexMap<mehari::annotate::seqvars::ann::Consequence, usize>,
}

/// Checks whether the variants pass through the query interpreter.
fn passes_for_gene(query: &CaseQuery, seqvars: &Vec<VariantRecord>) -> Result<bool, anyhow::Error> {
    // Short-circuit in case of disabled recessive mode.
    if query.genotype.recessive_mode == RecessiveMode::Disabled {
        return Ok(true);
    }

    // Extract family information for recessive mode.
    let (index, parents) = {
        let mut index = String::new();
        let mut parents = Vec::new();
        for (sample_name, SampleGenotypeChoice { genotype, .. }) in
            query.genotype.sample_genotypes.iter()
        {
            match genotype {
                GenotypeChoice::RecessiveIndex => {
                    index.clone_from(sample_name);
                }
                GenotypeChoice::RecessiveFather | GenotypeChoice::RecessiveMother => {
                    parents.push(sample_name.clone());
                }
                _ => (),
            }
        }
        (index, parents)
    };
    tracing::debug!("index = {}, parents ={:?}", &index, &parents);

    // All parents must have been seen as het. and hom. ref. at least once for compound
    // heterozygous mode.
    let mut seen_het_parents = BTreeSet::new();
    let mut seen_ref_parents = BTreeSet::new();
    let mut seen_index_het: usize = 0;

    // Go over all variants and try to find single variant compatible with hom. recessive
    // mode or at least two variants compatible with compound heterozygous mode.
    for seqvar in seqvars {
        // Get parsed index genotype.
        let index_gt: common::Genotype = seqvar
            .call_infos
            .get(&index)
            .expect("no call info for index")
            .genotype
            .as_ref()
            .expect("no GT for index")
            .parse()
            .map_err(|e| anyhow::anyhow!("could not parse index genotype: {}", e))?;

        tracing::debug!("seqvar = {:?}, index_gt = {:?}", &seqvar, &index_gt);

        // Get parent genotypes and count hom. alt parents and het. parents.
        let parent_gts = parents
            .iter()
            .map(|parent_name| {
                seqvar
                    .call_infos
                    .get(parent_name)
                    .expect("no call info for parent")
                    .genotype
                    .as_ref()
                    .expect("no GT for parent")
                    .parse::<common::Genotype>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let homalt_parents = parents
            .iter()
            .zip(parent_gts.iter())
            .filter(|(_, gt)| **gt == common::Genotype::HomAlt)
            .map(|(name, _)| name.clone())
            .collect::<Vec<_>>();
        let het_parents = parents
            .iter()
            .zip(parent_gts.iter())
            .filter(|(_, gt)| **gt == common::Genotype::Het)
            .map(|(name, _)| name.clone())
            .collect::<Vec<_>>();
        let ref_parents = parents
            .iter()
            .zip(parent_gts.iter())
            .filter(|(_, gt)| **gt == common::Genotype::HomRef)
            .map(|(name, _)| name.clone())
            .collect::<Vec<_>>();
        tracing::debug!(
            "seqvar = {:?}, homalt_parents = {:?}, het_parents = {:?}, ref_parents = {:?}",
            &seqvar,
            &homalt_parents,
            &het_parents,
            &ref_parents
        );
        if !homalt_parents.is_empty() {
            // Skip this variant, found homozygous parent.
            continue;
        }

        // We can pass in two cases:
        //
        // 1. index hom. alt, both parents het.
        // 2. index het, one parent het., other is ref.

        if index_gt == common::Genotype::HomAlt {
            if matches!(
                query.genotype.recessive_mode,
                RecessiveMode::Homozygous | RecessiveMode::Any
            ) {
                // Case 1: index hom. alt, any given parent must be het.
                if het_parents.len() != parent_gts.len() {
                    // Skip this variant, any given parent must be het.
                    continue;
                } else {
                    // All good, this variant supports the recessive mode for the gene.
                    return Ok(true);
                }
            }
        } else if index_gt == common::Genotype::Het {
            if matches!(
                query.genotype.recessive_mode,
                RecessiveMode::CompoundHeterozygous | RecessiveMode::Any
            ) {
                // Case 2: index het, one parent het./other. ref.?
                match parent_gts.len() {
                    0 => {
                        // No parents, all good.
                    }
                    1 => {
                        // Single parent, must be het. or hom. ref.
                        if het_parents.len() == 1 {
                            seen_het_parents
                                .insert(het_parents.into_iter().next().expect("checked above"));
                        } else if ref_parents.len() == 1 {
                            seen_ref_parents
                                .insert(ref_parents.into_iter().next().expect("checked above"));
                        } else {
                            // Skip this variant, single parent not het. or hom. ref.
                            continue;
                        }
                    }
                    2 => {
                        // Two parents, one must be het. and the other hom. ref.
                        if het_parents.len() == 1 && ref_parents.len() == 1 {
                            seen_het_parents
                                .insert(het_parents.into_iter().next().expect("checked above"));
                            seen_ref_parents
                                .insert(ref_parents.into_iter().next().expect("checked above"));
                        } else {
                            // Skip this variant, no comp. het. pattern.
                            continue;
                        }
                    }
                    _ => unreachable!("More than two parents?"),
                }
                seen_index_het += 1;
            }
        } else {
            // Skip this variant, index is ref.
            continue;
        }
    }

    Ok(
        if matches!(
            query.genotype.recessive_mode,
            RecessiveMode::CompoundHeterozygous | RecessiveMode::Any
        ) {
            // Check recessive condition.  We need to have at least two variants and all parents must
            // have been seen as het. and hom. ref.
            seen_index_het >= 2
                && seen_het_parents.len() == parents.len()
                && seen_ref_parents.len() == parents.len()
        } else {
            false
        },
    )
}

/// Run the `args.path_input` VCF file and run through the given `interpreter` writing to
/// `args.path_output`.
async fn run_query(
    interpreter: &interpreter::QueryInterpreter,
    pb_query: &pbs_query::CaseQuery,
    args: &Args,
    annotator: &annonars::Annotator,
    rng: &mut rand::rngs::StdRng,
) -> Result<QueryStats, anyhow::Error> {
    let start_time = common::now_as_pbjson_timestamp();
    let tmp_dir = tempfile::TempDir::new()?;

    let chrom_to_chrom_no = &CHROM_TO_CHROM_NO;
    let mut stats = QueryStats::default();

    // Buffer for generating UUIDs.
    let mut uuid_buf = [0u8; 16];

    // Open VCF file, create reader, and read header.
    let mut input_reader = common::noodles::open_vcf_reader(&args.path_input)
        .await
        .map_err(|e| {
            anyhow::anyhow!("could not open file {} for reading: {}", args.path_input, e)
        })?;
    let input_header = input_reader.read_header().await?;

    let path_unsorted = tmp_dir.path().join("unsorted.jsonl");
    let path_by_hgnc = tmp_dir.path().join("by_hgnc_filtered.jsonl");
    let path_by_coord = tmp_dir.path().join("by_coord.jsonl");
    let path_noheader = tmp_dir.path().join("noheader.jsonl");

    // Read through input records using the query interpreter as a filter and write to
    // temporary file for unsorted records.
    {
        // Create temporary output file.
        let mut tmp_unsorted = std::fs::File::create(&path_unsorted)
            .map(std::io::BufWriter::new)
            .map_err(|e| anyhow::anyhow!("could not create temporary unsorted file: {}", e))?;

        let mut records = input_reader.records(&input_header).await;
        while let Some(record_buf) = records.try_next().await? {
            stats.count_total += 1;
            let record_seqvar = VariantRecord::try_from_vcf(&record_buf, &input_header)
                .map_err(|e| anyhow::anyhow!("could not parse VCF record: {}", e))?;
            tracing::trace!("processing record {:?}", record_seqvar);

            if interpreter.passes(&record_seqvar, annotator)?.pass_all {
                stats.count_passed += 1;
                if let Some(ann) = record_seqvar.ann_fields.first() {
                    ann.consequences.iter().for_each(|csq| {
                        stats
                            .passed_by_consequences
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
            .chunk_by(|by_hgnc_id| by_hgnc_id.hgnc_id.clone())
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

    // Perform the annotation and write into file without header.
    {
        tracing::debug!("writing noheader file {}", path_noheader.display());
        let writer = tokio::fs::OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&path_noheader)
            .await
            .map_err(|e| anyhow::anyhow!("could not open output file: {}", e))?;
        let mut writer = tokio::io::BufWriter::new(writer);
        // Open reader for temporary by-coordinate file.
        let tmp_by_coord = std::fs::File::open(&path_by_coord)
            .map(std::io::BufReader::new)
            .map_err(|e| anyhow::anyhow!("could not open temporary by_coord file: {}", e))?;
        // Iterate through the temporary by-coordinate file, generate and write output records.
        for line in tmp_by_coord.lines() {
            // get next line into a String
            let line = if let Ok(line) = line {
                line
            } else {
                anyhow::bail!("error reading line from input file")
            };
            let seqvar: VariantRecord = serde_json::from_str(&line).map_err(|e| {
                anyhow::anyhow!(
                    "error parsing line from input file: {:?} (line: {:?})",
                    e,
                    &line
                )
            })?;

            create_and_write_record(
                seqvar,
                annotator,
                chrom_to_chrom_no,
                &mut writer,
                args,
                rng,
                &mut uuid_buf,
            )
            .await?;
        }

        // Properly flush the output file, so upload to S3 can be done if necessary.
        writer
            .flush()
            .await
            .map_err(|e| anyhow::anyhow!("could not flush output file before closing: {}", e))?;
    }

    // Finally, write out records in JSONL format in JSONL format.  The first line will contain the
    // header, the rest the records.
    //
    // Use output helper for semi-transparent upload to S3.
    let out_path_helper = crate::common::s3::OutputPathHelper::new(&args.path_output)?;
    {
        tracing::debug!("writing file {}", out_path_helper.path_out());
        // Open output file for writing (potentially temporary, then uploaded to S3 via helper).
        let file = std::fs::OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(out_path_helper.path_out())
            .map_err(|e| anyhow::anyhow!("could not open output file: {}", e))?;
        let mut writer = std::io::BufWriter::new(file);
        write_header(args, pb_query, &stats, start_time, &mut writer)?;
        // Open reader for file without header.
        let mut reader = std::fs::File::open(&path_noheader)
            .map(std::io::BufReader::new)
            .map_err(|e| anyhow::anyhow!("could not open temporary no_header file: {}", e))?;
        // Append the temporary file to the output file.
        std::io::copy(&mut reader, &mut writer)
            .map_err(|e| anyhow::anyhow!("could not copy temporary file to output file: {}", e))?;
        // Properly flush the output file, so upload to S3 can be done if necessary.
        writer
            .flush()
            .map_err(|e| anyhow::anyhow!("could not flush output file before closing: {}", e))?;
    }
    // Potentially upload the output file to S3.
    out_path_helper
        .upload_for_s3()
        .await
        .map_err(|e| anyhow::anyhow!("could not upload output file to S3: {}", e))?;

    Ok(stats)
}

/// Write the header to the output file.
fn write_header(
    args: &Args,
    pb_query: &pbs_query::CaseQuery,
    stats: &QueryStats,
    start_time: pbjson_types::Timestamp,
    writer: &mut std::io::BufWriter<std::fs::File>,
) -> Result<(), anyhow::Error> {
    let header = pbs_output::OutputHeader {
        genome_release: Into::<pbs_output::GenomeRelease>::into(args.genome_release) as i32,
        versions: vec![pbs_output::VersionEntry {
            name: "varfish-worker".to_string(),
            version: common::worker_version().to_string(),
        }],
        query: Some(pb_query.clone()),
        case_uuid: args.case_uuid.unwrap_or_default().to_string(),
        statistics: Some(pbs_output::OutputStatistics {
            count_total: stats.count_total as u64,
            count_passed: stats.count_passed as u64,
            passed_by_consequences: stats
                .passed_by_consequences
                .iter()
                .filter_map(|(csq, count)| -> Option<pbs_output::ConsequenceCount> {
                    // We ignore consequences that don't have a mapping into the protobuf.
                    if let Ok(csq) = TryInto::<pbs_query::Consequence>::try_into(*csq) {
                        Some(pbs_output::ConsequenceCount {
                            consequence: csq as i32,
                            count: *count as u32,
                        })
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        }),
        resources: Some(pbs_output::ResourcesUsed {
            start_time: Some(start_time),
            end_time: Some(common::now_as_pbjson_timestamp()),
            memory_used: common::rss_size()?,
        }),
        variant_score_columns: variant_related_annotation::score_columns(),
    };
    writeln!(
        writer,
        "{}",
        serde_json::to_string(&header)
            .map_err(|e| anyhow::anyhow!("could not convert header to JSON: {}", e))?
    )?;
    Ok(())
}

/// Trait for records that can be constructed from a `VariantRecord` and an `Annotator`.
trait WithSeqvarAndAnnotator: Sized {
    /// The error type to use.
    type Error;

    /// Construct the record from the given `seqvar` and `annotator`.
    ///
    /// # Arguments
    ///
    /// * `seqvar` - The variant record to use for construction.
    /// * `annotator` - The annotator to use for construction.
    ///
    /// # Error
    ///
    /// Returns an error if the record could not be constructed.
    fn with_seqvar_and_annotator(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Self, anyhow::Error>;
}

impl TryInto<crate::pbs::varfish::v1::seqvars::query::Consequence>
    for mehari::annotate::seqvars::ann::Consequence
{
    type Error = anyhow::Error;

    fn try_into(self) -> Result<crate::pbs::varfish::v1::seqvars::query::Consequence, Self::Error> {
        use crate::pbs::varfish::v1::seqvars::query;

        Ok(match self {
            mehari::annotate::seqvars::ann::Consequence::TranscriptAblation => {
                query::Consequence::TranscriptAblation
            }
            mehari::annotate::seqvars::ann::Consequence::ExonLossVariant => {
                query::Consequence::ExonLossVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SpliceAcceptorVariant => {
                query::Consequence::SpliceAcceptorVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SpliceDonorVariant => {
                query::Consequence::SpliceDonorVariant
            }
            mehari::annotate::seqvars::ann::Consequence::StopGained => {
                query::Consequence::StopGained
            }
            mehari::annotate::seqvars::ann::Consequence::FrameshiftVariant => {
                query::Consequence::FrameshiftVariant
            }
            mehari::annotate::seqvars::ann::Consequence::StopLost => query::Consequence::StopLost,
            mehari::annotate::seqvars::ann::Consequence::StartLost => query::Consequence::StartLost,
            mehari::annotate::seqvars::ann::Consequence::TranscriptAmplification => {
                query::Consequence::TranscriptAmplification
            }
            mehari::annotate::seqvars::ann::Consequence::DisruptiveInframeInsertion => {
                query::Consequence::DisruptiveInframeInsertion
            }
            mehari::annotate::seqvars::ann::Consequence::DisruptiveInframeDeletion => {
                query::Consequence::DisruptiveInframeDeletion
            }
            mehari::annotate::seqvars::ann::Consequence::ConservativeInframeInsertion => {
                query::Consequence::ConservativeInframeInsertion
            }
            mehari::annotate::seqvars::ann::Consequence::ConservativeInframeDeletion => {
                query::Consequence::ConservativeInframeDeletion
            }
            mehari::annotate::seqvars::ann::Consequence::MissenseVariant => {
                query::Consequence::MissenseVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SpliceDonorFifthBaseVariant => {
                query::Consequence::SpliceDonorFifthBaseVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SpliceRegionVariant => {
                query::Consequence::SpliceRegionVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SpliceDonorRegionVariant => {
                query::Consequence::SpliceDonorRegionVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SplicePolypyrimidineTractVariant => {
                query::Consequence::SplicePolypyrimidineTractVariant
            }
            mehari::annotate::seqvars::ann::Consequence::StartRetainedVariant => {
                query::Consequence::StartRetainedVariant
            }
            mehari::annotate::seqvars::ann::Consequence::StopRetainedVariant => {
                query::Consequence::StopRetainedVariant
            }
            mehari::annotate::seqvars::ann::Consequence::SynonymousVariant => {
                query::Consequence::SynonymousVariant
            }
            mehari::annotate::seqvars::ann::Consequence::CodingSequenceVariant => {
                query::Consequence::CodingSequenceVariant
            }
            mehari::annotate::seqvars::ann::Consequence::FivePrimeUtrExonVariant => {
                query::Consequence::FivePrimeUtrExonVariant
            }
            mehari::annotate::seqvars::ann::Consequence::FivePrimeUtrIntronVariant => {
                query::Consequence::FivePrimeUtrIntronVariant
            }
            mehari::annotate::seqvars::ann::Consequence::ThreePrimeUtrExonVariant => {
                query::Consequence::ThreePrimeUtrExonVariant
            }
            mehari::annotate::seqvars::ann::Consequence::ThreePrimeUtrIntronVariant => {
                query::Consequence::ThreePrimeUtrIntronVariant
            }
            mehari::annotate::seqvars::ann::Consequence::NonCodingTranscriptExonVariant => {
                query::Consequence::NonCodingTranscriptExonVariant
            }
            mehari::annotate::seqvars::ann::Consequence::NonCodingTranscriptIntronVariant => {
                query::Consequence::NonCodingTranscriptIntronVariant
            }
            mehari::annotate::seqvars::ann::Consequence::UpstreamGeneVariant => {
                query::Consequence::UpstreamGeneVariant
            }
            mehari::annotate::seqvars::ann::Consequence::DownstreamGeneVariant => {
                query::Consequence::DownstreamGeneVariant
            }
            mehari::annotate::seqvars::ann::Consequence::IntergenicVariant => {
                query::Consequence::IntergenicVariant
            }
            mehari::annotate::seqvars::ann::Consequence::IntronVariant => {
                query::Consequence::IntronVariant
            }
            mehari::annotate::seqvars::ann::Consequence::FeatureElongation
            | mehari::annotate::seqvars::ann::Consequence::FeatureTruncation
            | mehari::annotate::seqvars::ann::Consequence::MatureMirnaVariant
            | mehari::annotate::seqvars::ann::Consequence::TfbsAblation
            | mehari::annotate::seqvars::ann::Consequence::TfbsAmplification
            | mehari::annotate::seqvars::ann::Consequence::TfBindingSiteVariant
            | mehari::annotate::seqvars::ann::Consequence::RegulatoryRegionAblation
            | mehari::annotate::seqvars::ann::Consequence::RegulatoryRegionAmplification
            | mehari::annotate::seqvars::ann::Consequence::RegulatoryRegionVariant
            | mehari::annotate::seqvars::ann::Consequence::GeneVariant => {
                return Err(anyhow::anyhow!("Unsupported consequence: {:?}", self))
            }
        })
    }
}

impl WithSeqvarAndAnnotator for pbs_output::GeneRelatedAnnotation {
    type Error = anyhow::Error;

    fn with_seqvar_and_annotator(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Self, Self::Error> {
        if let Some(ann) = seqvar.ann_fields.first() {
            if !ann.gene_id.is_empty() && !ann.gene_symbol.is_empty() {
                let hgnc_id = ann.gene_id.clone();

                let gene_record = annotator
                    .query_genes(&hgnc_id)
                    .map_err(|e| anyhow::anyhow!("problem querying genes database: {}", e))?;
                let mois = annotator.hgnc_to_moi.get(&hgnc_id);

                return Ok(Self {
                    identity: Some(pbs_output::GeneIdentity {
                        hgnc_id: hgnc_id.clone(),
                        gene_symbol: ann.gene_symbol.clone(),
                    }),
                    consequences: gene_related_annotation::consequences(ann)?,
                    phenotypes: gene_related_annotation::phenotypes(&gene_record, mois),
                    constraints: gene_related_annotation::constraints(&gene_record)?,
                });
            }
        }

        Ok(Default::default())
    }
}

/// Supporting code for pbs_output::GeneRelatedAnnotation.
mod gene_related_annotation {
    use mehari::annotate::seqvars::ann;

    use super::*;

    pub(crate) fn consequences(
        ann: &ann::AnnField,
    ) -> Result<Option<pbs_output::GeneRelatedConsequences>, anyhow::Error> {
        let location = if ann.distance.is_none() {
            pbs_output::VariantLocation::Exon
        } else if ann
            .consequences
            .contains(&ann::Consequence::UpstreamGeneVariant)
        {
            pbs_output::VariantLocation::Upstream
        } else if ann
            .consequences
            .contains(&ann::Consequence::DownstreamGeneVariant)
        {
            pbs_output::VariantLocation::Downstream
        } else {
            pbs_output::VariantLocation::Intron
        };

        let (tx_accession, tx_version) = if ann.feature_id.is_empty() {
            (None, None)
        } else {
            let (accession, version) = ann
                .feature_id
                .split_once('.')
                .unwrap_or((&ann.feature_id, ""));
            let tx_accession = Some(accession.to_string());
            let tx_version = if version.is_empty() {
                None
            } else {
                version.parse::<i32>().ok()
            };
            (tx_accession, tx_version)
        };

        let (rank_ord, rank_total) = if let Some(rank) = ann.rank.as_ref() {
            (Some(rank.ord), Some(rank.total))
        } else {
            (None, None)
        };

        Ok(Some(pbs_output::GeneRelatedConsequences {
            hgvs_t: ann.hgvs_t.clone(),
            hgvs_p: ann.hgvs_p.clone(),
            consequences: ann
                .consequences
                .iter()
                .filter_map(|csq| -> Option<i32> {
                    // We ignore consequences that don't have a mapping into the protobuf.
                    if let Ok(csq) = TryInto::<pbs_query::Consequence>::try_into(*csq) {
                        Some(csq as i32)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
            tx_accession,
            tx_version,
            location: location as i32,
            rank_ord,
            rank_total,
        }))
    }

    pub(crate) fn phenotypes(
        gene_record: &Option<::annonars::pbs::genes::base::Record>,
        mois: Option<&indexmap::IndexSet<hpo::ModeOfInheritance>>,
    ) -> Option<pbs_output::GeneRelatedPhenotypes> {
        gene_record
            .as_ref()
            .map(|gene_record| pbs_output::GeneRelatedPhenotypes {
                is_acmg_sf: gene_record.acmg_sf.is_some(),
                is_disease_gene: gene_record.omim.is_some() || gene_record.orpha.is_some(),
                mode_of_inheritances: mois
                    .cloned()
                    .unwrap_or_default()
                    .into_iter()
                    .map(|moi| Into::<pbs_output::ModeOfInheritance>::into(moi) as i32)
                    .collect::<Vec<_>>(),
            })
    }

    pub(crate) fn constraints(
        gene_record: &Option<::annonars::pbs::genes::base::Record>,
    ) -> Result<Option<pbs_output::GeneRelatedConstraints>, anyhow::Error> {
        gene_record
            .as_ref()
            .map(
                |gene_record| -> Result<pbs_output::GeneRelatedConstraints, anyhow::Error> {
                    let gnomad =
                        gene_record
                            .gnomad_constraints
                            .as_ref()
                            .map(|gnomad_constraints| pbs_output::GnomadConstraints {
                                mis_z: gnomad_constraints
                                    .mis_z
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_lof: gnomad_constraints
                                    .oe_lof
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_lof_lower: gnomad_constraints
                                    .oe_lof_lower
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_lof_upper: gnomad_constraints
                                    .oe_lof_upper
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_mis: gnomad_constraints
                                    .oe_mis
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_mis_lower: gnomad_constraints
                                    .oe_mis_lower
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                oe_mis_upper: gnomad_constraints
                                    .oe_mis_upper
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                                pli: gnomad_constraints.pli.map(|x| x as f32).unwrap_or_default(),
                                syn_z: gnomad_constraints
                                    .syn_z
                                    .map(|x| x as f32)
                                    .unwrap_or_default(),
                            });
                    let decipher = gene_record.decipher_hi.as_ref().map(|decipher_hi| {
                        pbs_output::DecipherConstraints {
                            p_hi: decipher_hi.p_hi as f32,
                            hi_index: decipher_hi.hi_index as f32,
                        }
                    });
                    let rcnv = gene_record
                        .rcnv
                        .as_ref()
                        .map(|rcnv| pbs_output::RcnvConstraints {
                            p_haplo: rcnv.p_haplo as f32,
                            p_triplo: rcnv.p_triplo as f32,
                        });
                    let shet = gene_record
                        .shet
                        .as_ref()
                        .map(|shet| pbs_output::ShetConstraints {
                            s_het: shet.s_het as f32,
                        });
                    let clingen = gene_record
                        .clingen
                        .as_ref()
                        .map(|clingen| -> Result<pbs_output::ClingenDosageAnnotation, anyhow::Error> {
                            Ok(pbs_output::ClingenDosageAnnotation {
                                haplo: pbs_output::ClingenDosageScore::try_from(
                                    clingen.haploinsufficiency_score,
                                )
                                .map_err(|e| {
                                    anyhow::anyhow!(
                                        "could not convert haploinsufficiency score: {}",
                                        e
                                    )
                                })
                                .map(|x| x as i32)?,
                                triplo: pbs_output::ClingenDosageScore::try_from(
                                    clingen.triplosensitivity_score,
                                )
                                .map_err(|e| {
                                    anyhow::anyhow!(
                                        "could not convert triplosensitivity score: {}",
                                        e
                                    )
                                })
                                .map(|x| x as i32)?,
                            })
                        })
                        .transpose()?;
                    Ok(pbs_output::GeneRelatedConstraints {
                        gnomad,
                        decipher,
                        rcnv,
                        shet,
                        clingen,
                    })
                },
            )
            .transpose()
    }
}

/// Helper code for pbs_output::VariantRelatedAnnotation.
mod variant_related_annotation {
    use crate::pbs::varfish::v1::seqvars::output as pbs_output;
    use schema::data::Af;

    use super::*;

    /// Helper modules for score collection.
    pub mod score_collection {
        /// Trait for score collection.
        pub trait Collector {
            /// Register one column value.
            fn register(&mut self, column_name: &str, value: &serde_json::Value);
            /// Write collected scores to the given `indexmap::IndexMap``.
            fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>);
        }

        /// Simple implementation for collecting a single score.
        #[derive(Debug, Clone)]
        pub struct SingleValueCollector {
            /// The column name to collect.
            pub input_column_name: String,
            /// The output column name.
            pub output_column_name: String,
            /// Separator to split multiple columns by.
            pub separator: Option<char>,
            /// Whether to compute maximum (is minimum if separator is some).
            pub is_max: bool,
            /// The collected value, if any.
            pub value: Option<serde_json::Value>,
        }

        impl SingleValueCollector {
            /// Construct given column name.
            pub fn new(
                input_column_name: &str,
                output_column_name: &str,
                separator: Option<char>,
                is_max: Option<bool>,
            ) -> Self {
                Self {
                    input_column_name: input_column_name.into(),
                    output_column_name: output_column_name.into(),
                    separator,
                    is_max: is_max.unwrap_or(true),
                    value: None,
                }
            }
        }

        impl Collector for SingleValueCollector {
            fn register(&mut self, column_name: &str, value: &serde_json::Value) {
                if column_name == self.input_column_name && !value.is_null() {
                    if let Some(separator) = self.separator {
                        if let serde_json::Value::String(value_str) = value {
                            // Split string value if we have a separator and a string.
                            let values_f64 = value_str
                                .split(separator)
                                .collect::<Vec<_>>()
                                .into_iter()
                                .flat_map(|s| s.parse::<f64>().ok())
                                .collect::<Vec<_>>();
                            let f64_value = if self.is_max {
                                values_f64.into_iter().max_by(|a, b| a.total_cmp(b))
                            } else {
                                values_f64.into_iter().min_by(|a, b| a.total_cmp(b))
                            };
                            // Write out the value if we could convert into float, silently ignore the `.`
                            // entries (and everything else from dbNSFP) that could not be converted.
                            if let Some(f64_value) = f64_value {
                                if let Some(f64_value) = serde_json::Number::from_f64(f64_value) {
                                    self.value = Some(serde_json::Value::Number(f64_value));
                                }
                            }
                        }
                    } else {
                        self.value = Some(value.clone());
                    }
                }
            }

            fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>) {
                if let Some(value) = self.value.as_ref() {
                    dict.insert(self.output_column_name.clone(), value.clone());
                }
            }
        }

        /// Simple implementation for collecting aggregated min/max score.
        #[derive(Debug, Clone)]
        pub struct ExtremalValueCollector {
            /// Whether is max (else is min).
            pub is_max: bool,
            /// The output column name.
            pub output_column_name: String,
            /// The column names to collect.
            pub column_names: Vec<String>,
            /// The collected values, if any.
            pub values: indexmap::IndexMap<String, serde_json::Value>,
        }

        impl ExtremalValueCollector {
            /// Construct given column name.
            pub fn new(column_names: &[&str], output_column_name: &str, is_max: bool) -> Self {
                Self {
                    is_max,
                    output_column_name: output_column_name.to_string(),
                    column_names: column_names.iter().map(|s| s.to_string()).collect(),
                    values: Default::default(),
                }
            }
        }

        impl Collector for ExtremalValueCollector {
            fn register(&mut self, column_name: &str, value: &serde_json::Value) {
                if self.column_names.iter().any(|s| s.as_str() == column_name) && !value.is_null() {
                    self.values.insert(column_name.to_string(), value.clone());
                }
            }

            fn write_to(&self, dict: &mut indexmap::IndexMap<String, serde_json::Value>) {
                let mut sel_name = None;
                let mut sel_value = serde_json::Value::Null;
                for (name, value) in self.values.iter() {
                    if value.is_null() || !value.is_number() {
                        // can only select numeric values that are not null
                        continue;
                    }

                    if sel_value.is_null()
                        || (self.is_max
                            && value.as_f64().unwrap_or_default()
                                > sel_value.as_f64().unwrap_or_default())
                        || (!self.is_max
                            && value.as_f64().unwrap_or_default()
                                < sel_value.as_f64().unwrap_or_default())
                    {
                        sel_name = Some(name.clone());
                        sel_value = value.clone();
                    }
                }

                if let Some(sel_name) = sel_name {
                    dict.insert(self.output_column_name.clone(), sel_value);
                    dict.insert(
                        format!("{}_argmax", self.output_column_name),
                        serde_json::Value::String(sel_name),
                    );
                }
            }
        }
    }

    pub(crate) fn with_seqvar_and_annotator(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<pbs_output::VariantRelatedAnnotation, anyhow::Error> {
        Ok(pbs_output::VariantRelatedAnnotation {
            dbids: dbids(seqvar, annotator)?,
            frequency: frequency(seqvar),
            clinvar: clinvar(seqvar, annotator)?,
            scores: scores(seqvar, annotator)?,
        })
    }

    fn dbids(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Option<pbs_output::DbIds>, anyhow::Error> {
        let dbsnp_id = annotator
            .query_dbsnp(seqvar)
            .map_err(|e| anyhow::anyhow!("problem querying dbSNP: {}", e))?
            .map(|record| format!("rs{}", record.rs_id));
        if let Some(dbsnp_id) = dbsnp_id {
            Ok(Some(pbs_output::DbIds {
                dbsnp_id: Some(dbsnp_id),
            }))
        } else {
            Ok(None)
        }
    }

    fn frequency(seqvar: &VariantRecord) -> Option<pbs_output::FrequencyAnnotation> {
        Some(pbs_output::FrequencyAnnotation {
            gnomad_exomes: Some(pbs_output::NuclearFrequency {
                an: seqvar.population_frequencies.gnomad_exomes.an,
                het: seqvar.population_frequencies.gnomad_exomes.het,
                homalt: seqvar.population_frequencies.gnomad_exomes.hom,
                hemialt: seqvar.population_frequencies.gnomad_exomes.hemi,
                af: seqvar.population_frequencies.gnomad_exomes.af(),
            }),
            gnomad_genomes: Some(pbs_output::NuclearFrequency {
                an: seqvar.population_frequencies.gnomad_genomes.an,
                het: seqvar.population_frequencies.gnomad_genomes.het,
                homalt: seqvar.population_frequencies.gnomad_genomes.hom,
                hemialt: seqvar.population_frequencies.gnomad_genomes.hemi,
                af: seqvar.population_frequencies.gnomad_genomes.af(),
            }),
            gnomad_mtdna: Some(pbs_output::GnomadMitochondrialFrequency {
                an: seqvar.population_frequencies.gnomad_mtdna.an,
                het: seqvar.population_frequencies.gnomad_mtdna.het,
                homalt: seqvar.population_frequencies.gnomad_mtdna.hom,
                af: seqvar.population_frequencies.gnomad_mtdna.af(),
            }),
            helixmtdb: Some(pbs_output::HelixMtDbFrequency {
                an: seqvar.population_frequencies.helixmtdb.an,
                het: seqvar.population_frequencies.helixmtdb.het,
                homalt: seqvar.population_frequencies.helixmtdb.hom,
                af: seqvar.population_frequencies.helixmtdb.af(),
            }),
            inhouse: Some(pbs_output::NuclearFrequency {
                an: seqvar.population_frequencies.inhouse.an,
                het: seqvar.population_frequencies.inhouse.het,
                homalt: seqvar.population_frequencies.inhouse.hom,
                hemialt: seqvar.population_frequencies.inhouse.hemi,
                af: seqvar.population_frequencies.inhouse.af(),
            }),
        })
    }

    fn clinvar(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Option<pbs_output::ClinvarAnnotation>, anyhow::Error> {
        let record = annotator
            .query_clinvar_minimal(seqvar)
            .map_err(|e| anyhow::anyhow!("problem querying clinvar-minimal: {}", e))?;
        if let Some(record) = record.as_ref() {
            if record.records.is_empty() {
                tracing::error!(
                    "variant {:?} found with empty list in ClinVar (should not happen)",
                    seqvar
                );
                return Ok(None);
            } else if record.records.len() > 1 {
                tracing::warn!(
                    "variant {:?} found list with {} entries, using first",
                    seqvar,
                    record.records.len()
                );
            }
            let vcv_record = &record.records[0];
            let accession = vcv_record.accession.as_ref().expect("no accession?");
            let vcv_accession = format!("{}.{}", &accession.accession, accession.version);

            if let Some(agc) = vcv_record
                .classifications
                .as_ref()
                .and_then(|c| c.germline_classification.as_ref())
            {
                let germline_significance_description = if let Some(description) =
                    agc.description.as_ref()
                {
                    description.clone()
                } else {
                    tracing::error!("variant {:?} has germline classification without description (should not happen)", &seqvar);
                    return Ok(None);
                };
                let germline_review_status = agc.review_status;
                // TODO: search through all submitted records and pick the most significant one.
                let effective_germline_significance_description =
                    germline_significance_description.clone();

                Ok(Some(pbs_output::ClinvarAnnotation {
                    vcv_accession,
                    germline_significance_description,
                    germline_review_status,
                    effective_germline_significance_description,
                }))
            } else {
                tracing::trace!(
                    "variant {:?} has no germline classification (likely somatic only)",
                    &seqvar
                );
                Ok(None)
            }
        } else {
            Ok(None)
        }
    }

    /// Return information about the scores entries.
    pub(crate) fn score_columns() -> Vec<pbs_output::VariantScoreColumn> {
        vec![
            // Scores obtained from CADD file.
            pbs_output::VariantScoreColumn {
                name: "cadd_phred".to_string(),
                label: "CADD".to_string(),
                description: "PHRED-scaled CADD score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "mmsplice".to_string(),
                label: "MMSplice".to_string(),
                description: "MMSplice score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "mmsplice_argmax".to_string(),
                label: "MMSplice (which)".to_string(),
                description: "Which MMSplice score is maximal".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "polyphen".to_string(),
                label: "PolyPhen".to_string(),
                description: "PolyPhen score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "sift".to_string(),
                label: "SIFT".to_string(),
                description: "SIFT score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "spliceai".to_string(),
                label: "SpliceAI".to_string(),
                description: "SpliceAI score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "spliceai_argmax".to_string(),
                label: "SpliceAI (which)".to_string(),
                description: "Which SpliceAI score is maximal".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            // Scores obtained from dbNSFP file.
            pbs_output::VariantScoreColumn {
                name: "alphamissense".to_string(),
                label: "AlphaMissense".to_string(),
                description: "AlphaMissense score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "bayesdel_addaf".to_string(),
                label: "BayesDel".to_string(),
                description: "BayesDel AddAF score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "fathmm".to_string(),
                label: "FATHMM".to_string(),
                description: "FATHMM score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "fitcons_integrated".to_string(),
                label: "fitCons".to_string(),
                description: "The integrated fitCons score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "lrt".to_string(),
                label: "LRT".to_string(),
                description: "LRT score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "metasvm".to_string(),
                label: "MetaSVM".to_string(),
                description: "MetaSVM score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "polyphen2_hdiv".to_string(),
                label: "Polyphen2 HDIV".to_string(),
                description: "Polyphen2 HDIV score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "polyphen2_hvar".to_string(),
                label: "Polyphen2 HVAR".to_string(),
                description: "Polyphen2 HVAR score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "primateai".to_string(),
                label: "PrimateAI".to_string(),
                description: "PrimateAI score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "provean".to_string(),
                label: "PROVEAN".to_string(),
                description: "PROVEAN score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
            pbs_output::VariantScoreColumn {
                name: "revel".to_string(),
                label: "REVEL".to_string(),
                description: "REVEL score".to_string(),
                r#type: pbs_output::VariantScoreColumnType::Number as i32,
            },
        ]
    }

    /// Query precomputed scores for `seqvar` from annonars `annotator`.
    pub fn scores(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Option<pbs_output::ScoreAnnotations>, anyhow::Error> {
        use score_collection::*;
        let mut result = indexmap::IndexMap::new();

        // Extract values from CADD.
        if let Some(cadd_values) = annotator
            .query_cadd(seqvar)
            .as_ref()
            .map_err(|e| anyhow::anyhow!("problem querying CADD: {}", e))?
        {
            let mut collectors: Vec<Box<dyn Collector>> = vec![
                Box::new(SingleValueCollector::new("PHRED", "cadd_phred", None, None)),
                Box::new(SingleValueCollector::new("SIFTval", "sift", None, None)),
                Box::new(SingleValueCollector::new(
                    "PolyPhenVal",
                    "polyphen",
                    None,
                    None,
                )),
                Box::new(ExtremalValueCollector::new(
                    &[
                        "SpliceAI-acc-gain",
                        "SpliceAI-acc-loss",
                        "SpliceAI-don-gain",
                        "SpliceAI-don-loss",
                    ],
                    "spliceai",
                    true,
                )),
                Box::new(ExtremalValueCollector::new(
                    &[
                        "MMSp_acceptorIntron",
                        "MMSp_acceptor",
                        "MMSp_exon",
                        "MMSp_donor",
                        "MMSp_donorIntron",
                    ],
                    "mmsplice",
                    true,
                )),
            ];

            for (column, value) in annotator
                .annonars_dbs
                .cadd_ctx
                .schema
                .columns
                .iter()
                .zip(cadd_values.iter())
            {
                for collector in collectors.iter_mut() {
                    collector.register(column.name.as_str(), value);
                }
            }

            collectors.iter_mut().for_each(|collector| {
                collector.write_to(&mut result);
            })
        }

        // Extract values from dbNSFP

        if let Some(dbnsfp_values) = annotator
            .query_dbnsfp(seqvar)
            .as_ref()
            .map_err(|e| anyhow::anyhow!("problem querying dbNSFP: {}", e))?
        {
            let mut collectors: Vec<Box<dyn Collector>> = vec![
                Box::new(SingleValueCollector::new(
                    "AlphaMissense_score",
                    "alphamissense",
                    Some(';'),
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "BayesDel_addAF_score",
                    "bayesdel_addaf",
                    None,
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "FATHMM_score",
                    "fathmm",
                    Some(';'),
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "integrated_fitCons_score",
                    "fitcons_integrated",
                    None,
                    None,
                )),
                Box::new(SingleValueCollector::new("LRT_score", "lrt", None, None)),
                Box::new(SingleValueCollector::new(
                    "MetaSVM_score",
                    "metasvm",
                    None,
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "Polyphen2_HDIV_score",
                    "polyphen2_hdiv",
                    Some(';'),
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "Polyphen2_HVAR_score",
                    "polyphen2_hvar",
                    Some(';'),
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "PrimateAI_score",
                    "primateai",
                    None,
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "PROVEAN_score",
                    "provean",
                    Some(';'),
                    None,
                )),
                Box::new(SingleValueCollector::new(
                    "REVEL_score",
                    "revel",
                    Some(';'),
                    None,
                )),
            ];

            for (column, value) in annotator
                .annonars_dbs
                .dbnsfp_ctx
                .schema
                .columns
                .iter()
                .zip(dbnsfp_values.iter())
            {
                for collector in collectors.iter_mut() {
                    collector.register(column.name.as_str(), value);
                }
            }

            collectors.iter_mut().for_each(|collector| {
                collector.write_to(&mut result);
            })
        }

        Ok(Some(pbs_output::ScoreAnnotations {
            entries: result
                .into_iter()
                .map(
                    |(key, value)| -> Result<pbs_output::ScoreEntry, anyhow::Error> {
                        Ok(pbs_output::ScoreEntry {
                            key,
                            value: serde_json::from_value(value)
                                .map_err(|e| anyhow::anyhow!("could not convert value: {}", e))?,
                        })
                    },
                )
                .collect::<Result<Vec<_>, _>>()?,
        }))
    }
}

impl WithSeqvarAndAnnotator for pbs_output::VariantRelatedAnnotation {
    type Error = anyhow::Error;

    fn with_seqvar_and_annotator(
        seqvar: &VariantRecord,
        annotator: &Annotator,
    ) -> Result<Self, Self::Error> {
        variant_related_annotation::with_seqvar_and_annotator(seqvar, annotator)
    }
}

impl WithSeqvarAndAnnotator for pbs_output::CallRelatedAnnotation {
    type Error = anyhow::Error;

    fn with_seqvar_and_annotator(
        seqvar: &VariantRecord,
        _annotator: &Annotator,
    ) -> Result<Self, Self::Error> {
        Ok(Self {
            call_infos: seqvar
                .call_infos
                .iter()
                .map(|(sample, call_info)| pbs_output::SampleCallInfo {
                    sample: sample.clone(),
                    genotype: call_info.genotype.clone(),
                    dp: call_info.dp,
                    ad: call_info.ad,
                    gq: call_info.gq,
                    ps: call_info.ps,
                })
                .collect(),
        })
    }
}

/// Create output payload and write the record to the output file.
async fn create_and_write_record(
    seqvar: VariantRecord,
    annotator: &Annotator,
    chrom_to_chrom_no: &std::collections::HashMap<String, u32>,
    writer: &mut tokio::io::BufWriter<tokio::fs::File>,
    args: &Args,
    rng: &mut rand::rngs::StdRng,
    uuid_buf: &mut [u8; 16],
) -> Result<(), anyhow::Error> {
    // Build the output record protobuf.
    let record = pbs_output::OutputRecord {
        uuid: Uuid::from_bytes({
            rng.fill_bytes(uuid_buf);
            *uuid_buf
        })
        .to_string(),
        case_uuid: args.case_uuid.unwrap_or_default().to_string(),
        vcf_variant: Some(pbs_output::VcfVariant {
            genome_release: Into::<pbs_output::GenomeRelease>::into(args.genome_release) as i32,
            chrom: seqvar.vcf_variant.chrom.clone(),
            chrom_no: chrom_to_chrom_no
                .get(&seqvar.vcf_variant.chrom)
                .cloned()
                .unwrap_or_default() as i32,
            pos: seqvar.vcf_variant.pos,
            ref_allele: seqvar.vcf_variant.ref_allele.clone(),
            alt_allele: seqvar.vcf_variant.alt_allele.clone(),
        }),
        variant_annotation: Some(pbs_output::VariantAnnotation {
            gene: Some(
                pbs_output::GeneRelatedAnnotation::with_seqvar_and_annotator(&seqvar, annotator)
                    .map_err(|e| {
                        anyhow::anyhow!("problem creating gene-related annotation: {}", e)
                    })?,
            ),
            variant: Some(
                pbs_output::VariantRelatedAnnotation::with_seqvar_and_annotator(&seqvar, annotator)
                    .map_err(|e| {
                        anyhow::anyhow!("problem creating variant-related annotation: {}", e)
                    })?,
            ),
            call: Some(
                pbs_output::CallRelatedAnnotation::with_seqvar_and_annotator(&seqvar, annotator)
                    .map_err(|e| {
                        anyhow::anyhow!("problem creating call-related annotation: {}", e)
                    })?,
            ),
        }),
    };

    // Write out the record to JSONL.

    let mut buf = Vec::<u8>::new();
    writeln!(
        &mut buf,
        "{}",
        serde_json::to_string(&record)
            .map_err(|e| anyhow::anyhow!("could not convert record to JSON: {}", e))?
    )?;
    writer
        .write_all(&buf)
        .await
        .map_err(|e| anyhow::anyhow!("could not write record to output file: {}", e))
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

    tracing::info!("Loading query... {}", args.path_query_json);
    let pb_query: pbs_query::CaseQuery =
        serde_json::from_reader(std::fs::File::open(&args.path_query_json)?)?;
    let query = CaseQuery::try_from(pb_query.clone())?;

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
    let hgnc_allowlist =
        crate::strucvars::query::translate_genes(&query.locus.genes, &in_memory_dbs);

    tracing::info!("Running queries...");
    let before_query = Instant::now();
    let query_stats = run_query(
        &interpreter::QueryInterpreter::new(query, hgnc_allowlist),
        &pb_query.clone(),
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
    for (effect, count) in query_stats.passed_by_consequences.iter() {
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

    use super::schema::data::{CallInfo, VariantRecord};
    use crate::seqvars::query::schema::query::{CaseQuery, GenotypeChoice, RecessiveMode};

    #[rstest]
    #[case::comphet_het_het_ref_fails(
        RecessiveMode::CompoundHeterozygous,
        vec!["0/1,0/1,0/0"],
        false,
    )]
    #[case::comphet_hom_ref_ref_fails(
        RecessiveMode::CompoundHeterozygous,
        vec!["1/1,0/0,0/0"],
        false,
    )]
    #[case::comphet_het_het_hom_and_het_hom_het_passes(
        RecessiveMode::CompoundHeterozygous,
        vec!["0/1,0/1,0/0","0/1,0/0,0/1"],
        true,
    )]
    #[case::any_hom_ref_ref_and_het_ref_het_fails(
        RecessiveMode::Any,
        vec!["1/1,0/0,0/0","0/1,0/0,0/1"],
        false,
    )]
    #[case::any_hom_het_het_and_het_ref_het_passes(
        RecessiveMode::Any,
        vec!["1/1,0/1,0/1","0/1,0/0,0/1"],
        true,
    )]
    #[case::any_het_het_ref_and_het_ref_het_passes(
        RecessiveMode::Any,
        vec!["1/0,1/0,0/0","0/1,0/0,0/1"],
        true,
    )]
    #[case::any_het_het_ref_and_het_het_ref_fails(
        RecessiveMode::Any,
        vec!["1/0,1/0,0/0","1/0,0/1,0/0"],
        false,
    )]
    fn passes_for_gene_full_trio(
        #[case] recessive_mode: RecessiveMode,
        #[case] trio_gts: Vec<&str>,
        #[case] passes: bool,
    ) -> Result<(), anyhow::Error> {
        use crate::seqvars::query::schema::query::{QuerySettingsGenotype, SampleGenotypeChoice};

        let query = CaseQuery {
            genotype: QuerySettingsGenotype {
                recessive_mode,
                sample_genotypes: indexmap::indexmap! {
                    String::from("index") => SampleGenotypeChoice { sample: String::from("index"), genotype: GenotypeChoice::RecessiveIndex, ..Default::default() },
                    String::from("father") => SampleGenotypeChoice { sample: String::from("father"), genotype: GenotypeChoice::RecessiveFather, ..Default::default() },
                    String::from("mother") => SampleGenotypeChoice { sample: String::from("mother"), genotype: GenotypeChoice::RecessiveMother, ..Default::default() },
                },
            },
            ..Default::default()
        };
        let seqvars = trio_gts
            .iter()
            .map(|gts| {
                let gts: Vec<&str> = gts.split(',').collect();
                VariantRecord {
                    call_infos: indexmap::indexmap! {
                        String::from("index") =>
                            CallInfo {
                                sample: String::from("index"),
                                genotype: Some(gts[0].into()),
                                ..Default::default()
                            },
                        String::from("father") =>
                            CallInfo {
                                genotype: Some(gts[1].into()),
                                ..Default::default()
                            },
                        String::from("mother") =>
                            CallInfo {
                                genotype: Some(gts[2].into()),
                                ..Default::default()
                            },
                    },
                    ..Default::default()
                }
            })
            .collect::<Vec<_>>();

        assert_eq!(super::passes_for_gene(&query, &seqvars)?, passes);

        Ok(())
    }

    // TODO: re-enable smoke test
    // #[tracing_test::traced_test]
    // #[rstest::rstest]
    // #[case::case_1_ingested_vcf("tests/seqvars/query/Case_1.ingested.vcf")]
    // #[case::dragen_ingested_vcf("tests/seqvars/query/dragen.ingested.vcf")]
    // #[tokio::test]
    // async fn smoke_test(#[case] path_input: &str) -> Result<(), anyhow::Error> {
    //     mehari::common::set_snapshot_suffix!("{}", path_input.split('/').last().unwrap());

    //     let tmpdir = temp_testdir::TempDir::default();
    //     let path_output = format!("{}/out.tsv", tmpdir.to_string_lossy());
    //     let path_input: String = path_input.into();
    //     let path_query_json = path_input.replace(".ingested.vcf", ".query.json");

    //     let args_common = Default::default();
    //     let args = super::Args {
    //         genome_release: crate::common::GenomeRelease::Grch37,
    //         path_db: "tests/seqvars/query/db".into(),
    //         path_query_json,
    //         path_input,
    //         path_output,
    //         max_results: None,
    //         rng_seed: Some(42),
    //         max_tad_distance: 10_000,
    //         result_set_id: None,
    //         case_uuid_id: None,
    //     };
    //     super::run(&args_common, &args).await?;

    //     insta::assert_snapshot!(std::fs::read_to_string(args.path_output.as_str())?);

    //     Ok(())
    // }
}
