//! Code implementing the "seqvars query" sub command.

mod annonars;
mod interpreter;
mod schema;
mod sorting;

use std::io::{BufRead, Write};
use std::time::Instant;

use clap::{command, Parser};
use ext_sort::LimitedBufferBuilder;
use ext_sort::{ExternalSorter, ExternalSorterBuilder};
use itertools::Itertools;
use noodles_vcf as vcf;

use mehari::{annotate::seqvars::CHROM_TO_CHROM_NO, common::open_read_maybe_gz};
use rand_core::{RngCore, SeedableRng};
use thousands::Separable;
use uuid::Uuid;

use crate::common;
use crate::seqvars::query::schema::GenotypeChoice;
use crate::{common::trace_rss_now, common::GenomeRelease};

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

/// Gene-related information for the gene.
pub mod gene_related {
    use mehari::annotate::seqvars::ann::Consequence;

    use super::schema::SequenceVariant;

    /// Gene-related information for a `ResultPayload`.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Record {
        /// Gene identity related (for display of gene symbol).
        pub identity: Identity,
        /// Gene-related consequences of a variant.
        pub consequences: Consequences,
        /// Phenotype-related information, if any.
        pub phenotype: Option<Phenotype>,
        /// Gene-wise constraints on the gene, if any.
        pub constraints: Option<Constraints>,
    }

    impl Record {
        /// Construct given a `SequenceVariant` if the information is given in the annotation.
        ///
        /// Note that we will only look at the first annotation record as the ingest creates
        /// one `SequenceVariant` record per gene.
        ///
        /// # Error
        ///
        /// Returns an error if `seqvar` does not contain all necessary information.
        pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Option<Self>, anyhow::Error> {
            if let Some(ann) = seqvar.ann_fields.first() {
                if !ann.gene_id.is_empty() && !ann.gene_symbol.is_empty() {
                    return Ok(Some(Self {
                        identity: Identity::new(ann.gene_id.clone(), ann.gene_symbol.clone()),
                        consequences: Consequences::new(
                            ann.hgvs_t
                                .clone()
                                .ok_or_else(|| anyhow::anyhow!("missing hgvs_t annotation"))?,
                            ann.hgvs_p.clone(),
                            ann.consequences.clone(),
                        ),
                        // TODO: phenotype
                        // TODO: constraints
                        ..Default::default()
                    }));
                }
            }
            Ok(None)
        }
    }

    /// Result information for gene identity.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Identity {
        /// HGNC gene ID.
        pub hgnc_id: String,
        /// HGNC gene symbol.
        pub hgnc_symbol: String,
    }

    /// Result information related to phenotype / disease.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Phenotype {
        /// Whether the gene is a known disease gene.
        pub is_disease_gene: bool,
        // TODO: modes of inheritance
        // TODO: disease/phenotype terms?
    }

    /// Consequences related to a gene.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Consequences {
        /// HGVS.{c,n} code of variant
        pub hgvs_t: String,
        /// HGVS.p code of variant
        pub hgvs_p: Option<String>,

        /// The predicted variant consequences.
        pub consequences: Vec<Consequence>,
    }

    /// Result gene constraint information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Constraints {
        /// gnomAD loeuf score
        pub gnomad_loeuf: f32,
        /// gnomAD mis_z score
        pub gnomad_mis_z: f32,
        /// gnomAD oe_lof score
        pub gnomad_oe_lof: f32,
        /// gnomAD oe_lof_lower score
        pub gnomad_oe_lof_lower: f32,
        /// gnomAD oe_lof_upper score
        pub gnomad_oe_lof_upper: f32,
        /// gnomAD oe_mis score
        pub gnomad_oe_mis: f32,
        /// gnomAD oe_mis_lower score
        pub gnomad_oe_mis_lower: f32,
        /// gnomAD oe_mis_upper score
        pub gnomad_oe_mis_upper: f32,
        /// gnomAD pLI score
        pub gnomad_pli: f32,
        /// gnomAD syn_z score
        pub gnomad_syn_z: f32,
    }
}

/// Variant-related information.
pub mod variant_related {
    use super::{annonars::Annotator, schema::SequenceVariant};

    /// Record for variant-related scores.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Record {
        /// Precomputed scores.
        pub precomputed_scores: indexmap::IndexMap<String, f32>,
        /// Database identifiers.
        pub db_ids: DbIds,
        /// Clinvar information.
        pub clinvar: Option<Clinvar>,
        /// Frequency information.
        pub frequency: Frequency,
    }

    impl Record {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Self, anyhow::Error> {
            Ok(Self {
                precomputed_scores: Self::query_precomputed_scores(
                    seqvar, annotator,
                )?,
                db_ids: DbIds::with_seqvar_and_annotator(seqvar, annotator)?,
                clinvar: Clinvar::with_seqvar_and_annotator(seqvar, annotator)?,
                frequency: Frequency::with_seqvar(seqvar)?,
            })
        }

        /// Query precomputed scores for `seqvar` from annonars `annotator`.
        pub fn query_precomputed_scores(
            _seqvar: &SequenceVariant,
            _annotator: &Annotator,
        ) -> Result<indexmap::IndexMap<String, f32>, anyhow::Error> {
            // TODO: implement me!
            Ok(indexmap::IndexMap::default())
        }
    }

    /// Database identifiers.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct DbIds {
        /// The variant's dbSNP identifier present.
        pub dbsnp_rs: Option<String>,
    }

    impl DbIds {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Self, anyhow::Error> {
            // TODO: need to properly separate error from no result in annonars.
            Ok(Self {
                dbsnp_rs: annotator
                    .query_dbsnp(seqvar)
                    .ok()
                    .map(|record| format!("rs{}", record.rs_id)),
            })
        }
    }

    /// ClinVar-related information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Clinvar {
        // TODO: switch to VCV summary or hierarchical?
        /// The RCV accession.
        pub rcv: String,
        /// The clinical significance.
        pub significance: String,
        /// The review status.
        pub review_status: String,
    }

    impl Clinvar {
        /// Construct given sequence variant and annonars annotator.
        pub fn with_seqvar_and_annotator(
            seqvar: &SequenceVariant,
            annotator: &Annotator,
        ) -> Result<Option<Self>, anyhow::Error> {
            // TODO: need to properly separate error from no result in annonars.
            annotator.query_clinvar_minimal(seqvar).ok().map(|record| {
                let annonars::clinvar_minimal::pbs::Record {
                    rcv,
                    clinical_significance,
                    review_status,
                    ..
                } = record;
                Ok(Self {
                    rcv,
                    significance: match clinical_significance {
                        0 => "Pathogenic",
                        1 => "Likely pathogenic",
                        2 => "Uncertain significance",
                        3 => "Likely benign",
                        4 => "Benign",
                        _ => anyhow::bail!("invalid clinical significance enum: {}", clinical_significance)
                    }.into(),
                    review_status: match review_status {
                        0 => "no assertion provided",
                        1 => "no assertion criteria provided",
                        2 => "criteria provided, conflicting interpretations",
                        3 => "criteria provided, single submitter",
                        4 => "criteria provided, multiple submitters, no conflicts",
                        5 => "reviewed by expert panel",
                        6 => "practice guideline",
                        _ => anyhow::bail!("invalid review status enum: {}", review_status)
                    }.into(),
                })
            }).transpose()
        }
    }

    /// Frequency information.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Frequency {
        /// gnomAD-genomes frequency
        pub gnomad_genomes: Option<NuclearFrequency>,
        /// gnomAD-exomes frequency
        pub gnomad_exomes: Option<NuclearFrequency>,
        /// gnomad-mtDNA frequency
        pub gnomad_mtdna: Option<MtdnaFrequency>,
        /// HelixMtDb frequency
        pub helixmtdb: Option<MtdnaFrequency>,
        /// in-house frequency
        pub inhouse: Option<NuclearFrequency>,
    }

    impl Frequency {
        /// Extract frequency information from `seqvar`
        pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Frequency, anyhow::Error> {
            let chrom = annonars::common::cli::canonicalize(&seqvar.chrom);
            let frequency = if chrom == "MT" {
                FrequencyBuilder::default()
                    .gnomad_genomes(Some(NuclearFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                        seqvar.gnomad_genomes_hemi,
                    )))
                    .gnomad_exomes(Some(NuclearFrequency::new(
                        seqvar.gnomad_exomes_af(),
                        seqvar.gnomad_exomes_an,
                        seqvar.gnomad_exomes_het,
                        seqvar.gnomad_exomes_hom,
                        seqvar.gnomad_exomes_hemi,
                    )))
                    .build()
            } else {
                FrequencyBuilder::default()
                    .gnomad_mtdna(Some(MtdnaFrequency::new(
                        seqvar.gnomad_genomes_af(),
                        seqvar.gnomad_genomes_an,
                        seqvar.gnomad_genomes_het,
                        seqvar.gnomad_genomes_hom,
                    )))
                    .helixmtdb(Some(MtdnaFrequency::new(
                        seqvar.helixmtdb_af(),
                        seqvar.helix_an,
                        seqvar.helix_het,
                        seqvar.helix_hom,
                    )))
                    .build()
            }
            .map_err(|e| anyhow::anyhow!("could not build frequency information: {}", e))?;
            Ok(frequency)
        }
    }

    /// Nuclear frequency information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct NuclearFrequency {
        /// Overall allele frequency.
        pub allele_freq: f32,
        /// Number of alleles.
        pub allele_count: i32,
        /// Number of heterozygous carriers.
        pub het_carriers: i32,
        /// Number of homozygous carriers.
        pub hom_carriers: i32,
        /// Number of hemizygous carriers.
        pub hemi_carriers: i32,
    }

    /// Mitochondrial frequency information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct MtdnaFrequency {
        /// Overall allele frequency.
        pub allele_freq: f32,
        /// Number of alleles.
        pub allele_count: i32,
        /// Number of heterplasmic carriers.
        pub het_carriers: i32,
        /// Number of homoplasmic carriers.
        pub hom_carriers: i32,
    }
}

/// Call-related information.
pub mod call_related {
    use super::schema::SequenceVariant;

    /// Call-related record.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct Record {
        /// The genotype information for each sample.
        pub call_info: indexmap::IndexMap<String, CallInfo>,
    }

    impl Record {
        /// Construct a new `Record` from a `SequenceVariant`.
        ///
        /// # Error
        ///
        /// Returns an error if the `SequenceVariant` does not contain all necessary information.
        pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Self, anyhow::Error> {
            Ok(Self {
                call_info: seqvar
                    .call_info
                    .iter()
                    .map(|(sample_name, call_info)| {
                        (
                            sample_name.clone(),
                            CallInfo::new(
                                call_info.dp,
                                call_info.ad,
                                call_info.quality.map(|q| q as i32),
                                call_info.genotype.clone(),
                            ),
                        )
                    })
                    .collect(),
            })
        }
    }

    /// Genotype information.
    #[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
    pub struct CallInfo {
        /// Depth of coverage.
        pub dp: Option<i32>,
        /// Alternate read depth.
        pub ad: Option<i32>,
        /// Genotype quality.
        pub gq: Option<i32>,
        /// Genotype.
        pub gt: Option<String>,
    }
}

/// Result record information.
pub mod result {
    /// A result record from the query.
    ///
    /// These records are written to TSV for import into the database.   They contain the
    /// bare necessary information for sorting etc. in the database.  The actual main
    /// data is in the payload, serialized as JSON.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Record {
        /// UUID for the record.
        pub sodar_uuid: uuid::Uuid,
        /// Genome release for the coordinate.
        pub release: String,
        /// Chromosome name.
        pub chromosome: String,
        /// Chromosome number.
        pub chromosome_no: i32,
        /// Reference allele sequence.
        pub reference: String,
        /// Alternative allele sequence.
        pub alternative: String,
        /// UCSC bin of the record.
        pub bin: u32,
        /// Start position of the record.
        pub start: i32,
        /// End position of the record.
        pub end: i32,
        /// The result set ID as specified on the command line.
        pub smallvariantqueryresultset_id: String,
        /// The JSON-serialized `ResultPayload`.
        pub payload: String,
    }

    /// The structured result information of the result record.
    #[derive(
        Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder,
    )]
    pub struct Payload {
        /// Case UUID as specified on the command line.
        pub case_uuid: uuid::Uuid,
        /// The affected gene and consequence, if any.
        pub gene_related: Option<super::gene_related::Record>,
        /// Variant-related information, always present.
        pub variant_related: super::variant_related::Record,
        /// Genotypes call related, always present.
        pub call_related: super::call_related::Record,
    }
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
            &Some(GenotypeChoice::ComphetIndex) => {
                index_name = sample_name.clone();
                mode = Mode::ComphetRecessive;
            }
            &Some(GenotypeChoice::RecessiveIndex) => {
                index_name = sample_name.clone();
                mode = Mode::Recessive;
            }
            &Some(GenotypeChoice::RecessiveParent) => {
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
fn run_query(
    interpreter: &interpreter::QueryInterpreter,
    args: &Args,
    annotator: &annonars::Annotator,
    rng: &mut rand::rngs::StdRng,
) -> Result<QueryStats, anyhow::Error> {
    let tmp_dir = tempdir::TempDir::new("vfw")?;

    let chrom_to_chrom_no = &CHROM_TO_CHROM_NO;
    let mut stats = QueryStats::default();

    // Buffer for generating UUIDs.
    let mut uuid_buf = [0u8; 16];

    // Open VCF file, create reader, and read header.
    let mut input_reader = open_read_maybe_gz(&args.path_input).map_err(|e| {
        anyhow::anyhow!("could not open file {} for reading: {}", args.path_input, e)
    })?;
    let mut input_reader = vcf::Reader::new(&mut input_reader);
    let input_header = input_reader.read_header()?;

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

        for input_record in input_reader.records(&input_header) {
            stats.count_total += 1;
            let record_seqvar = SequenceVariant::from_vcf(&input_record?, &input_header)
                .map_err(|e| anyhow::anyhow!("could not parse VCF record: {}", e))?;
            tracing::debug!("processing record {:?}", record_seqvar);

            if interpreter.passes(&record_seqvar, &annotator)?.pass_all {
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
        tmp_unsorted.flush().map_err(|e| {
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
    let result_payload = result::PayloadBuilder::default()
        .case_uuid(args.case_uuid_id.unwrap_or_default().clone())
        .gene_related(
            gene_related::Record::with_seqvar(&seqvar)
                .map_err(|e| anyhow::anyhow!("problem creating gene-related payload: {}", e))?,
        )
        .variant_related(
            variant_related::Record::with_seqvar_and_annotator(&seqvar, annotator)
                .map_err(|e| anyhow::anyhow!("problem creating variant-related payload: {}", e))?,
        )
        .call_related(
            call_related::Record::with_seqvar(&seqvar)
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
            &result::RecordBuilder::default()
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
    )?;
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
                let gts: Vec<&str> = gts.split(",").map(|gt| gt).collect();
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
            case_uuid_id: None,
        };
        super::run(&args_common, &args)?;

        insta::assert_snapshot!(std::fs::read_to_string(args.path_output.as_str())?);

        Ok(())
    }
}
