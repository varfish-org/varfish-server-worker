//! Code implementing the "sv query" sub command.

pub mod bgdbs;
pub mod clinvar;
pub mod genes;
pub mod interpreter;
pub mod masked;
pub mod pathogenic;
pub mod records;
pub mod schema;
pub mod tads;

use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::File,
    time::Instant,
};

use clap::{command, Parser};
use csv::QuoteStyle;
use hgvs::static_data::ASSEMBLY_INFOS;
use indexmap::IndexMap;
use log::warn;
use mehari::{
    annotate::seqvars::provider::TxIntervalTrees,
    db::create::txs::data::{Strand, Transcript, TxSeqDatabase},
};
use serde::Serialize;
use thousands::Separable;
use uuid::Uuid;

use crate::{
    common::{build_chrom_map, numeric_gene_id, open_read_maybe_gz, trace_rss_now},
    db::conf::{Database, GenomeRelease, TadSet as TadSetChoice},
    sv::query::{
        interpreter::QueryInterpreter, pathogenic::Record as KnownPathogenicRecord,
        records::StructuralVariant as RecordSv, schema::CaseQuery,
        schema::StructuralVariant as SchemaSv,
    },
};

use self::{
    bgdbs::{load_bg_dbs, BgDbBundle, BgDbOverlaps},
    clinvar::{load_clinvar_sv, ClinvarSv},
    genes::{load_gene_db, GeneDb},
    masked::{load_masked_dbs, MaskedBreakpointCount, MaskedDbBundle},
    pathogenic::{load_patho_dbs, PathoDbBundle},
    schema::{CallInfo, StrandOrientation, SvSubType, SvType, TranscriptEffect},
    tads::{load_tads, TadSetBundle},
};

/// Length of the upstream/downstream region.
static X_STREAM: i32 = 5000;

/// Command line arguments for `sv query` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for SVs", long_about = None)]
pub struct Args {
    /// Genome release to assume.
    #[arg(long, value_enum, default_value_t = GenomeRelease::Grch37)]
    pub genome_release: GenomeRelease,
    /// Path to worker database to use for querying.
    #[arg(long, required = true)]
    pub path_db: String,
    /// Path to query JSON file.
    #[arg(long, required = true)]
    pub path_query_json: String,
    /// Path to input TSV file.
    #[arg(long, required = true)]
    pub path_input_svs: String,
    /// Path to the output TSV file.
    #[arg(long, required = true)]
    pub path_output_svs: String,
    /// Value to write to the result set id column.
    #[arg(long, default_value_t = 0)]
    pub result_set_id: u64,
    /// Optional maximal number of total records to write out.
    #[arg(long)]
    pub max_results: Option<usize>,

    /// Radius around BND sites used when building the database.
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: i32,
    /// Radius around INS sites used when building the database.
    #[arg(long, default_value_t = 50)]
    pub slack_ins: i32,
    /// Minimal reciprocal overlap for SVs of the same type, used when building
    /// the database.
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Maximal distance to TAD to consider.
    #[arg(long, default_value_t = 10_000)]
    pub max_tad_distance: i32,
}

/// Gene information.
#[derive(Debug, Default, Serialize)]
struct Gene {
    /// Gene symbol
    symbol: Option<String>,
    /// ENSEMBL gene ID
    ensembl_id: Option<String>,
    /// Entrez gene ID
    entrez_id: Option<u32>,
    /// Whether the gene is in the ACMG list for incidental findings.
    is_acmg: bool,
    /// Whether the gene is linked to an OMIM disease.
    is_disease_gene: bool,
}

/// Explanation of transcript effect per individual gene.
#[derive(Debug, Default, Serialize)]
struct GeneTranscriptEffects {
    /// Identifying gene.
    gene: Gene,
    /// Transcript effects for the gene.
    transcript_effects: Vec<TranscriptEffect>,
}

/// The structured result information of the result record.
#[derive(Debug, Default, Serialize)]
struct ResultPayload {
    /// The name of the calling tool.
    callers: Vec<String>,
    /// The overlapping VCVs
    clinvar_ovl_vcvs: Vec<String>,
    /// The directly overlapping genes.
    ovl_genes: Vec<Gene>,
    /// Genes that are not directly overlapping but contained in overlapping
    /// TADs.
    tad_genes: Vec<Gene>,
    /// Overlapping known pathogenic SV records.
    known_pathogenic: Vec<KnownPathogenicRecord>,
    /// Information about the call support from the structural variant.
    call_info: IndexMap<String, CallInfo>,
    /// Whether there is an overlap with a disease gene in the overlap.
    ovl_disease_gene: bool,
    /// Whether there is an overlap with a disease gene in the overlapping TADs.
    tad_disease_gene: bool,
    /// The size of the SV, None for ins and BND
    sv_length: Option<u32>,
    /// Overlap counts with background databases.
    overlap_counts: BgDbOverlaps,
    /// Overlap counts with masked sequenced.
    masked_breakpoints: MaskedBreakpointCount,
    /// Distance to next TAD boundary.
    tad_boundary_distance: Option<u32>,
    /// Effects on the transcripts per gene.
    tx_effects: Vec<GeneTranscriptEffects>,
}

/// A result record from the query.
#[derive(Debug, Default, Serialize)]
struct ResultRecord {
    sodar_uuid: Uuid,
    svqueryresultset: u64,
    release: String,
    chromosome: String,
    chromosome_no: i32,
    bin: u32,
    chromosome2: String,
    chromosome_no2: i32,
    bin2: u32,
    start: i32,
    end: i32,
    pe_orientation: StrandOrientation,
    sv_type: SvType,
    sv_sub_type: SvSubType,
    payload: String,
}

fn resolve_gene_id(database: Database, gene_db: &GeneDb, gene_id: u32) -> Vec<Gene> {
    let record_idxs = match database {
        Database::RefSeq => gene_db.xlink.from_entrez.get_vec(&gene_id),
        Database::Ensembl => gene_db.xlink.from_ensembl.get_vec(&gene_id),
    };
    if let Some(record_idxs) = record_idxs {
        record_idxs
            .iter()
            .map(|record_idx| {
                let record = &gene_db.xlink.records[*record_idx as usize];
                Gene {
                    symbol: Some(record.symbol.clone()),
                    ensembl_id: Some(format!("ENSG{:011}", record.ensembl_gene_id)),
                    entrez_id: Some(record.entrez_id),
                    is_acmg: gene_db.acmg.contains(record.entrez_id),
                    is_disease_gene: gene_db.omim.contains(record.entrez_id),
                }
            })
            .collect()
    } else {
        match database {
            Database::RefSeq => vec![Gene {
                entrez_id: Some(gene_id),
                symbol: None,
                ensembl_id: None,
                is_acmg: gene_db.acmg.contains(gene_id),
                is_disease_gene: gene_db.omim.contains(gene_id),
            }],
            Database::Ensembl => vec![Gene {
                ensembl_id: Some(format!("ENSG{gene_id:011}")),
                symbol: None,
                entrez_id: None,
                is_acmg: false,
                is_disease_gene: false,
            }],
        }
    }
}

/// Utility struct to store statistics about counts.
#[derive(Debug, Default)]
struct QueryStats {
    pub count_passed: usize,
    pub count_total: usize,
    pub by_sv_type: BTreeMap<SvType, usize>,
}

/// Open the SV file at `path_sv_tsv` and run through the given `interpreter`.
fn run_query(
    interpreter: &QueryInterpreter,
    args: &Args,
    dbs: &Databases,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    chrom_to_acc: &HashMap<String, String>,
) -> Result<QueryStats, anyhow::Error> {
    let chrom_map = build_chrom_map();
    let mut stats = QueryStats::default();

    // Construct reader and writer for CSV records
    let mut csv_reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(open_read_maybe_gz(&args.path_input_svs)?);
    let mut csv_writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .quote_style(QuoteStyle::Never)
        .from_path(&args.path_output_svs)?;

    // Read through input records using the query interpreter as a filter
    for record in csv_reader.deserialize() {
        stats.count_total += 1;
        let record_sv: RecordSv = record?;
        let schema_sv: SchemaSv = record_sv.clone().into();

        let mut result_payload = ResultPayload {
            call_info: schema_sv.call_info.clone(),
            callers: record_sv.callers.clone(),
            ..ResultPayload::default()
        };

        let mut ovl_gene_ids = Vec::new();

        let sv_query = if matches!(schema_sv.sv_type, SvType::Ins | SvType::Bnd) {
            schema_sv.pos.saturating_sub(1)..schema_sv.pos
        } else {
            schema_sv.pos.saturating_sub(1)..schema_sv.end
        };

        let passes = interpreter.passes(
            &schema_sv,
            &mut |sv: &SchemaSv| {
                result_payload.overlap_counts = dbs.bg_dbs.count_overlaps(
                    sv,
                    &interpreter.query,
                    &chrom_map,
                    args.slack_ins,
                    args.slack_bnd,
                );
                result_payload.overlap_counts.clone()
            },
            &mut |sv: &SchemaSv| {
                result_payload.masked_breakpoints =
                    dbs.masked.masked_breakpoint_count(sv, &chrom_map);
                result_payload.masked_breakpoints.clone()
            },
            &mut |sv: &SchemaSv| {
                ovl_gene_ids = dbs.genes.overlapping_gene_ids(
                    interpreter.query.database,
                    *chrom_map
                        .get(&sv.chrom)
                        .expect("cannot translate chromosome") as u32,
                    sv_query.clone(),
                );
                ovl_gene_ids.sort();
                ovl_gene_ids.clone()
            },
            &mut |sv: &SchemaSv| {
                result_payload.tx_effects =
                    compute_tx_effects(sv, mehari_tx_db, mehari_tx_idx, &dbs.genes, chrom_to_acc);
                let mut res = Vec::new();
                for tx_effect in &result_payload.tx_effects {
                    res.extend(tx_effect.transcript_effects.iter())
                }
                res.sort();
                res.dedup();
                res
            },
        )?;

        if passes.pass_all {
            if schema_sv.sv_type != SvType::Ins && schema_sv.sv_type != SvType::Bnd {
                result_payload.sv_length = Some((schema_sv.end - schema_sv.pos + 1) as u32);
            }

            // Copy effective and compatible genotypes to output.
            for (sample, compatible) in passes.compatible.iter() {
                let call_info = result_payload
                    .call_info
                    .get_mut(sample)
                    .expect("must exist");
                call_info.effective_genotype = *passes.effective.get(sample).expect("must exist");
                call_info.matched_gt_criteria = Some(compatible.clone());
            }

            // Count passing record in statistics
            stats.count_passed += 1;
            *stats.by_sv_type.entry(schema_sv.sv_type).or_default() += 1;

            // Get overlaps with known pathogenic SVs and ClinVar SVs
            result_payload.known_pathogenic =
                dbs.patho_dbs.overlapping_records(&schema_sv, &chrom_map);
            result_payload.clinvar_ovl_vcvs = dbs
                .clinvar_sv
                .overlapping_vcvs(
                    &schema_sv,
                    &chrom_map,
                    interpreter.query.clinvar_sv_min_pathogenicity,
                    interpreter.query.clinvar_sv_min_overlap,
                )
                .into_iter()
                .map(|vcv| format!("VCV{vcv:09}"))
                .collect();

            // Get genes in overlapping TADs
            let tad_gene_ids = {
                let gene_ids: HashSet<_> = HashSet::from_iter(ovl_gene_ids.iter());
                let tads =
                    dbs.tad_sets
                        .overlapping_tads(TadSetChoice::Hesc, &schema_sv, &chrom_map);
                let mut tad_gene_ids = Vec::new();
                tads.iter()
                    .map(|tad| {
                        dbs.genes.overlapping_gene_ids(
                            interpreter.query.database,
                            tad.chrom_no,
                            (tad.begin - 1)..tad.end,
                        )
                    })
                    .for_each(|mut v| tad_gene_ids.append(&mut v));
                let tad_gene_ids: HashSet<_> = HashSet::from_iter(tad_gene_ids.into_iter());
                let mut tad_gene_ids = Vec::from_iter(tad_gene_ids);
                tad_gene_ids.retain(|gene_id| !gene_ids.contains(gene_id));
                tad_gene_ids.sort();
                tad_gene_ids
            };
            result_payload.tad_boundary_distance =
                dbs.tad_sets
                    .boundary_dist(TadSetChoice::Hesc, &schema_sv, &chrom_map);

            // Convert the genes into more verbose records and put them into the result
            ovl_gene_ids.iter().for_each(|gene_id| {
                result_payload.ovl_genes.append(&mut resolve_gene_id(
                    interpreter.query.database,
                    &dbs.genes,
                    *gene_id,
                ))
            });
            result_payload.ovl_disease_gene = result_payload
                .ovl_genes
                .iter()
                .any(|gene| gene.is_disease_gene);
            tad_gene_ids.iter().for_each(|gene_id| {
                result_payload.tad_genes.append(&mut resolve_gene_id(
                    interpreter.query.database,
                    &dbs.genes,
                    *gene_id,
                ))
            });
            result_payload.tad_disease_gene = result_payload
                .tad_genes
                .iter()
                .any(|gene| gene.is_disease_gene);

            if let Some(max_results) = args.max_results {
                if stats.count_total > max_results {
                    warn!(
                        "stopping writing {} records but there are more results!",
                        stats.count_total
                    );
                }
            }

            // Finally, write out the record.
            csv_writer.serialize(&ResultRecord {
                sodar_uuid: uuid::Uuid::new_v4(),
                svqueryresultset: args.result_set_id,
                release: record_sv.release.clone(),
                chromosome: record_sv.chromosome.clone(),
                chromosome_no: record_sv.chromosome_no,
                start: record_sv.start,
                bin: record_sv.bin,
                chromosome2: record_sv
                    .chromosome2
                    .unwrap_or(record_sv.chromosome)
                    .clone(),
                chromosome_no2: record_sv.chromosome_no2,
                bin2: record_sv.bin2,
                end: record_sv.end,
                pe_orientation: record_sv.pe_orientation,
                sv_type: record_sv.sv_type,
                sv_sub_type: record_sv.sv_sub_type,
                payload: serde_json::to_string(&result_payload)?,
            })?;
        }
    }

    Ok(stats)
}

/// Generate `Gene` record for a given Entrez gene ID.
fn construct_gene(entrez_id: u32, gene_db: &GeneDb) -> Gene {
    let idx = gene_db
        .xlink
        .from_entrez
        .get(&entrez_id)
        .expect("gene must exist at this point");
    let record = &gene_db.xlink.records[*idx as usize];
    Gene {
        entrez_id: Some(entrez_id),
        symbol: Some(record.symbol.clone()),
        ensembl_id: Some(format!("ENSG{:011}", record.ensembl_gene_id)),
        is_acmg: gene_db.acmg.contains(record.entrez_id),
        is_disease_gene: gene_db.omim.contains(record.entrez_id),
    }
}

/// Ad-hoc data structure for `tx_regions`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct TxRegion {
    // 0-based begin position
    begin: i32,
    // 0-based end position
    end: i32,
    // "arbitrary" number
    no: usize,
    // effect of the transcript (encodes region type)
    effect: TranscriptEffect,
}

/// Return list of half-open intervals for a given transcript.
fn tx_regions(tx: &Transcript) -> Vec<TxRegion> {
    assert_eq!(
        tx.genome_alignments.len(),
        1,
        "only one alignment supported"
    );
    if tx.genome_alignments[0].exons.is_empty() {
        // no exons? skip!
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut tx_start = None;
    let mut tx_end = None;
    let genome_alignment = &tx.genome_alignments[0];

    // Loop over all exons to determine leftmost and rightmost genome position.
    for exon_alignment in &genome_alignment.exons {
        if let Some(value) = tx_start {
            if exon_alignment.alt_start_i < value {
                tx_start = Some(exon_alignment.alt_start_i);
            }
        } else {
            tx_start = Some(exon_alignment.alt_start_i);
        }
        if let Some(value) = tx_end {
            if exon_alignment.alt_end_i > value {
                tx_end = Some(exon_alignment.alt_end_i);
            }
        } else {
            tx_end = Some(exon_alignment.alt_end_i);
        }
    }

    let tx_start = tx_start.expect("must have been set");
    let tx_end = tx_end.expect("must have been set");

    // Perform the actual region extraction.
    let mut prev_alt_end_i = 0;
    for (no, exon_alignment) in genome_alignment.exons.iter().enumerate() {
        if exon_alignment.alt_start_i == tx_start {
            // is first, register upstream/downstream
            result.push(TxRegion {
                begin: exon_alignment.alt_start_i - X_STREAM,
                end: exon_alignment.alt_start_i - 1,
                no,
                effect: if genome_alignment.strand == Strand::Plus as i32 {
                    TranscriptEffect::UpstreamVariant
                } else {
                    TranscriptEffect::DownstreamVariant
                },
            });
        } else {
            // is not first, register splice region on the left boundary
            result.push(TxRegion {
                begin: (exon_alignment.alt_start_i - 1) - 8,
                end: (exon_alignment.alt_start_i - 1) + 3,
                no,
                effect: TranscriptEffect::SpliceRegionVariant,
            })
        }

        if exon_alignment.alt_end_i == tx_end {
            // is last, register upstream/downstream
            result.push(TxRegion {
                begin: exon_alignment.alt_end_i,
                end: exon_alignment.alt_end_i + X_STREAM,
                no,
                effect: if genome_alignment.strand == Strand::Plus as i32 {
                    TranscriptEffect::DownstreamVariant
                } else {
                    TranscriptEffect::UpstreamVariant
                },
            });
        } else {
            // is not last, register splice region on the right boundary
            result.push(TxRegion {
                begin: exon_alignment.alt_end_i - 3,
                end: exon_alignment.alt_end_i + 8,
                no,
                effect: TranscriptEffect::SpliceRegionVariant,
            })
        }

        // register the exon
        result.push(TxRegion {
            begin: exon_alignment.alt_start_i - 1,
            end: exon_alignment.alt_end_i,
            no,
            effect: TranscriptEffect::ExonVariant,
        });

        if exon_alignment.alt_start_i != tx_start {
            // is not first exon, register intron "right" of it
            result.push(TxRegion {
                begin: prev_alt_end_i,
                end: exon_alignment.alt_start_i - 1,
                no,
                effect: TranscriptEffect::IntronVariant,
            });
        }

        // store end of prev exon for next intron's start
        prev_alt_end_i = exon_alignment.alt_end_i;
    }

    result
}

/// Return the transcript region / effect for the given breakpoint.
fn gene_tx_effects_for_bp(tx: &Transcript, pos: i32) -> Vec<TranscriptEffect> {
    // Obtain list of regions for transcript.
    let regions = tx_regions(tx);

    // Determine how this relates to the breakpoint.
    let pos = pos - 1; // 1-based to 0-based
    let mut result = regions
        .iter()
        .filter(|r| r.begin <= pos && pos < r.end)
        .map(|r| r.effect)
        .collect::<Vec<_>>();
    if result.is_empty() {
        result.push(TranscriptEffect::IntergenicVariant);
    } else {
        result.sort();
        result.dedup();
    }
    result
}

/// Return the transcript region / effect for the given range.
fn gene_tx_effect_for_range(tx: &Transcript, pos: i32, end: i32) -> Vec<TranscriptEffect> {
    // Obtain list of regions for transcript.
    let regions = tx_regions(tx);

    // Determine how this relates to the left and right breakpoints.
    let pos = pos - 1; // 1-based to 0-based
    let mut result = regions
        .iter()
        .filter(|region| pos < region.end && region.begin < end)
        .map(|region| region.effect)
        .collect::<Vec<_>>();

    // Remove any duplicates.
    result.sort();
    result.dedup();

    // If we have both upstream and downstream then the full transcript is affected.
    if result.contains(&TranscriptEffect::UpstreamVariant)
        && result.contains(&TranscriptEffect::DownstreamVariant)
    {
        result.push(TranscriptEffect::TranscriptVariant);
    }

    result
}

/// Helper that computes effects on transcripts for a single breakend, e.g., one side of BND or INS.
fn compute_tx_effects_for_breakpoint(
    sv: &SchemaSv,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    gene_db: &GeneDb,
    chrom_to_acc: &HashMap<String, String>,
) -> Vec<GeneTranscriptEffects> {
    // Shortcut to the `TranscriptDb`.
    let tx_db = mehari_tx_db
        .tx_db
        .as_ref()
        .expect("transcripts must be present");
    // Compute canonical chromosome name and map to accession.
    let chrom = chrom_to_acc.get(&annonars::common::cli::canonicalize(&sv.chrom));
    if chrom.is_none() {
        return Default::default();
    }
    let chrom = chrom.expect("chromosome must be known at this point");
    // Create range to query the interval trees for.
    let query = (sv.pos - X_STREAM)..(sv.pos + X_STREAM);

    if let Some(idx) = mehari_tx_idx.contig_to_idx.get(chrom) {
        let mut effects_by_gene: HashMap<_, Vec<_>> = HashMap::new();

        // Collect all transcripts that overlap the INS and compute the effect of the INS on
        // the transcript.
        let tree = &mehari_tx_idx.trees[*idx];
        for it in tree.find(query) {
            let tx = &tx_db.transcripts[*it.data() as usize];
            if let Some(&idx) = gene_db
                .xlink
                .from_hgnc
                .get(&format!("HGNC:{}", &tx.gene_id))
            {
                let entrez_id = gene_db.xlink.records[idx as usize].entrez_id;
                effects_by_gene
                    .entry(entrez_id)
                    .or_default()
                    .extend(gene_tx_effects_for_bp(tx, sv.pos));
            } else {
                tracing::warn!("could not resolve HGNC gene ID {:?}", tx.gene_id)
            }
        }

        // Deduplicate results.
        effects_by_gene.iter_mut().for_each(|(_, v)| {
            v.sort();
            v.dedup()
        });

        // Convert the results into the final format.
        effects_by_gene
            .into_iter()
            .map(|(entrez_id, transcript_effects)| GeneTranscriptEffects {
                gene: construct_gene(entrez_id, gene_db),
                transcript_effects,
            })
            .collect()
    } else {
        // We do not have any transcripts for this chromosome.
        Default::default()
    }
}

/// Compute effect for linear SVs.
fn compute_tx_effects_for_linear(
    sv: &SchemaSv,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    gene_db: &GeneDb,
    chrom_to_acc: &HashMap<String, String>,
) -> Vec<GeneTranscriptEffects> {
    // Shortcut to the `TranscriptDb`.
    let tx_db = mehari_tx_db
        .tx_db
        .as_ref()
        .expect("transcripts must be present");
    // Compute canonical chromosome name and map to accession.
    let chrom = chrom_to_acc.get(&annonars::common::cli::canonicalize(&sv.chrom));
    if chrom.is_none() {
        return Default::default();
    }
    let chrom = chrom.expect("chromosome must be known at this point");
    // Create range to query the interval trees for.
    let query = (sv.pos - X_STREAM)..(sv.end + X_STREAM);

    if let Some(idx) = mehari_tx_idx.contig_to_idx.get(chrom) {
        let mut effects_by_gene: HashMap<_, Vec<_>> = HashMap::new();

        // Collect all transcripts that overlap the SV and compute the effect of the SV on
        // the transcript.
        let tree = &mehari_tx_idx.trees[*idx];
        for it in tree.find(query) {
            let tx = &tx_db.transcripts[*it.data() as usize];
            if let Some(&idx) = gene_db
                .xlink
                .from_hgnc
                .get(&format!("HGNC:{}", &tx.gene_id))
            {
                let entrez_id = gene_db.xlink.records[idx as usize].entrez_id;
                effects_by_gene
                    .entry(entrez_id)
                    .or_default()
                    .extend(gene_tx_effect_for_range(tx, sv.pos, sv.end));
            } else {
                tracing::warn!("could not resolve HGNC gene ID {:?}", tx.gene_id)
            }
        }

        // Deduplicate results.
        effects_by_gene.iter_mut().for_each(|(_, v)| {
            v.sort();
            v.dedup()
        });

        // Convert the results into the final format.
        effects_by_gene
            .into_iter()
            .map(|(entrez_id, transcript_effects)| GeneTranscriptEffects {
                gene: construct_gene(entrez_id, gene_db),
                transcript_effects,
            })
            .collect()
    } else {
        // We do not have any transcripts for this chromosome.
        Default::default()
    }
}

/// Compute effect(s) of `sv` on transcript of genes.
fn compute_tx_effects(
    sv: &SchemaSv,
    mehari_tx_db: &TxSeqDatabase,
    mehari_tx_idx: &TxIntervalTrees,
    gene_db: &GeneDb,
    chrom_to_acc: &HashMap<String, String>,
) -> Vec<GeneTranscriptEffects> {
    match sv.sv_type {
        SvType::Ins | SvType::Bnd => compute_tx_effects_for_breakpoint(
            sv,
            mehari_tx_db,
            mehari_tx_idx,
            gene_db,
            chrom_to_acc,
        ),
        SvType::Del | SvType::Dup | SvType::Inv | SvType::Cnv => {
            compute_tx_effects_for_linear(sv, mehari_tx_db, mehari_tx_idx, gene_db, chrom_to_acc)
        }
    }
}

/// Bundle the used database to reduce argument count.
#[derive(Default, Debug)]
pub struct Databases {
    pub bg_dbs: BgDbBundle,
    pub patho_dbs: PathoDbBundle,
    pub tad_sets: TadSetBundle,
    pub masked: MaskedDbBundle,
    pub genes: GeneDb,
    pub clinvar_sv: ClinvarSv,
}

/// Translate gene allow list to gene identifier sfrom
fn translate_gene_allowlist(
    gene_allowlist: &Vec<String>,
    dbs: &Databases,
    query: &CaseQuery,
) -> HashSet<u32> {
    let mut result = HashSet::new();

    let re_entrez = regex::Regex::new(r"^\d+").expect("invalid regex in source code");
    let re_ensembl = regex::Regex::new(r"ENSG\d+").expect("invalid regex in source code");

    let symbol_to_id: HashMap<_, _> =
        HashMap::from_iter(dbs.genes.xlink.records.iter().map(|record| {
            (
                record.symbol.clone(),
                match query.database {
                    Database::RefSeq => record.entrez_id,
                    Database::Ensembl => record.ensembl_gene_id,
                },
            )
        }));

    for gene in gene_allowlist {
        let gene = gene.trim();
        if re_entrez.is_match(gene) {
            if let Ok(gene_id) = numeric_gene_id(gene) {
                match query.database {
                    Database::RefSeq => {
                        result.insert(gene_id);
                    }
                    Database::Ensembl => {
                        if let Some(record_ids) = dbs.genes.xlink.from_ensembl.get_vec(&gene_id) {
                            for record_id in record_ids {
                                result
                                    .insert(dbs.genes.xlink.records[*record_id as usize].entrez_id);
                            }
                        };
                    }
                }
            } else {
                warn!("Cannot map candidate Entrez gene identifier {}", &gene);
                continue;
            }
        } else if re_ensembl.is_match(gene) {
            if let Ok(gene_id) = numeric_gene_id(gene) {
                match query.database {
                    Database::RefSeq => {
                        if let Some(record_ids) = dbs.genes.xlink.from_entrez.get_vec(&gene_id) {
                            for record_id in record_ids {
                                result.insert(
                                    dbs.genes.xlink.records[*record_id as usize].ensembl_gene_id,
                                );
                            }
                        };
                    }
                    Database::Ensembl => {
                        result.insert(gene_id);
                    }
                }
            } else {
                warn!("Cannot map candidate ENSEMBL gene identifier {}", &gene);
                continue;
            }
        } else if let Some(gene_id) = symbol_to_id.get(gene) {
            result.insert(*gene_id);
        } else {
            warn!("Could not map candidate gene symbol {}", &gene);
        }
    }

    result
}

/// Load database from the given path with the given genome release.
fn load_databases(
    path_worker_db: &str,
    genome_release: GenomeRelease,
    max_tad_distance: i32,
) -> Result<Databases, anyhow::Error> {
    Ok(Databases {
        bg_dbs: load_bg_dbs(path_worker_db, genome_release)?,
        patho_dbs: load_patho_dbs(path_worker_db, genome_release)?,
        tad_sets: load_tads(path_worker_db, genome_release, max_tad_distance)?,
        masked: load_masked_dbs(path_worker_db, genome_release)?,
        genes: load_gene_db(path_worker_db, genome_release)?,
        clinvar_sv: load_clinvar_sv(path_worker_db, genome_release)?,
    })
}

/// Main entry point for `sv query` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = Instant::now();
    tracing::info!("args_common = {:?}", &args_common);
    tracing::info!("args = {:?}", &args);

    tracing::info!("Loading query...");
    let query: CaseQuery = serde_json::from_reader(File::open(&args.path_query_json)?)?;
    tracing::info!(
        "... done loading query = {}",
        &serde_json::to_string(&query)?
    );

    tracing::info!("Loading worker databases...");
    let before_loading = Instant::now();
    let path_worker_db = format!("{}/worker", &args.path_db);
    let dbs = load_databases(&path_worker_db, args.genome_release, args.max_tad_distance)?;
    tracing::info!(
        "...done loading databases in {:?}",
        before_loading.elapsed()
    );

    trace_rss_now();

    tracing::info!("Loading mehari tx database...");
    let before_loading = Instant::now();
    let path_mehari_tx_db = format!(
        "{}/mehari/{}/txs.bin.zst",
        &args.path_db,
        &args.genome_release.to_string()
    );
    let mehari_tx_db = mehari::annotate::seqvars::load_tx_db(&path_mehari_tx_db)?;
    tracing::info!(
        "...done loading mehari tx database in {:?}",
        before_loading.elapsed()
    );
    tracing::info!("Building mehari index data structures...");
    let before_building = Instant::now();
    let mehari_tx_idx = TxIntervalTrees::new(&mehari_tx_db, args.genome_release.into());
    let chrom_to_acc = ASSEMBLY_INFOS[args.genome_release.into()]
        .sequences
        .iter()
        .map(|record| {
            (
                annonars::common::cli::canonicalize(&record.name),
                record.refseq_ac.clone(),
            )
        })
        .collect::<HashMap<_, _>>();
    tracing::info!(
        "...done building mehari index data structures in {:?}",
        before_building.elapsed()
    );

    trace_rss_now();

    tracing::info!("Translating gene allow list...");
    let gene_allowlist = if let Some(gene_allowlist) = &query.gene_allowlist {
        if gene_allowlist.is_empty() {
            None
        } else {
            Some(translate_gene_allowlist(gene_allowlist, &dbs, &query))
        }
    } else {
        None
    };

    tracing::info!("Running queries...");
    let before_query = Instant::now();
    let query_stats = run_query(
        &QueryInterpreter::new(query, gene_allowlist),
        args,
        &dbs,
        &mehari_tx_db,
        &mehari_tx_idx,
        &chrom_to_acc,
    )?;
    tracing::info!("... done running query in {:?}", before_query.elapsed());
    tracing::info!(
        "summary: {} records passed out of {}",
        query_stats.count_passed.separate_with_commas(),
        query_stats.count_total.separate_with_commas()
    );
    tracing::info!("passing records by SV type");
    for (sv_type, count) in query_stats.by_sv_type.iter() {
        tracing::info!("{:?} -- {}", sv_type, count);
    }

    trace_rss_now();

    tracing::info!(
        "All of `sv query` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}
