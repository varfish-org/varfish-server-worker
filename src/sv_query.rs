use clap::Parser;

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use std::time::Instant;
use std::{collections::HashMap, fs::File, path::Path};

use flate2::read::GzDecoder;
use serde::Deserialize;
use thousands::Separable;

use crate::common::Args as CommonArgs;

/// Background SV database record to be kept in memory.
#[derive(Debug)]
pub struct BgSvRecord {
    /// The 0-based begin position.
    pub begin: i32,
    /// The 0-based end position.
    pub end: i32,

    /// Total number of carriers.
    pub carriers: i32,
    /// Number of het. carriers.
    pub carriers_het: i32,
    /// Number of hom. carriers.
    pub carriers_hom: i32,
    /// Number of hemi. carriers.
    pub carriers_hemi: i32,
}

/// Background SV database record as read from TSV file.
#[derive(Debug, Deserialize)]
pub struct BgSvFileRecord {
    pub id: i32,
    pub release: String,
    pub chromosome: String,
    pub chromosome_no: i32,
    /// start position, 1-based
    pub start: i32,
    pub chromosome2: String,
    pub chromosome_no2: i32,
    /// end position, 1-based
    pub end: i32,
    pub pe_orientation: String,
    pub sv_type: String,
    pub bin: i32,
    pub src_count: i32,
    pub carriers: i32,
    pub carriers_het: i32,
    pub carriers_hom: i32,
    pub carriers_hemi: i32,
    pub bg_sv_set_id: i32,
}

impl BgSvFileRecord {
    /// Convert on-disk record from TSV to in-memory record.
    pub fn to_in_memory(&self) -> BgSvRecord {
        BgSvRecord {
            begin: self.start - 1,
            end: self.end,
            carriers: self.carriers,
            carriers_het: self.carriers_het,
            carriers_hom: self.carriers_hom,
            carriers_hemi: self.carriers_hemi,
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Base directory path for datbases
    #[arg(long)]
    pub db_base_dir: String,
    /// Integer to write as result set ID
    #[arg(long)]
    pub result_set_id: i32,
    /// Query JSON file
    #[arg(long)]
    pub query_json: String,
    /// Path to input VCF file
    #[arg(long)]
    pub input_vcf: String,
    /// Path to output VCF file
    #[arg(long)]
    pub output_vcf: String,
}

const CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M",
];

pub fn build_chrom_map() -> HashMap<String, usize> {
    let mut result = HashMap::new();
    for (i, &chrom_name) in CHROMS.iter().enumerate() {
        result.insert(chrom_name.to_owned(), i);
        result.insert(format!("chr{}", chrom_name).to_owned(), i);
    }
    result.insert("MT".to_owned(), 24);
    result.insert("chrMT".to_owned(), 24);
    result
}

fn load_bg_sv_records(
    args: &Args,
    term: &console::Term,
    chrom_map: HashMap<String, usize>,
) -> Result<Vec<Vec<BgSvRecord>>, anyhow::Error> {
    let mut rec_by_contig: Vec<Vec<BgSvRecord>> = Vec::new();
    for _i in 0..25 {
        rec_by_contig.push(Vec::new());
    }
    let before_parsing = Instant::now();
    let path = Path::new(&args.db_base_dir)
        .join("bg-inhouse")
        .join("varfish-sv.tsv.gz");
    term.write_line(&format!("Parsing {}", path.display()))?;
    let file = File::open(&path)?;
    let decoder = GzDecoder::new(&file);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(decoder);
    for result in rdr.deserialize() {
        let record: BgSvFileRecord = result?;
        let idx = *chrom_map
            .get(&record.chromosome)
            .unwrap_or_else(|| panic!("unknown chromosome {}", &record.chromosome));
        rec_by_contig[idx].push(record.to_in_memory());
    }
    term.write_line(&format!(
        "-- time spent parsing: {:?}",
        before_parsing.elapsed()
    ))?;
    Ok(rec_by_contig)
}

fn build_bg_sv_tree(
    term: &console::Term,
    rec_by_contig: &[Vec<BgSvRecord>],
) -> Result<Vec<ArrayBackedIntervalTree<i32, usize>>, anyhow::Error> {
    term.write_line("Building trees...")?;
    let mut trees: Vec<ArrayBackedIntervalTree<i32, usize>> = Vec::new();
    let before_building = Instant::now();
    for (i, records) in rec_by_contig.iter().enumerate() {
        let before_tree = Instant::now();
        let mut tree = ArrayBackedIntervalTree::new();
        for (i, record) in records.iter().enumerate() {
            tree.insert((record.begin)..(record.end), i);
        }
        trees.push(tree);
        term.write_line(&format!(
            "  time for {} ({} intervals): {:?}",
            CHROMS[i],
            records.len().separate_with_commas(),
            before_tree.elapsed()
        ))?;
    }
    term.write_line(&format!(
        "-- total time spent building trees: {:?}",
        before_building.elapsed()
    ))?;
    Ok(trees)
}

pub(crate) fn run(
    term: &console::Term,
    common: &CommonArgs,
    args: &Args,
) -> Result<(), anyhow::Error> {
    term.write_line("Starting sv-query")?;
    term.write_line(&format!("common = {:?}", &common))?;
    term.write_line(&format!("args = {:?}", &args))?;
    term.write_line("")?;

    let chrom_map = build_chrom_map();

    let bg_sv_records_by_config = load_bg_sv_records(args, term, chrom_map)?;
    let _bg_sv_trees = build_bg_sv_tree(term, &bg_sv_records_by_config)?;

    Ok(())
}
