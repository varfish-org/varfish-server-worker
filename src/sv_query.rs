pub mod dbrecords;

use std::time::Instant;
use std::{collections::HashMap, fs::File, path::Path};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use byte_unit::Byte;
use clap::Parser;
use flate2::read::GzDecoder;
use serde::de::DeserializeOwned;
use thousands::Separable;

use crate::common::Args as CommonArgs;
use crate::sv_query::dbrecords::{BeginEnd, ChromosomeCoordinate, ToInMemory};

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
    result.insert("x".to_owned(), 22);
    result.insert("y".to_owned(), 23);
    result.insert("chrx".to_owned(), 22);
    result.insert("chry".to_owned(), 23);
    result.insert("mt".to_owned(), 24);
    result.insert("m".to_owned(), 24);
    result.insert("chrmt".to_owned(), 24);
    result.insert("chrm".to_owned(), 24);
    result.insert("MT".to_owned(), 24);
    result.insert("chrMT".to_owned(), 24);
    result
}

fn load_bg_sv_records<
    ResultRecord,
    FileRecord: ChromosomeCoordinate + DeserializeOwned + ToInMemory<ResultRecord>,
>(
    bg_sv_path: &Path,
    term: &console::Term,
    chrom_map: &HashMap<String, usize>,
) -> Result<Vec<Vec<ResultRecord>>, anyhow::Error> {
    let mut rec_by_contig: Vec<Vec<ResultRecord>> = Vec::new();
    for _i in 0..25 {
        rec_by_contig.push(Vec::new());
    }
    let before_parsing = Instant::now();
    term.write_line(&format!("Parsing {}", bg_sv_path.display()))?;
    let file = File::open(bg_sv_path)?;
    let decoder = GzDecoder::new(&file);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(decoder);
    for result in rdr.deserialize() {
        let record: FileRecord = result?;
        let idx = *chrom_map
            .get(record.chromosome())
            .unwrap_or_else(|| panic!("unknown chromosome {}", record.chromosome()));
        rec_by_contig[idx].push(record.to_in_memory());
    }
    term.write_line(&format!(
        "-- time spent parsing: {:?}",
        before_parsing.elapsed()
    ))?;
    Ok(rec_by_contig)
}

fn build_bg_sv_tree<Record: BeginEnd>(
    term: &console::Term,
    rec_by_contig: &[Vec<Record>],
) -> Result<Vec<ArrayBackedIntervalTree<i32, usize>>, anyhow::Error> {
    term.write_line("Building trees...")?;
    let mut trees: Vec<ArrayBackedIntervalTree<i32, usize>> = Vec::new();
    let before_building = Instant::now();
    for (i, contig_records) in rec_by_contig.iter().enumerate() {
        let before_tree = Instant::now();
        let mut tree = ArrayBackedIntervalTree::new();
        for (i, record) in contig_records.iter().enumerate() {
            tree.insert((record.begin())..(record.end()), i);
        }
        trees.push(tree);
        term.write_line(&format!(
            "  time for {} ({} intervals): {:?}",
            CHROMS[i],
            contig_records.len().separate_with_commas(),
            before_tree.elapsed()
        ))?;
    }
    term.write_line(&format!(
        "-- total time spent building trees: {:?}",
        before_building.elapsed()
    ))?;
    Ok(trees)
}

fn print_rss_now(term: &console::Term) -> Result<(), anyhow::Error> {
    let me = procfs::process::Process::myself().unwrap();
    let page_size = procfs::page_size().unwrap();
    term.write_line(&format!(
        "RSS now: {}",
        Byte::from_bytes((me.stat().unwrap().rss * page_size) as u128).get_appropriate_unit(true)
    ))?;
    Ok(())
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

    print_rss_now(term)?;

    let bg_sv_path = Path::new(&args.db_base_dir)
        .join("bg-inhouse")
        .join("varfish-sv.tsv.gz");
    let bg_sv_records_by_contig = load_bg_sv_records::<
        dbrecords::bg_sv::Record,
        dbrecords::bg_sv::FileRecord,
    >(&bg_sv_path, term, &chrom_map)?;
    let _bg_sv_trees = build_bg_sv_tree(term, &bg_sv_records_by_contig)?;
    print_rss_now(term)?;

    let gnomad_sv_path = Path::new(&args.db_base_dir)
        .join("bg-public")
        .join("gnomad-sv.tsv.gz");
    let gnomad_sv_records_by_contig = load_bg_sv_records::<
        dbrecords::gnomad_sv::Record,
        dbrecords::gnomad_sv::FileRecord,
    >(&gnomad_sv_path, term, &chrom_map)?;
    let _gnomad_sv_trees = build_bg_sv_tree(term, &gnomad_sv_records_by_contig)?;
    print_rss_now(term)?;

    let dbvar_sv_path = Path::new(&args.db_base_dir)
        .join("bg-public")
        .join("dbvar-sv.tsv.gz");
    let dbvar_sv_records_by_contig = load_bg_sv_records::<
        dbrecords::db_var::Record,
        dbrecords::db_var::FileRecord,
    >(&dbvar_sv_path, term, &chrom_map)?;
    let _dbvar_sv_trees = build_bg_sv_tree(term, &dbvar_sv_records_by_contig)?;
    print_rss_now(term)?;

    Ok(())
}
