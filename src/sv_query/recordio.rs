use std::{collections::HashMap, fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use byte_unit::Byte;
use flate2::read::GzDecoder;
use serde::de::DeserializeOwned;
use thousands::Separable;

use super::dbrecords::{self, BeginEnd, ChromosomeCoordinate, ToInMemory};

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

/// This struct bundles the background database records by chromosome.
pub struct SvRecordsByChrom {
    pub bg_sv_records: Vec<Vec<dbrecords::bg_sv::Record>>,
    pub bg_sv_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
    pub gnomad_sv_records: Vec<Vec<dbrecords::gnomad_sv::Record>>,
    pub gnomad_sv_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
    pub dbvar_records: Vec<Vec<dbrecords::dbvar::Record>>,
    pub dbvar_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
    pub dgv_records: Vec<Vec<dbrecords::dgv::Record>>,
    pub dgv_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
    pub dgv_gs_records: Vec<Vec<dbrecords::dgv_gs::Record>>,
    pub dgv_gs_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
    pub exac_cnv_records: Vec<Vec<dbrecords::exac_cnv::Record>>,
    pub exac_cnv_trees: Vec<ArrayBackedIntervalTree<i32, usize>>,
}

pub fn load_sv_records(
    term: &console::Term,
    db_base_dir: &str,
) -> Result<SvRecordsByChrom, anyhow::Error> {
    let chrom_map = build_chrom_map();

    let before_parsing = Instant::now();
    print_rss_now(term)?;

    let bg_sv_path = Path::new(&db_base_dir)
        .join("bg-inhouse")
        .join("varfish-sv.tsv.gz");
    let bg_sv_records = load_bg_sv_records::<dbrecords::bg_sv::Record, dbrecords::bg_sv::FileRecord>(
        &bg_sv_path,
        term,
        &chrom_map,
    )?;
    let bg_sv_trees = build_bg_sv_tree(term, &bg_sv_records)?;
    print_rss_now(term)?;

    let gnomad_sv_path = Path::new(&db_base_dir)
        .join("bg-public")
        .join("gnomad-sv.tsv.gz");
    let gnomad_sv_records = load_bg_sv_records::<
        dbrecords::gnomad_sv::Record,
        dbrecords::gnomad_sv::FileRecord,
    >(&gnomad_sv_path, term, &chrom_map)?;
    let gnomad_sv_trees = build_bg_sv_tree(term, &gnomad_sv_records)?;
    print_rss_now(term)?;

    let dbvar_path = Path::new(&db_base_dir)
        .join("bg-public")
        .join("dbvar-sv.tsv.gz");
    let dbvar_records = load_bg_sv_records::<dbrecords::dbvar::Record, dbrecords::dbvar::FileRecord>(
        &dbvar_path,
        term,
        &chrom_map,
    )?;
    let dbvar_trees = build_bg_sv_tree(term, &dbvar_records)?;
    print_rss_now(term)?;

    let dgv_path = Path::new(&db_base_dir)
        .join("bg-public")
        .join("dgv-sv.tsv.gz");
    let dgv_records = load_bg_sv_records::<dbrecords::dgv::Record, dbrecords::dgv::FileRecord>(
        &dgv_path, term, &chrom_map,
    )?;
    let dgv_trees = build_bg_sv_tree(term, &dgv_records)?;
    print_rss_now(term)?;

    let dgv_gs_path = Path::new(&db_base_dir)
        .join("bg-public")
        .join("dgv-gs-sv.tsv.gz");
    let dgv_gs_records = load_bg_sv_records::<
        dbrecords::dgv_gs::Record,
        dbrecords::dgv_gs::FileRecord,
    >(&dgv_gs_path, term, &chrom_map)?;
    let dgv_gs_trees = build_bg_sv_tree(term, &dgv_gs_records)?;
    print_rss_now(term)?;

    let exac_cnv_path = Path::new(&db_base_dir)
        .join("bg-public")
        .join("exac-cnv.tsv.gz");
    let exac_cnv_records = load_bg_sv_records::<
        dbrecords::exac_cnv::Record,
        dbrecords::exac_cnv::FileRecord,
    >(&exac_cnv_path, term, &chrom_map)?;
    let exac_cnv_trees = build_bg_sv_tree(term, &exac_cnv_records)?;
    print_rss_now(term)?;

    term.write_line(&format!(
        "Total time spent parsing: {:?}",
        before_parsing.elapsed()
    ))?;

    Ok(SvRecordsByChrom {
        bg_sv_records,
        bg_sv_trees,
        gnomad_sv_records,
        gnomad_sv_trees,
        dbvar_records,
        dbvar_trees,
        dgv_records,
        dgv_trees,
        dgv_gs_records,
        dgv_gs_trees,
        exac_cnv_records,
        exac_cnv_trees,
    })
}
