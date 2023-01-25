//! Implementation of `sv build-bgdb` command

pub mod recordio;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
    time::Instant,
};

use bio::data_structures::interval_tree::IntervalTree;
use clap::Parser;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use serde_json::to_writer;
use serde_jsonlines::JsonLinesReader;
use strum::IntoEnumIterator;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, print_rss_now, Args as CommonArgs, CHROMS},
    sv_query::schema::{StrandOrientation, SvType},
};

use self::recordio::FileRecord;

/// Representation of the fields for the in-house background database.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct BackgroundDbRecord {
    /// genome build
    pub release: String,
    /// chromosome name
    pub chromosome: String,
    /// start position, 1-based
    pub start: i32,
    /// chromosome2 name
    pub chromosome2: String,
    /// end position, 1-based
    pub end: i32,
    /// paired-end orientation
    pub pe_orientation: StrandOrientation,
    /// type of the SV
    pub sv_type: SvType,
    /// number of overall carriers
    pub carriers: u32,
    /// number of het. carriers
    pub carriers_het: u32,
    /// number of hom. carriers
    pub carriers_hom: u32,
    /// number of hemi. carriers
    pub carriers_hemi: u32,
}

impl BackgroundDbRecord {
    /// Compute reciprocal overlap between `self` and `other`.
    pub fn overlap(&self, other: &BackgroundDbRecord) -> f32 {
        let s1 = if self.start > 0 { self.start - 1 } else { 0 };
        let e1 = self.end + 1;
        let s2 = if other.start > 0 { other.start - 1 } else { 0 };
        let e2 = other.end;

        let ovl_s = std::cmp::max(s1, s2);
        let ovl_e = std::cmp::min(e1, e2);
        if ovl_e <= ovl_s {
            0.0
        } else {
            let len1 = (e1 - s1) as f32;
            let len2 = (e2 - s2) as f32;
            let ovl_len = (ovl_e - ovl_s) as f32;
            (ovl_len / len1).min(ovl_len / len2)
        }
    }

    pub fn merge_into(&mut self, other: &BackgroundDbRecord) {
        self.carriers += other.carriers;
        self.carriers_het += other.carriers_het;
        self.carriers_hom += other.carriers_hom;
        self.carriers_hemi += other.carriers_hemi;
    }

    fn from_db_record(record: FileRecord) -> Self {
        BackgroundDbRecord {
            release: record.release,
            chromosome: record.chromosome,
            start: record.start,
            chromosome2: record.chromosome2,
            end: record.end,
            pe_orientation: record.pe_orientation,
            sv_type: record.sv_type,
            carriers: record.num_het + record.num_hom_alt + record.num_hemi_alt,
            carriers_het: record.num_het,
            carriers_hom: record.num_hom_alt,
            carriers_hemi: record.num_hemi_alt,
        }
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

/// Struct to parse command line arguments into.
#[derive(Parser, Debug)]
#[command(author, version, about = "Build SV background database", long_about = None)]
pub struct Args {
    /// Output TSV file to write with resulting SV background db.
    #[arg(long)]
    pub output_tsv: String,
    /// The value to write to the output background SV set ID column.
    #[arg(long)]
    pub output_bg_sv_set_id: i32,
    /// Input files to cluster, prefix with `@` to file with line-wise paths.
    #[arg(required = true)]
    pub input_tsv: Vec<String>,
    /// Minimal reciprocal overlap to use (slightly more strict that the normal query value of 0.75).
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Padding to use for BNDs
    #[arg(long, default_value_t = 50)]
    pub bnd_slack: u32,
    /// Padding to use for INS
    #[arg(long, default_value_t = 50)]
    pub ins_slack: u32,
}

/// Create one file with records for each chromosome and SV type.
fn create_tmp_files(
    tmp_dir: &tempdir::TempDir,
) -> Result<HashMap<(usize, SvType), BufWriter<File>>, anyhow::Error> {
    let mut files = HashMap::new();

    for (chrom_no, chrom) in CHROMS.iter().enumerate() {
        for sv_type in SvType::iter() {
            let path = tmp_dir
                .path()
                .join(format!("records.chr{}.{:?}.tsv", *chrom, sv_type));
            files.insert((chrom_no, sv_type), BufWriter::new(File::create(path)?));
        }
    }

    Ok(files)
}

/// Split the input into one file in `tmp_dir` for each chromosome and SV type.
fn split_input_by_chrom_and_sv_type(
    tmp_dir: &tempdir::TempDir,
    input_tsv_paths: Vec<String>,
    term: &console::Term,
) -> Result<(), anyhow::Error> {
    let mut tmp_files = create_tmp_files(tmp_dir)?;
    let chrom_map = build_chrom_map();
    let before_parsing = Instant::now();
    let mut count_files = 0;
    for path in &input_tsv_paths {
        term.write_line(&format!("- parsing {:?}", &path))?;

        let file = File::open(path)?;
        let decoder = GzDecoder::new(&file);
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(decoder);
        let before_parsing = Instant::now();
        let mut count_records = 0;
        for result in rdr.deserialize() {
            let record: FileRecord = result?;

            let chrom_no = *chrom_map
                .get(&record.chromosome)
                .expect("unknown chromosome");
            let sv_type = record.sv_type;
            let mut tmp_file = tmp_files
                .get_mut(&(chrom_no, sv_type))
                .expect("no file for chrom/sv_type");
            to_writer(&mut tmp_file, &record)?;
            tmp_file.write_all(&[b'\n'])?;

            count_records += 1;
        }
        print_rss_now(term)?;
        term.write_line(&format!(
            "-- total time spent parsing {} records: {:?}",
            count_records.separate_with_commas(),
            before_parsing.elapsed()
        ))?;

        count_files += 1;
    }
    term.write_line(&format!(
        "== total time spent parsing {} files: {:?}",
        count_files.separate_with_commas(),
        before_parsing.elapsed()
    ))?;
    print_rss_now(term)?;
    for (_, mut f) in tmp_files.drain() {
        f.get_mut().sync_all()?
    }
    Ok(())
}

/// Read in all records from `reader`, merge overlapping ones.
///
/// The idea to merge here is to get rid of large stacks of SVs with a reciprocal overlap
/// that is more strict than the 0.75 that is generally used for querying.  We merge with
/// existing clusters with the reciprocal overlap is >=0.8 for all members.
fn merge_to_out(
    args: &Args,
    reader: &mut BufReader<File>,
    writer: &mut csv::Writer<File>,
    term: &console::Term,
) -> Result<(), anyhow::Error> {
    let mut clusters: Vec<Vec<usize>> = vec![];
    let mut tree: IntervalTree<i32, usize> = IntervalTree::new();
    let mut records: Vec<BackgroundDbRecord> = Vec::new();

    // Read in all records and perform the "merge compression"
    let mut reader = JsonLinesReader::new(reader);
    while let Ok(Some(record)) = reader.read::<FileRecord>() {
        let record = BackgroundDbRecord::from_db_record(record);
        let slack = match record.sv_type {
            SvType::Bnd => args.bnd_slack as i32,
            SvType::Ins => args.ins_slack as i32,
            _ => 0,
        };
        let query = (record.start - 1 - slack)..(record.end + slack);
        let mut found_any_cluster = false;
        for mut it_tree in tree.find_mut(&query) {
            let cluster_idx = *it_tree.data();
            let mut match_all_in_cluster = true;
            for it_cluster in &clusters[cluster_idx] {
                let record_id = it_cluster;
                let match_this = match record.sv_type {
                    SvType::Bnd | SvType::Ins => true,
                    _ => {
                        let ovl = record.overlap(&records[*record_id]);
                        assert!(ovl >= 0f32);
                        ovl >= args.min_overlap
                    }
                };
                match_all_in_cluster = match_all_in_cluster && match_this;
            }
            if match_all_in_cluster {
                // extend cluster
                // println!("extending cluster {:?}", cluster_idx);
                clusters[cluster_idx].push(records.len());
                found_any_cluster = true;
                break;
            }
        }
        if !found_any_cluster {
            // create new cluster
            tree.insert((record.start - 1)..(record.end), clusters.len());
            clusters.push(vec![records.len()]);
        }
        // always register the record
        records.push(record);
    }

    print_rss_now(term)?;

    // Sort the cluster representatives by start coordinate.
    let mut sorted_idxs = vec![0; clusters.len()];
    for i in 0..clusters.len() {
        sorted_idxs[i] = i;
    }
    sorted_idxs.sort_by(|a, b| {
        (records[clusters[*a][0]].start, records[clusters[*a][0]].end)
            .partial_cmp(&(records[clusters[*b][0]].start, records[clusters[*b][0]].end))
            .unwrap()
    });

    // Finally, write out all records in sorted order
    for cluster in clusters {
        let mut out_record = records[cluster[0]].clone();
        for record_id in &cluster[1..] {
            out_record.merge_into(&records[*record_id]);
        }
        writer.serialize(&out_record)?;
    }

    Ok(())
}

/// Perform (chrom, sv_type) wise merging of records in temporary files.
fn merge_split_files(
    tmp_dir: &tempdir::TempDir,
    args: &Args,
    term: &console::Term,
) -> Result<(), anyhow::Error> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&args.output_tsv)?;

    for chrom in CHROMS {
        for sv_type in SvType::iter() {
            let filename = format!("records.chr{}.{:?}.tsv", *chrom, sv_type);
            let path = tmp_dir.path().join(&filename);
            term.write_line(&format!("reading from {}", &filename))?;
            let mut reader = BufReader::new(File::open(path)?);
            merge_to_out(args, &mut reader, &mut writer, term)?;
        }
    }

    writer.flush()?;

    Ok(())
}

/// Main entry point for the `sv build-bgdb` command.
pub fn run(term: &console::Term, common: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line("Starting sv build-bgdb")?;
    term.write_line(&format!("common = {:?}", &common))?;
    term.write_line(&format!("args = {:?}", &args))?;
    term.write_line("")?;

    // Create final list of input paths (expand `@file.tsv`)
    let mut input_tsv_paths = Vec::new();
    for path in &args.input_tsv {
        if let Some(path) = path.strip_prefix('@') {
            let path = shellexpand::tilde(&path);
            let lines = read_lines(path.into_owned())?;
            for line in lines {
                input_tsv_paths.push(line?.clone());
            }
        } else {
            let path = shellexpand::tilde(&path);
            input_tsv_paths.push(path.into_owned())
        }
    }
    term.write_line(&format!(
        "final input TSV file list (#: {}): {:?}",
        input_tsv_paths.len(),
        &input_tsv_paths
    ))?;

    print_rss_now(term)?;

    // Read all input files and write all records by chromosome and SV type
    let tmp_dir = tempdir::TempDir::new("tmp.vsw_sv_bgdb")?;
    term.write_line(&format!("using tmpdir={:?}", &tmp_dir))?;
    split_input_by_chrom_and_sv_type(&tmp_dir, input_tsv_paths, term)?;

    // Read the output of the previous step by chromosome and SV type, perform overlapping
    // and merge such "compressed" data set to the final output file.
    merge_split_files(&tmp_dir, args, term)?;

    Ok(())
}
