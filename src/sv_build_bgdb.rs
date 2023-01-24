//! Implementation of `sv build-bgdb` command

pub mod recordio;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufWriter, Write},
    path::Path,
    time::Instant,
};

use clap::Parser;
use flate2::read::GzDecoder;
use serde_json::to_writer;
use strum::IntoEnumIterator;
use thousands::Separable;

use crate::{
    common::{build_chrom_map, print_rss_now, Args as CommonArgs, CHROMS},
    sv_query::schema::SvType,
};

use self::recordio::FileRecord;

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

    // Open one temporary file per chrom and SV type.
    let tmp_dir = tempdir::TempDir::new("tmp.vsw_sv_bgdb.")?;
    term.write_line(&format!("using tmpdir={:?}", &tmp_dir))?;
    let mut tmp_files = create_tmp_files(&tmp_dir)?;
    let chrom_map = build_chrom_map();

    // Read all input files, split by chromosome and SV type and write to temporary files.
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

    // Flush and close all temporary files.
    for (_, mut f) in tmp_files.drain() {
        f.get_mut().sync_all()?
    }

    Ok(())
}
