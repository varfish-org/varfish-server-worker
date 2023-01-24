//! Implementation of `sv build-bgdb` command

pub mod recordio;

use std::{
    fs::File,
    io::{self, BufRead},
    path::Path,
    time::Instant,
};

use clap::Parser;
use flate2::read::GzDecoder;
use thousands::Separable;

use crate::common::Args as CommonArgs;

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
            let _record: FileRecord = result?;
            count_records += 1;
        }
        term.write_line(&format!(
            "-- total time spent parsing {} record: {:?}",
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

    Ok(())
}
