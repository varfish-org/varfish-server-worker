//! Implementation of `sv build-bgdb` command

use std::{
    fs::File,
    io::{self, BufRead},
    path::Path,
};

use clap::Parser;

use crate::common::Args as CommonArgs;

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

    Ok(())
}
