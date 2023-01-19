pub mod dbrecords;
pub mod recordio;

use clap::Parser;

use self::recordio::load_sv_records;
use crate::common::Args as CommonArgs;

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

pub(crate) fn run(
    term: &console::Term,
    common: &CommonArgs,
    args: &Args,
) -> Result<(), anyhow::Error> {
    term.write_line("Starting sv-query")?;
    term.write_line(&format!("common = {:?}", &common))?;
    term.write_line(&format!("args = {:?}", &args))?;
    term.write_line("")?;

    let _sv_records = load_sv_records(term, &args.db_base_dir)?;

    Ok(())
}
