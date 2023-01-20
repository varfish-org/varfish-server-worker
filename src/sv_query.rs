pub mod dbrecords;
pub mod recordio;
pub mod schema;

use std::{collections::HashSet, time::Instant};

use anyhow::anyhow;
use clap::Parser;
use rust_htslib::bcf::{self, Read};
use thousands::Separable;

use self::recordio::{build_chrom_map, load_sv_records};
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

pub fn run(term: &console::Term, common: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line("Starting sv-query")?;
    term.write_line(&format!("common = {:?}", &common))?;
    term.write_line(&format!("args = {:?}", &args))?;
    term.write_line("")?;

    let sv_records = load_sv_records(term, &args.db_base_dir)?;
    let chrom_map = build_chrom_map();

    let mut bcf_reader = bcf::Reader::from_path(&args.input_vcf)?;

    let chrom_names = {
        let mut chrom_names = Vec::new();
        for i in 0..bcf_reader.header().contig_count() {
            let contig_name = std::str::from_utf8(bcf_reader.header().rid2name(i)?)?;
            chrom_names.push(contig_name.to_owned());
        }
        chrom_names
    };

    let mut unexpected_contigs: HashSet<String> = HashSet::new();

    term.write_line("Starting queries...")?;
    let mut count = 0;
    let before_query = Instant::now();
    for record_result in bcf_reader.records() {
        let record = record_result?;
        if record.rid().is_none() {
            continue;
        }

        let rid = record.rid().unwrap() as usize;
        let chrom_name = chrom_names.get(rid).unwrap();
        let chrom_idx = if let Some(chrom_idx) = chrom_map.get(chrom_name) {
            *chrom_idx
        } else {
            if !unexpected_contigs.contains(chrom_name) {
                term.write_line(&format!(
                    "Will skip records on contig {} (this is only printed once)",
                    &chrom_name
                ))?;
                unexpected_contigs.insert(chrom_name.to_owned());
            }
            chrom_names.len()
        };
        if chrom_idx == chrom_names.len() {
            continue; // skip
        }

        let begin = record.pos() as u32;
        let end = if let Some(end) = record.info(b"END").integer()? {
            if end.len() < 1 {
                Err(anyhow!("Empty END tag"))
            } else {
                Ok(end[0] as u32)
            }
        } else {
            Err(anyhow!("Could not decode INFO/end"))
        }?;

        sv_records.bg_sv_trees[chrom_idx].find(begin..end);
        sv_records.gnomad_sv_trees[chrom_idx].find(begin..end);
        sv_records.dbvar_trees[chrom_idx].find(begin..end);
        sv_records.dgv_trees[chrom_idx].find(begin..end);
        sv_records.dgv_gs_trees[chrom_idx].find(begin..end);
        sv_records.exac_cnv_trees[chrom_idx].find(begin..end);
        count += 1;
    }
    term.write_line(&format!(
        "-- total time spent querying for {} records: {:?}",
        count.separate_with_commas(),
        before_query.elapsed()
    ))?;

    Ok(())
}
