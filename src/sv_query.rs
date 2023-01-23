pub mod dbrecords;
pub mod interpreter;
pub mod recordio;
pub mod schema;

use std::{
    collections::{HashMap, HashSet},
    time::Instant,
};

use anyhow::anyhow;
use clap::Parser;
use rust_htslib::bcf::{self, Read};
use thousands::Separable;

use self::{
    recordio::{build_chrom_map, BgRecordsByChrom},
    schema::{CallInfo, StructuralVariant, SvSubType, SvType},
};
use crate::common::Args as CommonArgs;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Base directory path for databases
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

    let sv_records = BgRecordsByChrom::load_from_dir(term, &args.db_base_dir)?;
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

    let samples = bcf_reader
        .header()
        .samples()
        .iter()
        .map(|sample| {
            std::str::from_utf8(sample)
                .expect("invalid UTF-8 in sample name")
                .to_owned()
        })
        .collect::<Vec<String>>();

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

        let (sv_type, sv_sub_type): (SvType, SvSubType) =
            if let Some(svtype) = record.info(b"SVTYPE").string()? {
                if svtype.len() < 1 {
                    Err(anyhow!("Empty SVTYPE tag"))
                } else {
                    let sv_sub_type_str = std::str::from_utf8(svtype[0])?;
                    let sv_type_str = sv_sub_type_str.split(':').next().unwrap();
                    Ok((
                        serde_json::from_str(&format!("\"{}\"", &sv_type_str))?,
                        serde_json::from_str(&format!("\"{}\"", &sv_sub_type_str))?,
                    ))
                }
            } else {
                Err(anyhow!("Could not decode INFO/SVTYPE"))
            }?;
        let strand_orientation = if let Some(ct) = record.info(b"CT").string()? {
            if ct.len() < 1 {
                Err(anyhow!("Empty INFO/CT tag"))
            } else {
                let ct_str = std::str::from_utf8(ct[0])?;
                Ok(serde_json::from_str(&format!("\"{}\"", &ct_str))?)
            }
        } else {
            Err(anyhow!("Could not decode INFO/CT"))
        }?;

        let call_info = {
            let mut call_infos = HashMap::new();
            for (i, sample) in samples.iter().enumerate() {
                let genotype = record
                    .genotypes()
                    .map_or(None, |genotypes| Some(format!("{}", genotypes.get(i))));
                let quality = record
                    .format(b"GQ")
                    .integer()
                    .map_or(None, |quality| Some(quality[0][0] as f32));
                let paired_end_ref = record
                    .format(b"DR")
                    .integer()
                    .map_or(None, |dr| Some(dr[0][0] as u32));
                let paired_end_var = record
                    .format(b"DV")
                    .integer()
                    .map_or(None, |dv| Some(dv[0][0] as u32));
                let split_read_ref = record
                    .format(b"RR")
                    .integer()
                    .map_or(None, |rr| Some(rr[0][0] as u32));
                let split_read_var = record
                    .format(b"RV")
                    .integer()
                    .map_or(None, |rv| Some(rv[0][0] as u32));
                call_infos.insert(
                    sample.clone(),
                    CallInfo {
                        genotype,
                        quality,
                        paired_end_cov: paired_end_ref.map_or(None, |paired_end_ref| {
                            paired_end_var.map_or(None, |paired_end_var| {
                                Some(paired_end_ref + paired_end_var)
                            })
                        }),
                        paired_end_var: paired_end_var,
                        split_read_cov: split_read_ref.map_or(None, |split_read_ref| {
                            split_read_var.map_or(None, |split_read_var| {
                                Some(split_read_ref + split_read_var)
                            })
                        }),
                        split_read_var: split_read_var,
                        copy_number: None,
                        average_normalized_cov: None,
                        point_count: None,
                        average_mapping_quality: None,
                    },
                );
            }
            call_infos
        };

        let sv = if sv_type == SvType::Bnd {
            let chr2 = if let Some(chr2) = record.info(b"CHR2").string()? {
                if chr2.len() < 1 {
                    Err(anyhow!("Empty INFO/CHR2"))
                } else {
                    Ok(std::str::from_utf8(chr2[0])?.to_owned())
                }
            } else {
                Err(anyhow!("Could not decode INFO/CHR2"))
            }?;
            let pos2 = if let Some(pos2) = record.info(b"POS2").integer()? {
                if pos2.len() < 1 {
                    Err(anyhow!("Empty INFO/POS2"))
                } else {
                    Ok(pos2[0] as u32)
                }
            } else {
                Err(anyhow!("Could not decode INFO/POS2"))
            }?;

            StructuralVariant {
                chrom: chrom_name.clone(),
                pos: record.pos() as u32,
                sv_type,
                sv_sub_type,
                chrom2: Some(chr2),
                end: pos2,
                strand_orientation: Some(strand_orientation),
                call_info: call_info,
            }
        } else {
            let end = if let Some(end) = record.info(b"END").integer()? {
                if end.len() < 1 {
                    Err(anyhow!("Empty INFO/END"))
                } else {
                    Ok(end[0] as u32)
                }
            } else {
                Err(anyhow!("Could not decode INFO/END"))
            }?;

            StructuralVariant {
                chrom: chrom_name.clone(),
                pos: record.pos() as u32,
                sv_type,
                sv_sub_type,
                chrom2: None,
                end: end,
                strand_orientation: Some(strand_orientation),
                call_info: call_info,
            }
        };

        let _overlaps = sv_records.count_overlaps(&chrom_map, &sv);
        // println!("{:?}", &sv);
        // println!("  `--> {:?}", &overlaps);
        // println!("--");

        count += 1;
    }
    term.write_line(&format!(
        "-- total time spent querying for {} records: {:?}",
        count.separate_with_commas(),
        before_query.elapsed()
    ))?;

    Ok(())
}
