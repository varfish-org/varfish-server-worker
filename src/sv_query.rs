//! Implementation of `sv query` command.

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
    interpreter::QueryInterpreter,
    recordio::BgRecordsByChrom,
    schema::{
        CallInfo, CaseQuery, Database, GenotypeChoice, GenotypeCriteria, StructuralVariant,
        SvSubType, SvType,
    },
};
use crate::common::{build_chrom_map, Args as CommonArgs};

#[derive(Parser, Debug)]
#[command(author, version, about = "Run query for SVs", long_about = None)]
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

pub fn build_query(samples: &[String]) -> CaseQuery {
    let genotype_criteria = vec![
        // CNVs -- variant
        GenotypeCriteria {
            select_sv_sub_type: SvSubType::all(),
            min_srpr_var: Some(6),
            ..GenotypeCriteria::new(GenotypeChoice::Variant)
        },
        // CNVs -- non-variant
        GenotypeCriteria {
            select_sv_sub_type: SvSubType::all(),
            max_srpr_var: Some(5),
            ..GenotypeCriteria::new(GenotypeChoice::NonVariant)
        },
    ];

    CaseQuery {
        svdb_gnomad_enabled: true,
        svdb_gnomad_min_overlap: Some(0.8),
        svdb_gnomad_max_carriers: Some(10),
        svdb_exac_enabled: true,
        svdb_exac_min_overlap: Some(0.8),
        svdb_exac_max_carriers: Some(10),
        svdb_dbvar_enabled: true,
        svdb_dbvar_min_overlap: Some(0.8),
        svdb_dbvar_max_carriers: Some(20),
        svdb_g1k_enabled: true,
        svdb_g1k_min_overlap: Some(0.8),
        svdb_g1k_max_alleles: Some(10),
        svdb_inhouse_enabled: true,
        svdb_inhouse_min_overlap: Some(0.8),
        svdb_inhouse_max_carriers: Some(10),

        sv_size_min: Some(500),
        sv_types: SvType::all(),
        sv_sub_types: SvSubType::all(),

        genotype: HashMap::from_iter(
            samples
                .iter()
                .map(|sample| (sample.clone(), GenotypeChoice::Variant)),
        ),
        genotype_criteria,

        ..CaseQuery::new(Database::Refseq)
    }
}

pub fn run(term: &console::Term, common: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line("Starting sv query")?;
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

    term.write_line("Buildling query...")?;
    let query = build_query(&samples);
    let interpreter = QueryInterpreter::new(query.clone());

    term.write_line("Starting queries...")?;
    let mut count_total = 0;
    let mut count_passes = 0;
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
                        paired_end_cov: paired_end_ref.and_then(|paired_end_ref| {
                            paired_end_var.map(|paired_end_var| paired_end_ref + paired_end_var)
                        }),
                        paired_end_var,
                        split_read_cov: split_read_ref.and_then(|split_read_ref| {
                            split_read_var.map(|split_read_var| split_read_ref + split_read_var)
                        }),
                        split_read_var,
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
                call_info,
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
                end,
                strand_orientation: Some(strand_orientation),
                call_info,
            }
        };

        count_total += 1;
        if interpreter.passes(&sv, |sv| sv_records.count_overlaps(&chrom_map, &query, sv))? {
            count_passes += 1;
            if count_passes <= 100 {
                println!(
                    "{:?}\n{:?}\n--",
                    &sv,
                    &sv_records.count_overlaps(&chrom_map, &query, &sv)
                );
            }
        }
    }
    term.write_line(&format!(
        "-- total time spent querying for {} (passing: {}) records: {:?}",
        count_total.separate_with_commas(),
        count_passes.separate_with_commas(),
        before_query.elapsed()
    ))?;

    Ok(())
}
