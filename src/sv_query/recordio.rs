use std::{collections::HashMap, fs::File, path::Path, time::Instant};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;

use flate2::read::GzDecoder;
use serde::de::DeserializeOwned;
use thousands::Separable;

use crate::common::{build_chrom_map, print_rss_now, CHROMS};

use super::{
    dbrecords::{
        self, reciprocal_overlap, BeginEnd, ChromosomeCoordinate, Count, SvOverlapCounts,
        ToInMemory,
    },
    interpreter::{BND_SLACK, INS_SLACK},
    schema::{CaseQuery, StructuralVariant, SvType},
};

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
        if let Some(mem_record) = record.to_in_memory()? {
            rec_by_contig[idx].push(mem_record);
        }
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
) -> Result<Vec<ArrayBackedIntervalTree<u32, u32>>, anyhow::Error> {
    term.write_line("Building trees...")?;
    let mut trees: Vec<ArrayBackedIntervalTree<u32, u32>> = Vec::new();
    let before_building = Instant::now();
    for (i, contig_records) in rec_by_contig.iter().enumerate() {
        let before_tree = Instant::now();
        let mut tree = ArrayBackedIntervalTree::new();
        for (i, record) in contig_records.iter().enumerate() {
            tree.insert(record.begin()..record.end(), i as u32);
        }
        tree.index();
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

/// This struct bundles the background database records by chromosome.
pub struct BgRecordsByChrom {
    pub bg_sv_records: Vec<Vec<dbrecords::bg_sv::Record>>,
    pub bg_sv_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub gnomad_sv_records: Vec<Vec<dbrecords::gnomad_sv::Record>>,
    pub gnomad_sv_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub g1k_sv_records: Vec<Vec<dbrecords::g1k_sv::Record>>,
    pub g1k_sv_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub dbvar_records: Vec<Vec<dbrecords::dbvar::Record>>,
    pub dbvar_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub dgv_records: Vec<Vec<dbrecords::dgv::Record>>,
    pub dgv_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub dgv_gs_records: Vec<Vec<dbrecords::dgv_gs::Record>>,
    pub dgv_gs_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
    pub exac_cnv_records: Vec<Vec<dbrecords::exac_cnv::Record>>,
    pub exac_cnv_trees: Vec<ArrayBackedIntervalTree<u32, u32>>,
}

impl BgRecordsByChrom {
    pub fn load_from_dir(
        term: &console::Term,
        db_base_dir: &str,
    ) -> Result<BgRecordsByChrom, anyhow::Error> {
        let chrom_map = build_chrom_map();

        let before_parsing = Instant::now();
        print_rss_now(term)?;

        let bg_sv_path = Path::new(&db_base_dir)
            .join("bg-inhouse")
            .join("varfish-sv.tsv.gz");
        let bg_sv_records = load_bg_sv_records::<
            dbrecords::bg_sv::Record,
            dbrecords::bg_sv::FileRecord,
        >(&bg_sv_path, term, &chrom_map)?;
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

        let g1k_sv_path = Path::new(&db_base_dir)
            .join("bg-public")
            .join("g1k-sv.tsv.gz");
        let g1k_sv_records = load_bg_sv_records::<
            dbrecords::g1k_sv::Record,
            dbrecords::g1k_sv::FileRecord,
        >(&g1k_sv_path, term, &chrom_map)?;
        let g1k_sv_trees = build_bg_sv_tree(term, &g1k_sv_records)?;
        print_rss_now(term)?;

        let dbvar_path = Path::new(&db_base_dir)
            .join("bg-public")
            .join("dbvar-sv.tsv.gz");
        let dbvar_records = load_bg_sv_records::<
            dbrecords::dbvar::Record,
            dbrecords::dbvar::FileRecord,
        >(&dbvar_path, term, &chrom_map)?;
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

        Ok(BgRecordsByChrom {
            bg_sv_records,
            bg_sv_trees,
            gnomad_sv_records,
            gnomad_sv_trees,
            g1k_sv_records,
            g1k_sv_trees,
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

    pub fn count_overlaps(
        &self,
        chrom_map: &HashMap<String, usize>,
        query: &CaseQuery,
        sv: &StructuralVariant,
    ) -> SvOverlapCounts {
        let chrom_idx = *chrom_map.get(&sv.chrom).unwrap();
        let range = if sv.sv_type == SvType::Ins {
            sv.pos.saturating_sub(INS_SLACK)..sv.pos.saturating_add(INS_SLACK)
        } else if sv.sv_type == SvType::Bnd {
            sv.pos.saturating_sub(BND_SLACK)..sv.pos.saturating_add(BND_SLACK)
        } else {
            sv.pos..sv.end
        };

        let inhouse_carriers: u32 = self.bg_sv_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.bg_sv_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_inhouse_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_inhouse_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let gnomad_carriers: u32 = self.gnomad_sv_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.gnomad_sv_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_gnomad_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_gnomad_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let g1k_alleles: u32 = self.g1k_sv_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.g1k_sv_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_g1k_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_g1k_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let dbvar_carriers: u32 = self.dbvar_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.dbvar_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_dbvar_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_dbvar_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let dgv_carriers: u32 = self.dgv_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.dgv_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_dgv_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_dgv_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let dgv_gs_carriers: u32 = self.dgv_gs_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.dgv_gs_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_dgv_gs_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_dgv_gs_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;
        let exac_carriers: u32 = self.exac_cnv_trees[chrom_idx]
            .find(range.clone())
            .iter()
            .map(|e| &self.exac_cnv_records[chrom_idx][*e.data() as usize])
            .filter(|record| record.sv_type.is_compatible(sv.sv_type))
            .filter(|record| {
                query.svdb_exac_enabled
                    && (record.sv_type == SvType::Ins
                        || record.sv_type == SvType::Bnd
                        || query.svdb_exac_min_overlap.map_or(true, |min_overlap| {
                            (reciprocal_overlap(*record, &range)) >= min_overlap
                        }))
            })
            .map(|record| record.count())
            .sum::<usize>() as u32;

        SvOverlapCounts {
            dgv_carriers,
            dgv_gs_carriers,
            gnomad_carriers,
            exac_carriers,
            dbvar_carriers,
            g1k_alleles,
            inhouse_carriers,
        }
    }
}
