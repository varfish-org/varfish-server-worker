//! Command line interface

use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use bio::data_structures::interval_tree::IntervalTree;
use clap::{command, Parser};
use mehari::common::open_read_maybe_gz;
use noodles_vcf as vcf;
use serde_json::to_writer;
use serde_jsonlines::JsonLinesReader;
use strum::IntoEnumIterator;
use thousands::Separable;

use crate::{
    common::{
        build_chrom_map, open_write_maybe_gz, read_lines, trace_rss_now, GenomeRelease, CHROMS,
    },
    strucvars::query::schema::SvType,
};

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
    input_vcf_paths: Vec<String>,
    genome_release: GenomeRelease,
) -> Result<(), anyhow::Error> {
    tracing::info!("parse all input files and split them up");
    let mut tmp_files = create_tmp_files(tmp_dir)?;
    let chrom_map = build_chrom_map();
    let before_parsing = Instant::now();
    let mut count_files = 0;
    for path_input in &input_vcf_paths {
        tracing::debug!("parsing {:?}", &path_input);
        let mut input_reader = open_read_maybe_gz(path_input).map_err(|e| {
            anyhow::anyhow!("could not open file {} for reading: {}", path_input, e)
        })?;
        let mut input_reader = vcf::Reader::new(&mut input_reader);
        let input_header = input_reader.read_header()?;

        let (pedigree, _) = crate::common::extract_pedigree_and_case_uuid(&input_header)?;
        let mut prev = std::time::Instant::now();
        let records = input_reader.records(&input_header);
        let before_parsing = Instant::now();
        let mut count_records = 0;
        for input_record in records {
            let input_record = super::output::Record::from_vcf(
                &input_record?,
                &input_header,
                genome_release,
                &pedigree,
            )?;

            let chrom_no = *chrom_map
                .get(&input_record.chromosome)
                .expect("unknown chromosome");
            let sv_type = input_record.sv_type;
            let mut tmp_file = tmp_files
                .get_mut(&(chrom_no, sv_type))
                .expect("no file for chrom/sv_type");
            to_writer(&mut tmp_file, &input_record)?;
            tmp_file.write_all(&[b'\n'])?;

            // Write out progress indicator every 60 seconds.
            if prev.elapsed().as_secs() >= 60 {
                tracing::info!("at {}:{}", &input_record.chromosome, input_record.begin + 1);
                prev = std::time::Instant::now();
            }

            count_records += 1;
        }

        trace_rss_now();
        tracing::debug!(
            "total time spent parsing {} records: {:?}",
            count_records.separate_with_commas(),
            before_parsing.elapsed()
        );

        count_files += 1;
    }
    tracing::info!(
        "total time spent parsing {} files: {:?}",
        count_files.separate_with_commas(),
        before_parsing.elapsed()
    );
    trace_rss_now();
    for (_, mut f) in tmp_files.drain() {
        f.get_mut().sync_all()?
    }
    Ok(())
}

/// Read in all records from `reader`, merge overlapping ones.
///
/// The idea to merge here is to get rid of large stacks of SVs with a
/// reciprocal overlap that is more strict than the 0.75 that is generally used
/// for querying.  We merge with existing clusters with the reciprocal overlap
/// is >=0.8 for all members.
fn merge_to_out(
    args: &Args,
    reader: &mut BufReader<File>,
    writer: &mut csv::Writer<impl Write>,
) -> Result<usize, anyhow::Error> {
    let mut clusters: Vec<Vec<usize>> = vec![];
    let mut tree: IntervalTree<i32, usize> = IntervalTree::new();
    let mut records: Vec<super::output::Record> = Vec::new();

    // Read in all records and perform the "merge compression"
    let mut reader = JsonLinesReader::new(reader);
    while let Ok(Some(record)) = reader.read::<super::output::Record>() {
        let begin = match record.sv_type {
            SvType::Bnd => record.begin - 1 - args.slack_bnd,
            SvType::Ins => record.begin - 1 - args.slack_ins,
            _ => record.begin,
        };
        let end = match record.sv_type {
            SvType::Bnd => record.begin + args.slack_bnd,
            SvType::Ins => record.begin + args.slack_ins,
            _ => record.end,
        };
        let query = begin..end;
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
                clusters[cluster_idx].push(records.len());
                found_any_cluster = true;
                break;
            }
        }
        if !found_any_cluster {
            // create new cluster
            tree.insert(
                match record.sv_type {
                    SvType::Bnd | SvType::Ins => (record.begin - 1)..record.begin,
                    _ => (record.begin - 1)..record.end,
                },
                clusters.len(),
            );
            clusters.push(vec![records.len()]);
        }
        // always register the record
        records.push(record);
    }

    trace_rss_now();

    // Sort the cluster representatives by start coordinate.
    let mut sorted_idxs = vec![0; clusters.len()];
    for (i, sorted_idx) in sorted_idxs.iter_mut().enumerate() {
        *sorted_idx = i;
    }
    sorted_idxs.sort_by(|a, b| {
        (records[clusters[*a][0]].begin, records[clusters[*a][0]].end)
            .partial_cmp(&(records[clusters[*b][0]].begin, records[clusters[*b][0]].end))
            .unwrap()
    });

    // Finally, write out all records in sorted order
    let mut out_records = 0;
    for cluster in clusters {
        let mut out_record = records[cluster[0]].clone();
        for record_id in &cluster[1..] {
            out_record.merge_into(&records[*record_id]);
        }
        out_records += 1;
        writer.serialize(&out_record)?;
    }

    Ok(out_records)
}

/// Perform (chrom, sv_type) wise merging of records in temporary files.
fn merge_split_files(
    tmp_dir: &tempdir::TempDir,
    args: &Args,
    path_output_tsv: &Path,
) -> Result<(), anyhow::Error> {
    tracing::info!("merge all files to {:?}...", &path_output_tsv);
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(open_write_maybe_gz(path_output_tsv).map_err(|e| {
            anyhow::anyhow!("Cannot open {:?} for writing: {:?}", &path_output_tsv, e)
        })?);

    // Write header as comment.
    writer.write_record([
        "#chromosome",
        "begin",
        "chromosome2",
        "end",
        "pe_orientation",
        "sv_type",
        "carriers",
        "carriers_het",
        "carriers_hom",
        "carriers_hemi",
    ])?;

    let mut out_records = 0;
    for chrom in CHROMS {
        for sv_type in SvType::iter() {
            let filename = format!("records.chr{}.{:?}.tsv", *chrom, sv_type);
            let path = tmp_dir.path().join(&filename);
            tracing::debug!("reading from {}", &filename);
            let mut reader = BufReader::new(File::open(path)?);
            out_records += merge_to_out(args, &mut reader, &mut writer)?;
        }
    }
    tracing::info!("wrote a total of {} records", out_records);

    writer.flush()?;

    Ok(())
}

/// Command line arguments for `db build` sub command.
#[derive(Parser, Debug)]
#[command(about = "Build inhouse database", long_about = None)]
pub struct Args {
    /// Genome build to use in the build.
    #[arg(long, value_enum, default_value_t = GenomeRelease::Grch37)]
    pub genome_release: GenomeRelease,
    /// Path to output TSV file.
    #[arg(long)]
    pub path_output: PathBuf,
    /// Input files to cluster, prefix with `@` to file with line-wise paths.
    #[arg(required = true)]
    pub path_input: Vec<String>,

    /// Minimal reciprocal overlap to use (slightly more strict that the normal
    /// query value of 0.75).
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Padding to use for BNDs
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: i32,
    /// Padding to use for INS
    #[arg(long, default_value_t = 50)]
    pub slack_ins: i32,
}

/// Main entry point for the `sv bg-db-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("Starting `strucvars aggregate`");
    tracing::info!("  common_args = {:?}", &common_args);
    tracing::info!("  args = {:?}", &args);

    // Create final list of input paths (expand `@file.tsv`)
    let mut input_vcf_paths = Vec::new();
    for input_vcf in &args.path_input {
        if let Some(path) = input_vcf.strip_prefix('@') {
            let path = shellexpand::tilde(&path);
            let lines = read_lines(path.into_owned())?;
            for line in lines {
                input_vcf_paths.push(line?.clone());
            }
        } else {
            let path = shellexpand::tilde(&input_vcf);
            input_vcf_paths.push(path.into_owned())
        }
    }
    tracing::debug!(
        "final input VCF file list is (#: {}): {:?}",
        input_vcf_paths.len(),
        &input_vcf_paths
    );

    trace_rss_now();

    // Read all input files and write all records by chromosome and SV type
    let tmp_dir = tempdir::TempDir::new("vfw")?;
    tracing::debug!("using tmpdir={:?}", &tmp_dir);
    split_input_by_chrom_and_sv_type(&tmp_dir, input_vcf_paths, args.genome_release)?;

    // Read the output of the previous step by chromosome and SV type, perform
    // overlapping and merge such "compressed" data set to the final output
    // file.
    tracing::info!("Merging to output TSV file...");
    merge_split_files(&tmp_dir, args, &args.path_output)?;
    tracing::info!("... done - don't forget to convert to binary");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{run, Args};
    use crate::common::{Args as CommonArgs, GenomeRelease};
    use clap_verbosity_flag::Verbosity;
    use temp_testdir::TempDir;

    #[test]
    fn run_smoke_gts_path() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            genome_release: GenomeRelease::Grch37,
            path_output: tmp_dir.join("out.tsv"),
            path_input: vec![String::from("tests/strucvars/aggregate/oneline.vcf")],
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
        };

        run(&common_args, &args)?;

        let output = std::fs::read_to_string(tmp_dir.join("out.tsv"))?;
        insta::assert_snapshot!(output);

        Ok(())
    }

    #[test]
    fn run_smoke_gts_twice() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            genome_release: GenomeRelease::Grch37,
            path_output: tmp_dir.join("out.tsv"),
            path_input: vec![
                String::from("tests/strucvars/aggregate/oneline.vcf"),
                String::from("tests/strucvars/aggregate/oneline.vcf"),
            ],
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
        };

        run(&common_args, &args)?;

        let output = std::fs::read_to_string(tmp_dir.join("out.tsv"))?;
        insta::assert_snapshot!(output);

        Ok(())
    }

    #[test]
    fn run_smoke_at_path() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            genome_release: GenomeRelease::Grch37,
            path_output: tmp_dir.join("out.tsv"),
            path_input: vec!["@tests/strucvars/aggregate/list.txt".to_string()],
            min_overlap: 0.8,
            slack_bnd: 50,
            slack_ins: 50,
        };

        run(&common_args, &args)?;

        let output = std::fs::read_to_string(tmp_dir.join("out.tsv"))?;
        insta::assert_snapshot!(output);

        Ok(())
    }
}
