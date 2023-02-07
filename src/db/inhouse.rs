//! Code supporting the `db mk-inhouse` sub command.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use bio::data_structures::interval_tree::IntervalTree;
use clap::{command, Parser};
use serde_json::to_writer;
use serde_jsonlines::JsonLinesReader;
use strum::IntoEnumIterator;
use thousands::Separable;
use tracing::{debug, info};

use crate::{
    common::{
        build_chrom_map, md5sum, open_read_maybe_gz, open_write_maybe_gz, read_lines, sha256sum,
        trace_rss_now, CHROMS,
    },
    db::{
        conf::{DbDef, GenomeRelease, Top},
        to_bin::{self, vardbs::InputFileType},
    },
    sv::query::schema::SvType,
};

use super::compile::ArgGenomeRelease;

/// Code for reading the input file.
pub mod input {
    use serde::{de::IntoDeserializer, Deserialize, Deserializer, Serialize};

    use crate::sv::query::schema::{StrandOrientation, SvType};

    /// Representation of the fields from the `StructuralVariant` table from VarFish Server
    /// that we need for building the background records.
    #[derive(Debug, Deserialize, Serialize, Clone)]
    pub struct Record {
        /// genome build
        pub release: String,
        /// chromosome name
        pub chromosome: String,
        /// UCSC bin
        pub bin: u32,
        /// start position, 1-based
        pub start: u32,
        /// chromosome2 name
        pub chromosome2: String,
        /// end position, 1-based
        pub end: u32,
        /// paired-end orientation
        #[serde(deserialize_with = "from_varfish_pe_orientation")]
        pub pe_orientation: StrandOrientation,
        /// SV type of the record
        #[serde(deserialize_with = "from_varfish_sv_type")]
        pub sv_type: SvType,
        /// number of hom. alt. carriers
        pub num_hom_alt: u32,
        /// number of hom. ref. carriers
        pub num_hom_ref: u32,
        /// number of het. carriers
        pub num_het: u32,
        /// number of hemi. alt. carriers
        pub num_hemi_alt: u32,
        /// number of hemi. ref. carriers
        pub num_hemi_ref: u32,
    }

    /// Deserialize "sv_type" from VarFish database.
    ///
    /// This function will strip everything after the first underscore.
    fn from_varfish_sv_type<'de, D>(deserializer: D) -> Result<SvType, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        let end = s.find('_').unwrap_or(s.len());
        SvType::deserialize(s[..end].into_deserializer())
    }

    /// Deserialize "pe_orientation" from VarFish database.
    ///
    /// This function will convert `"."` to `StrandOrientation::NotApplicable`
    fn from_varfish_pe_orientation<'de, D>(deserializer: D) -> Result<StrandOrientation, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        if s.eq(".") {
            Ok(StrandOrientation::NotApplicable)
        } else {
            StrandOrientation::deserialize(s.into_deserializer())
        }
    }
}

/// Code for writing the output file.
pub mod output {
    use serde::{Deserialize, Serialize};

    use crate::sv::query::schema::{StrandOrientation, SvType};

    use super::input::Record as InputRecord;

    /// Representation of the fields for the in-house background database.
    #[derive(Debug, Deserialize, Serialize, Clone)]
    pub struct Record {
        /// chromosome name
        pub chromosome: String,
        /// start position, 1-based
        pub begin: u32,
        /// chromosome2 name
        pub chromosome2: String,
        /// end position, 1-based
        pub end: u32,
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

    impl Record {
        /// Compute reciprocal overlap between `self` and `other`.
        pub fn overlap(&self, other: &Record) -> f32 {
            let s1 = if self.begin > 0 { self.begin - 1 } else { 0 };
            let e1 = self.end + 1;
            let s2 = if other.begin > 0 { other.begin - 1 } else { 0 };
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

        pub fn merge_into(&mut self, other: &Record) {
            self.carriers += other.carriers;
            self.carriers_het += other.carriers_het;
            self.carriers_hom += other.carriers_hom;
            self.carriers_hemi += other.carriers_hemi;
        }

        pub fn from_db_record(record: InputRecord) -> Self {
            Record {
                chromosome: record.chromosome,
                begin: record.start,
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
    genome_release: ArgGenomeRelease,
) -> Result<(), anyhow::Error> {
    info!("parse all input files and split them up");
    let genome_release = genome_release.to_string().to_lowercase();
    let mut tmp_files = create_tmp_files(tmp_dir)?;
    let chrom_map = build_chrom_map();
    let before_parsing = Instant::now();
    let mut count_files = 0;
    for path in &input_tsv_paths {
        debug!("parsing {:?}", &path);

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(open_read_maybe_gz(path)?);
        let before_parsing = Instant::now();
        let mut count_records = 0;
        for result in rdr.deserialize() {
            let record: input::Record = result?;

            if record.release.to_lowercase() != genome_release {
                return Err(anyhow::anyhow!(
                    "Unexpected release {}, expected {}",
                    record.release.to_lowercase(),
                    &genome_release
                ));
            }

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
        trace_rss_now();
        debug!(
            "total time spent parsing {} records: {:?}",
            count_records.separate_with_commas(),
            before_parsing.elapsed()
        );

        count_files += 1;
    }
    info!(
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
/// The idea to merge here is to get rid of large stacks of SVs with a reciprocal overlap
/// that is more strict than the 0.75 that is generally used for querying.  We merge with
/// existing clusters with the reciprocal overlap is >=0.8 for all members.
fn merge_to_out(
    args: &Args,
    reader: &mut BufReader<File>,
    writer: &mut csv::Writer<impl Write>,
) -> Result<usize, anyhow::Error> {
    let mut clusters: Vec<Vec<usize>> = vec![];
    let mut tree: IntervalTree<i32, usize> = IntervalTree::new();
    let mut records: Vec<output::Record> = Vec::new();

    // Read in all records and perform the "merge compression"
    let mut reader = JsonLinesReader::new(reader);
    while let Ok(Some(record)) = reader.read::<input::Record>() {
        let record = output::Record::from_db_record(record);
        let slack = match record.sv_type {
            SvType::Bnd => args.slack_bnd,
            SvType::Ins => args.slack_ins,
            _ => 0,
        };
        let query = ((record.begin - 1 - slack) as i32)..((record.end + slack) as i32);
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
                ((record.begin - 1) as i32)..(record.end as i32),
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
    info!("merge all files to {:?}...", &path_output_tsv);
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(open_write_maybe_gz(path_output_tsv)?);

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
            debug!("reading from {}", &filename);
            let mut reader = BufReader::new(File::open(path)?);
            out_records += merge_to_out(args, &mut reader, &mut writer)?;
        }
    }
    info!("wrote a total of {} records", out_records);

    writer.flush()?;

    Ok(())
}

/// Command line arguments for `db build` sub command.
#[derive(Parser, Debug)]
#[command(about = "Build inhouse database", long_about = None)]
pub struct Args {
    /// Genome build to use in the build.
    #[arg(long, value_enum, default_value_t = ArgGenomeRelease::Grch37)]
    pub genome_release: ArgGenomeRelease,
    /// Path to `varfish-server-worker` directory (output).
    #[arg(long, required = true)]
    pub path_worker_db: PathBuf,
    /// Input files to cluster, prefix with `@` to file with line-wise paths.
    #[arg(required = true)]
    pub path_input_tsvs: Vec<String>,
    /// Minimal reciprocal overlap to use (slightly more strict that the normal query value of 0.75).
    #[arg(long, default_value_t = 0.8)]
    pub min_overlap: f32,
    /// Padding to use for BNDs
    #[arg(long, default_value_t = 50)]
    pub slack_bnd: u32,
    /// Padding to use for INS
    #[arg(long, default_value_t = 50)]
    pub slack_ins: u32,
}

/// Main entry point for the `sv bg-db-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting `db mk-inhouse`");
    info!("  common_args = {:?}", &common_args);
    info!("  args = {:?}", &args);

    // Create final list of input paths (expand `@file.tsv`)
    let mut input_tsv_paths = Vec::new();
    for input_tsv in &args.path_input_tsvs {
        if let Some(path) = input_tsv.strip_prefix('@') {
            let path = shellexpand::tilde(&path);
            let lines = read_lines(path.into_owned())?;
            for line in lines {
                input_tsv_paths.push(line?.clone());
            }
        } else {
            let path = shellexpand::tilde(&input_tsv);
            input_tsv_paths.push(path.into_owned())
        }
    }
    debug!(
        "final input TSV file list is (#: {}): {:?}",
        input_tsv_paths.len(),
        &input_tsv_paths
    );

    trace_rss_now();

    // Read all input files and write all records by chromosome and SV type
    let tmp_dir = tempdir::TempDir::new("vfw")?;
    debug!("using tmpdir={:?}", &tmp_dir);
    split_input_by_chrom_and_sv_type(&tmp_dir, input_tsv_paths, args.genome_release)?;

    let base_path = args
        .path_worker_db
        .join("vardbs")
        .join(args.genome_release.to_string())
        .join("strucvar");
    let path_out_tsv = base_path.join("inhouse.tsv.gz");

    // Read the output of the previous step by chromosome and SV type, perform overlapping
    // and merge such "compressed" data set to the final output file.
    merge_split_files(&tmp_dir, args, &path_out_tsv)?;

    // Update database configuration.
    info!("Reading database configuration...");
    let genome_release = match args.genome_release {
        ArgGenomeRelease::Grch37 => GenomeRelease::Grch37,
        ArgGenomeRelease::Grch38 => GenomeRelease::Grch38,
    };
    let conf_toml = std::fs::read_to_string(args.path_worker_db.join("conf.toml"))?;
    let mut conf: Top = toml::from_str(&conf_toml)?;
    conf.vardbs[genome_release].strucvar.inhouse = Some(DbDef {
        path: path_out_tsv
            .strip_prefix(&args.path_worker_db)
            .unwrap()
            .to_str()
            .unwrap()
            .to_owned(),
        md5: Some(md5sum(&path_out_tsv)?),
        sha256: Some(sha256sum(&path_out_tsv)?),
        ..Default::default()
    });

    // Convert in-house database to binary.
    info!("Converting to binary...");
    let path_out_bin = base_path.join("inhouse.bin");
    if let Some(inhouse) = &mut conf.vardbs[genome_release].strucvar.inhouse {
        to_bin::vardbs::convert_to_bin(
            &path_out_tsv,
            &path_out_bin,
            &args.path_worker_db,
            InputFileType::InhouseDb,
            inhouse,
        )?;
    }

    info!("Writing database configuration...");
    std::fs::write(
        args.path_worker_db.join("conf.toml"),
        toml::to_string_pretty(&conf)?,
    )?;

    Ok(())
}
