//! Code related to the "db genes build" subcommand.

use std::{collections::HashMap, io::BufRead, time::Instant};

use clap::Parser;
use indicatif::ProgressIterator;
use tracing::info;

use crate::{
    db::genes::data,
    pheno::prepare::{indicatif_style, VERSION},
};

use super::data::{acmg_sf, gnomad_constraints, hgnc, ncbi};

/// Command line arguments for `db genes build` sub command.
#[derive(Parser, Debug)]
#[command(author, version, about = "Build genes database", long_about = None)]
pub struct Args {
    /// Path to the directory with the output of the download directory.
    #[arg(long, required = true)]
    pub path_in_download: String,
    /// Path to output RocksDB.
    #[arg(long, required = true)]
    pub path_out_rocksdb: String,
}

/// Load ACMG SF list.
///
/// # Result
///
/// A map from HGNC ID to ACMG SF record.
fn load_acmg(path: String) -> Result<HashMap<String, acmg_sf::Record>, anyhow::Error> {
    info!("  loading ACMG SF list from {}", path);
    let mut result = HashMap::new();

    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    for record in reader.deserialize::<acmg_sf::Record>() {
        let record = record?;
        result.insert(record.hgnc_id.clone(), record);
    }

    Ok(result)
}

/// Load gnomAD constraints.
///
/// # Result
///
/// A map from ENSEMBL gene ID to gnomAD constraints record.
fn load_gnomad_constraints(
    path: String,
) -> Result<HashMap<String, gnomad_constraints::Record>, anyhow::Error> {
    info!("  loading gnomAD constraints from {}", path);
    let mut result = HashMap::new();

    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    for record in reader.deserialize::<gnomad_constraints::Record>() {
        let record = record?;
        result.insert(record.ensembl_gene_id.clone(), record);
    }

    Ok(result)
}

/// Load HGNC information.
///
/// # Result
///
/// A map from HGNC ID to HGNC record.
fn load_hgnc(path: String) -> Result<HashMap<String, hgnc::Record>, anyhow::Error> {
    info!("  loading HGNC information from {}", path);
    let mut result = HashMap::new();

    let reader = std::fs::File::open(path).map(|f| std::io::BufReader::new(f))?;
    for line in reader.lines() {
        let line = line?;
        let record = serde_json::from_str::<hgnc::Record>(&line)?;
        result.insert(record.hgnc_id.clone(), record);
    }

    Ok(result)
}

/// Load NCBI information.
///
/// # Result
///
/// A map from NCBI gene ID to NCBI record.
fn load_ncbi(path: String) -> Result<HashMap<String, ncbi::Record>, anyhow::Error> {
    info!("  loading NCBI information from {}", path);
    let mut result = HashMap::new();

    let reader = std::fs::File::open(path).map(|f| std::io::BufReader::new(f))?;
    for line in reader.lines() {
        let line = line?;
        let record = serde_json::from_str::<ncbi::Record>(&line)?;
        result.insert(record.gene_id.clone(), record);
    }

    Ok(result)
}

/// Construct tuned RocksDB options.
fn build_rocksdb_options() -> rocksdb::Options {
    let mut options = rocksdb::Options::default();

    options.create_if_missing(true);
    options.create_missing_column_families(true);
    options.prepare_for_bulk_load();
    options.set_disable_auto_compactions(true);

    // Compress all files with zstd.
    options.set_compression_per_level(&[]);
    options.set_compression_type(rocksdb::DBCompressionType::Zstd);
    // We only want to set level to 2 but have to set the rest as well using the Rust interface.
    // The (default) values for the other levels were taken from the output of a RocksDB
    // output folder created with default settings.
    options.set_compression_options(-14, 2, 0, 0);

    options
}

/// Write gene database to a RocksDB.
fn write_rocksdb(
    acmg_by_hgnc_id: HashMap<String, acmg_sf::Record>,
    constraints_by_ensembl_id: HashMap<String, gnomad_constraints::Record>,
    hgnc: HashMap<String, hgnc::Record>,
    ncbi_by_ncbi_id: HashMap<String, ncbi::Record>,
    args: &&Args,
) -> Result<(), anyhow::Error> {
    // Construct RocksDB options and open file for writing.
    let options = build_rocksdb_options();
    let db = rocksdb::DB::open_cf(&options, &args.path_out_rocksdb, ["meta", "genes"])?;

    let cf_meta = db.cf_handle("meta").unwrap();
    let cf_genes = db.cf_handle("genes").unwrap();

    tracing::info!("  writing meta data to database");
    db.put_cf(&cf_meta, "builder-version", format!("{VERSION:?}"))?;
    // TODO: read meta information about input data and write out

    tracing::info!("  compose genes data into database");
    let style = indicatif_style();
    for hgnc_record in hgnc.values().into_iter().progress_with_style(style) {
        let hgnc_id = hgnc_record.hgnc_id.clone();
        let record = data::Record {
            acmg_sf: acmg_by_hgnc_id.get(&hgnc_id).cloned(),
            gnomad_constraints: hgnc_record
                .ensembl_gene_id
                .as_ref()
                .map(|ensembl_gene_id| constraints_by_ensembl_id.get(ensembl_gene_id).cloned())
                .unwrap_or_default(),
            hgnc: hgnc_record.clone(),
            ncbi: hgnc_record
                .entrez_id
                .as_ref()
                .map(|entrez_id| ncbi_by_ncbi_id.get(entrez_id).cloned())
                .unwrap_or_default(),
        };
        db.put_cf(&cf_genes, hgnc_id, serde_json::to_string(&record)?)?;
    }

    // Finally, compact manually.
    tracing::info!("  enforce manual compaction");
    db.compact_range_cf(&cf_meta, None::<&[u8]>, None::<&[u8]>);
    db.compact_range_cf(&cf_genes, None::<&[u8]>, None::<&[u8]>);

    let compaction_start = Instant::now();
    let mut last_printed = compaction_start;
    while db
        .property_int_value(rocksdb::properties::COMPACTION_PENDING)?
        .unwrap()
        > 0
        || db
            .property_int_value(rocksdb::properties::NUM_RUNNING_COMPACTIONS)?
            .unwrap()
            > 0
    {
        std::thread::sleep(std::time::Duration::from_millis(100));
        if last_printed.elapsed() > std::time::Duration::from_millis(1000) {
            log::info!(
                "  ... waiting for compaction for {:?}",
                compaction_start.elapsed()
            );
            last_printed = Instant::now();
        }
    }

    info!("  RocksDB file is closed now.  Do not forget to remove the zero-byte WAL .log file.");

    Ok(())
}

/// Main entry point for the `sv bg-db-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting `db gene build`");
    info!("  common_args = {:?}", &common_args);
    info!("  args = {:?}", &args);

    let before_loading = Instant::now();
    info!("Loading genes data files...");
    let acmg_by_hgnc_id = load_acmg(format!("{}/genes/acmg/acmg.tsv", args.path_in_download))?;
    let constraints_by_ensembl_id = load_gnomad_constraints(format!(
        "{}/genes/gnomad_constraints/gnomad_constraints.tsv",
        args.path_in_download
    ))?;
    let hgnc = load_hgnc(format!(
        "{}/genes/hgnc/hgnc_info.jsonl",
        args.path_in_download
    ))?;
    let ncbi_by_ncbi_id = load_ncbi(format!(
        "{}/genes/ncbi/gene_info.jsonl",
        args.path_in_download
    ))?;
    info!(
        "... done loadin genes data files in {:?}",
        before_loading.elapsed()
    );

    let before_writing = Instant::now();
    info!("Writing genes database...");
    write_rocksdb(
        acmg_by_hgnc_id,
        constraints_by_ensembl_id,
        hgnc,
        ncbi_by_ncbi_id,
        &args,
    )?;
    info!(
        "... done writing genes database in {:?}",
        before_writing.elapsed()
    );

    Ok(())
}

#[cfg(test)]
pub mod test {
    use super::*;

    use crate::common::Args as CommonArgs;
    use clap_verbosity_flag::Verbosity;
    use temp_testdir::TempDir;

    #[test]
    fn smoke_test() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(1, 0),
        };
        let args = Args {
            path_in_download: String::from("tests/db/genes"),
            path_out_rocksdb: tmp_dir
                .to_path_buf()
                .into_os_string()
                .into_string()
                .unwrap(),
        };

        run(&common_args, &args)?;

        Ok(())
    }
}
