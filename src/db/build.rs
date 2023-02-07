//! Code supporting the `db build` sub command.

use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use enum_map::{enum_map, EnumMap};
use strum_macros::Display;
use thousands::Separable;
use tracing::{debug, info, trace};

use crate::{
    common::{md5sum, read_md5_file, sha256sum},
    db::conf::{GenomeRelease, Top},
};

use super::conf::{Database, DbDef, TadSet};

/// Copy one database definition file set.
fn copy_db_def(
    path_out: &PathBuf,
    base_path: &PathBuf,
    file_name: &str,
    file_ext: &str,
) -> Result<DbDef, anyhow::Error> {
    debug!(
        "Copying database definition file {}/{}{}",
        base_path.to_str().unwrap(),
        file_name,
        file_ext
    );

    debug!("Compute and verify MD5 checksum");
    // First create path to `.md5` file, then read in the MD5 string.
    let md5_path = base_path.join(&format!("{}{}.md5", file_name, file_ext));
    let md5_str = read_md5_file(&md5_path)?.to_lowercase();

    let file_path = base_path.join(&format!("{}{}", file_name, file_ext));
    let chk_md5 = md5sum(&file_path)?.to_lowercase();
    if md5_str != chk_md5 {
        return Err(anyhow::anyhow!(
            "MD5 mismatch, expected={} vs actual={}",
            md5_str,
            chk_md5
        ));
    }

    debug!("Compute and verify SHA256 checksum");
    let chk_sha256 = sha256sum(&file_path)?.to_lowercase();

    debug!("Create directory and copy file");
    trace!("Create {:?}", &path_out);
    std::fs::create_dir_all(&path_out)?;
    for ext in &[file_ext, ".spec.json"] {
        trace!("Copy {}{}", &file_name, ext);
        let name = format!("{}{}", file_name, ext);
        let bytes_copied = std::fs::copy(base_path.join(&name), path_out.join(&name))?;
        trace!("  copied {} bytes", bytes_copied.separate_with_commas());
    }

    Ok(DbDef {
        path: format!("{:?}", file_path),
        md5: Some(chk_md5),
        sha256: Some(chk_sha256),
        ..Default::default()
    })
}

/// Copy over the gene regions.
fn copy_gene_regions(args: &Args) -> Result<EnumMap<Database, DbDef>, anyhow::Error> {
    let path_out = args
        .path_worker_db
        .join("features")
        .join(args.genome_release.to_string())
        .join("gene_regions");
    let base_path = args
        .path_db_downloader
        .join("features")
        .join(args.genome_release.to_string())
        .join("gene_regions");

    Ok(enum_map! {
        Database::Ensembl => copy_db_def(&path_out, &base_path, "ensembl", ".bed.gz")?,
        Database::RefSeq => copy_db_def(&path_out, &base_path, "refseq", ".bed.gz")?,
    })
}

/// Copy over the tads.
fn copy_tads(args: &Args) -> Result<EnumMap<TadSet, DbDef>, anyhow::Error> {
    let path_out = args
        .path_worker_db
        .join("features")
        .join(args.genome_release.to_string())
        .join("tads");
    let base_path = args
        .path_db_downloader
        .join("features")
        .join(args.genome_release.to_string())
        .join("tads");

    Ok(enum_map! {
        TadSet::Hesc => copy_db_def(&path_out, &base_path, "hesc", ".bed")?,
        TadSet::Imr90 => copy_db_def(&path_out, &base_path, "imr90", ".bed")?,
    })
}

/// Local genome release for command line arguments.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, ValueEnum, Display)]
#[strum(serialize_all = "lowercase")]
pub enum ArgGenomeRelease {
    Grch37,
    Grch38,
}

/// Command line arguments for `db build` sub command.
#[derive(Parser, Debug)]
#[command(about = "Build worker database from downloader directory", long_about = None)]
pub struct Args {
    /// Genome build to use in the build.
    #[arg(long, value_enum, default_value_t = ArgGenomeRelease::Grch37)]
    pub genome_release: ArgGenomeRelease,
    /// Path to `varfish-db-downloader` directory (input).
    #[arg(long, required = true)]
    pub path_db_downloader: PathBuf,
    /// Path to `varfish-server-worker` directory (output).
    #[arg(long, required = true)]
    pub path_worker_db: PathBuf,
}

/// Main entry point for the `sv bg-db-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    info!("Starting `db build`");
    info!("  common_args = {:?}", &common_args);
    info!("  args = {:?}", &args);

    let genome_release = match args.genome_release {
        ArgGenomeRelease::Grch37 => GenomeRelease::Grch37,
        ArgGenomeRelease::Grch38 => GenomeRelease::Grch38,
    };
    let mut conf = Top {
        release_enabled: enum_map! {
            GenomeRelease::Grch37 => genome_release == GenomeRelease::Grch37,
            GenomeRelease::Grch38 => genome_release == GenomeRelease::Grch38,
        },
        ..Default::default()
    };

    info!("Copying features/{:?}/gene_regions", genome_release);
    conf.features[genome_release].gene_regions = copy_gene_regions(&args)?;

    info!("Copying features/{:?}/tads", genome_release);
    conf.features[genome_release].tads = copy_tads(&args)?;

    info!("Copying genes/acmg");
    info!("Copying genes/gnomad_constraints");
    info!("Copying genes/mim2gene");
    info!("Copying genes/xlink");
    info!("Copying vardbs/{:?}/strucvars", genome_release);

    info!(
        "Converting features/{:?}/gene_regions to binary",
        genome_release
    );
    info!("Converting vardbs/genes/xlink to binary");
    info!(
        "Converting vardbs/{:?}/strucvars/inhouse to binary",
        genome_release
    );
    info!(
        "Converting vardbs/{:?}/strucvars/clinvar to binary",
        genome_release
    );

    info!(
        "Writing out configuration file to {:?}",
        args.path_worker_db.join("conf.toml")
    );
    std::fs::write(
        args.path_worker_db.join("conf.toml"),
        toml::to_string_pretty(&conf)?,
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{run, Args};
    use crate::{common::Args as CommonArgs, db::build::ArgGenomeRelease};
    use clap_verbosity_flag::Verbosity;
    use temp_testdir::TempDir;

    #[test]
    fn run_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::default();
        let common_args = CommonArgs {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            genome_release: ArgGenomeRelease::Grch37,
            path_db_downloader: "tests/db/build/varfish-db-downloader".into(),
            path_worker_db: tmp_dir.to_path_buf(),
        };

        run(&common_args, &args)?;

        Ok(())
    }
}
