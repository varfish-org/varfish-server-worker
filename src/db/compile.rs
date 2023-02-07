//! Code supporting the `db build` sub command.

use std::path::{Path, PathBuf};

use clap::{Parser, ValueEnum};
use enum_map::{enum_map, EnumMap};
use strum_macros::Display;
use thousands::Separable;
use tracing::{debug, info, trace};

use crate::{
    common::{md5sum, read_md5_file, sha256sum},
    db::{
        conf::{GeneXlink, GenomeRelease, Top},
        to_bin,
    },
};

use super::conf::{
    default_min_overlap, default_slack_bnd, default_slack_ins, Database, DbDef, StrucVarDbs, TadSet,
};

/// Copy one database definition file set.
fn copy_db_def(
    db_path: &Path,
    path_out: &PathBuf,
    base_path: &Path,
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
    let md5_path = base_path.join(&format!("{file_name}{file_ext}.md5"));
    let md5_str = read_md5_file(&md5_path)?.to_lowercase();

    let file_path = base_path.join(&format!("{file_name}{file_ext}"));
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
    std::fs::create_dir_all(path_out)?;
    for ext in &[file_ext, ".spec.json"] {
        trace!("Copy {}{}", &file_name, ext);
        let name = format!("{file_name}{ext}");
        let bytes_copied = std::fs::copy(base_path.join(&name), path_out.join(&name))?;
        trace!("  copied {} bytes", bytes_copied.separate_with_commas());
    }

    Ok(DbDef {
        path: path_out
            .join(format!("{}{}", &file_name, &file_ext))
            .strip_prefix(db_path)?
            .to_str()
            .unwrap()
            .to_owned(),
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
        Database::Ensembl => copy_db_def(&args.path_worker_db, &path_out, &base_path, "ensembl", ".bed.gz")?,
        Database::RefSeq => copy_db_def(&args.path_worker_db, &path_out, &base_path, "refseq", ".bed.gz")?,
    })
}

/// Convert gene regions to binary.
fn convert_gene_regions(
    args: &Args,
    dbdefs: &mut EnumMap<Database, DbDef>,
) -> Result<(), anyhow::Error> {
    debug!("Converting gene regions");

    let base_path = args
        .path_worker_db
        .join("features")
        .join(args.genome_release.to_string())
        .join("gene_regions");

    to_bin::gene_region::convert_to_bin(
        base_path.join("ensembl.bed.gz"),
        base_path.join("ensembl.bin"),
        &args.path_worker_db,
        &mut dbdefs[Database::Ensembl],
    )?;
    to_bin::gene_region::convert_to_bin(
        base_path.join("refseq.bed.gz"),
        base_path.join("refseq.bin"),
        &args.path_worker_db,
        &mut dbdefs[Database::RefSeq],
    )?;

    Ok(())
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
        TadSet::Hesc => copy_db_def(&args.path_worker_db, &path_out, &base_path, "hesc", ".bed")?,
        TadSet::Imr90 => copy_db_def(&args.path_worker_db, &path_out, &base_path, "imr90", ".bed")?,
    })
}

/// Copy over the ACMG secondary findings list.
fn copy_acmg(args: &Args) -> Result<DbDef, anyhow::Error> {
    let path_out = args.path_worker_db.join("genes").join("acmg");
    let base_path = args.path_db_downloader.join("genes").join("acmg");

    copy_db_def(&args.path_worker_db, &path_out, &base_path, "acmg", ".tsv")
}

/// Copy over the gnomAD gene constraints file.
fn copy_gnomad_constraints(args: &Args) -> Result<DbDef, anyhow::Error> {
    let path_out = args.path_worker_db.join("genes").join("gnomad_constraints");
    let base_path = args
        .path_db_downloader
        .join("genes")
        .join("gnomad_constraints");

    copy_db_def(
        &args.path_worker_db,
        &path_out,
        &base_path,
        "gnomad_constraints",
        ".tsv",
    )
}

/// Copy over the mim2gene mapping file.
fn copy_mim2gene(args: &Args) -> Result<DbDef, anyhow::Error> {
    let path_out = args.path_worker_db.join("genes").join("mim2gene");
    let base_path = args.path_db_downloader.join("genes").join("mim2gene");

    copy_db_def(
        &args.path_worker_db,
        &path_out,
        &base_path,
        "mim2gene",
        ".tsv",
    )
}

/// Copy over the xlink mapping files.
fn copy_xlink(args: &Args) -> Result<EnumMap<GeneXlink, DbDef>, anyhow::Error> {
    let path_out = args.path_worker_db.join("genes").join("xlink");
    let base_path = args.path_db_downloader.join("genes").join("xlink");

    Ok(enum_map! {
        GeneXlink::Hgnc => copy_db_def(&args.path_worker_db, &path_out, &base_path, "hgnc", ".tsv")?,
        GeneXlink::Ensembl => copy_db_def(&args.path_worker_db, &path_out, &base_path, "ensembl", ".tsv")?,
    })
}

/// Convert gene ID xlink to binary.
fn convert_xlink(args: &Args, dbdefs: &mut EnumMap<GeneXlink, DbDef>) -> Result<(), anyhow::Error> {
    debug!("Converting xlink");

    let base_path = args.path_worker_db.join("genes").join("xlink");

    to_bin::xlink::convert_to_bin(
        base_path.join("ensembl.tsv"),
        base_path.join("ensembl.bin"),
        &args.path_worker_db,
        &mut dbdefs[GeneXlink::Ensembl],
    )?;
    to_bin::xlink::convert_to_bin(
        base_path.join("hgnc.tsv"),
        base_path.join("refseq.bin"),
        &args.path_worker_db,
        &mut dbdefs[GeneXlink::Hgnc],
    )?;

    Ok(())
}

/// Copy over the strucvar files.
fn copy_strucvar(args: &Args) -> Result<StrucVarDbs, anyhow::Error> {
    let path_out = args
        .path_worker_db
        .join("vardbs")
        .join(args.genome_release.to_string())
        .join("strucvar");
    let base_path = args
        .path_db_downloader
        .join("vardbs")
        .join(args.genome_release.to_string())
        .join("strucvar");

    Ok(StrucVarDbs {
        slack_bnd: default_slack_bnd(),
        slack_ins: default_slack_ins(),
        min_overlap: default_min_overlap(),

        gnomad_sv: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "gnomad_sv",
            ".bed.gz",
        )?,
        dbvar: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "dbvar",
            ".bed.gz",
        )?,
        dgv: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "dgv",
            ".bed.gz",
        )?,
        dgv_gs: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "dgv_gs",
            ".bed.gz",
        )?,
        exac: Some(copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "exac",
            ".bed.gz",
        )?),
        g1k: Some(copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "g1k",
            ".bed.gz",
        )?),

        inhouse: None,

        patho_mms: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "patho_mms",
            ".bed",
        )?,
        clinvar: copy_db_def(
            &args.path_worker_db,
            &path_out,
            &base_path,
            "clinvar",
            ".bed.gz",
        )?,
    })
}

/// Convert structural variants database.
fn convert_strucvar(args: &Args, dbdefs: &mut StrucVarDbs) -> Result<(), anyhow::Error> {
    use super::to_bin::vardbs::InputFileType;

    debug!("Converting strucvar");

    let base_path = args
        .path_worker_db
        .join("vardbs")
        .join(args.genome_release.to_string())
        .join("strucvar");

    to_bin::vardbs::convert_to_bin(
        base_path.join("gnomad_sv.bed.gz"),
        base_path.join("gnomad_sv.bin"),
        &args.path_worker_db,
        InputFileType::Gnomad,
        &mut dbdefs.gnomad_sv,
    )?;
    to_bin::vardbs::convert_to_bin(
        base_path.join("dbvar.bed.gz"),
        base_path.join("dbvar.bin"),
        &args.path_worker_db,
        InputFileType::Dbvar,
        &mut dbdefs.dbvar,
    )?;
    to_bin::vardbs::convert_to_bin(
        base_path.join("dgv.bed.gz"),
        base_path.join("dgv.bin"),
        &args.path_worker_db,
        InputFileType::Dgv,
        &mut dbdefs.dgv,
    )?;
    to_bin::vardbs::convert_to_bin(
        base_path.join("dgv_gs.bed.gz"),
        base_path.join("dgv_gs.bin"),
        &args.path_worker_db,
        InputFileType::DgvGs,
        &mut dbdefs.dgv_gs,
    )?;
    if let Some(exac) = &mut dbdefs.exac {
        to_bin::vardbs::convert_to_bin(
            base_path.join("exac.bed.gz"),
            base_path.join("exac.bin"),
            &args.path_worker_db,
            InputFileType::Exac,
            exac,
        )?;
    }
    if let Some(g1k) = &mut dbdefs.g1k {
        to_bin::vardbs::convert_to_bin(
            base_path.join("g1k.bed.gz"),
            base_path.join("g1k.bin"),
            &args.path_worker_db,
            InputFileType::G1k,
            g1k,
        )?;
    }

    to_bin::clinvar::convert_to_bin(
        base_path.join("clinvar.bed.gz"),
        base_path.join("clinvar.bin"),
        &args.path_worker_db,
        &mut dbdefs.clinvar,
    )?;

    Ok(())
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

    let path_conf_toml = args.path_worker_db.join("conf.toml");
    let mut conf = if path_conf_toml.exists() {
        info!("Reading database configuration...");
        let conf_toml = std::fs::read_to_string(&path_conf_toml)?;
        toml::from_str(&conf_toml)?
    } else {
        info!("Creating fresh database configuration...");
        Top {
            release_enabled: enum_map! {
                GenomeRelease::Grch37 => genome_release == GenomeRelease::Grch37,
                GenomeRelease::Grch38 => genome_release == GenomeRelease::Grch38,
            },
            ..Default::default()
        }
    };

    info!("Copying features/{:?}/gene_regions", genome_release);
    conf.features[genome_release].gene_regions = copy_gene_regions(args)?;

    info!("Copying features/{:?}/tads", genome_release);
    conf.features[genome_release].tads = copy_tads(args)?;

    info!("Copying genes/acmg");
    conf.genes.acmg = copy_acmg(args)?;

    info!("Copying genes/gnomad_constraints");
    conf.genes.gnomad_constraints = copy_gnomad_constraints(args)?;

    info!("Copying genes/mim2gene");
    conf.genes.mim2gene = copy_mim2gene(args)?;

    info!("Copying genes/xlink");
    conf.genes.xlink = copy_xlink(args)?;

    info!("Copying vardbs/{:?}/strucvar", genome_release);
    conf.vardbs[genome_release].strucvar = copy_strucvar(args)?;

    info!(
        "Converting features/{:?}/gene_regions to binary",
        genome_release
    );
    convert_gene_regions(args, &mut conf.features[genome_release].gene_regions)?;

    info!("Converting vardbs/genes/xlink to binary");
    convert_xlink(args, &mut conf.genes.xlink)?;

    info!("Converting vardbs/{:?}/strucvar to binary", genome_release);
    convert_strucvar(args, &mut conf.vardbs[genome_release].strucvar)?;

    info!(
        "Writing out configuration file to {:?}",
        args.path_worker_db.join("conf.toml")
    );
    std::fs::write(&path_conf_toml, toml::to_string_pretty(&conf)?)?;

    Ok(())
}

// #[cfg(test)]
// mod tests {
//     use super::{run, Args};
//     use crate::{common::Args as CommonArgs, db::compile::ArgGenomeRelease};
//     use clap_verbosity_flag::Verbosity;
//     use temp_testdir::TempDir;

//     #[test]
//     fn run_smoke() -> Result<(), anyhow::Error> {
//         let tmp_dir = TempDir::default();
//         let common_args = CommonArgs {
//             verbose: Verbosity::new(0, 0),
//         };
//         let args = Args {
//             genome_release: ArgGenomeRelease::Grch37,
//             path_db_downloader: "tests/db/build/varfish-db-downloader".into(),
//             path_worker_db: tmp_dir.to_path_buf(),
//         };

//         run(&common_args, &args)?;

//         Ok(())
//     }
// }
