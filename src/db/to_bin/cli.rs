//! Command line interface for "db to-bin".

use std::path::PathBuf;

use clap::Parser;

use crate::{
    common::trace_rss_now,
    db::to_bin::{
        clinvar, gene_region, masked,
        vardbs::{self, InputFileType},
        xlink,
    },
};

/// Select input/conversion type.
#[derive(
    clap::ValueEnum,
    Clone,
    Copy,
    Debug,
    strum::Display,
    PartialEq,
    Eq,
    enum_map::Enum,
    PartialOrd,
    Ord,
    Hash,
)]
pub enum InputType {
    /// Convert ClinVar Structural Variation to binary.
    ClinvarSv,
    /// Convert In-house SV to binary.
    StrucvarInhouse,
    /// Convert DbVar to binary.
    StrucvarDbVar,
    /// Convert DGV to binary.
    StrucvarDgv,
    /// Convert DGV gold standard to binary.
    StrucvarDgvGs,
    /// Convert ExAC CNV to binary.
    StrucvarExacCnv,
    /// Convert Thousand Genomes to binary.
    StrucvarG1k,
    /// Convert gnomAD SV to binary.
    StrucvarGnomadSv,
    /// Convert gene region to binary.
    GeneRegion,
    /// Convert masked region to binary.
    MaskedRegion,
    /// Convert cross-link to binary.
    Xlink,
}

/// Command line arguments for `db build` sub command.
#[derive(Parser, Debug)]
#[command(about = "Convert to binary protobuf files", long_about = None)]
pub struct Args {
    /// Input type to convert to binary.
    #[arg(long, value_enum)]
    pub input_type: InputType,
    /// Input file to convert to binary.
    #[arg(long)]
    pub path_input: String,
    /// Path to output BIN file.
    #[arg(long)]
    pub path_output_bin: PathBuf,
}

/// Main entry point for the `sv bg-db-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("Starting `db to-bin`");
    tracing::info!("  common_args = {:?}", &common_args);
    tracing::info!("  args = {:?}", &args);

    trace_rss_now();

    tracing::info!("Starting conversion...");
    match args.input_type {
        InputType::ClinvarSv => clinvar::convert_to_bin(&args.path_input, &args.path_output_bin)?,
        InputType::StrucvarInhouse => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output_bin,
            InputFileType::InhouseDb,
        )?,
        InputType::StrucvarDbVar => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output_bin,
            InputFileType::Dbvar,
        )?,
        InputType::StrucvarDgv => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output_bin, InputFileType::Dgv)?
        }
        InputType::StrucvarDgvGs => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output_bin,
            InputFileType::DgvGs,
        )?,
        InputType::StrucvarExacCnv => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output_bin, InputFileType::Exac)?
        }
        InputType::StrucvarG1k => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output_bin, InputFileType::G1k)?
        }
        InputType::StrucvarGnomadSv => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output_bin,
            InputFileType::Gnomad,
        )?,
        InputType::GeneRegion => {
            gene_region::convert_to_bin(&args.path_input, &args.path_output_bin)?
        }
        InputType::MaskedRegion => masked::convert_to_bin(&args.path_input, &args.path_output_bin)?,
        InputType::Xlink => xlink::convert_to_bin(&args.path_input, &args.path_output_bin)?,
    }
    tracing::info!("... done with conversion");

    trace_rss_now();

    Ok(())
}

#[cfg(test)]
mod test {
    use clap_verbosity_flag::Verbosity;

    use crate::common;

    use super::{Args, InputType};

    #[test]
    fn run_clinvar_sv_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::ClinvarSv,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/clinvar.bed.gz",
            ),
            path_output_bin: tmp_dir.join("clinvar.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_inhouse_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarInhouse,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/inhouse.tsv",
            ),
            path_output_bin: tmp_dir.join("strucvar_inhouse.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_dbvar_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarDbVar,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/dbvar.bed.gz",
            ),
            path_output_bin: tmp_dir.join("strucvar_dbvar.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_dgv_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarDgv,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/dgv.bed.gz",
            ),
            path_output_bin: tmp_dir.join("strucvar_dgv.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_dgv_gs_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarDgvGs,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/dgv_gs.bed.gz",
            ),
            path_output_bin: tmp_dir.join("strucvar_dgv_gs.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_exac_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarExacCnv,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/exac.bed.gz",
            ),
            path_output_bin: tmp_dir.join("exac.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_g1k_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarG1k,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/g1k.bed.gz",
            ),
            path_output_bin: tmp_dir.join("g1k.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_gnomad_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarGnomadSv,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/vardbs/grch37/strucvar/gnomad_sv.bed.gz",
            ),
            path_output_bin: tmp_dir.join("gnomad.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_gene_region_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::GeneRegion,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/features/grch37/gene_regions/refseq.bed.gz",
            ),
            path_output_bin: tmp_dir.join("refseq.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_masked_region_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::MaskedRegion,
            path_input: String::from(
                "tests/db/compile/varfish-db-downloader/features/grch37/masked/repeat.bed.gz",
            ),
            path_output_bin: tmp_dir.join("masked.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_xlink_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::Xlink,
            path_input: String::from("tests/db/compile/varfish-db-downloader/genes/xlink/hgnc.tsv"),
            path_output_bin: tmp_dir.join("xlink.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }
}
