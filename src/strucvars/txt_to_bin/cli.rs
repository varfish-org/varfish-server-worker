//! Command line interface for "db to-bin".

use std::path::PathBuf;

use clap::Parser;

use crate::{
    common::trace_rss_now,
    strucvars::txt_to_bin::{
        clinvar, masked,
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
    /// Convert gnomAD SV v2 to binary.
    StrucvarGnomadSv2,
    /// Convert gnomAD CNV v4 to binary.
    StrucvarGnomadCnv4,
    /// Convert gnomAD SV v4 to binary.
    StrucvarGnomadSv4,
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
    pub path_output: PathBuf,
}

/// Main entry point for the `strucvars txt-to-bin` command.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("Starting `strucvars txt-to-bin`");
    tracing::info!("  common_args = {:?}", &common_args);
    tracing::info!("  args = {:?}", &args);

    trace_rss_now();

    tracing::info!("Starting conversion...");
    match args.input_type {
        InputType::ClinvarSv => clinvar::convert_to_bin(&args.path_input, &args.path_output)?,
        InputType::StrucvarInhouse => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output,
            InputFileType::InhouseDb,
        )?,
        InputType::StrucvarDbVar => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output, InputFileType::Dbvar)?
        }
        InputType::StrucvarDgv => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output, InputFileType::Dgv)?
        }
        InputType::StrucvarDgvGs => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output, InputFileType::DgvGs)?
        }
        InputType::StrucvarExacCnv => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output, InputFileType::Exac)?
        }
        InputType::StrucvarG1k => {
            vardbs::convert_to_bin(&args.path_input, &args.path_output, InputFileType::G1k)?
        }
        InputType::StrucvarGnomadSv2 => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output,
            InputFileType::GnomadSv2,
        )?,
        InputType::StrucvarGnomadCnv4 => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output,
            InputFileType::GnomadCnv4,
        )?,
        InputType::StrucvarGnomadSv4 => vardbs::convert_to_bin(
            &args.path_input,
            &args.path_output,
            InputFileType::GnomadSv4,
        )?,
        InputType::MaskedRegion => masked::convert_to_bin(&args.path_input, &args.path_output)?,
        InputType::Xlink => xlink::convert_to_bin(&args.path_input, &args.path_output)?,
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

    #[rstest::rstest]
    fn run_clinvar_sv_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::ClinvarSv,
            path_input: String::from(
                "tests/db/to-bin/varfish-db-downloader/vardbs/clinvar/clinvar-svs.jsonl.gz",
            ),
            path_output: tmp_dir.join("clinvar.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/inhouse.tsv",
            ),
            path_output: tmp_dir.join("strucvar_inhouse.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dbvar.bed.gz",
            ),
            path_output: tmp_dir.join("strucvar_dbvar.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dgv.bed.gz",
            ),
            path_output: tmp_dir.join("strucvar_dgv.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/dgv_gs.bed.gz",
            ),
            path_output: tmp_dir.join("strucvar_dgv_gs.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_gnomad_exomes_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarExacCnv,
            path_input: String::from(
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/exac.bed.gz",
            ),
            path_output: tmp_dir.join("exac.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/g1k.bed.gz",
            ),
            path_output: tmp_dir.join("g1k.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_gnomad_sv2_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarGnomadSv2,
            path_input: String::from(
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch37/strucvar/gnomad_sv.bed.gz",
            ),
            path_output: tmp_dir.join("gnomad.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_gnomad_cnv4_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarGnomadCnv4,
            path_input: String::from(
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch38/strucvar/gnomad-cnv.bed.gz",
            ),
            path_output: tmp_dir.join("gnomad-cnv.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }

    #[test]
    fn run_strucvar_gnomad_sv4_smoke() -> Result<(), anyhow::Error> {
        let tmp_dir = temp_testdir::TempDir::default();
        let common_args = common::Args {
            verbose: Verbosity::new(0, 0),
        };
        let args = Args {
            input_type: InputType::StrucvarGnomadSv4,
            path_input: String::from(
                "tests/db/to-bin/varfish-db-downloader/vardbs/grch38/strucvar/gnomad-sv.bed.gz",
            ),
            path_output: tmp_dir.join("gnomad-sv.bin"),
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
                "tests/db/to-bin/varfish-db-downloader/features/grch37/masked/repeat.bed.gz",
            ),
            path_output: tmp_dir.join("masked.bin"),
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
            path_input: String::from("tests/db/to-bin/varfish-db-downloader/genes/xlink/hgnc.tsv"),
            path_output: tmp_dir.join("xlink.bin"),
        };

        super::run(&common_args, &args)?;

        Ok(())
    }
}
