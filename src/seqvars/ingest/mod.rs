//! Implementation of `seqvars ingest` subcommand.

use crate::common::{self, GenomeRelease};
use noodles_vcf as vcf;

pub mod header;

/// Command line arguments for `seqvars ingest` subcommand.
#[derive(Debug, clap::Parser)]
#[command(author, version, about = "ingest sequence variant VCF", long_about = None)]
pub struct Args {
    /// The assumed genome build.
    #[clap(long)]
    pub genomebuild: GenomeRelease,
    /// Path to input file.
    #[clap(long)]
    pub path_in: String,
    /// Path to output file.
    #[clap(long)]
    pub path_out: String,
}

/// Return the version of the `varfish-server-worker` crate and `x.y.z` in tests.
fn worker_version() -> &'static str {
    if cfg!(test) {
        "x.y.z"
    } else {
        env!("CARGO_PKG_VERSION")
    }
}

/// Main entry point for `seqvars ingest` sub command.
pub fn run(args_common: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    let before_anything = std::time::Instant::now();
    tracing::info!("args_common = {:#?}", &args_common);
    tracing::info!("args = {:#?}", &args);

    common::trace_rss_now();

    tracing::info!("opening input file...");
    let mut input_reader = {
        let file = std::fs::File::open(&args.path_in)
            .map_err(|e| anyhow::anyhow!("could not open input file {}: {}", &args.path_in, e))
            .map(std::io::BufReader::new)?;
        vcf::reader::Builder::default()
            .build_from_reader(file)
            .map_err(|e| anyhow::anyhow!("could not build VCF reader: {}", e))?
    };

    tracing::info!("analyzing header...");
    let input_header = input_reader
        .read_header()
        .map_err(|e| anyhow::anyhow!("problem reading VCF header: {}", e))?;
    let output_header =
        header::build_output_header(&input_header, args.genomebuild, worker_version())?;

    let mut output_writer = {
        let writer = std::fs::File::create(&args.path_out).map_err(|e| {
            anyhow::anyhow!(
                "could not output file for writing {}: {}",
                &args.path_out,
                e
            )
        })?;
        let writer = std::io::BufWriter::new(writer);
        vcf::writer::Writer::new(writer)
    };
    output_writer.write_header(&output_header)?;

    tracing::info!(
        "All of `seqvars ingest` completed in {:?}",
        before_anything.elapsed()
    );
    Ok(())
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use crate::common::GenomeRelease;

    macro_rules! set_snapshot_suffix {
        ($($expr:expr),*) => {
            let mut settings = insta::Settings::clone_current();
            settings.set_snapshot_suffix(format!($($expr,)*));
            let _guard = settings.bind_to_scope();
        }
    }

    #[rstest]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn smoke_test_run(#[case] path: &str) {
        set_snapshot_suffix!("{:?}", path.split('/').last().unwrap().replace(".", "_"));

        let tmpdir = temp_testdir::TempDir::default();

        let args_common = Default::default();
        let args = super::Args {
            genomebuild: GenomeRelease::Grch37,
            path_in: path.into(),
            path_out: tmpdir
                .join("out.vcf")
                .to_str()
                .expect("invalid path")
                .into(),
        };
        super::run(&args_common, &args).unwrap();
    }
}
