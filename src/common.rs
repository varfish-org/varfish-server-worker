//! Common functionality.

use byte_unit::Byte;
use clap_verbosity_flag::Verbosity;

use clap::Parser;

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity,
}

/// Helper to print the current memory resident set size to a `Term`.
pub fn print_rss_now(term: &console::Term) -> Result<(), anyhow::Error> {
    let me = procfs::process::Process::myself().unwrap();
    let page_size = procfs::page_size().unwrap();
    term.write_line(&format!(
        "RSS now: {}",
        Byte::from_bytes((me.stat().unwrap().rss * page_size) as u128).get_appropriate_unit(true)
    ))?;
    Ok(())
}
