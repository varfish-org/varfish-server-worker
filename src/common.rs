use clap_verbosity_flag::Verbosity;

use clap::Parser;

#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity,
}
