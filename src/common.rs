//! Common functionality.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    ops::Range,
    path::Path,
};

use byte_unit::Byte;
use clap_verbosity_flag::{InfoLevel, Verbosity};

use clap::Parser;
use flate2::{bufread::MultiGzDecoder, write::GzEncoder, Compression};
use hgvs::static_data::Assembly;
use indexmap::IndexMap;

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}

/// Helper to print the current memory resident set size via `tracing`.
pub fn trace_rss_now() {
    let me = procfs::process::Process::myself().unwrap();
    let page_size = procfs::page_size();
    tracing::debug!(
        "RSS now: {}",
        Byte::from_bytes((me.stat().unwrap().rss * page_size) as u128).get_appropriate_unit(true)
    );
}

/// Definition of canonical chromosome names.
pub const CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M",
];

/// Build mapping of chromosome names to chromosome counts.
pub fn build_chrom_map() -> IndexMap<String, usize> {
    let mut result = IndexMap::new();
    for (i, &chrom_name) in CHROMS.iter().enumerate() {
        result.insert(chrom_name.to_owned(), i);
        result.insert(format!("chr{chrom_name}").to_owned(), i);
    }
    result.insert("x".to_owned(), 22);
    result.insert("y".to_owned(), 23);
    result.insert("chrx".to_owned(), 22);
    result.insert("chry".to_owned(), 23);
    result.insert("mt".to_owned(), 24);
    result.insert("m".to_owned(), 24);
    result.insert("chrmt".to_owned(), 24);
    result.insert("chrm".to_owned(), 24);
    result.insert("MT".to_owned(), 24);
    result.insert("chrMT".to_owned(), 24);
    result
}

/// Transparently open a file with gzip decoder.
pub fn open_read_maybe_gz<P>(path: P) -> Result<Box<dyn BufRead>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        tracing::trace!("Opening {:?} as gzip for reading", path.as_ref());
        let file = File::open(path)?;
        let bufreader = BufReader::new(file);
        let decoder = MultiGzDecoder::new(bufreader);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        tracing::trace!("Opening {:?} as plain text for reading", path.as_ref());
        let file = File::open(path).map(BufReader::new)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Transparently opena  file with gzip encoder.
pub fn open_write_maybe_gz<P>(path: P) -> Result<Box<dyn Write>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        tracing::trace!("Opening {:?} as gzip for writing", path.as_ref());
        let file = File::create(path)?;
        let bufwriter = BufWriter::new(file);
        let encoder = GzEncoder::new(bufwriter, Compression::default());
        Ok(Box::new(encoder))
    } else {
        tracing::trace!("Opening {:?} as plain text for writing", path.as_ref());
        let file = File::create(path)?;
        Ok(Box::new(file))
    }
}

// Compute reciprocal overlap between two ranges.
pub fn reciprocal_overlap(lhs: Range<i32>, rhs: Range<i32>) -> f32 {
    let lhs_b = lhs.start;
    let lhs_e = lhs.end;
    let rhs_b = rhs.start;
    let rhs_e = rhs.end;
    let ovl_b = std::cmp::max(lhs_b, rhs_b);
    let ovl_e = std::cmp::min(lhs_e, rhs_e);
    if ovl_b >= ovl_e {
        0f32
    } else {
        let ovl_len = (ovl_e - ovl_b) as f32;
        let x1 = ovl_len / (lhs_e - lhs_b) as f32;
        let x2 = ovl_len / (rhs_e - rhs_b) as f32;
        x1.min(x2)
    }
}

/// Helper to convert ENSEMBL and RefSeq gene ID to u32.
pub fn numeric_gene_id(raw_id: &str) -> Result<u32, anyhow::Error> {
    let clean_id = if raw_id.starts_with("ENSG") {
        // Strip "ENSG" prefix and as many zeroes as follow
        raw_id
            .chars()
            .skip("ENSG".len())
            .skip_while(|c| *c == '0')
            .collect()
    } else {
        raw_id.to_owned()
    };

    clean_id
        .parse::<u32>()
        .map_err(|e| anyhow::anyhow!("could not parse gene id {:?}: {}", &clean_id, &e))
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P: AsRef<Path>>(
    filename: P,
) -> std::io::Result<std::io::Lines<std::io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(std::io::BufReader::new(file).lines())
}

/// Select the genome release to use.
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
pub enum GenomeRelease {
    // GRCh37 / hg19
    #[strum(serialize = "GRCh37")]
    Grch37,
    /// GRCh38 / hg38
    #[strum(serialize = "GRCh38")]
    Grch38,
}

impl GenomeRelease {
    pub fn name(&self) -> String {
        match self {
            GenomeRelease::Grch37 => String::from("GRCh37"),
            GenomeRelease::Grch38 => String::from("GRCh38"),
        }
    }
}

impl From<GenomeRelease> for Assembly {
    fn from(val: GenomeRelease) -> Self {
        match val {
            GenomeRelease::Grch37 => Assembly::Grch37p10,
            GenomeRelease::Grch38 => Assembly::Grch38,
        }
    }
}

impl From<Assembly> for GenomeRelease {
    fn from(assembly: Assembly) -> Self {
        match assembly {
            Assembly::Grch37 | Assembly::Grch37p10 => GenomeRelease::Grch37,
            Assembly::Grch38 => GenomeRelease::Grch38,
        }
    }
}

impl std::str::FromStr for GenomeRelease {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_ascii_lowercase();
        if s.starts_with("grch37") {
            Ok(GenomeRelease::Grch37)
        } else if s.starts_with("grch38") {
            Ok(GenomeRelease::Grch38)
        } else {
            Err(anyhow::anyhow!("Unknown genome release: {}", s))
        }
    }
}

/// The version of `varfish-server-worker` package.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn numeric_gene_id_simple() -> Result<(), anyhow::Error> {
        assert_eq!(1, numeric_gene_id("ENSG0000000001")?);
        assert_eq!(1, numeric_gene_id("ENSG1")?);
        assert_eq!(1, numeric_gene_id("1")?);

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::io::Read;

    #[test]
    fn trace_rss_now_smoke() {
        super::trace_rss_now();
    }

    #[test]
    fn build_chrom_map_snapshot() {
        let map = super::build_chrom_map();
        insta::assert_yaml_snapshot!(map);
    }

    #[rstest::rstest]
    #[case(true)]
    #[case(false)]
    fn open_write_maybe_gz(#[case] is_gzip: bool) -> Result<(), anyhow::Error> {
        let filename = if is_gzip { "test.txt" } else { "test.txt.gz" };
        let tmp_dir = temp_testdir::TempDir::default();

        {
            let mut f = super::open_write_maybe_gz(tmp_dir.join(filename))?;
            f.flush()?;
        }

        let mut f = std::fs::File::open(tmp_dir.join(filename)).map(std::io::BufReader::new)?;
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;

        insta::assert_snapshot!(format!("{:x?}", &buf));

        Ok(())
    }

    #[rstest::rstest]
    #[case(true)]
    #[case(false)]
    fn open_read_maybe_gz(#[case] is_gzip: bool) -> Result<(), anyhow::Error> {
        let mut f = super::open_read_maybe_gz(if is_gzip {
            "tests/common/test.txt.gz"
        } else {
            "tests/common/test.txt"
        })?;

        let mut buf = String::new();
        f.read_to_string(&mut buf)?;

        insta::assert_snapshot!(&buf);

        Ok(())
    }

    #[rstest::rstest]
    #[case(0..10, 0..10, 1.0)]
    #[case(0..10, 5..15, 0.5)]
    #[case(5..15, 0..10, 0.5)]
    #[case(0..10, 10..20, 0.0)]
    #[case(0..2, 0..10, 0.2)]
    #[case(0..10, 0..2, 0.2)]
    fn reciprocal_overlap(
        #[case] lhs: std::ops::Range<i32>,
        #[case] rhs: std::ops::Range<i32>,
        #[case] expected: f32,
    ) {
        let actual = super::reciprocal_overlap(lhs, rhs);
        assert!(float_cmp::approx_eq!(f32, expected, actual, ulps = 2))
    }

    #[rstest::rstest]
    #[case("ENSG0000000142", 142)]
    #[case("42", 42)]
    fn numeric_gene_id(#[case] raw_id: &str, #[case] expected: u32) {
        let actual = super::numeric_gene_id(raw_id).unwrap();
        assert_eq!(expected, actual);
    }

    #[test]
    fn read_lines() -> Result<(), anyhow::Error> {
        let lines = super::read_lines("tests/common/lines.txt")?.collect::<Result<Vec<_>, _>>()?;

        insta::assert_yaml_snapshot!(lines);

        Ok(())
    }

    #[rstest::rstest]
    #[case(crate::common::GenomeRelease::Grch37, "GRCh37")]
    #[case(crate::common::GenomeRelease::Grch38, "GRCh38")]
    fn genome_release_name(#[case] release: super::GenomeRelease, #[case] expected: &str) {
        assert_eq!(expected, release.name());
    }

    #[rstest::rstest]
    #[case(
        crate::common::GenomeRelease::Grch37,
        hgvs::static_data::Assembly::Grch37p10
    )]
    #[case(
        crate::common::GenomeRelease::Grch38,
        hgvs::static_data::Assembly::Grch38
    )]
    fn assembly_from_genome_release(
        #[case] release: super::GenomeRelease,
        #[case] assembly: hgvs::static_data::Assembly,
    ) -> Result<(), anyhow::Error> {
        let res: hgvs::static_data::Assembly = release.into();

        assert_eq!(res, assembly);

        Ok(())
    }

    #[rstest::rstest]
    #[case(
        crate::common::GenomeRelease::Grch37,
        hgvs::static_data::Assembly::Grch37
    )]
    #[case(
        crate::common::GenomeRelease::Grch37,
        hgvs::static_data::Assembly::Grch37p10
    )]
    #[case(
        crate::common::GenomeRelease::Grch38,
        hgvs::static_data::Assembly::Grch38
    )]
    fn genome_release_from_assembly(
        #[case] release: super::GenomeRelease,
        #[case] assembly: hgvs::static_data::Assembly,
    ) -> Result<(), anyhow::Error> {
        let res: super::GenomeRelease = assembly.into();

        assert_eq!(res, release);

        Ok(())
    }

    #[rstest::rstest]
    #[case(crate::common::GenomeRelease::Grch37, "grch37")]
    #[case(crate::common::GenomeRelease::Grch38, "grch38")]
    fn genome_relese_from_str(
        #[case] release: super::GenomeRelease,
        #[case] s: &str,
    ) -> Result<(), anyhow::Error> {
        let res: super::GenomeRelease = s.parse()?;

        assert_eq!(res, release);

        Ok(())
    }
}
