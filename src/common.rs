//! Common functionality.

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    ops::Range,
    path::Path,
};

use byte_unit::Byte;
use clap_verbosity_flag::{InfoLevel, Verbosity};

use clap::Parser;
use flate2::{bufread::MultiGzDecoder, write::GzEncoder, Compression};
use md5::{Digest, Md5};
use sha2::Sha256;
use tracing::{debug, trace};

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
    let page_size = procfs::page_size().unwrap();
    debug!(
        "RSS now: {}",
        Byte::from_bytes((me.stat().unwrap().rss * page_size) as u128).get_appropriate_unit(true)
    );
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

/// Definition of canonical chromosome names.
pub const CHROMS: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y", "M",
];

/// Build mapping of chromosome names to chromosome counts.
pub fn build_chrom_map() -> HashMap<String, usize> {
    let mut result = HashMap::new();
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
pub fn open_read_maybe_gz<P>(path: P) -> Result<Box<dyn Read>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        trace!("Opening {:?} as gzip for reading", path.as_ref());
        let file = File::open(path)?;
        let bufreader = BufReader::new(file);
        let decoder = MultiGzDecoder::new(bufreader);
        Ok(Box::new(decoder))
    } else {
        trace!("Opening {:?} as plain text for reading", path.as_ref());
        let file = File::open(path)?;
        Ok(Box::new(file))
    }
}

/// Transparently opena  file with gzip encoder.
pub fn open_write_maybe_gz<P>(path: P) -> Result<Box<dyn Write>, anyhow::Error>
where
    P: AsRef<Path>,
{
    if path.as_ref().extension().map(|s| s.to_str()) == Some(Some("gz")) {
        trace!("Opening {:?} as gzip for writing", path.as_ref());
        let file = File::create(path)?;
        let bufwriter = BufWriter::new(file);
        let encoder = GzEncoder::new(bufwriter, Compression::default());
        Ok(Box::new(encoder))
    } else {
        trace!("Opening {:?} as plain text for writing", path.as_ref());
        let file = File::open(path)?;
        Ok(Box::new(file))
    }
}

// Compute reciprocal overlap between two ranges.
pub fn reciprocal_overlap(lhs: Range<u32>, rhs: Range<u32>) -> f32 {
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

/// Compute SHA256 sum for file at `path`.
pub fn sha256sum<P>(path: P) -> Result<String, anyhow::Error>
where
    P: AsRef<Path> + std::fmt::Debug,
{
    debug!("Computing SHA256 checksum for {:?}", &path);

    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let _n = std::io::copy(&mut file, &mut hasher)?;
    let hash = hasher.finalize();
    let mut buf = [0u8; 64];
    let checksum = base16ct::lower::encode_str(&hash, &mut buf).unwrap();
    debug!(" SHA256 = {}", &checksum);
    Ok(checksum.to_owned())
}

/// Compute MD5 sum for file at `path`.
pub fn md5sum<P>(path: P) -> Result<String, anyhow::Error>
where
    P: AsRef<Path> + std::fmt::Debug,
{
    debug!("Computing MD5 checksum for {:?}", &path);
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Md5::new();
    const BUF_LEN: usize = 65_536;
    let mut _bytes_read = 0;
    let mut buffer = [0; BUF_LEN];

    loop {
        let n = file.read(&mut buffer)?;
        hasher.update(&buffer[..n]);
        _bytes_read += n;
        if n != BUF_LEN {
            break;
        }
    }
    let hash = hasher.finalize();
    let mut buf = [0u8; 64];
    let checksum = base16ct::lower::encode_str(&hash, &mut buf).unwrap();

    debug!(" MD5 = {}", &checksum);
    Ok(checksum.to_owned())
}

/// Read .md5 file and return String representation of result.
pub fn read_md5_file<P>(path: P) -> Result<String, anyhow::Error>
where
    P: AsRef<Path> + std::fmt::Debug,
{
    let fcontents = std::fs::read_to_string(&path)
        .map_err(|e| anyhow::anyhow!("Could not open .md5 file {:?}: {}", &path, e))?;
    let md5_str = fcontents
        .split_whitespace()
        .into_iter()
        .next()
        .ok_or(anyhow::anyhow!("Could not get MD5 sum from {:?}", path))?;
    Ok(md5_str.to_owned())
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

#[cfg(test)]
mod tests {
    use crate::common::{md5sum, sha256sum};

    use super::{numeric_gene_id, read_md5_file};
    use pretty_assertions::assert_eq;

    #[test]
    fn load_test_md5() -> Result<(), anyhow::Error> {
        let expected = "d41d8cd98f00b204e9800998ecf8427e";
        let actual = read_md5_file("tests/common/test.md5")?;

        assert_eq!(expected, actual);

        Ok(())
    }
    #[test]
    fn sha256sum_on_payload_txt() -> Result<(), anyhow::Error> {
        let expected = "da8868ebc063390856b50f147c8ed8f55e9115df1b21ce97973deec5646eb217";
        let actual = sha256sum("tests/common/payload.txt")?;

        assert_eq!(expected, actual);

        Ok(())
    }

    #[test]
    fn md5sum_on_payload_txt() -> Result<(), anyhow::Error> {
        let expected = "310b5a47c12b03c6e71684d170b287a6";
        let actual = md5sum("tests/common/payload.txt")?;

        assert_eq!(expected, actual);

        Ok(())
    }

    #[test]
    fn numeric_gene_id_simple() -> Result<(), anyhow::Error> {
        assert_eq!(1, numeric_gene_id("ENSG0000000001")?);
        assert_eq!(1, numeric_gene_id("ENSG1")?);
        assert_eq!(1, numeric_gene_id("1")?);

        Ok(())
    }
}
