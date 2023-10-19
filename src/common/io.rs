//! Common, IO-related code.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use flate2::{bufread::MultiGzDecoder, write::GzEncoder, Compression};

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

/// Transparently open afile with gzip encoder.
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
mod test {
    use std::io::Read;

    #[rstest::rstest]
    #[case(true)]
    #[case(false)]
    fn open_write_maybe_gz(#[case] is_gzip: bool) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{:?}", is_gzip);

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

    #[test]
    fn read_lines() -> Result<(), anyhow::Error> {
        let lines = super::read_lines("tests/common/lines.txt")?.collect::<Result<Vec<_>, _>>()?;

        insta::assert_yaml_snapshot!(lines);

        Ok(())
    }
}
