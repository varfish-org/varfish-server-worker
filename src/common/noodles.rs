//! Helper code for using noodles with S3.
//!
//! At the moment, we only support reading but not yet writing.
//!
//! The main reason is that the awslabs SDK for Rust currently supports easily
//! getting an `AsyncRead` from an S3 object but not an `AsyncWrite`.  For
//! this, we would have to create a wrapper that writes to a multipart upload
//! or similar.

use async_compression::tokio::bufread::GzipDecoder;
use mehari::common::io::{std::is_gz, tokio::open_read_maybe_gz};
use mehari::common::noodles::{AsyncVcfReader, VariantReader};
use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles::tabix;
use noodles::vcf;
use noodles::vcf::variant::Record;
use std::{path::Path, pin::Pin};
use tokio::io::{AsyncBufRead, BufReader};

/// Build TBI for file at `path_src` and write to `path_dst`.
pub async fn build_tbi<S, D>(path_src: S, path_dst: D) -> Result<(), anyhow::Error>
where
    S: AsRef<std::path::Path>,
    D: AsRef<std::path::Path>,
{
    let mut reader = tokio::fs::File::open(path_src.as_ref())
        .await
        .map(bgzf::AsyncReader::new)
        .map(vcf::AsyncReader::new)
        .map_err(|e| anyhow::anyhow!("error input file for tbi creation: {}", e))?;

    let header = reader
        .read_header()
        .await
        .map_err(|e| anyhow::anyhow!("error reading header: {}", e))?;

    let mut record = vcf::Record::default();

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut start_position = reader.get_ref().virtual_position();

    while reader
        .read_record(&mut record)
        .await
        .map_err(|e| anyhow::anyhow!("problem reading record: {}", e))?
        != 0
    {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record.reference_sequence_name().to_string();
        let start = record
            .variant_start()
            .expect("no variant_start?")
            .map_err(|e| anyhow::anyhow!("error converting start position: {}", e))?;
        let end = record
            .variant_end(&header)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
            .and_then(|position| {
                Position::try_from(usize::from(position))
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
            })
            .map_err(|e| anyhow::anyhow!("error converting end position: {}", e))?;

        indexer
            .add_record(&reference_sequence_name, start, end, chunk)
            .map_err(|e| anyhow::anyhow!("error adding record to tabix index: {}", e))?;

        start_position = end_position;
    }

    let index = indexer.build();

    {
        let mut writer = tokio::fs::File::create(path_dst.as_ref())
            .await
            .map(tokio::io::BufWriter::new)
            .map(tabix::AsyncWriter::new)
            .map_err(|e| anyhow::anyhow!("error output file for tbi creation: {}", e))?;
        writer
            .write_index(&index)
            .await
            .map_err(|e| anyhow::anyhow!("error writing tabix index: {}", e))?;
        writer
            .shutdown()
            .await
            .map_err(|e| anyhow::anyhow!("error flushing tabix index: {}", e))?;

        tokio::time::sleep(std::time::Duration::from_millis(200)).await;
    }

    Ok(())
}

/// Open plain text or gzip reader via S3.
pub async fn s3_open_read_maybe_gz<P>(path: P) -> Result<Pin<Box<dyn AsyncBufRead>>, anyhow::Error>
where
    P: AsRef<Path>,
{
    // Get configuration from environment variables.
    let access_key = if let Ok(access_key) = std::env::var("AWS_ACCESS_KEY_ID") {
        access_key
    } else {
        anyhow::bail!("could not access key from env AWS_ACCESS_KEY_ID")
    };
    let secret_key = if let Ok(secret_key) = std::env::var("AWS_SECRET_ACCESS_KEY") {
        secret_key
    } else {
        anyhow::bail!("could not get secret key from env AWS_SECRET_ACCESS_KEY")
    };
    let endpoint_url = if let Ok(endpoint_url) = std::env::var("AWS_ENDPOINT_URL") {
        endpoint_url
    } else {
        anyhow::bail!("could not get endpoint url from env AWS_ENDPOINT_URL")
    };
    let region = if let Ok(region) = std::env::var("AWS_REGION") {
        region
    } else {
        anyhow::bail!("could not AWS region from env AWS_REGION")
    };

    let cred =
        aws_sdk_s3::config::Credentials::new(access_key, secret_key, None, None, "loaded-from-env");
    let s3_config = aws_sdk_s3::config::Builder::new()
        .endpoint_url(&endpoint_url)
        .credentials_provider(cred)
        .region(aws_config::Region::new(region))
        .force_path_style(true) // apply bucketname as path param instead of pre-domain
        .build();

    // Split bucket and path from input path.
    let path_string = format!("{}", path.as_ref().display());
    let (bucket, key) = if let Some((bucket, key)) = path_string.split_once('/') {
        (bucket.to_string(), key.to_string())
    } else {
        anyhow::bail!("invalid S3 path: {}", path.as_ref().display());
    };

    // Setup S3 client and access object.
    let client = aws_sdk_s3::Client::from_conf(s3_config);
    let object = client.get_object().bucket(&bucket).key(&key).send().await?;

    let path_is_gzip = is_gz(path.as_ref());
    tracing::debug!(
        "Opening S3 object {} as {} for reading (async)",
        path.as_ref().display(),
        if path_is_gzip {
            "gzip (allow multi-member)"
        } else {
            "plain text"
        }
    );

    if path_is_gzip {
        let bufreader = BufReader::new(object.body.into_async_read());
        let decoder = {
            let mut decoder = GzipDecoder::new(bufreader);
            decoder.multiple_members(true);
            decoder
        };
        Ok(Box::pin(BufReader::new(decoder)))
    } else {
        Ok(Box::pin(BufReader::new(object.body.into_async_read())))
    }
}

/// Helper function that opens a list of paths as VCF readers.
pub async fn open_vcf_readers(paths: &[String]) -> Result<Vec<VariantReader>, anyhow::Error> {
    let mut result = Vec::new();
    for path in paths.iter() {
        let buf_read = if super::s3::s3_mode() && !path.starts_with('/') {
            s3_open_read_maybe_gz(path).await?
        } else {
            open_read_maybe_gz(path).await?
        };
        result.push(VariantReader::Vcf(AsyncVcfReader::new(buf_read)));
    }
    Ok(result)
}

/// Helper function that opens one VCF reader at the given path.
///
/// The behaviour is as follows:
///
/// - If `path_in` is "-" then open stdin and read as plain text.
/// - If environment variable `AWS_ACCESS_KEY_ID` is set then enable S3 mode.
/// - If `path_in` is absolute or S3 mode is disabled then open `path_in` as local file
/// - Otherwise, attempt to open `path_in` as S3 object.
pub async fn open_vcf_reader(path_in: &str) -> Result<VariantReader, anyhow::Error> {
    if super::s3::s3_mode() && path_in != "-" && !path_in.starts_with('/') {
        tracing::debug!("Opening S3 object {} for reading (async)", path_in);
        Ok(VariantReader::Vcf(vcf::AsyncReader::new(
            s3_open_read_maybe_gz(path_in)
                .await
                .map_err(|e| anyhow::anyhow!("could not build VCF reader from S3 file: {}", e))?,
        )))
    } else {
        tracing::debug!("Opening local file {} for reading (async)", path_in);
        Ok(VariantReader::Vcf(vcf::AsyncReader::new(
            open_read_maybe_gz(path_in).await.map_err(|e| {
                anyhow::anyhow!("could not build VCF reader from local file: {}", e)
            })?,
        )))
    }
}

#[cfg(test)]
mod test {
    #[tokio::test]
    async fn build_tbi() -> Result<(), anyhow::Error> {
        let tmpdir = temp_testdir::TempDir::default();

        let path_in: String = "tests/seqvars/ingest/NA12878_dragen.vcf.gz".into();
        let path_out: String = tmpdir
            .join("out.vcf.gz.tbi")
            .to_str()
            .expect("invalid path")
            .into();
        super::build_tbi(&path_in, &path_out).await?;

        // This must be fixed first https://github.com/zaeleus/noodles/issues/213
        // let mut buffer: Vec<u8> = Vec::new();
        // hxdmp::hexdump(&crate::common::read_to_bytes(&path_out)?, &mut buffer)?;
        // insta::assert_snapshot!(String::from_utf8_lossy(&buffer));

        Ok(())
    }
}
