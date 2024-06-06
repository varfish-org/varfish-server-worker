//! Helper code for working with S3.

use mehari::common::io::std::is_gz;

/// Helper that returns whether S3 mode has been enabled via `AWS_ACCESS_KEY_ID`.
pub fn s3_mode() -> bool {
    let result = std::env::var("AWS_ACCESS_KEY_ID").is_ok();
    tracing::trace!("S3 mode is {}", if result { "enabled" } else { "disabled" });
    result
}

/// Return the S3 configuration from environment variables.
pub async fn config_from_env() -> Result<aws_sdk_s3::config::Config, anyhow::Error> {
    let config = aws_config::from_env().load().await;
    let endpoint_url = std::env::var("AWS_ENDPOINT_URL").map_err(|e| {
        anyhow::anyhow!(
            "Could not get endpoint url from env AWS_ENDPOINT_URL: {}",
            e
        )
    });
    match endpoint_url {
        Ok(endpoint_url) => {
            tracing::trace!("will use endpoint url {:?}", &endpoint_url);
            Ok(aws_sdk_s3::config::Builder::from(&config)
                .endpoint_url(endpoint_url)
                .force_path_style(true)
                .build())
        }
        Err(e) => Err(e),
    }
}

pub async fn upload_file(src: &str, dst: &str) -> Result<(), anyhow::Error> {
    let client = aws_sdk_s3::Client::from_conf(config_from_env().await?);

    let (bucket, key) = if let Some((bucket, key)) = dst.split_once('/') {
        (bucket.to_string(), key.to_string())
    } else {
        anyhow::bail!("invalid S3 path: {}", dst);
    };

    tracing::debug!("will upload to bucket {:?} and key {:?}", &bucket, &key);

    // // Create bucket if it does not exist yet.
    // let exists = client
    //     .bucket_exists(&minio::s3::args::BucketExistsArgs::new(&bucket).unwrap())
    //     .await
    //     .unwrap();
    // tracing::debug!("!!");
    // if !exists {
    //     client
    //         .make_bucket(&minio::s3::args::MakeBucketArgs::new(&bucket).unwrap())
    //         .await
    //         .unwrap();
    // }
    // tracing::debug!("!!");

    let body = aws_sdk_s3::primitives::ByteStream::from_path(std::path::Path::new(src))
        .await
        .map_err(|e| anyhow::anyhow!("could not open file {:?}: {}", src, e))?;
    client
        .put_object()
        .bucket(bucket)
        .key(key)
        .body(body)
        .send()
        .await
        .map_err(|e| anyhow::anyhow!("could not upload file {:?}: {}", src, e))?;

    Ok(())
}

/// Helper struct to encapsulate VCF S3 file upload and TBI creation.
pub struct OutputPathHelper {
    /// Temporary directory to use.
    #[allow(dead_code)] // keep around for RAII
    tmpdir: tempfile::TempDir,
    /// Original output path.
    path_out_orig: String,
    /// Effective output path.
    path_out_effective: String,
}

impl OutputPathHelper {
    pub fn new(path_out: &str) -> Result<Self, anyhow::Error> {
        let tmpdir = tempfile::tempdir().map_err(|e| {
            anyhow::anyhow!("could not create temporary directory for S3 upload: {}", e)
        })?;
        Ok(Self {
            path_out_orig: path_out.to_string(),
            path_out_effective: if s3_mode() {
                tracing::debug!(
                    "S3 mode, using temporary directory: {}",
                    tmpdir.path().display()
                );
                let p = std::path::Path::new(path_out);
                format!(
                    "{}",
                    tmpdir
                        .path()
                        .join(p.file_name().expect("no file name"))
                        .display()
                )
            } else {
                path_out.to_string()
            },
            tmpdir,
        })
    }

    /// Return output path.
    pub fn path_out(&self) -> &str {
        &self.path_out_effective
    }

    /// Create TBI file if necessary.
    pub async fn create_tbi_for_bgzf(&self) -> Result<(), anyhow::Error> {
        if is_gz(&self.path_out_orig) {
            tracing::info!("Creating TBI index for BGZF VCF file...");
            crate::common::noodles::build_tbi(
                &self.path_out_effective,
                &format!("{}.tbi", &self.path_out_effective),
            )
            .await
            .map_err(|e| anyhow::anyhow!("problem building TBI: {}", e))?;
            tracing::info!("... done writing TBI index");
        } else {
            tracing::info!("(not building TBI index for plain text VCF file");
        }

        Ok(())
    }

    /// Upload to S3 if necessary.
    pub async fn upload_for_s3(&self) -> Result<(), anyhow::Error> {
        if s3_mode() {
            tracing::info!("Uploading to S3...");
            upload_file(&self.path_out_effective, &self.path_out_orig).await?;
            if is_gz(&self.path_out_orig) {
                upload_file(
                    &format!("{}.tbi", &self.path_out_effective),
                    &format!("{}.tbi", &self.path_out_orig),
                )
                .await?;
            }
            tracing::info!("... done uploading to S3");
        }

        Ok(())
    }
}
