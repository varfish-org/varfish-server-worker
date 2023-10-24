//! Helper code for working with S3.

/// Helper that returns whether S3 mode has been enabled via `AWS_ACCESS_KEY_ID`.
pub fn s3_mode() -> bool {
    std::env::var("AWS_ACCESS_KEY_ID").is_ok()
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

    let (bucket, key) = if let Some((bucket, key)) = dst.split_once("/") {
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
