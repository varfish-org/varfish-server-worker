//! Common utility code for noodles.

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use noodles_vcf as vcf;

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
    indexer.set_header(csi::index::header::Builder::vcf().build());

    let mut start_position = reader.get_ref().virtual_position();

    while reader
        .read_record(&header, &mut record)
        .await
        .map_err(|e| anyhow::anyhow!("problem reading record: {}", e))?
        != 0
    {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record.chromosome().to_string();
        let start = Position::try_from(usize::from(record.position()))
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
            .map_err(|e| anyhow::anyhow!("error converting start position: {}", e))?;
        let end = record
            .end()
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
