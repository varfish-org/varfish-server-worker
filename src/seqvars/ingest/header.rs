use noodles_vcf as vcf;

use crate::common::GenomeRelease;

/// Enumeration for the known variant callers.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub enum VariantCaller {
    GatkHaplotypeCaller { version: String },
    GatkUnifiedGenotyper { version: String },
    Dragen { version: String },
    Other,
}

impl VariantCaller {
    pub fn guess(header: &vcf::Header) -> Option<Self> {
        for (other, collection) in header.other_records() {
            if other.as_ref().starts_with("GATKCommandLine")
                || other.as_ref().starts_with("DRAGENCommandLine")
            {
                use vcf::header::record::value::collection::Collection;
                if let Collection::Structured(map) = collection {
                    for (key, values) in map.iter() {
                        if let ("HaplotypeCaller", Some(version)) =
                            (key.as_str(), values.other_fields().get("Version").cloned())
                        {
                            return Some(VariantCaller::GatkHaplotypeCaller { version });
                        }
                        if let ("UnifiedGenotyper", Some(version)) =
                            (key.as_str(), values.other_fields().get("Version").cloned())
                        {
                            return Some(VariantCaller::GatkUnifiedGenotyper { version });
                        }
                        if let ("dragen", Some(version)) =
                            (key.as_str(), values.other_fields().get("Version").cloned())
                        {
                            return Some(VariantCaller::Dragen { version });
                        }
                    }
                }
            }
        }
        None
    }
}

/// Generate the output header from the input header.
pub fn build_output_header(
    input_header: &vcf::Header,
    genomebuild: GenomeRelease,
) -> Result<vcf::Header, anyhow::Error> {
    let variant_caller = VariantCaller::guess(input_header)
        .ok_or_else(|| anyhow::anyhow!("Unable to guess variant caller"))?;
    todo!()
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::VariantCaller;

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
    fn variant_caller_guess(#[case] path: &str) -> Result<(), anyhow::Error> {
        set_snapshot_suffix!("{:?}", path.split('/').last().unwrap());

        let vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;

        insta::assert_yaml_snapshot!(VariantCaller::guess(&vcf_header));

        Ok(())
    }

    #[rstest]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn build_output_header(#[case] path: &str) -> Result<(), anyhow::Error> {
        set_snapshot_suffix!("{:?}", path.split('/').last().unwrap());
        let tmpdir = temp_testdir::TempDir::default();

        let input_vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;
        let output_vcf_header =
            super::build_output_header(&input_vcf_header, crate::common::GenomeRelease::Grch37)?;

        let out_path = tmpdir.join("out.vcf");
        let out_path_str = out_path.to_str().expect("invalid path");
        {
            noodles_vcf::writer::Writer::new(std::fs::File::create(out_path_str)?)
                .write_header(&input_vcf_header)?;
        }

        insta::assert_snapshot!(std::fs::read_to_string(out_path_str)?);

        Ok(())
    }
}
