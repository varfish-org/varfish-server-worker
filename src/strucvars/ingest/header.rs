use std::collections::HashSet;

use noodles_vcf as vcf;
use vcf::header::{record::value::map::AlternativeAllele, SampleNames};

use crate::common::{add_contigs_37, add_contigs_38, GenomeRelease};

/// Return token for caller name.
fn caller_name(sv_caller: &mehari::annotate::strucvars::SvCaller) -> &'static str {
    match sv_caller {
        mehari::annotate::strucvars::SvCaller::Delly { .. } => "Delly",
        mehari::annotate::strucvars::SvCaller::DragenSv { .. } => "DragenSv",
        mehari::annotate::strucvars::SvCaller::DragenCnv { .. } => "DragenCnv",
        mehari::annotate::strucvars::SvCaller::Gcnv { .. } => "Gcnv",
        mehari::annotate::strucvars::SvCaller::Manta { .. } => "Manta",
        mehari::annotate::strucvars::SvCaller::Melt { .. } => "Melt",
        mehari::annotate::strucvars::SvCaller::Popdel { .. } => "Popdel",
        mehari::annotate::strucvars::SvCaller::ClinCnv { .. } => "ClinCnv",
        mehari::annotate::strucvars::SvCaller::Sniffles2 { .. } => "Sniffles2",
    }
}

/// Return caller version.
fn caller_version(sv_caller: &mehari::annotate::strucvars::SvCaller) -> String {
    match sv_caller {
        mehari::annotate::strucvars::SvCaller::Delly { version }
        | mehari::annotate::strucvars::SvCaller::DragenSv { version }
        | mehari::annotate::strucvars::SvCaller::DragenCnv { version }
        | mehari::annotate::strucvars::SvCaller::Gcnv { version }
        | mehari::annotate::strucvars::SvCaller::Manta { version }
        | mehari::annotate::strucvars::SvCaller::Melt { version }
        | mehari::annotate::strucvars::SvCaller::Popdel { version }
        | mehari::annotate::strucvars::SvCaller::ClinCnv { version }
        | mehari::annotate::strucvars::SvCaller::Sniffles2 { version } => version.clone(),
    }
}

/// Generate the output header from the input header.
pub fn build_output_header(
    input_sample_names: &SampleNames,
    input_sv_callers: &[&mehari::annotate::strucvars::SvCaller],
    pedigree: Option<&mehari::ped::PedigreeByName>,
    genomebuild: GenomeRelease,
    file_date: &str,
    worker_version: &str,
    case_uuid: &str,
) -> Result<vcf::Header, anyhow::Error> {
    use vcf::header::record::value::{
        map::{format, info, Filter, Format, Info},
        Map,
    };
    use vcf::header::Number;
    use vcf::record::genotypes::keys::key;

    let builder = vcf::Header::builder()
        .insert(
            "fileDate".parse()?,
            vcf::header::record::Value::from(file_date),
        )?
        .add_filter("PASS", Map::<Filter>::new("All filters passed"))
        .add_info(
            vcf::record::info::field::key::IS_IMPRECISE,
            Map::<Info>::from(&vcf::record::info::field::key::IS_IMPRECISE),
        )
        .add_info(
            vcf::record::info::field::key::END_POSITION,
            Map::<Info>::from(&vcf::record::info::field::key::END_POSITION),
        )
        .add_info(
            vcf::record::info::field::key::SV_TYPE,
            Map::<Info>::from(&vcf::record::info::field::key::SV_TYPE),
        )
        .add_info(
            vcf::record::info::field::key::SV_LENGTHS,
            Map::<Info>::from(&vcf::record::info::field::key::SV_LENGTHS),
        )
        .add_info(
            vcf::record::info::field::key::SV_CLAIM,
            Map::<Info>::from(&vcf::record::info::field::key::SV_CLAIM),
        )
        .add_info(
            "callers".parse()?,
            Map::<Info>::new(
                Number::Unknown,
                info::Type::String,
                "Callers that called the variant",
            ),
        )
        .add_info(
            "chr2".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                info::Type::String,
                "Second chromosome, if not equal to CHROM",
            ),
        )
        .add_info(
            "annsv".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                info::Type::String,
                "Effect annotations: 'Allele | Annotation | Gene_Name | Gene_ID'",
            ),
        )
        .add_format(
            key::CONDITIONAL_GENOTYPE_QUALITY,
            Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
        )
        .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
        .add_format(
            "pec".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Total coverage with paired-end reads",
            ),
        )
        .add_format(
            "pev".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Paired-end reads supporting the variant",
            ),
        )
        .add_format(
            "src".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Total coverage with split reads",
            ),
        )
        .add_format(
            "srv".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Split reads supporting the variant",
            ),
        )
        .add_format(
            "amq".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Float,
                "Average mapping quality over the variant",
            ),
        )
        .add_format(
            "cn".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Copy number of the variant in the sample",
            ),
        )
        .add_format(
            "anc".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Float,
                "Average normalized coverage over the variant in the sample",
            ),
        )
        .add_format(
            "pc".parse()?,
            Map::<Format>::new(
                Number::Count(1),
                format::Type::Integer,
                "Point count (windows/targets/probes)",
            ),
        )
        .add_alternative_allele("DEL".parse()?, Map::<AlternativeAllele>::new("Deletion"))
        .add_alternative_allele("DUP".parse()?, Map::<AlternativeAllele>::new("Duplication"))
        .add_alternative_allele("INS".parse()?, Map::<AlternativeAllele>::new("Insertion"))
        .add_alternative_allele(
            "CNV".parse()?,
            Map::<AlternativeAllele>::new("Copy Number Variation"),
        )
        .add_alternative_allele("INV".parse()?, Map::<AlternativeAllele>::new("Inversion"));

    let mut builder = match genomebuild {
        GenomeRelease::Grch37 => add_contigs_37(builder),
        GenomeRelease::Grch38 => add_contigs_38(builder),
    }
    .map_err(|e| anyhow::anyhow!("problem adding contigs: {}", e))?;

    if let Some(pedigree) = pedigree {
        let ped_idv = pedigree
            .individuals
            .iter()
            .map(|(name, _)| name.clone())
            .collect::<HashSet<_>>();
        let input_idv = input_sample_names.iter().cloned().collect::<HashSet<_>>();
        if !ped_idv.eq(&input_idv) {
            anyhow::bail!(
                "pedigree individuals = {:?} != input individuals: {:?}",
                &ped_idv,
                &input_idv
            )
        }

        for name in input_sample_names {
            let i = pedigree
                .individuals
                .get(name)
                .expect("checked equality above");
            if input_sample_names.contains(&i.name) {
                builder = builder.add_sample_name(i.name.clone());
            }

            // Add SAMPLE entry.
            builder = builder.insert(
                "SAMPLE".parse()?,
                noodles_vcf::header::record::Value::Map(
                    i.name.clone(),
                    Map::<Other>::builder()
                        .insert(
                            "Sex".parse()?,
                            mehari::annotate::strucvars::vcf_header::sex_str(i.sex),
                        )
                        .insert(
                            "Disease".parse()?,
                            mehari::annotate::strucvars::vcf_header::disease_str(i.disease),
                        )
                        .build()?,
                ),
            )?;

            // Add PEDIGREE entry.
            let mut map_builder = Map::<Other>::builder();
            if let Some(father) = i.father.as_ref() {
                map_builder = map_builder.insert("Father".parse()?, father.clone());
            }
            if let Some(mother) = i.mother.as_ref() {
                map_builder = map_builder.insert("Mother".parse()?, mother.clone());
            }
            builder = builder.insert(
                "PEDIGREE".parse()?,
                noodles_vcf::header::record::Value::Map(i.name.clone(), map_builder.build()?),
            )?;
        }
    } else {
        for name in input_sample_names {
            builder = builder.add_sample_name(name.clone());
        }
    }

    use vcf::header::record::value::map::Other;

    let mut builder = builder
        .insert(
            "x-varfish-case-uuid".parse()?,
            vcf::header::record::Value::from(case_uuid),
        )?
        .insert(
            "x-varfish-version".parse()?,
            vcf::header::record::Value::Map(
                String::from("varfish-server-worker"),
                Map::<Other>::builder()
                    .insert("Version".parse()?, worker_version)
                    .build()?,
            ),
        )?;

    for sv_caller in input_sv_callers.iter() {
        builder = builder.insert(
            "x-varfish-version".parse()?,
            vcf::header::record::Value::Map(
                caller_name(sv_caller).into(),
                Map::<Other>::builder()
                    .insert("Name".parse()?, caller_name(sv_caller))
                    .insert("Version".parse()?, caller_version(sv_caller))
                    .build()?,
            ),
        )?;
    }

    Ok(builder.build())
}

#[cfg(test)]
mod test {
    use mehari::ped::PedigreeByName;
    use rstest::rstest;

    #[rstest]
    #[case("tests/strucvars/ingest/delly2-min.vcf")]
    #[case("tests/strucvars/ingest/dragen-cnv-min.vcf")]
    #[case("tests/strucvars/ingest/dragen-sv-min.vcf")]
    #[case("tests/strucvars/ingest/gcnv-min.vcf")]
    #[case("tests/strucvars/ingest/manta-min.vcf")]
    #[case("tests/strucvars/ingest/melt-min.vcf")]
    #[case("tests/strucvars/ingest/popdel-min.vcf")]
    #[case("tests/strucvars/ingest/sniffles2-min.vcf")]
    #[tokio::test]
    async fn build_output_header_37(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path.split('/').last().unwrap());
        let tmpdir = temp_testdir::TempDir::default();

        let pedigree = PedigreeByName::from_path(path.replace(".vcf", ".ped")).unwrap();

        let input_vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;
        let sv_callers = {
            let mut reader = mehari::common::noodles::open_vcf_reader(path).await?;
            vec![mehari::annotate::strucvars::guess_sv_caller(&mut reader).await?]
        };
        let sv_caller_refs = sv_callers.iter().collect::<Vec<_>>();
        let output_vcf_header = super::build_output_header(
            input_vcf_header.sample_names(),
            &sv_caller_refs,
            Some(&pedigree),
            crate::common::GenomeRelease::Grch37,
            "20230421",
            "x.y.z",
            "d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c",
        )?;

        let out_path = tmpdir.join("out.vcf");
        let out_path_str = out_path.to_str().expect("invalid path");
        {
            noodles_vcf::writer::Writer::new(std::fs::File::create(out_path_str)?)
                .write_header(&output_vcf_header)?;
        }

        insta::assert_snapshot!(std::fs::read_to_string(out_path_str)?);

        Ok(())
    }

    #[rstest]
    #[case("tests/strucvars/ingest/delly2-min.vcf")]
    #[case("tests/strucvars/ingest/dragen-cnv-min.vcf")]
    #[case("tests/strucvars/ingest/dragen-sv-min.vcf")]
    #[case("tests/strucvars/ingest/gcnv-min.vcf")]
    #[case("tests/strucvars/ingest/manta-min.vcf")]
    #[case("tests/strucvars/ingest/melt-min.vcf")]
    #[case("tests/strucvars/ingest/popdel-min.vcf")]
    #[case("tests/strucvars/ingest/sniffles2-min.vcf")]
    #[tokio::test]
    async fn build_output_header_38(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path.split('/').last().unwrap());
        let tmpdir = temp_testdir::TempDir::default();

        let pedigree = PedigreeByName::from_path(path.replace(".vcf", ".ped")).unwrap();

        let input_vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;
        let sv_callers = {
            let mut reader = mehari::common::noodles::open_vcf_reader(path).await?;
            vec![mehari::annotate::strucvars::guess_sv_caller(&mut reader).await?]
        };
        let sv_caller_refs = sv_callers.iter().collect::<Vec<_>>();
        let output_vcf_header = super::build_output_header(
            input_vcf_header.sample_names(),
            &sv_caller_refs,
            Some(&pedigree),
            crate::common::GenomeRelease::Grch38,
            "20230421",
            "x.y.z",
            "d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c",
        )?;

        let out_path = tmpdir.join("out.vcf");
        let out_path_str = out_path.to_str().expect("invalid path");
        {
            noodles_vcf::writer::Writer::new(std::fs::File::create(out_path_str)?)
                .write_header(&output_vcf_header)?;
        }

        insta::assert_snapshot!(std::fs::read_to_string(out_path_str)?);

        Ok(())
    }
}
