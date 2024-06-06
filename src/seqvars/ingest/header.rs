use std::collections::HashSet;

use noodles_vcf as vcf;

use crate::common::GenomeRelease;

/// Enumeration for the known variant callers.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub enum VariantCaller {
    GatkHaplotypeCaller {
        version: String,
    },
    GatkUnifiedGenotyper {
        version: String,
    },
    Glnexus {
        version: String,
        config_name: Option<String>,
    },
    Dragen {
        version: String,
    },
    Other,
}

impl VariantCaller {
    /// Return a string with the name of the variant caller for the VCF header string.
    fn name(&self) -> &'static str {
        match self {
            VariantCaller::GatkHaplotypeCaller { .. } => "GatkHaplotypeCaller",
            VariantCaller::GatkUnifiedGenotyper { .. } => "GatkUnifiedGenotyper",
            VariantCaller::Dragen { .. } => "Dragen",
            VariantCaller::Glnexus { .. } => "Glnexus",
            VariantCaller::Other => "Other",
        }
    }
}

impl VariantCaller {
    pub fn guess(header: &vcf::Header) -> Option<Self> {
        use vcf::header::record::value::collection::Collection;

        let mut glnexus_version: Option<String> = None;
        let mut glnexus_config_name: Option<String> = None;

        for (other, collection) in header.other_records() {
            if ["GATKCommandLine", "DRAGENCommandLine"]
                .iter()
                .any(|k| other.as_ref().starts_with(k))
            {
                if let Collection::Structured(map) = collection {
                    for (key, values) in map.iter() {
                        match (key.as_str(), values.other_fields().get("Version").cloned()) {
                            ("HaplotypeCaller", Some(version)) => {
                                return Some(VariantCaller::GatkHaplotypeCaller { version })
                            }
                            ("UnifiedGenotyper", Some(version)) => {
                                return Some(VariantCaller::GatkUnifiedGenotyper { version })
                            }
                            ("dragen", Some(version)) => {
                                return Some(VariantCaller::Dragen { version })
                            }
                            _ => (),
                        }
                    }
                }
            } else if other.as_ref().starts_with("GLnexusVersion") {
                if let Collection::Unstructured(values) = collection {
                    glnexus_version = Some(values[0].clone());
                }
            } else if other.as_ref().starts_with("GLnexusConfigName") {
                if let Collection::Unstructured(values) = collection {
                    glnexus_config_name = Some(values[0].clone());
                }
            }
        }

        if let Some(version) = glnexus_version {
            return Some(VariantCaller::Glnexus {
                version,
                config_name: glnexus_config_name,
            });
        }

        None
    }
}

/// Add contigs for GRCh37.
fn add_contigs_37(builder: vcf::header::Builder) -> Result<vcf::header::Builder, anyhow::Error> {
    use vcf::header::record::value::map::Contig;
    use vcf::header::record::value::Map;

    let mut builder = builder;

    let specs: &[(&str, usize); 25] = &[
        ("1", 249250621),
        ("2", 243199373),
        ("3", 198022430),
        ("4", 191154276),
        ("5", 180915260),
        ("6", 171115067),
        ("7", 159138663),
        ("8", 146364022),
        ("9", 141213431),
        ("10", 135534747),
        ("11", 135006516),
        ("12", 133851895),
        ("13", 115169878),
        ("14", 107349540),
        ("15", 102531392),
        ("16", 90354753),
        ("17", 81195210),
        ("18", 78077248),
        ("19", 59128983),
        ("20", 63025520),
        ("21", 48129895),
        ("22", 51304566),
        ("X", 155270560),
        ("Y", 59373566),
        ("MT", 16569),
    ];

    for (contig, length) in specs {
        builder = builder.add_contig(
            contig
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid contig: {}", contig))?,
            Map::<Contig>::builder()
                .set_length(*length)
                .insert(
                    "assembly"
                        .parse()
                        .map_err(|_| anyhow::anyhow!("invalid key: assembly"))?,
                    "GRCh37",
                )
                .insert(
                    "species"
                        .parse()
                        .map_err(|_| anyhow::anyhow!("invalid key: species"))?,
                    "Homo sapiens",
                )
                .build()?,
        );
    }

    Ok(builder)
}

/// Add contigs for GRCh38.
fn add_contigs_38(builder: vcf::header::Builder) -> Result<vcf::header::Builder, anyhow::Error> {
    use vcf::header::record::value::map::Contig;
    use vcf::header::record::value::Map;

    let mut builder = builder;

    let specs: &[(&str, usize); 25] = &[
        ("chr1", 248956422),
        ("chr2", 242193529),
        ("chr3", 198295559),
        ("chr4", 190214555),
        ("chr5", 181538259),
        ("chr6", 170805979),
        ("chr7", 159345973),
        ("chr8", 145138636),
        ("chr9", 138394717),
        ("chr10", 133797422),
        ("chr11", 135086622),
        ("chr12", 133275309),
        ("chr13", 114364328),
        ("chr14", 107043718),
        ("chr15", 101991189),
        ("chr16", 90338345),
        ("chr17", 83257441),
        ("chr18", 80373285),
        ("chr19", 58617616),
        ("chr20", 64444167),
        ("chr21", 46709983),
        ("chr22", 50818468),
        ("chrX", 156040895),
        ("chrY", 57227415),
        ("chrM", 16569),
    ];

    for (contig, length) in specs {
        builder = builder.add_contig(
            contig
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid contig: {}", contig))?,
            Map::<Contig>::builder()
                .set_length(*length)
                .insert(
                    "assembly"
                        .parse()
                        .map_err(|_| anyhow::anyhow!("invalid key: assembly"))?,
                    "GRCh38",
                )
                .insert(
                    "species"
                        .parse()
                        .map_err(|_| anyhow::anyhow!("invalid key: species"))?,
                    "Homo sapiens",
                )
                .build()?,
        );
    }

    Ok(builder)
}

/// Generate the output header from the input header.
pub fn build_output_header(
    input_header: &vcf::Header,
    pedigree: &Option<mehari::ped::PedigreeByName>,
    id_mapping: &Option<indexmap::IndexMap<String, String>>,
    genomebuild: GenomeRelease,
    file_date: &str,
    case_uuid: &uuid::Uuid,
    worker_version: &str,
) -> Result<vcf::Header, anyhow::Error> {
    use vcf::header::record::value::{
        map::{info::Type, Filter, Format, Info},
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
            "gnomad_exomes_an".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of alleles in gnomAD exomes",
            ),
        )
        .add_info(
            "gnomad_exomes_hom".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of hom. alt. carriers in gnomAD exomes",
            ),
        )
        .add_info(
            "gnomad_exomes_het".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of het. alt. carriers in gnomAD exomes",
            ),
        )
        .add_info(
            "gnomad_exomes_hemi".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of hemi. alt. carriers in gnomAD exomes",
            ),
        )
        .add_info(
            "gnomad_genomes_an".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of alleles in gnomAD genomes",
            ),
        )
        .add_info(
            "gnomad_genomes_hom".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of hom. alt. carriers in gnomAD genomes",
            ),
        )
        .add_info(
            "gnomad_genomes_het".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of het. alt. carriers in gnomAD genomes",
            ),
        )
        .add_info(
            "gnomad_genomes_hemi".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of hemi. alt. carriers in gnomAD genomes",
            ),
        )
        .add_info(
            "helix_an".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of alleles in HelixMtDb",
            ),
        )
        .add_info(
            "helix_hom".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of hom. alt. carriers in HelixMtDb",
            ),
        )
        .add_info(
            "helix_het".parse()?,
            Map::<Info>::new(
                Number::Count(1),
                Type::Integer,
                "Number of het. alt. carriers in HelixMtDb",
            ),
        )
        .add_info(
            "ANN".parse()?,
            Map::<Info>::new(
                Number::Unknown,
                Type::String,
                "Functional annotations: 'Allele | Annotation | Annotation_Impact | \
                Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | \
                HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | \
                AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'",
            ),
        )
        .add_format(key::READ_DEPTHS, Map::<Format>::from(&key::READ_DEPTHS))
        .add_format(key::READ_DEPTH, Map::<Format>::from(&key::READ_DEPTH))
        .add_format(
            key::CONDITIONAL_GENOTYPE_QUALITY,
            Map::<Format>::from(&key::CONDITIONAL_GENOTYPE_QUALITY),
        )
        .add_format(key::GENOTYPE, Map::<Format>::from(&key::GENOTYPE))
        .add_format(key::PHASE_SET, Map::<Format>::from(&key::PHASE_SET));

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
        let input_idv = input_header
            .sample_names()
            .iter()
            .cloned()
            .collect::<HashSet<_>>();
        if let Some(id_mapping) = id_mapping {
            let src_idv = id_mapping.keys().cloned().collect::<HashSet<_>>();
            if !src_idv.eq(&input_idv) {
                anyhow::bail!(
                    "mapping source individuals = {:?} != input individuals: {:?}",
                    &src_idv,
                    &input_idv
                )
            }
            let dst_idv = id_mapping.values().cloned().collect::<HashSet<_>>();
            if !dst_idv.eq(&ped_idv) {
                anyhow::bail!(
                    "mapping destination individuals = {:?} != pedigree individuals: {:?}",
                    &dst_idv,
                    &ped_idv
                )
            }
        } else if !ped_idv.eq(&input_idv) {
            anyhow::bail!(
                "pedigree individuals = {:?} != input individuals: {:?}",
                &ped_idv,
                &input_idv
            )
        }

        let mut sample_names = Vec::new();
        for src_name in input_header.sample_names() {
            let dst_name = if let Some(id_mapping) = id_mapping {
                id_mapping.get(src_name).expect("checked above")
            } else {
                src_name
            };

            let i = pedigree.individuals.get(dst_name).expect("checked above");
            if input_header.sample_names().contains(src_name) {
                sample_names.push(dst_name.clone());
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

        sample_names.sort();
        builder = builder.set_sample_names(sample_names.into_iter().collect());
    } else {
        for name in input_header.sample_names() {
            let name = if let Some(id_mapping) = id_mapping {
                id_mapping
                    .get(name)
                    .ok_or_else(|| anyhow::anyhow!("mapping given but sample {} not found", name))?
                    .clone()
            } else {
                name.clone()
            };
            builder = builder.add_sample_name(name);
        }
    }

    use vcf::header::record::value::map::Other;

    let orig_caller = VariantCaller::guess(input_header)
        .ok_or_else(|| anyhow::anyhow!("unable to guess original variant caller"))?;

    let builder = builder
        .insert(
            "x-varfish-case-uuid".parse()?,
            vcf::header::record::Value::String(case_uuid.to_string()),
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

    let builder = match &orig_caller {
        VariantCaller::GatkHaplotypeCaller { version }
        | VariantCaller::GatkUnifiedGenotyper { version }
        | VariantCaller::Dragen { version } => builder.insert(
            "x-varfish-version".parse()?,
            vcf::header::record::Value::Map(
                String::from("orig-caller"),
                Map::<Other>::builder()
                    .insert("Name".parse()?, orig_caller.name())
                    .insert("Version".parse()?, version)
                    .build()?,
            ),
        )?,
        VariantCaller::Glnexus {
            version,
            config_name,
        } => builder.insert(
            "x-varfish-version".parse()?,
            vcf::header::record::Value::Map(
                String::from("orig-caller"),
                Map::<Other>::builder()
                    .insert("Name".parse()?, orig_caller.name())
                    .insert("Version".parse()?, version)
                    .insert(
                        "ConfigName".parse()?,
                        config_name.clone().unwrap_or_default(),
                    )
                    .build()?,
            ),
        )?,
        VariantCaller::Other => builder.insert(
            "x-varfish-version".parse()?,
            vcf::header::record::Value::Map(
                String::from("orig-caller"),
                Map::<Other>::builder()
                    .insert("Name".parse()?, "Other")
                    .build()?,
            ),
        )?,
    };

    Ok(builder.build())
}

#[cfg(test)]
mod test {
    use mehari::ped::PedigreeByName;
    use noodles_vcf as vcf;
    use rstest::rstest;

    use super::VariantCaller;

    #[rstest]
    #[case("tests/seqvars/ingest/clair3_glnexus.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn variant_caller_guess(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path.split('/').last().unwrap());

        let vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;

        insta::assert_yaml_snapshot!(VariantCaller::guess(&vcf_header));

        Ok(())
    }

    #[rstest]
    #[case("tests/seqvars/ingest/clair3_glnexus.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn build_output_header_37(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path.split('/').last().unwrap());
        let tmpdir = temp_testdir::TempDir::default();

        let pedigree = PedigreeByName::from_path(path.replace(".vcf", ".ped")).unwrap();

        let mut input_vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;
        let output_vcf_header = super::build_output_header(
            &input_vcf_header,
            &Some(pedigree),
            &None,
            crate::common::GenomeRelease::Grch37,
            "20230421",
            &uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
            "x.y.z",
        )?;

        // Work around glnexus issue with RNC.
        if let Some(format) = input_vcf_header.formats_mut().get_mut("RNC") {
            *format.number_mut() = vcf::header::Number::Count(1);
            *format.type_mut() = vcf::header::record::value::map::format::Type::String;
        }

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
    #[case("tests/seqvars/ingest/clair3_glnexus.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.4.vcf")]
    #[case("tests/seqvars/ingest/example_dragen.07.021.624.3.10.9.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.3.7-0.vcf")]
    #[case("tests/seqvars/ingest/example_gatk_hc.4.4.0.0.vcf")]
    fn build_output_header_38(#[case] path: &str) -> Result<(), anyhow::Error> {
        mehari::common::set_snapshot_suffix!("{}", path.split('/').last().unwrap());
        let tmpdir = temp_testdir::TempDir::default();

        let pedigree = PedigreeByName::from_path(path.replace(".vcf", ".ped")).unwrap();

        let mut input_vcf_header = noodles_vcf::reader::Builder::default()
            .build_from_path(path)?
            .read_header()?;
        let output_vcf_header = super::build_output_header(
            &input_vcf_header,
            &Some(pedigree),
            &None,
            crate::common::GenomeRelease::Grch38,
            "20230421",
            &uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
            "x.y.z",
        )?;

        // Work around glnexus issue with RNC.
        if let Some(format) = input_vcf_header.formats_mut().get_mut("RNC") {
            *format.number_mut() = vcf::header::Number::Count(1);
            *format.type_mut() = vcf::header::record::value::map::format::Type::String;
        }

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
