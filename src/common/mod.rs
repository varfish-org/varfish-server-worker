//! Common functionality.

use std::ops::Range;

use biocommons_bioutils::assemblies::Assembly;
use byte_unit::Byte;
use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use indexmap::IndexMap;
use noodles_vcf as vcf;

pub mod noodles;
pub mod s3;

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            verbose: Verbosity::new(0, 0),
        }
    }
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
    #[strum(serialize = "grch37")]
    Grch37,
    /// GRCh38 / hg38
    #[strum(serialize = "grch38")]
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

/// Helper type for encoding genotypes in parsing.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Genotype {
    /// hom. ref.
    HomRef,
    /// het.
    Het,
    /// hom. alt.
    HomAlt,
    /// other, includes no-call
    WithNoCall,
}

impl std::str::FromStr for Genotype {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "0/0" | "0|0" | "0" => Genotype::HomRef,
            "0/1" | "1/0" | "0|1" | "1|0" => Genotype::Het,
            "1/1" | "1|1" | "1" => Genotype::HomAlt,
            "./." | "./0" | "./1" | "0/." | "1/." => Genotype::WithNoCall,
            _ => anyhow::bail!("invalid genotype value: {:?}", s),
        })
    }
}

#[derive(Copy, Clone)]
pub enum Chrom {
    Auto, // or chrMT, but does not matter for carrier computation
    X,
    Y,
}

impl std::str::FromStr for Chrom {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "X" => Chrom::X,
            "Y" => Chrom::Y,
            _ => Chrom::Auto,
        })
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

/// Return the version of the `varfish-server-worker` crate and `x.y.z` in tests.
pub fn worker_version() -> &'static str {
    if cfg!(test) {
        "x.y.z"
    } else {
        env!("CARGO_PKG_VERSION")
    }
}

/// Add contigs for GRCh37.
pub fn add_contigs_37(
    builder: vcf::header::Builder,
) -> Result<vcf::header::Builder, anyhow::Error> {
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

    builder = builder.insert(
        "x-varfish-genome-build".parse()?,
        vcf::header::record::Value::from("GRCh37"),
    )?;

    Ok(builder)
}

/// Extract a PedigreeByName from the VCF header.
pub fn extract_pedigree_and_case_uuid(
    header: &vcf::Header,
) -> Result<(mehari::ped::PedigreeByName, uuid::Uuid), anyhow::Error> {
    let mut case_uuid = uuid::Uuid::nil();
    let mut pedigree = mehari::ped::PedigreeByName::default();

    if let vcf::header::record::value::Collection::Structured(sample_map) = header
        .other_records()
        .get("SAMPLE")
        .ok_or_else(|| anyhow::anyhow!("no SAMPLE record in VCF header"))?
    {
        for (sample_name, sample_values) in sample_map.iter() {
            let sex_value = sample_values
                .other_fields()
                .get("Sex")
                .ok_or_else(|| anyhow::anyhow!("no Sex field in SAMPLE header?"))?
                .as_ref();
            let sex_value = match sex_value {
                "Male" => mehari::ped::Sex::Male,
                "Female" => mehari::ped::Sex::Female,
                "Unknown" => mehari::ped::Sex::Unknown,
                _ => anyhow::bail!("invalid value for Sex: {}", sex_value),
            };

            let disease_value = sample_values
                .other_fields()
                .get("Disease")
                .ok_or_else(|| anyhow::anyhow!("no Disease field in SAMPLE header?"))?
                .as_ref();
            let disease_value = match disease_value {
                "Affected" => mehari::ped::Disease::Affected,
                "Unaffected" => mehari::ped::Disease::Unaffected,
                "Unknown" => mehari::ped::Disease::Unknown,
                _ => anyhow::bail!("invalid value for Disease: {}", disease_value),
            };

            pedigree.individuals.insert(
                sample_name.clone(),
                mehari::ped::Individual {
                    family: "FAM".into(),
                    name: sample_name.clone(),
                    sex: sex_value,
                    disease: disease_value,
                    ..Default::default()
                },
            );
        }
    }

    if let vcf::header::record::value::Collection::Structured(sample_map) = header
        .other_records()
        .get("PEDIGREE")
        .ok_or_else(|| anyhow::anyhow!("no PEDIGREE record in VCF header"))?
    {
        for (sample_name, pedigree_values) in sample_map.iter() {
            let father_value = pedigree_values.other_fields().get("Father");
            let mother_value = pedigree_values.other_fields().get("Mother");

            let individual = pedigree.individuals.get_mut(sample_name).ok_or_else(|| {
                anyhow::anyhow!("individual {} not found in SAMPLE header", sample_name)
            })?;
            if let Some(father_value) = father_value {
                individual.father = Some(father_value.clone());
            }
            if let Some(mother_value) = mother_value {
                individual.mother = Some(mother_value.clone());
            }
        }
    }

    if let vcf::header::record::value::Collection::Unstructured(lines) = header
        .other_records()
        .get("x-varfish-case-uuid")
        .ok_or_else(|| anyhow::anyhow!("no x-varfish-case-uuid record in VCF header"))?
    {
        case_uuid = lines
            .first()
            .ok_or_else(|| {
                anyhow::anyhow!("no x-varfish-case-uuid record in VCF header, but expected one")
            })?
            .parse()
            .map_err(|e| {
                anyhow::anyhow!(
                    "could not parse x-varfish-case-uuid record in VCF header: {}",
                    e
                )
            })?;
    }

    Ok((pedigree, case_uuid))
}

/// Add contigs for GRCh38.
pub fn add_contigs_38(
    builder: vcf::header::Builder,
) -> Result<vcf::header::Builder, anyhow::Error> {
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

    builder = builder.insert(
        "x-varfish-genome-build".parse()?,
        vcf::header::record::Value::from("GRCh38"),
    )?;

    Ok(builder)
}

#[cfg(test)]
pub(crate) fn read_to_bytes<P>(path: P) -> Result<Vec<u8>, anyhow::Error>
where
    P: AsRef<std::path::Path>,
{
    use std::io::Read;

    let mut f = std::fs::File::open(&path).expect("no file found");
    let metadata = std::fs::metadata(&path).expect("unable to read metadata");
    let mut buffer = vec![0; metadata.len() as usize];
    f.read_exact(&mut buffer).expect("buffer overflow");
    Ok(buffer)
}

/// Enum for the supported TADs.
#[derive(
    serde::Serialize,
    serde::Deserialize,
    enum_map::Enum,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Debug,
    Default,
    strum::EnumString,
    strum::Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum TadSet {
    /// hESC
    #[default]
    Hesc,
}

#[macro_export]
macro_rules! flush_and_shutdown {
    ($output_writer:expr) => {
        let mut write = $output_writer.into_inner();
        // Flush all buffers.
        write
            .as_mut()
            .flush()
            .await
            .map_err(|e| anyhow::anyhow!("problem flushing output VCF file: {}", e))?;
        // Shutdown the writer.
        write
            .shutdown()
            .await
            .map_err(|e| anyhow::anyhow!("problem shutting down output VCF file: {}", e))?;
        // Give the OS some time to actually close the file.
        tokio::time::sleep(std::time::Duration::from_millis(200)).await;
    };
}

#[cfg(test)]
mod test {
    use noodles_vcf as vcf;

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

    #[rstest::rstest]
    #[case(crate::common::GenomeRelease::Grch37, "GRCh37")]
    #[case(crate::common::GenomeRelease::Grch38, "GRCh38")]
    fn genome_release_name(#[case] release: super::GenomeRelease, #[case] expected: &str) {
        assert_eq!(expected, release.name());
    }

    #[rstest::rstest]
    #[case(
        crate::common::GenomeRelease::Grch37,
        biocommons_bioutils::assemblies::Assembly::Grch37p10
    )]
    #[case(
        crate::common::GenomeRelease::Grch38,
        biocommons_bioutils::assemblies::Assembly::Grch38
    )]
    fn assembly_from_genome_release(
        #[case] release: super::GenomeRelease,
        #[case] assembly: biocommons_bioutils::assemblies::Assembly,
    ) -> Result<(), anyhow::Error> {
        let res: biocommons_bioutils::assemblies::Assembly = release.into();

        assert_eq!(res, assembly);

        Ok(())
    }

    #[rstest::rstest]
    #[case(
        crate::common::GenomeRelease::Grch37,
        biocommons_bioutils::assemblies::Assembly::Grch37
    )]
    #[case(
        crate::common::GenomeRelease::Grch37,
        biocommons_bioutils::assemblies::Assembly::Grch37p10
    )]
    #[case(
        crate::common::GenomeRelease::Grch38,
        biocommons_bioutils::assemblies::Assembly::Grch38
    )]
    fn genome_release_from_assembly(
        #[case] release: super::GenomeRelease,
        #[case] assembly: biocommons_bioutils::assemblies::Assembly,
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

    #[test]
    fn extract_pedigree_snapshot() {
        let path = "tests/seqvars/aggregate/ingest.vcf";
        let mut vcf_reader = vcf::reader::Builder.build_from_path(path).unwrap();
        let header = vcf_reader.read_header().unwrap();

        let (pedigree, case_uuid) = super::extract_pedigree_and_case_uuid(&header).unwrap();
        insta::assert_debug_snapshot!(pedigree);
        insta::assert_debug_snapshot!(case_uuid);
    }
}
