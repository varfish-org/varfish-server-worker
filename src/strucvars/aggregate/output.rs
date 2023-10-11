//! Code for writing the output file.

use mehari::annotate::strucvars::{
    bnd::Breakend, csq::interface::StrandOrientation, PeOrientation,
};
use noodles_vcf as vcf;

use crate::{
    common::{Chrom, Genotype},
    strucvars::query::schema::SvType,
};

/// Representation of the fields for the in-house background database.
#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct Record {
    /// chromosome name
    pub chromosome: String,
    /// start position, 1-based
    pub begin: i32,
    /// chromosome2 name
    pub chromosome2: String,
    /// end position, 1-based
    pub end: i32,
    /// paired-end orientation
    pub pe_orientation: StrandOrientation,
    /// type of the SV
    pub sv_type: SvType,
    /// number of overall carriers
    pub carriers: u32,
    /// number of het. carriers
    pub carriers_het: u32,
    /// number of hom. carriers
    pub carriers_hom: u32,
    /// number of hemi. carriers
    pub carriers_hemi: u32,
}

impl Record {
    /// Compute reciprocal overlap between `self` and `other`.
    pub fn overlap(&self, other: &Record) -> f32 {
        let s1 = if self.begin > 0 { self.begin - 1 } else { 0 };
        let e1 = self.end + 1;
        let s2 = if other.begin > 0 { other.begin - 1 } else { 0 };
        let e2 = other.end;

        let ovl_s = std::cmp::max(s1, s2);
        let ovl_e = std::cmp::min(e1, e2);
        if ovl_e <= ovl_s {
            0.0
        } else {
            let len1 = (e1 - s1) as f32;
            let len2 = (e2 - s2) as f32;
            let ovl_len = (ovl_e - ovl_s) as f32;
            (ovl_len / len1).min(ovl_len / len2)
        }
    }

    pub fn merge_into(&mut self, other: &Record) {
        self.carriers += other.carriers;
        self.carriers_het += other.carriers_het;
        self.carriers_hom += other.carriers_hom;
        self.carriers_hemi += other.carriers_hemi;
    }

    /// Convert VCF record into a `Record`.
    pub fn from_vcf(
        record: &vcf::Record,
        header: &vcf::Header,
        _genome_release: crate::common::GenomeRelease,
        pedigree: &mehari::ped::PedigreeByName,
    ) -> Result<Self, anyhow::Error> {
        let chromosome = record.chromosome().to_string();
        let begin = {
            let position: usize = record.position().into();
            position.saturating_sub(1) as i32
        };
        let chromosome2 = if let Some(Some(vcf::record::info::field::Value::String(chromosome2))) =
            record.info().get(
                &"chr2"
                    .parse::<vcf::record::info::field::Key>()
                    .expect("invalid key chr2?"),
            ) {
            chromosome2.clone()
        } else {
            chromosome.clone()
        };
        let end = if let Some(Some(vcf::record::info::field::Value::Integer(end))) = record
            .info()
            .get(&vcf::record::info::field::key::END_POSITION)
        {
            *end
        } else {
            anyhow::bail!("missing INFO/END")
        };

        let sv_type = if let Some(Some(vcf::record::info::field::Value::String(sv_type))) =
            record.info().get(&vcf::record::info::field::key::SV_TYPE)
        {
            sv_type.parse()?
        } else {
            anyhow::bail!("missing INFO/SVTYPE")
        };

        let pe_orientation = match sv_type {
            SvType::Del => StrandOrientation::ThreeToFive,
            SvType::Dup => StrandOrientation::FiveToThree,
            SvType::Inv => StrandOrientation::FiveToFive,
            SvType::Bnd => {
                let pe_orientation = Breakend::from_ref_alt_str(
                    &record.reference_bases().to_string(),
                    &record
                        .alternate_bases()
                        .first()
                        .ok_or_else(|| anyhow::anyhow!("no alternate allele?"))?
                        .to_string(),
                )?
                .pe_orientation;
                match pe_orientation {
                    PeOrientation::ThreeToThree => StrandOrientation::ThreeToThree,
                    PeOrientation::FiveToFive => StrandOrientation::FiveToFive,
                    PeOrientation::ThreeToFive => StrandOrientation::ThreeToFive,
                    PeOrientation::FiveToThree => StrandOrientation::FiveToThree,
                    PeOrientation::Other => StrandOrientation::NotApplicable,
                }
            }
            SvType::Ins | SvType::Cnv => StrandOrientation::NotApplicable,
        };

        let chrom: Chrom =
            annonars::common::cli::canonicalize(record.chromosome().to_string().as_str())
                .as_str()
                .parse()?;

        let (mut carriers_het, mut carriers_hom, mut carriers_hemi) = (0, 0, 0);
        for (name, sample) in header
            .sample_names()
            .iter()
            .zip(record.genotypes().values())
        {
            let individual = pedigree
                .individuals
                .get(name)
                .ok_or_else(|| anyhow::anyhow!("individual {} not found in pedigree", name))?;

            let genotype: Genotype =
                if let Some(Some(gt)) = sample.get(&vcf::record::genotypes::keys::key::GENOTYPE) {
                    if let vcf::record::genotypes::sample::Value::String(gt) = gt {
                        gt.as_str().parse()?
                    } else {
                        anyhow::bail!("invalid genotype value: {:?}", gt);
                    }
                } else {
                    continue; // skip, no-call or empty
                };

            match (chrom, individual.sex, genotype) {
                (_, _, Genotype::WithNoCall) => continue,
                // on the autosomes, male/female count the same
                (Chrom::Auto, _, Genotype::HomRef) => (),
                (Chrom::Auto, _, Genotype::Het) => {
                    carriers_het += 1;
                }
                (Chrom::Auto, _, Genotype::HomAlt) => {
                    carriers_hom += 1;
                }
                // on the gonomosomes, we handle call male variant calls as hemizygous
                (Chrom::X, mehari::ped::Sex::Male, Genotype::HomRef)
                | (Chrom::Y, mehari::ped::Sex::Male, Genotype::HomRef) => (),
                (Chrom::X, mehari::ped::Sex::Male, Genotype::HomAlt)
                | (Chrom::X, mehari::ped::Sex::Male, Genotype::Het)
                | (Chrom::Y, mehari::ped::Sex::Male, Genotype::Het)
                | (Chrom::Y, mehari::ped::Sex::Male, Genotype::HomAlt) => {
                    carriers_hemi += 1;
                }
                // for female samples, we handle chrX as biallelic
                (Chrom::X, mehari::ped::Sex::Female, Genotype::HomRef)
                | (Chrom::X, mehari::ped::Sex::Female, Genotype::Het) => {
                    carriers_het += 1;
                }
                (Chrom::X, mehari::ped::Sex::Female, Genotype::HomAlt) => {
                    carriers_hom += 1;
                }
                // we ignore calls to chrY for female samples
                (Chrom::Y, mehari::ped::Sex::Female, Genotype::HomRef)
                | (Chrom::Y, mehari::ped::Sex::Female, Genotype::Het)
                | (Chrom::Y, mehari::ped::Sex::Female, Genotype::HomAlt) => (),
                // do not count samples with unknown sex on gonomosomes
                (Chrom::X, mehari::ped::Sex::Unknown, _)
                | (Chrom::Y, mehari::ped::Sex::Unknown, _) => (),
            };
        }

        Ok(Self {
            chromosome,
            begin,
            chromosome2,
            end,
            pe_orientation,
            sv_type,
            carriers_het,
            carriers_hom,
            carriers_hemi,
            carriers: carriers_het + carriers_hom + carriers_hemi,
        })
    }
}
