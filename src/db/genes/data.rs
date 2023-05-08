//! Code for storing gene-related data.

use serde::{Deserialize, Serialize};

/// Entry in the genes RocksDB database.
///
/// Note that the HGNC ID is used for the keys, e.g., `"HGNC:5"`.
#[derive(Debug, Serialize, Deserialize)]
pub struct Record {}

/// Code for data from the ACMG secondary findings list.
pub mod acmg_sf {
    use serde::{Deserialize, Serialize};

    /// A record from the ACMG secondary findings list.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Record {
        /// The HGNC ID.
        pub hgnc_id: String,
        /// The Ensembl gene ID.
        pub ensembl_gene_id: String,
        /// The NCBI gene ID.
        pub ncbi_gene_id: String,
        /// The HGNC gene symbol.
        pub gene_symbol: String,
        /// The MIM gene ID.
        pub mim_gene_id: String,
        /// The disease phenotype.
        pub disease_phenotype: String,
        /// The disease MIM id.
        pub disorder_mim: String,
        /// The phenotype category.
        pub phenotype_category: String,
        /// The mode of inheritance.
        pub inheritance: String,
        /// The version of the ACMG SF list of first appearance.
        pub sf_list_version: String,
        /// The variants to report according to ACMG SF.
        pub variants_to_report: String,
    }
}

/// Code for data from the NCBI gene database (aka "Entrez").
pub mod ncbi {
    use serde::{Deserialize, Serialize};

    /// A record from the NCBI gene database.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Record {}
}

/// Code for data from the HGNC database.
pub mod hgnc {
    use std::{fmt::Display, str::FromStr};

    use chrono::naive::NaiveDate;
    use serde::{Deserialize, Serialize};
    use serde_with::DisplayFromStr;

    /// Status of the symbol report, which can be either "Approved" or "Entry Withdrawn".
    #[derive(Debug, Serialize, Deserialize)]
    pub enum Status {
        #[serde(rename = "Approved")]
        Approve,
        #[serde(rename = "Entry Withdrawn")]
        Withdrawn,
    }

    /// Information from the locus-specific dabase.
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Lsdb {
        /// The name of the Locus Specific Mutation Database.
        pub name: String,
        /// The URL for the gene.
        pub url: String,
    }

    impl Display for Lsdb {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}|{}", &self.name, &self.url)
        }
    }

    impl FromStr for Lsdb {
        type Err = anyhow::Error;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let mut vals: Vec<String> = s.splitn(2, '|').map(|s| s.to_string()).collect();
            if vals.len() != 2 {
                anyhow::bail!("invalid LSDB string: {}", s);
            } else {
                let name = vals.pop().unwrap();
                let url = vals.pop().unwrap();
                Ok(Lsdb { name, url })
            }
        }
    }

    /// A record from the HGNC database.
    ///
    /// Also see the [HGNC website](https://www.genenames.org/download/archive/).
    #[serde_with::skip_serializing_none]
    #[serde_with::serde_as]
    #[derive(Debug, Serialize, Deserialize)]
    pub struct Record {
        /// HGNC ID. A unique ID created by the HGNC for every approved symbol.
        pub hgnc_id: String,
        /// The HGNC approved gene symbol.
        pub symbol: String,
        /// HGNC approved name for the gene.
        pub name: String,
        /// A group name for a set of related locus types as defined by the HGNC (e.g. non-coding
        /// RNA).
        pub locus_group: Option<String>,
        /// The locus type as defined by the HGNC (e.g. RNA, transfer).
        pub locus_type: Option<String>,
        /// Status of the symbol report.
        pub status: Status,
        /// Cytogenetic location of the gene (e.g. 2q34).
        pub location: Option<String>,
        /// Sortable cytogenic location of the gene (e.g. 02q34).
        pub location_sortable: Option<String>,
        /// Other symbols used to refer to this gene.
        pub alias_symbol: Option<Vec<String>>,
        /// Other names used to refer to this gene.
        pub alias_name: Option<Vec<String>>,
        /// Prevous symbols used to refer to this gene.
        pub prev_symbol: Option<Vec<String>>,
        /// Previous names used to refer to this gene.
        pub prev_name: Option<Vec<String>>,
        /// Name given to a gene group.
        pub gene_group: Option<Vec<String>>,
        /// ID used to designate a gene group.
        pub gene_group_id: Option<Vec<u32>>,
        /// The date the entry was first approved.
        pub date_approved_reserved: Option<NaiveDate>,
        /// The date the gene symbol was last changed.
        pub date_symbol_changed: Option<NaiveDate>,
        /// The date the gene name was last changed.
        pub date_name_changed: Option<NaiveDate>,
        /// Date the entry was last modified.
        pub date_modified: Option<NaiveDate>,
        /// Entrez gene id.
        pub entrez_id: Option<String>,
        /// Ensembl gene id.
        pub ensembl_gene_id: Option<String>,
        /// Vega gene id.
        pub vega_id: Option<String>,
        /// UCSC gene id.
        pub ucsc_id: Option<String>,
        /// ENA accession number(s).
        pub ena: Option<Vec<String>>,
        /// RefSeq nucleotide accession(s).
        pub refseq_accession: Option<Vec<String>>,
        /// Consensus CDS ID(ds).
        pub ccds_id: Option<Vec<String>>,
        /// Uniprot IDs.
        pub uniprot_ids: Option<Vec<String>>,
        // Pubmed IDs.
        pub pubmed_id: Option<Vec<u32>>,
        /// Mouse genome informatics database ID(s).
        pub mgd_id: Option<Vec<String>>,
        /// Rat genome database gene ID(s).
        pub rgd_id: Option<Vec<String>>,
        /// The name of the Locus Specific Mutation Database and URL for the gene.
        #[serde_as(as = "Option<Vec<DisplayFromStr>>")]
        pub lsdb: Option<Vec<Lsdb>>,
        /// Symbol used within COSMIC.
        pub cosmic: Option<String>,
        /// OMIM ID(s).
        pub omim_id: Option<Vec<String>>,
        /// miRBase ID.
        pub mirbase: Option<String>,
        /// Homeobox Database ID.
        pub homeodb: Option<u32>,
        /// snoRNABase ID.
        pub snornabase: Option<String>,
        /// Symbol used to link to the SLC tables database at bioparadigms.org for the gene.
        pub bioparadigms_slc: Option<String>,
        /// Orphanet ID.
        pub orphanet: Option<u32>,
        /// Pseudogene.org.
        #[serde(rename = "pseudogene.org")]
        pub pseudogene_org: Option<String>,
        /// Symbol used within HORDE for the gene.
        pub horde_id: Option<String>,
        /// ID used to link to the MEROPS peptidase database.
        pub merops: Option<String>,
        /// Symbol used within international ImMunoGeneTics information system.
        pub imgt: Option<String>,
        /// The objectId used to link to the IUPHAR/BPS Guide to PHARMACOLOGY database.
        pub iuphar: Option<String>,
        /// ID to link to the Mamit-tRNA database
        #[serde(rename = "mamit-trnadb")]
        pub mamit_trnadb: Option<String>,
        /// Symbol used within the Human Cell Differentiation Molecule database.
        pub cd: Option<String>,
        /// lncRNA Database ID.
        pub lncrnadb: Option<String>,
        /// ENZYME EC accession number.
        pub enzyme_id: Option<Vec<String>>,
        /// ID used to link to the Human Intermediate Filament Database.
        pub intermediate_filament_db: Option<String>,
        /// The HGNC ID that the Alliance of Genome Resources (AGR) use.
        pub agr: Option<String>,
        /// NCBI and Ensembl transcript IDs/acessions including the version number.
        pub mane_select: Option<Vec<String>>,
    }
}

#[cfg(test)]
mod tests {
    use super::{acmg_sf, hgnc};

    #[test]
    fn deserialize_acmg_sf_record() {}

    #[test]
    fn deserialize_ncbi_record() -> Result<(), anyhow::Error> {
        let path_tsv = "tests/db/genes/acmg_sf.example.tsv";
        let str_tsv = std::fs::read_to_string(path_tsv).unwrap();
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(str_tsv.as_bytes());
        let records = rdr
            .deserialize()
            .collect::<Result<Vec<acmg_sf::Record>, csv::Error>>()?;
        insta::assert_yaml_snapshot!(records);

        Ok(())
    }

    #[test]
    fn deserialize_hgnc_record() -> Result<(), anyhow::Error> {
        let path_json = "tests/db/genes/hgnc_info.example.jsonl";
        let str_json = std::fs::read_to_string(path_json)?;
        let records = str_json
            .lines()
            .map(|s| serde_json::from_str::<hgnc::Record>(s).unwrap())
            .collect::<Vec<_>>();

        insta::assert_yaml_snapshot!(records);

        Ok(())
    }
}
