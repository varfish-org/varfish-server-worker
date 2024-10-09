//! Code for accessing HPO-related information.

use crate::pbs::varfish::v1::seqvars::output as pbs_output;

/// Enumeration for modes of inheritance.
#[derive(
    serde::Serialize,
    serde::Deserialize,
    enum_map::Enum,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Clone,
    Copy,
    Debug,
    strum::EnumString,
    strum::Display,
)]
pub enum ModeOfInheritance {
    /// Autosomal dominant inheritance (HP:0000006).
    AutosomalDominant,
    /// Autosomal recessive inheritance (HP:0000007).
    AutosomalRecessive,
    /// X-linked dominant inheritance (HP:0001419).
    XLinkedDominant,
    /// X-linked recessive inheritance (HP:0001423).
    XLinkedRecessive,
    /// Y-linked inheritance (HP:0001450).
    YLinked,
    /// Mitochondrial inheritance (HP:0001427).
    Mitochondrial,
}

impl ModeOfInheritance {
    /// Allow parsing of `ModeOfInheritance` from HPO ID.
    pub fn from_hpo_id(hpo_id: &str) -> Option<Self> {
        match hpo_id {
            "HP:0000006" | "HP:0012275" => Some(Self::AutosomalDominant),
            "HP:0000007" => Some(Self::AutosomalRecessive),
            "HP:0001419" => Some(Self::XLinkedDominant),
            "HP:0001423" => Some(Self::XLinkedRecessive),
            "HP:0001450" => Some(Self::YLinked),
            "HP:0001427" => Some(Self::Mitochondrial),
            _ => None,
        }
    }
}

impl From<ModeOfInheritance> for pbs_output::ModeOfInheritance {
    fn from(val: ModeOfInheritance) -> Self {
        match val {
            ModeOfInheritance::AutosomalDominant => {
                pbs_output::ModeOfInheritance::AutosomalDominant
            }
            ModeOfInheritance::AutosomalRecessive => {
                pbs_output::ModeOfInheritance::AutosomalRecessive
            }
            ModeOfInheritance::XLinkedDominant => pbs_output::ModeOfInheritance::XLinkedDominant,
            ModeOfInheritance::XLinkedRecessive => pbs_output::ModeOfInheritance::XLinkedRecessive,
            ModeOfInheritance::YLinked => pbs_output::ModeOfInheritance::YLinked,
            ModeOfInheritance::Mitochondrial => pbs_output::ModeOfInheritance::Mitochondrial,
        }
    }
}

pub type HgncToMoiMap = indexmap::IndexMap<String, indexmap::IndexSet<ModeOfInheritance>>;

/// Load the the `phenotype_to_genes.tsv` and `hgnc_xlink.tsv` files from the `hpo`
/// directory and build a map from HGNC gene ID to set of `ModeOfInheritance`
/// values.
pub fn load_hgnc_to_inheritance_map<P: AsRef<std::path::Path>>(
    path: &P,
) -> Result<HgncToMoiMap, anyhow::Error> {
    let phenotypes_to_genes =
        phenotype_to_genes::load_entries(&path.as_ref().join("phenotype_to_genes.txt"))
            .map_err(|e| anyhow::anyhow!("error loading phenotype_to_genes.txt: {}", e))?;
    let ncbi_to_hgnc = hgnc_xlink::load_ncbi_to_hgnc(path.as_ref().join("hgnc_xlink.tsv"))
        .map_err(|e| anyhow::anyhow!("error loading hgnc_xlink.tsv: {}", e))?;

    let mut result = indexmap::IndexMap::new();

    for ptg_entry in phenotypes_to_genes {
        if let Some(ncbi_gene_id) = ptg_entry.ncbi_gene_id {
            if let Some(hgnc_ids) = ncbi_to_hgnc.get(&ncbi_gene_id) {
                for hgnc_id in hgnc_ids {
                    if let Some(mode_of_inheritance) =
                        ModeOfInheritance::from_hpo_id(&ptg_entry.hpo_id)
                    {
                        let modes_of_inheritance = result
                            .entry(hgnc_id.clone())
                            .or_insert_with(indexmap::IndexSet::new);
                        modes_of_inheritance.insert_sorted(mode_of_inheritance);
                    }
                }
            }
        }
    }

    Ok(result)
}

/// Code for accessing the `phenotype_to_genes.tsv` file.
pub(super) mod phenotype_to_genes {
    /// Data structure for representing an entry of the table.
    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    #[serde_with::skip_serializing_none]
    pub struct Entry {
        /// Entrez gene ID.
        #[serde(alias = "entrez_id")]
        pub ncbi_gene_id: Option<u32>,
        /// Gene symbol.
        pub gene_symbol: String,
        /// HPO ID.
        pub hpo_id: String,
        // HPO Name.
        pub hpo_name: String,
    }

    /// Read the `phenotype_to_genes.tsv` file using the `csv` crate via serde.
    ///
    /// # Errors
    ///
    /// In the case that the file could not be read.
    pub fn load_entries<P: AsRef<std::path::Path>>(path: &P) -> Result<Vec<Entry>, anyhow::Error> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(path.as_ref())?;
        let mut entries = Vec::new();
        for result in rdr.deserialize() {
            let entry: Entry = result?;
            entries.push(entry);
        }
        Ok(entries)
    }
}

/// Code for accessing the `hgnc_xlink.tsv` file.
pub(super) mod hgnc_xlink {
    use std::collections::HashMap;

    /// Data structure for representing an entry of the table.
    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    #[serde_with::skip_serializing_none]
    pub struct Entry {
        /// HGNC gene ID.
        pub hgnc_id: String,
        /// Ensembl gene ID.
        pub ensembl_gene_id: Option<String>,
        /// Entrez gene ID.
        #[serde(alias = "entrez_id")]
        pub ncbi_gene_id: Option<u32>,
        /// Gene symbol.
        pub gene_symbol: String,
    }

    /// Read the `hgnc_xlink.tsv` file using the `csv` crate via serde.
    ///
    /// # Errors
    ///
    /// In the case that the file could not be read.
    pub fn load_entries<P: AsRef<std::path::Path>>(path: &P) -> Result<Vec<Entry>, anyhow::Error> {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(path.as_ref())?;
        let mut entries = Vec::new();
        for result in rdr.deserialize() {
            let entry: Entry = result?;
            entries.push(entry);
        }
        Ok(entries)
    }

    /// Read the `hgnc_xlink.tsv` into a map from NCBI gene ID to HGNC gene ID.
    ///
    /// # Errors
    ///
    /// In the case that the file could not be read.
    pub fn load_ncbi_to_hgnc<P: AsRef<std::path::Path>>(
        path: P,
    ) -> Result<HashMap<u32, Vec<String>>, anyhow::Error> {
        let mut map = HashMap::new();
        for entry in load_entries(&path)? {
            if let Some(ncbi_gene_id) = entry.ncbi_gene_id {
                let hgnc_id = entry.hgnc_id;
                let entry = map.entry(ncbi_gene_id).or_insert_with(Vec::new);
                entry.push(hgnc_id);
            }
        }
        Ok(map)
    }
}

#[cfg(test)]
mod test {
    #[test]
    pub fn test_mode_of_inheritance_from_hpo_id() {
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0000006"),
            Some(super::ModeOfInheritance::AutosomalDominant)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0012275"),
            Some(super::ModeOfInheritance::AutosomalDominant)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0000007"),
            Some(super::ModeOfInheritance::AutosomalRecessive)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0001419"),
            Some(super::ModeOfInheritance::XLinkedDominant)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0001423"),
            Some(super::ModeOfInheritance::XLinkedRecessive)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0001450"),
            Some(super::ModeOfInheritance::YLinked)
        );
        assert_eq!(
            super::ModeOfInheritance::from_hpo_id("HP:0001427"),
            Some(super::ModeOfInheritance::Mitochondrial)
        );

        assert_eq!(super::ModeOfInheritance::from_hpo_id("HP:0000000"), None);
    }

    #[test]
    pub fn test_mode_of_inheritance_into_pbs() {
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(
                super::ModeOfInheritance::AutosomalDominant
            ),
            super::pbs_output::ModeOfInheritance::AutosomalDominant
        );
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(
                super::ModeOfInheritance::AutosomalRecessive
            ),
            super::pbs_output::ModeOfInheritance::AutosomalRecessive
        );
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(
                super::ModeOfInheritance::XLinkedDominant
            ),
            super::pbs_output::ModeOfInheritance::XLinkedDominant
        );
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(
                super::ModeOfInheritance::XLinkedRecessive
            ),
            super::pbs_output::ModeOfInheritance::XLinkedRecessive
        );
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(super::ModeOfInheritance::YLinked),
            super::pbs_output::ModeOfInheritance::YLinked
        );
        assert_eq!(
            Into::<super::pbs_output::ModeOfInheritance>::into(
                super::ModeOfInheritance::Mitochondrial
            ),
            super::pbs_output::ModeOfInheritance::Mitochondrial
        );
    }

    #[test]
    fn load_hgnc_to_inheritance_map() -> Result<(), anyhow::Error> {
        let path = std::path::Path::new("tests/seqvars/query/db/hpo");
        let map = super::load_hgnc_to_inheritance_map(&path)?;

        assert_eq!(map.len(), 4599);
        insta::assert_yaml_snapshot!(&map[0..5]);

        Ok(())
    }

    #[test]
    fn test_phenotype_to_genes_load_entries() -> Result<(), anyhow::Error> {
        let path = std::path::Path::new("tests/seqvars/query/db/hpo/phenotype_to_genes.txt");
        let entries = super::phenotype_to_genes::load_entries(&path)?;

        assert_eq!(entries.len(), 988720);
        insta::assert_yaml_snapshot!(&entries[0..5]);

        Ok(())
    }

    #[test]
    fn test_hgnc_to_xlink_load_entries() -> Result<(), anyhow::Error> {
        let path = std::path::Path::new("tests/seqvars/query/db/hpo/hgnc_xlink.tsv");
        let map = super::hgnc_xlink::load_entries(&path)?;

        assert_eq!(map.len(), 43839);
        insta::assert_yaml_snapshot!(&map[0..5]);

        Ok(())
    }

    #[test]
    fn test_load_ncbi_to_hgnc() -> Result<(), anyhow::Error> {
        let path = std::path::Path::new("tests/seqvars/query/db/hpo/hgnc_xlink.tsv");
        let map = super::hgnc_xlink::load_ncbi_to_hgnc(&path)?;

        assert_eq!(map.len(), 43792);
        assert_eq!(map.get(&1), Some(&vec!["HGNC:5".to_string()]));
        assert_eq!(map.get(&2), Some(&vec!["HGNC:7".to_string()]));

        Ok(())
    }
}
