//! Code for supporting the database configuration (description) file.

use clap::ValueEnum;
use enum_map::Enum;
use hgvs::static_data::Assembly;
use serde::{Deserialize, Serialize};
use strum_macros::EnumString;

/// Enum for the supported genome releases.
#[derive(
    Serialize,
    Deserialize,
    Enum,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Debug,
    Default,
    EnumString,
    ValueEnum,
    strum::Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum GenomeRelease {
    /// GRCh37
    #[default]
    Grch37,
    /// GRCh38
    Grch38,
}

impl From<GenomeRelease> for Assembly {
    fn from(val: GenomeRelease) -> Self {
        match val {
            GenomeRelease::Grch37 => Assembly::Grch37p10,
            GenomeRelease::Grch38 => Assembly::Grch38,
        }
    }
}

/// Enum for the supported gene/transcript databases.
#[derive(
    Serialize,
    Deserialize,
    Enum,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Debug,
    Default,
    EnumString,
    strum::Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum Database {
    /// RefSeq
    #[default]
    RefSeq,
    /// ENSEMBL
    Ensembl,
}

/// Enum for the supported TADs.
#[derive(
    Serialize,
    Deserialize,
    Enum,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Debug,
    Default,
    EnumString,
    strum::Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum TadSet {
    /// hESC
    #[default]
    Hesc,
}

/// Enum for the supported gene ID cross-link tables.
#[derive(
    Serialize,
    Deserialize,
    Enum,
    PartialEq,
    Eq,
    Clone,
    Copy,
    Debug,
    Default,
    EnumString,
    strum::Display,
)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum GeneXlink {
    /// HGNC (complete for Entrez/NCBI gene IDs)
    #[default]
    Hgnc,
    /// ENSEMBL (complete for ENSEMBL gene IDs)
    Ensembl,
}

/// A relative path to a file with its checksum and specification.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct DbDef {
    /// The (relative) path to the (possibly compressed) text-based file.
    pub path: String,
    /// Optional MD5  checksum of file at `path`.
    pub md5: Option<String>,
    /// Optional SHA256 checksum of file at `path`.
    pub sha256: Option<String>,

    // The (optional, relative) path to the binary representation (for faster loading).
    pub bin_path: Option<String>,
    /// Optional MD5  checksum of file at `path_bin`.
    pub bin_md5: Option<String>,
    /// Optional SHA256 checksum of file at `path_bin`.
    pub bin_sha256: Option<String>,
}
