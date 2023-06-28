//! Code for supporting the database configuration (description) file.

use clap::ValueEnum;
use enum_map::{Enum, EnumMap};
use serde::{Deserialize, Serialize};
use strum_macros::EnumString;

/// Enum for the supported genome releases.
#[derive(
    Serialize, Deserialize, Enum, PartialEq, Eq, Clone, Copy, Debug, Default, EnumString, ValueEnum,
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

/// Enum for the supported gene/transcript databases.
#[derive(Serialize, Deserialize, Enum, PartialEq, Eq, Clone, Copy, Debug, Default, EnumString)]
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
#[derive(Serialize, Deserialize, Enum, PartialEq, Eq, Clone, Copy, Debug, Default, EnumString)]
#[serde(rename_all = "lowercase")]
#[strum(serialize_all = "lowercase")]
pub enum TadSet {
    /// hESC
    #[default]
    Hesc,
}

/// Enum for the supported gene ID cross-link tables.
#[derive(Serialize, Deserialize, Enum, PartialEq, Eq, Clone, Copy, Debug, Default, EnumString)]
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

/// Top level configuration file.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct Top {
    /// Whether the given genome release is enabled.
    pub release_enabled: EnumMap<GenomeRelease, bool>,
    /// Gene-centric annotation.
    pub genes: GeneDbs,
    /// Genome features annotation.
    pub features: EnumMap<GenomeRelease, FeatureDbs>,
    /// Configuration of structural variant databases for each genome build.
    pub vardbs: EnumMap<GenomeRelease, VarDbs>,
}

/// Configuration of masked region databases.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct MaskedRegionDbs {
    /// Relative path to the file with repeat masker database with checksum.
    pub repeat: DbDef,
    /// Relative path to the file with segmental duplication database with
    /// checksum.
    pub segdup: DbDef,
}

/// Configuration of genomic features.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct FeatureDbs {
    /// Maximal distance to feature to consider.
    #[serde(default = "default_max_dist")]
    pub max_dist: i32,

    /// Gene regions for the supported databases.
    pub gene_regions: EnumMap<Database, DbDef>,
    /// TAD definitions (by domain, not boundary).
    pub tads: EnumMap<TadSet, DbDef>,
    /// Definition of masked regions for repeats.
    pub masked: MaskedRegionDbs,
}

pub fn default_max_dist() -> i32 {
    10_000
}
/// Configuration of gene databases.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct GeneDbs {
    /// ACMG secondary findings list.
    pub acmg: DbDef,
    /// gnomAD gene constraints.
    pub gnomad_constraints: DbDef,
    /// Mapping from Entrez/NCBI gene ID to OMIM disease ID.
    pub mim2gene: DbDef,
    /// Gene ID interlink tables.
    pub xlink: EnumMap<GeneXlink, DbDef>,
}

/// Configuration of variant databases.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct VarDbs {
    /// Sequence variant databases
    pub seqvar: SeqVarDbs,
    /// Structural variant databases
    pub strucvar: StrucVarDbs,
}

/// Configuration of sequence variant databases.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct SeqVarDbs {
    // TODO
}

/// Configuration of structural variant databases.
#[derive(Serialize, Deserialize, PartialEq, Debug, Default)]
pub struct StrucVarDbs {
    /// Radius around BND sites used when building the database.
    #[serde(default = "default_slack_bnd")]
    pub slack_bnd: i32,
    /// Radius around INS sites used when building the database.
    #[serde(default = "default_slack_ins")]
    pub slack_ins: i32,
    /// Minimal reciprocal overlap for SVs of the same type, used when building
    /// the database.
    #[serde(default = "default_min_overlap")]
    pub min_overlap: f32,

    /// Relative path to the file with gnomAD-SV database with checksum.
    pub gnomad_sv: DbDef,
    /// Relative path to the file with dbVar SV database with checksum.
    pub dbvar: DbDef,
    /// Relative path to the file with DGV database with checksum.
    pub dgv: DbDef,
    /// Relative path to the file with DGV GS database with checksum.
    pub dgv_gs: DbDef,
    /// Relative path to the file with ExAC database with checksum.
    pub exac: Option<DbDef>,
    /// Relative path to the file with Thousan Genomes database with checksum.
    pub g1k: Option<DbDef>,

    /// Relative path to the file with the in-house variants with checksum.
    pub inhouse: Option<DbDef>,

    /// Relative path to the file with the known pathogenic variants.
    pub patho_mms: DbDef,
    /// Relative path to the ClinVar SV file.
    pub clinvar: DbDef,
}

pub fn default_slack_bnd() -> i32 {
    50
}

pub fn default_slack_ins() -> i32 {
    50
}

pub fn default_min_overlap() -> f32 {
    0.8
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use toml::to_string_pretty;

    use crate::db::conf::Top;

    #[test]
    fn run_smoke() -> Result<(), anyhow::Error> {
        let expected_str = std::fs::read_to_string("tests/db/conf/default.toml")?;

        let toml_str = to_string_pretty(&Top::default())?;

        assert_eq!(expected_str, toml_str);

        Ok(())
    }
}
