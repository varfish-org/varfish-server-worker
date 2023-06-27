//! Module with code supporting the parsing.

use std::collections::HashMap;

use crate::sv::query::schema::{Pathogenicity, VariationType};
use serde::{de, Deserialize, Deserializer};
use tracing::warn;

lazy_static::lazy_static! {
    static ref VARIATION_TYPE_LABELS: HashMap<&'static str, VariationType> = {
        let mut m = HashMap::new();
        m.insert("Complex", VariationType::Complex);
        m.insert("copy number gain", VariationType::Dup);
        m.insert("copy number loss", VariationType::Del);
        m.insert("Deletion", VariationType::Del);
        m.insert("Duplication", VariationType::Dup);
        m.insert("fusion", VariationType::Bnd);
        m.insert("Indel", VariationType::Cnv);
        m.insert("Insertion", VariationType::Ins);
        m.insert("Inversion", VariationType::Inv);
        m.insert("Microsatellite", VariationType::Microsatellite);
        m.insert("Tandem duplication", VariationType::Dup);
        m.insert("Translocation", VariationType::Bnd);
        m
    };
}

impl VariationType {
    pub fn from_label(label: &str) -> Result<VariationType, anyhow::Error> {
        if let Some(result) = VARIATION_TYPE_LABELS.get(label) {
            Ok(*result)
        } else {
            Err(anyhow::anyhow!("Invalid VariationType label: {}", label))
        }
    }
}

lazy_static::lazy_static! {
    static ref PATHOGENICITY_LABELS: HashMap<&'static str, Pathogenicity> = {
        let mut m = HashMap::new();
        m.insert("{\"benign\"}", Pathogenicity::Benign);
        m.insert("{\"benign\",\"likely benign\"}", Pathogenicity::LikelyBenign);
        m.insert("{\"likely benign\"}", Pathogenicity::LikelyBenign);
        m.insert("{\"likely pathogenic\"}", Pathogenicity::LikelyPathogenic);
        m.insert("{\"likely pathogenic\",\"pathogenic\"}", Pathogenicity::LikelyPathogenic);
        m.insert("{\"pathogenic\"}", Pathogenicity::Pathogenic);
        m.insert("{\"uncertain significance\"}", Pathogenicity::Uncertain);

        m
    };
}

impl Pathogenicity {
    pub fn from_label(label: &str) -> Result<Self, anyhow::Error> {
        if let Some(pathogenicity) = PATHOGENICITY_LABELS.get(label) {
            Ok(*pathogenicity)
        } else {
            warn!("Cannot decode pathogenicity from {}", label);
            Ok(Pathogenicity::Uncertain)
        }
    }
}

/// Record as created by VarFish DB Downloader.
#[derive(Debug, Deserialize)]
pub struct Record {
    /// Chromosome name
    pub chromosome: String,
    /// 1-based start position
    pub begin: i32,
    /// 1-based end position
    pub end: i32,
    /// unused
    #[allow(dead_code)]
    bin: u32,
    /// unused
    #[allow(dead_code)]
    reference: String,
    /// unused
    #[allow(dead_code)]
    alternative: String,
    /// unused
    #[allow(dead_code)]
    clinvar_version: String,
    /// unused
    #[allow(dead_code)]
    set_type: String,
    /// ClinVar SV variation type
    #[serde(deserialize_with = "from_variation_type_label")]
    pub variation_type: VariationType,
    /// unused
    #[allow(dead_code)]
    symbols: String,
    /// unused
    #[allow(dead_code)]
    hgnc_ids: String,
    /// The ClinVar VCV identifier
    pub vcv: String,
    /// unused
    #[allow(dead_code)]
    summary_clinvar_review_status_label: String,
    /// unused
    #[allow(dead_code)]
    summary_clinvar_pathogenicity_label: String,
    /// Pathogenicity
    #[serde(
        alias = "summary_clinvar_pathogenicity",
        deserialize_with = "from_pathogenicity_summary"
    )]
    pub pathogenicity: Pathogenicity,
    /// unused
    #[allow(dead_code)]
    summary_clinvar_gold_stars: String,
    /// unused
    #[allow(dead_code)]
    summary_paranoid_review_status_label: String,
    /// unused
    #[allow(dead_code)]
    summary_paranoid_pathogenicity_label: String,
    /// unused
    #[allow(dead_code)]
    summary_paranoid_pathogenicity: String,
    /// unused
    #[allow(dead_code)]
    summary_paranoid_gold_stars: String,
    /// unused
    #[allow(dead_code)]
    details: String,
}

/// Deserialize "VariationType" from ClinVar TSV file
///
/// This function will strip everything after the first underscore.
fn from_variation_type_label<'de, D>(deserializer: D) -> Result<VariationType, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    VariationType::from_label(s).map_err(de::Error::custom)
}

/// Deserialize "Pathogenicity" from ClinVar TSV file
///
/// This function will strip everything after the first underscore.
fn from_pathogenicity_summary<'de, D>(deserializer: D) -> Result<Pathogenicity, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;
    Pathogenicity::from_label(s).map_err(de::Error::custom)
}
