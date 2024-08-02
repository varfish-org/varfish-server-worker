//! Code for sorting `VariantRecord` records by HGNC ID or coordinate.

use super::schema::data::VariantRecord;

/// Helper wrapper that allows to sort `VariantRecord` by HGNC ID.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct ByHgncId {
    pub hgnc_id: String,
    pub seqvar: VariantRecord,
}

impl From<VariantRecord> for ByHgncId {
    fn from(val: VariantRecord) -> Self {
        Self {
            hgnc_id: if !val.ann_fields.is_empty() {
                val.ann_fields[0].gene_id.clone()
            } else {
                String::new()
            },
            seqvar: val,
        }
    }
}

impl PartialEq for ByHgncId {
    fn eq(&self, other: &Self) -> bool {
        self.hgnc_id == other.hgnc_id
    }
}

impl Eq for ByHgncId {}

impl PartialOrd for ByHgncId {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ByHgncId {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hgnc_id.cmp(&other.hgnc_id)
    }
}

/// Helper wrapper that allows to sort `VariantRecord` by coordinate.
#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct ByCoordinate {
    pub coordinate: (String, i32),
    pub seqvar: VariantRecord,
}

impl From<VariantRecord> for ByCoordinate {
    fn from(val: VariantRecord) -> Self {
        Self {
            coordinate: (val.vcf_variant.chrom.clone(), val.vcf_variant.pos),
            seqvar: val,
        }
    }
}

impl PartialEq for ByCoordinate {
    fn eq(&self, other: &Self) -> bool {
        self.coordinate == other.coordinate
    }
}

impl Eq for ByCoordinate {}

impl PartialOrd for ByCoordinate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ByCoordinate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.coordinate.cmp(&other.coordinate)
    }
}
