//! Call-related information.

use crate::seqvars::query::schema::SequenceVariant;

/// Call-related record.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct Record {
    /// The genotype information for each sample.
    pub call_info: indexmap::IndexMap<String, CallInfo>,
}

impl Record {
    /// Construct a new `Record` from a `SequenceVariant`.
    ///
    /// # Error
    ///
    /// Returns an error if the `SequenceVariant` does not contain all necessary information.
    pub fn with_seqvar(seqvar: &SequenceVariant) -> Result<Self, anyhow::Error> {
        Ok(Self {
            call_info: seqvar
                .call_info
                .iter()
                .map(|(sample_name, call_info)| {
                    (
                        sample_name.clone(),
                        CallInfo::new(
                            call_info.dp,
                            call_info.ad,
                            call_info.quality.map(|q| q as i32),
                            call_info.genotype.as_ref().map(|s| {
                                if s.starts_with('|') || s.starts_with('/') {
                                    s[1..].to_string()
                                } else {
                                    s.to_string()
                                }
                            }),
                        ),
                    )
                })
                .collect(),
        })
    }
}

/// Genotype information.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_new::new)]
pub struct CallInfo {
    /// Depth of coverage.
    pub dp: Option<i32>,
    /// Alternate read depth.
    pub ad: Option<i32>,
    /// Genotype quality.
    pub gq: Option<i32>,
    /// Genotype.
    pub gt: Option<String>,
}
