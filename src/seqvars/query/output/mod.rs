//! Data structureds for writing the output.

pub mod gene_related;
pub mod variant_related;

pub mod call_related;

/// A result record from the query.
///
/// These records are written to TSV for import into the database.   They contain the
/// bare necessary information for sorting etc. in the database.  The actual main
/// data is in the payload, serialized as JSON.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder)]
pub struct Record {
    /// UUID for the record.
    pub sodar_uuid: uuid::Uuid,
    /// Genome release for the coordinate.
    pub release: String,
    /// Chromosome name.
    pub chromosome: String,
    /// Chromosome number.
    pub chromosome_no: i32,
    /// Reference allele sequence.
    pub reference: String,
    /// Alternative allele sequence.
    pub alternative: String,
    /// UCSC bin of the record.
    pub bin: u32,
    /// Start position of the record.
    pub start: i32,
    /// End position of the record.
    pub end: i32,
    /// The result set ID as specified on the command line.
    pub smallvariantqueryresultset_id: String,
    /// The JSON-serialized `ResultPayload`.
    pub payload: String,
}

/// The structured result information of the result record.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize, derive_builder::Builder)]
pub struct Payload {
    /// Case UUID as specified on the command line.
    pub case_uuid: uuid::Uuid,
    /// The affected gene and consequence, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gene_related: Option<gene_related::Record>,
    /// Variant-related information, always present.
    #[serde(skip_serializing_if = "variant_related::Record::is_empty")]
    pub variant_related: variant_related::Record,
    /// Genotypes call related, always present.
    pub call_related: call_related::Record,
}
