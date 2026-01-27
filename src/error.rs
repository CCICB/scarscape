//! Error types for scarscape.
//!
//! At the top of each rust script we can add `use crate::error::*``;
//! to gain access to our error enums AND a convenient Result alias that forces use of this
//! internal error type

#[derive(Debug, thiserror::Error, Clone, PartialEq, Eq)]
pub enum Error {
    // Manifest Parsing Errors
    #[error("missing required column '{0}'")]
    ManifestMissingColumn(String),

    #[error("missing required value at record {record}, column '{col}'")]
    ManifestMissingValue { record: usize, col: String },

    #[error("invalid sample_id at record {record}: {reason}")]
    ManifestInvalidSampleId { record: usize, reason: &'static str },

    #[error("duplicate sample identifier '{sample_id}'")]
    ManifestDuplicateSampleId { sample_id: String },

    #[error("failed to parse manifest")]
    ManifestParse,

    #[error("Sample Id is empty at record {record}")]
    SampleIdEmpty { record: usize },

    #[error("Failed to convert string slice [{0}] to PathBuf")]
    StringSliceToPathBuf(String),
}

// A custom Result type that forces use of scarscape's internal Error type
pub type Result<T> = std::result::Result<T, crate::error::Error>;
