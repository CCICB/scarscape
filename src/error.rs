//! Error types for scarscape.
//!
//! At the top of each rust script we can add `use crate::error::*``;
//! to gain access to our error enums AND a convenient Result alias that forces use of this
//! internal error type

use std::path::PathBuf;

use crate::model::SampleId;

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

    #[error("The {col} file for sample {sample} [record {record}] does NOT exist [{}]", path.display())]
    ManifestEntryFileNotFound {
        col: String,
        sample: SampleId,
        record: usize,
        path: PathBuf,
    },

    #[error("File not found: {}", path.display())]
    FileNotFound { path: PathBuf },

    #[error("Unexpected file type for path [{}]. Expected {expected}", path.display())]
    WrongFileType { path: PathBuf, expected: String },

    #[error("Sample Id is empty at record {record}")]
    SampleIdEmpty { record: usize },

    #[error("Failed to convert string slice [{0}] to PathBuf")]
    StringSliceToPathBuf(String),

    #[error("path is not valid UTF-8: {path}")]
    NonUtf8Path { path: PathBuf },

    #[error("unsupported file extension for {path}. Expected one of: {expected}")]
    UnsupportedExtension { path: PathBuf, expected: String },

    #[error("failed to open file for BGZF check: {path}")]
    OpenFailed { path: PathBuf },

    #[error("failed to read file header for BGZF check: {path}")]
    ReadFailed { path: PathBuf },

    // --- VCF I/O / parsing ---
    #[error("failed to read VCF header")]
    VcfHeaderRead,

    #[error("failed to read VCF record at index {record}")]
    VcfRecordRead {
        record: usize, // 1-based
    },

    // --- Record content errors ---
    #[error("VCF record {record} has unsupported ALT allele: {alt}")]
    UnsupportedVcfAlt { record: usize, alt: String },

    #[error("VCF record {record} has invalid {which} allele sequence: {allele}")]
    InvalidAlleleSequence {
        record: usize,
        which: &'static str, // "REF" or "ALT"
        allele: String,
    },
}

// A custom Result type that forces use of scarscape's internal Error type
pub type Result<T> = std::result::Result<T, crate::error::Error>;
