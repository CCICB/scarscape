//! Central error types for scarscape.

use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ManifestError {
    #[error("failed to read manifest: {path}")]
    ReadFailed {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("failed to parse manifest TSV: {path}")]
    ParseFailed { path: PathBuf, source: csv::Error },

    #[error("manifest is missing required column: {col}")]
    MissingColumn { col: &'static str },

    #[error("manifest contains duplicate sample_id: {sample_id}")]
    DuplicateSampleId { sample_id: String },

    #[error("invalid sample_id at row {row}: {reason}")]
    InvalidSampleId { row: usize, reason: &'static str },

    #[error("missing required value at row {row}, column '{col}'")]
    MissingValue { row: usize, col: &'static str },

    #[error("input file does not exist for sample '{sample_id}': {path}")]
    MissingFile { sample_id: String, path: PathBuf },

    #[error("input path is not a file for sample '{sample_id}': {path}")]
    NotAFile { sample_id: String, path: PathBuf },
}

#[derive(Debug, Error)]
pub enum ScarscapeError {
    #[error(transparent)]
    Manifest(#[from] ManifestError),
    // Later:
    // #[error(transparent)]
    // Vcf(#[from] VcfError),
    // #[error(transparent)]
    // Cnv(#[from] CnvError),
}
