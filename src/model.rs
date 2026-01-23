//! Domain model
//!
//! This module defines the shared vocabulary used across scarscape:
//! sample identity, genome builds, coordinates, and normalized event types.
//!
//! Domain types encode invariants and semantics but do not perform I/O.

use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SampleId(String);

impl SampleId {
    /// Create a new SampleId with minimal validation.
    ///
    /// Keep this conservative and user-friendly: enforce non-empty and trim whitespace.
    /// Add stricter character rules later if you want.
    pub fn new(raw: &str) -> Result<Self, &'static str> {
        let s = raw.trim();
        if s.is_empty() {
            return Err("sample_id is empty");
        }
        Ok(Self(s.to_owned()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

#[derive(Debug, Clone)]
pub struct SampleInputs {
    pub sample: SampleId,
    pub snv: PathBuf,
    pub sv: Option<PathBuf>,
    pub cnv: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct Manifest {
    pub samples: Vec<SampleInputs>,
}
