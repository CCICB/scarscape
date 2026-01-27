//! Domain model
//!
//! This module defines the shared vocabulary used across scarscape:
//! sample identity, genome builds, coordinates, and normalized event types.
//!
//! Domain types encode invariants and semantics but do not perform I/O.

use crate::error::*;
use core::fmt;
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SampleId(String);

impl SampleId {
    /// Create a new SampleId with minimal validation.
    ///
    /// Keep this conservative and user-friendly: enforce non-empty and trim whitespace.
    /// Add stricter character rules later if you want.
    pub fn new(raw: &str, record_index: usize) -> Result<Self> {
        let s = raw.trim();
        if s.is_empty() {
            return Err(Error::SampleIdEmpty {
                record: record_index,
            });
        }
        Ok(Self(s.to_owned()))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for SampleId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
/// The normalized, domain-level record produced by the Manifest parser.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleInputs {
    pub sample: SampleId,
    pub snv_vcf: PathBuf,
    pub sv_vcf: Option<PathBuf>,
    pub cnv_segments: Option<PathBuf>,
}

impl fmt::Display for SampleInputs {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "-------------------")?;
        writeln!(f, "Sample ID: {}", self.sample)?;
        writeln!(f, "-------------------")?;
        writeln!(f, "SNV path: {}", self.snv_vcf.display())?;

        match &self.sv_vcf {
            Some(path) => writeln!(f, "SV path: {}", path.display())?,
            None => writeln!(f, "SV path: not supplied")?,
        }

        match &self.cnv_segments {
            Some(path) => writeln!(f, "CNV path: {}", path.display())?,
            None => writeln!(f, "CNV path: not supplied")?,
        }

        Ok(())
    }
}
