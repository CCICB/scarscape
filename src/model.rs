//! Domain model
//!
//! This module defines the shared vocabulary used across scarscape:
//! sample identity, genome builds, coordinates, and normalized event types.
//!
//! Domain types encode invariants and semantics but do not perform I/O.

use crate::error::*;
use crate::io::utils::*;
use bedrs::Segment;
use core::fmt;
use noodles_vcf::record::samples::Sample;
use std::os::unix::fs::FileTypeExt;
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SampleId(String);

/// A stable identifier for a sample.
///
/// `SampleId` is used as a primary key across the pipeline and in output reports.
/// Construction enforces basic normalization so downstream components can assume:
/// - leading/trailing whitespace has been removed
/// - the identifier is not empty
///
/// # Invariants
/// - The inner string is non-empty after trimming.
/// - The inner string is the trimmed form of the input.
///
/// # Errors
/// Returns [`Error::SampleIdEmpty`] if the input is empty after trimming.
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

/// SegmentFile type
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct SegmentFile {
    path: PathBuf,
}

impl SegmentFile {
    /// Create a new segment file. Make sure it exists
    pub fn new(path: PathBuf) -> Result<Self> {
        if !path.exists() {
            return Err(Error::FileNotFound { path });
        }

        Ok(Self { path })
    }
}

impl std::fmt::Display for SegmentFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.path.display())
    }
}

pub struct VcfFile {
    path: PathBuf,
    filetype: VcfFileType,
}

impl std::fmt::Display for VcfFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Path: {}\nFiletype: {}",
            self.path.display(),
            self.filetype
        )
    }
}

impl VcfFile {
    pub fn new(path: PathBuf) -> Result<Self> {
        // Check File Exists
        if !path.exists() {
            return Err(Error::FileNotFound { path: path.clone() });
        }

        // Check Extension Makes Sense
        if !has_extension(&path, ["vcf.gz", "vcf", "bcf"])? {
            return Err(Error::WrongFileType {
                path: path.clone(),
                expected: "VCF/BCF (.vcf, .vcf.gz, .bcf)".to_string(),
            });
        }
        // Guess Compression Algorith
        let filetype = path_to_filetype(&path)?;

        Ok(Self { path, filetype })
    }
}

/// Parsed inputs for a single sample as specified by the manifest.
///
/// This is the *unvalidated* form produced by the manifest reader. It captures:
/// - the sample identifier
/// - required file paths (e.g., SNV VCF)
/// - optional file paths (e.g., SV/CNV)
///
/// Paths may be absolute or relative depending on how this entry was constructed.
/// When produced by [`ManifestReader::from_path`], relative paths are resolved
/// against the manifest file's directory.
///
/// # Design
/// Keep this type free of filesystem I/O. Validation is performed by
/// [`ManifestEntry::validate`] which produces a [`ValidatedManifestEntry`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ManifestEntry {
    pub sample: SampleId,
    pub snv_vcf: PathBuf,
    pub sv_vcf: Option<PathBuf>,
    pub cnv_segments: Option<PathBuf>,
    pub record_index: usize,
}

impl fmt::Display for ManifestEntry {
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
/// A manifest entry that has passed filesystem and type validation.
///
/// This type guarantees that required inputs exist on disk and that each
/// referenced file conforms to the expected kind (e.g., VCF vs FASTA).
///
/// Downstream pipeline stages should prefer this type so they can avoid
/// repeated precondition checks.
///
/// # Invariants
/// - All required files exist.
/// - Any optional file paths present also exist.
/// - Each file matches the expected kind checks configured by the validator.
///
/// # Notes
/// Validation is intentionally lightweight: it may check extensions, magic bytes,
/// and/or small header probes, but should not parse entire files.
pub struct ValidManifestEntry {
    pub sample: SampleId,
    pub snv_vcf: VcfFile,
    pub sv_vcf: Option<VcfFile>,
    pub cnv_segments: Option<SegmentFile>,
    pub record_index: usize,
}
impl fmt::Display for ValidManifestEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "-------------------")?;
        writeln!(f, "Sample ID: {}", self.sample)?;
        writeln!(f, "-------------------")?;
        writeln!(f, "SNV path: {}", self.snv_vcf)?;

        match &self.sv_vcf {
            Some(vcf) => writeln!(f, "SV File: {vcf}")?,
            None => writeln!(f, "SV file: not supplied")?,
        }

        match &self.cnv_segments {
            Some(segment) => writeln!(f, "CNV path: {segment}")?,
            None => writeln!(f, "CNV path: not supplied")?,
        }

        Ok(())
    }
}

impl ManifestEntry {
    pub fn validate(self) -> Result<ValidManifestEntry> {
        // Create VcfFiles Objects (this includes performing checks the file exists)
        let snv_vcf = VcfFile::new(self.snv_vcf)?;
        let sv_vcf = match self.sv_vcf {
            Some(path) => Some(VcfFile::new(path)?),
            None => None,
        };

        // Create SegmentFile Objects
        let cnv_segments = match self.cnv_segments {
            Some(path) => Some(SegmentFile::new(path)?),
            None => None,
        };

        //  SampleId does not need to be checked as its constructor runs validation

        Ok(ValidManifestEntry {
            sample: self.sample,
            snv_vcf,
            sv_vcf,
            cnv_segments,
            record_index: self.record_index,
        })
    }
}
