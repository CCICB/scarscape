//! I/O adapters
//!
//! This module owns all parsing, validation, and normalization boundaries
//! for external data sources (VCFs, FASTA, manifests, reference bundles).
//!
//! ## Responsibilities
//! - Parse external file formats
//! - Validate structural correctness
//! - Enforce strictness rules where required
//! - Emit normalized event streams
//!
//! ## Non-responsibilities
//! - Statistical computation
//! - Output formatting
//! - Business logic

pub mod cnv;
pub mod fasta;
pub mod manifest;
pub mod references;
pub mod vcf;
