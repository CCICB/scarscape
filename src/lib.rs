//! scarscape
//!
//! `scarscape` is a self-contained, streaming CLI for computing summary
//! statistics from cancer genome data.
//!
//! ## Architecture overview
//!
//! The codebase is organized around four concerns:
//!
//! 1. I/O adapters (`io::*`) — parsing, validation, normalization boundaries
//! 2. Domain model (`model`) — shared types and invariants
//! 3. Stats engines (`stats`) — pure computations over normalized streams
//! 4. Reporting (`report`) — stable schemas and serialization
//!
//! Each layer has a narrow responsibility and communicates via explicit
//! handoff types.

pub mod cli;
pub mod pipeline;

pub mod model;

pub mod io;

pub mod report;
pub mod stats;

pub mod error;
