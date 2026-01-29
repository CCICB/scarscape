//! VCF small-variant adapter (VCF -> DnaSmallMutation)
//!
//! This module provides a **streaming adapter** that reads a VCF and yields
//! normalized small-variant domain events (`DnaSmallMutation`).
//!
//! Design goals:
//! - filesystem-free construction via [`SmallVariantReader::from_reader`] for unit tests
//! - convenience construction via [`SmallVariantReader::from_path`] for CLI/pipelines
//! - streaming iteration (`Iterator<Item = Result<DnaSmallMutation>>`) with minimal buffering
//! - multi-allelic records are expanded into one mutation per ALT allele
//! - symbolic / breakend ALT alleles are rejected here (SV belongs in its own adapter)
//!
//! # Output contract
//! Each yielded [`DnaSmallMutation`] is constructed directly from VCF record fields:
//! - `chromosome` from `CHROM`
//! - `position` is **1-based** from `POS`
//! - `reference` from `REF`
//! - `alternative` from a single `ALT` allele
//! - `multiallelic` is `true` when the source record had >1 ALT allele
//! - `pass` is derived from the VCF `FILTER` field using noodles' `Filters::is_pass`
//! - `context` is always `None` here (computed later, e.g. via reference lookup)
//!
//! # Error semantics
//! This adapter is intentionally conservative: if a record contains unsupported ALT
//! encodings (symbolic/breakend) or invalid DNA alleles, an error is yielded.
//!
//! For multi-allelic records, one ALT allele may yield `Ok(...)` while another yields
//! `Err(...)`; the iterator preserves streamability by buffering per-record expansions.
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use noodles_vcf as vcf;
use noodles_vcf::variant::record::{AlternateBases, Filters};

use crate::error::{Error, Result};
use seqlib::mutations::DnaSmallMutation;
use seqlib::sequences::DnaSeq;

/// A streaming VCF reader that yields normalized [`DnaSmallMutation`] events.
///
/// This is an **adapter**: it performs lightweight parsing + normalization and emits
/// domain events suitable for downstream classification/tallying.
///
/// ## Streaming behavior
/// - Internally, one VCF record may expand into multiple mutations (one per ALT).
/// - These per-record results are stored in a small buffer and drained by [`Iterator::next`].
/// - The iterator yields `Result<DnaSmallMutation>` to allow parse/content errors to be
///   surfaced without panicking.
///
/// ## File provenance
/// `vcf_path` is optional and is only intended for debugging / improved error messages.
/// It is not required for correct parsing.
pub struct SmallVariantReader<R> {
    rdr: vcf::io::reader::Reader<R>,
    header: vcf::Header,
    vcf_path: Option<PathBuf>,
    record_index: usize, // 1-based VCF records (excluding header)
    done: bool,

    // When a record has multiple ALT alleles we expand them into separate mutations.
    buffer: VecDeque<Result<DnaSmallMutation>>,
}

impl SmallVariantReader<BufReader<File>> {
    /// Construct a [`SmallVariantReader`] from a file path.
    ///
    /// This is a convenience wrapper that:
    /// - checks file existence
    /// - opens the file
    /// - wraps it in a [`BufReader`]
    /// - delegates to [`SmallVariantReader::from_reader`]
    ///
    /// Prefer this in CLI/pipeline code. Prefer [`SmallVariantReader::from_reader`]
    /// for unit tests.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();

        // If you want to enforce existence *here*, do it here (not in from_reader).
        if !path.exists() {
            return Err(Error::FileNotFound { path });
        }

        let f = File::open(&path).map_err(|_| Error::OpenFailed { path: path.clone() })?;
        let reader = BufReader::new(f);
        Self::from_reader(reader, Some(path))
    }
}

impl<R: BufRead> SmallVariantReader<R> {
    /// Construct a [`SmallVariantReader`] from any [`BufRead`] input.
    ///
    /// This is the primary constructor and the preferred entrypoint for unit tests
    /// because it avoids filesystem concerns.
    ///
    /// ## Parameters
    /// - `reader`: any buffered input implementing [`BufRead`] (e.g. `BufReader<File>`,
    ///   `Cursor<Vec<u8>>`, stdin wrapper, etc.)
    /// - `vcf_path`: optional path used for provenance/debugging (may be `None`)
    ///
    /// ## Errors
    /// Returns [`Error::VcfHeaderRead`] if the VCF header cannot be read.
    pub fn from_reader(reader: R, vcf_path: Option<PathBuf>) -> Result<Self> {
        let mut rdr = vcf::io::reader::Reader::new(reader);

        let header = rdr.read_header().map_err(|_| Error::VcfHeaderRead)?;

        Ok(Self {
            rdr,
            header,
            vcf_path,
            record_index: 0,
            done: false,
            buffer: VecDeque::new(),
        })
    }

    fn refill_buffer_from_next_record(&mut self) -> Result<()> {
        let mut record = vcf::Record::default();

        match self.rdr.read_record(&mut record) {
            Ok(0) => {
                self.done = true;
                Ok(())
            }
            Ok(_) => {
                self.record_index += 1;

                let expanded = expand_record_to_mutations(&self.header, &record, self.record_index);

                self.buffer.extend(expanded);

                Ok(())
            }
            Err(_) => Err(Error::VcfRecordRead {
                record: self.record_index + 1, // next record (1-based) we attempted to read
            }),
        }
    }
}

impl<R: BufRead> Iterator for SmallVariantReader<R> {
    type Item = Result<DnaSmallMutation>;

    /// Yield the next parsed mutation (or error) from the underlying VCF stream.
    ///
    /// This iterator:
    /// - is forward-only and streaming
    /// - expands multi-allelic records into separate yielded items
    /// - stops after the first unrecoverable I/O record read error
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(item) = self.buffer.pop_front() {
                return Some(item);
            }

            if self.done {
                return None;
            }

            if let Err(e) = self.refill_buffer_from_next_record() {
                self.done = true;
                return Some(Err(e));
            }
        }
    }
}

/// Expand one VCF record into 0..N [`DnaSmallMutation`] (one per ALT).
///
/// ## `pass` handling
/// The `pass` field is derived using noodles' [`Filters::is_pass`] which interprets
/// the VCF `FILTER` field relative to the header.
///
/// ## Notes
/// - ALT="." is skipped
/// - symbolic ALTs (`<DEL>`) and breakends (`[...` / `...]`) are rejected
/// - allele strings are validated as DNA via [`DnaSeq::new`]
fn expand_record_to_mutations(
    header: &vcf::Header,
    record: &vcf::Record,
    record_index: usize,
) -> Vec<Result<DnaSmallMutation>> {
    let chrom = record.reference_sequence_name().to_string();

    // POS: Option<Result<Position, _>>
    let pos_1based: u64 = match record.variant_start() {
        Some(Ok(p)) => p.get() as u64,
        Some(Err(_)) | None => {
            return vec![Err(Error::VcfInvalidPos {
                record: record_index,
            })];
        }
    };

    // Compute pass once per record from FILTER
    //
    // Convention:
    // - FILTER="PASS" => pass=true
    // - FILTER="."    => pass=true (no filters applied / not filtered)
    // - otherwise     => pass=false
    //
    // If FILTER parsing itself errors, we surface a content error.
    let pass: bool = match record.filters().is_pass(header) {
        Ok(b) => b,
        Err(_) => {
            return vec![Err(Error::VcfInvalidFilter {
                record: record_index,
            })];
        }
    };

    let ref_bases = record.reference_bases().to_string();
    let ref_seq = match DnaSeq::new(&ref_bases) {
        Ok(s) => s,
        Err(_) => {
            return vec![Err(Error::InvalidAlleleSequence {
                record: record_index,
                which: "REF",
                allele: ref_bases,
            })];
        }
    };

    // ALT: iterator yields Result<&str, _>
    let mut alts: Vec<String> = Vec::new();
    for a in record.alternate_bases().iter() {
        match a {
            Ok(s) => alts.push(s.to_string()),
            Err(_) => {
                return vec![Err(Error::VcfInvalidAlt {
                    record: record_index,
                })];
            }
        }
    }

    let multiallelic = alts.len() > 1;
    let mut out = Vec::new();

    for alt in alts {
        if alt == "." {
            continue;
        }

        // Reject symbolic/breakend alts in the small-variant reader.
        if alt.starts_with('<') || alt.contains('[') || alt.contains(']') {
            out.push(Err(Error::UnsupportedVcfAlt {
                record: record_index,
                alt,
            }));
            continue;
        }

        let alt_seq = match DnaSeq::new(&alt) {
            Ok(s) => s,
            Err(_) => {
                out.push(Err(Error::InvalidAlleleSequence {
                    record: record_index,
                    which: "ALT",
                    allele: alt,
                }));
                continue;
            }
        };
        out.push(Ok(DnaSmallMutation::new(
            chrom.clone(),
            pos_1based,
            ref_seq.clone(),
            alt_seq,
            multiallelic,
            pass,
            None,
        )));
    }

    out
}
