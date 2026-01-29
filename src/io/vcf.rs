//! VCF small-variant adapter (VCF -> DnaSmallMutation)
//!
//! Aligned with DEV_PHILOSOPHY / ManifestReader:
//! - `from_reader` is filesystem-free and unit-testable
//! - `from_path` is convenience for CLI/pipeline
//! - streaming Iterator<Item = Result<DnaSmallMutation>>
//! - expands multi-allelic records into one mutation per ALT allele
//! - rejects symbolic/breakend ALTs here (SV should have its own adapter)

use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use noodles_vcf as vcf;
use noodles_vcf::variant::record::{AlternateBases, Filters};

use crate::error::{Error, Result};
use seqlib::mutations::DnaSmallMutation;
use seqlib::sequences::DnaSeq;

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
    /// Convenience for real files. Keeps filesystem concerns here.
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
    /// Primary constructor for library + tests (filesystem-free).
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

/// Expand one VCF record into 0..N DnaSmallMutation (one per ALT).
///
/// Note: we keep this intentionally conservative for the "small variant" adapter:
/// - ALT="." is skipped
/// - symbolic ALTs (<DEL>) and breakends ([...], ...]) are rejected
///
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
