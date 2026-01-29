use crate::error::{Error, Result};
use noodles_vcf as vcf;
use seqlib::base::DnaBase;
use seqlib::mutations::DnaSmallMutation;
use seqlib::sequences::Seq;
use std::{
    collections::VecDeque,
    io::{self, BufRead},
    path::PathBuf,
};
pub struct SmallVariantReader<R> {
    rdr: vcf::io::reader::Reader<R>,
    header: vcf::Header,
    vcf_path: Option<PathBuf>,
    record_index: usize, // 1-based VCF records (excluding header)
    done: bool,

    // When a record has multiple ALT alleles we expand them into separate mutations.
    buffer: VecDeque<Result<DnaSmallMutation>>,
}

impl<R: io::Read> SmallVariantReader<R> {
    pub fn new(mut rdr: vcf::io::reader::Reader<R>) -> Result<Self> {
        let header = rdr
            .read_header()
            .map_err(|e| Error::VcfParse { source: e })?;
        Ok(Self {
            rdr,
            header,
            vcf_path: None,
            record_index: 0,
            done: false,
            buffer: VecDeque::new(),
        })
    }
}

impl<R: io::Read> Iterator for SmallVariantReader<R> {
    type Item = Result<DnaSmallMutation>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.buffer.pop_front() {
            return Some(item);
        }

        if self.done {
            return None;
        }

        let mut record = vcf::Record::default();
        match self.rdr.read_record(&self.header, &mut record) {
            Ok(0) => {
                self.done = true;
                None
            }
            Ok(_) => {
                self.record_index += 1;

                // Expand this record into per-ALT mutations.
                // Any parsing errors become buffered Err(...) items.
                expand_record_to_mutations(&record, self.record_index, self.vcf_path.as_deref())
                    .into_iter()
                    .for_each(|x| self.buffer.push_back(x));

                self.buffer.pop_front()
            }
            Err(e) => {
                self.done = true;
                Some(Err(Error::VcfParse { source: e }))
            }
        }
    }
}

use noodles_vcf::variant::record::AlternateBases;

fn expand_record_to_mutations(
    record: &noodles_vcf::Record,
    record_index: usize,
    vcf_path: Option<&std::path::Path>,
) -> Vec<Result<DnaSmallMutation>> {
    let chrom = record.reference_sequence_name().to_string();
    let pos = record.variant_start().map(|p| p.get() as u64).unwrap_or(0); // you likely want a real error if missing

    let ref_bases = record.reference_bases().to_string();

    // Collect ALT alleles
    let alts: Vec<String> = record
        .alternate_bases()
        .iter()
        .map(|a| a.to_string())
        .collect();

    let multiallelic = alts.len() > 1;

    // If ALT is "." => no variant; skip
    // (noodles represents this as AlternateBases::Missing / "." depending on version;
    // easiest is to string-check and skip.)
    let mut out = Vec::new();

    for alt in alts {
        if alt == "." {
            continue;
        }

        // Reject symbolic / breakend in the "small mutation" reader
        if alt.starts_with('<') || alt.contains('[') || alt.contains(']') {
            out.push(Err(Error::UnsupportedVcfAlt {
                record: record_index,
                alt,
            }));
            continue;
        }

        let reference = match Seq::<DnaBase>::new(&ref_bases) {
            Ok(s) => s,
            Err(e) => {
                out.push(Err(Error::InvalidAlleleSequence {
                    record: record_index,
                    which: "REF",
                    allele: ref_bases.clone(),
                    source: e,
                }));
                continue;
            }
        };

        let alternative = match Seq::<DnaBase>::new(&alt) {
            Ok(s) => s,
            Err(e) => {
                out.push(Err(Error::InvalidAlleleSequence {
                    record: record_index,
                    which: "ALT",
                    allele: alt.clone(),
                    source: e,
                }));
                continue;
            }
        };

        out.push(Ok(DnaSmallMutation::new(
            chrom.clone(),
            pos,
            reference,
            alternative,
            multiallelic,
            None, // context computed later
        )));
    }

    out
}
// /// Implement Iterator for ManifestReader so we can iterate through records
// /// This requires implementation of the functions: `next`
// impl<R: BufRead> Iterator for SmallVariantReader<R> {
//     type Item = Result<seqlib::mutations::DnaSmallMutation>;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         // Here Item is a Result<SampleInput, Error> result
//         // as described above
//         if self.done {
//             return None;
//         }
//
//         let mut rec = StringRecord::new();
//         match self.rdr.read_record(&mut rec) {
//             // if there's no more records to read, return None and stop iterating
//             Ok(false) => {
//                 self.done = true;
//                 None
//             }
//             // If there are more records to read
//             Ok(true) => {
//                 self.record_index += 1;
//
//                 // Parse the string record into SampleInputs object
//                 let result_sample_input = parse_record(
//                     &rec,
//                     &self.columns,
//                     self.record_index,
//                     self.base_dir.as_deref(),
//                 )
//                 .and_then(|m| m.validate());
//
//                 Some(result_sample_input)
//             }
//             Err(_) => {
//                 self.done = true;
//                 Some(Err(Error::ManifestParse))
//             }
//         }
//     }
// }
