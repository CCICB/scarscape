use crate::error::{Error, Result};
use crate::model::{SampleId, SampleInputs};
use csv::StringRecord;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RequiredCol {
    name: String,
    idx: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OptionalCol {
    name: String,
    idx: Option<usize>,
}

// impl ColInfo {
//     pub fn validate(&self) -> Result<()> {
//         if self.required && self.idx.is_none() {
//             return Err(Error::ManifestMissingColumn(self.name.clone()));
//         } else {
//             return Ok(());
//         }
//     }
// }

// For each column that could be described in manifest,
// Describe their name, idx, and whether they were found in the manifest
pub struct ColumnIndexes {
    sample: RequiredCol, // Column index describing sample identifier
    snv: RequiredCol,    // Column index describing snv VCF (DNA)
    sv: OptionalCol,     // Column index describing structural variant VCF (DNA)
    cnv: OptionalCol,    // Column index describing copy number segment file
}

/// Opinionated manifest reader.
///
/// This is a stateful object that remembers where it is in the file, what columns mean what, and
/// how to turn a row into one SampleInput. By implementing an Iterator for this reader we can
/// then iterate over every sample in our manifest
///
/// - TSV with headers
/// - `#` comment lines allowed
/// - required columns: `sample`, `snv`
/// - optional columns: `sv`, `cnv`
///
/// Output: iterator of `Result<SampleInput, ManifestError>`.
pub struct ManifestReader<R: BufRead> {
    rdr: csv::Reader<R>,
    columns: ColumnIndexes,
    record_index: usize, // 1-based data records (excluding header)
    manifest_path: Option<PathBuf>,
    base_dir: Option<PathBuf>,
    done: bool,
}

impl ManifestReader<BufReader<File>> {
    /// Production constructor: open file + set `base_dir` for resolving relative paths.
    pub fn from_path(path: impl AsRef<Path>) -> io::Result<Self> {
        // Grab file path information
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let base_dir = path.parent().map(|p| p.to_path_buf());

        // Delegate to common constructor
        Self::from_reader_with_context(BufReader::new(file), Some(path), base_dir)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
    }
}

impl<R: Read> ManifestReader<BufReader<R>> {
    /// Test-friendly constructor: accepts any reader.
    ///
    /// Note: relative paths are NOT resolved unless the adapter has a base_dir
    /// (base_dir is only set by `from_path`).
    pub fn from_reader(reader: R) -> Result<Self> {
        Self::from_reader_with_context(BufReader::new(reader), None, None)
    }
}

// Implement behaviour for ManifestReader
impl<R: BufRead> ManifestReader<R> {
    /// The core constructor for ManifestReader. Does lots of validation to make sure header is
    /// correct and fileformat looks sensible
    fn from_reader_with_context(
        reader: R,
        manifest_path: Option<PathBuf>,
        base_dir: Option<PathBuf>,
    ) -> Result<Self> {
        // Leverage csv::ReaderBuilder to build our generic CSV reader object for us (saves us a bunch of work)
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b',')
            .has_headers(true)
            .comment(Some(b'#'))
            .flexible(false)
            .from_reader(reader);

        // Grab headers from our reader
        let headers = rdr.headers().map_err(|_| Error::ManifestParse)?.clone();

        // Make sure it has the headings we expect (sample, snv, cnv, sv) ( get column index or
        // error if it doesn't exist)
        // else return error)
        let mut columns = ColumnIndexes {
            sample: RequiredCol {
                name: "sample".to_string(),
                idx: 0,
            },
            snv: RequiredCol {
                name: "sample".to_string(),
                idx: 0,
            },
            sv: OptionalCol {
                name: "sv".to_string(),
                idx: None,
            },
            cnv: OptionalCol {
                name: "cnv".to_string(),
                idx: None,
            },
        };

        // Set column index for each field
        // Required Columns:
        columns.sample.idx = col_index(&headers, &columns.sample.name)
            .ok_or_else(|| Error::ManifestMissingColumn(columns.sample.name.clone()))?;
        columns.snv.idx = col_index(&headers, &columns.snv.name)
            .ok_or_else(|| Error::ManifestMissingColumn(columns.snv.name.clone()))?;
        // Optional Columns
        columns.sv.idx = col_index(&headers, &columns.sv.name);
        columns.cnv.idx = col_index(&headers, &columns.cnv.name);

        // Create the new reader
        Ok(Self {
            rdr,
            columns,
            record_index: 0,
            manifest_path,
            base_dir,
            done: false,
        })
    }
}

/// Implement Iterator for ManifestReader so we can iterate through records
/// This requires implementation of the functions: `next`
impl<R: BufRead> Iterator for ManifestReader<R> {
    type Item = Result<SampleInputs>;

    fn next(&mut self) -> Option<Self::Item> {
        // Here Item is a Result<SampleInput, Error> result
        // as described above
        if self.done {
            return None;
        }

        let mut rec = StringRecord::new();
        match self.rdr.read_record(&mut rec) {
            // if there's no more records to read, return None and stop iterating
            Ok(false) => {
                self.done = true;
                None
            }
            // If there are more records to read
            Ok(true) => {
                self.record_index += 1;

                // Parse the string record into SampleInputs object
                let result_sample_input = parse_record(
                    &rec,
                    &self.columns,
                    self.record_index,
                    self.base_dir.as_deref(),
                );

                Some(result_sample_input)
            }
            Err(_) => {
                self.done = true;
                Some(Err(Error::ManifestParse))
            }
        }
    }
}

fn parse_record(
    record: &StringRecord,
    columns: &ColumnIndexes,
    record_index: usize,
    base_dir: Option<&Path>,
) -> Result<SampleInputs> {
    let sample_raw = get_required_field(record, &columns.sample, record_index)?;
    let snv_raw = get_required_field(record, &columns.snv, record_index)?;

    let sample_id = SampleId::new(sample_raw, record_index)?;

    let snv_vcf = resolve_path(base_dir, snv_raw);

    let sv_vcf = get_optional_field(record, &columns.sv).map(|s| resolve_path(base_dir, s));

    let cnv_segments = get_optional_field(record, &columns.cnv).map(|s| resolve_path(base_dir, s));

    Ok(SampleInputs {
        sample: sample_id,
        snv_vcf,
        sv_vcf,
        cnv_segments,
    })
}

/// Get a RequiredCol value from a string record
/// Returns a string borrowed from the record
fn get_required_field<'a>(
    record: &'a StringRecord,
    col: &RequiredCol,
    record_index: usize,
) -> Result<&'a str> {
    let v = record.get(col.idx).unwrap_or("").trim();
    if v.is_empty() {
        return Err(Error::ManifestMissingValue {
            record: record_index,
            col: col.name.clone(),
        });
    }
    Ok(v)
}

fn get_optional_field<'a>(record: &'a StringRecord, col: &OptionalCol) -> Option<&'a str> {
    let idx: usize = col.idx?;

    // Return Record
    record.get(idx)
}

/// Fetch the index of the column named 'col'
fn col_index(headers: &StringRecord, col: &str) -> Option<usize> {
    headers.iter().position(|h| h.trim() == col)
}

fn resolve_path(base_dir: Option<&Path>, raw: &str) -> PathBuf {
    let p = PathBuf::from(raw);
    if p.is_relative() {
        if let Some(base) = base_dir {
            return base.join(p);
        }
    }
    p
}
