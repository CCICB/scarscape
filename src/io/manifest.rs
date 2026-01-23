use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

use csv::StringRecord;

use crate::{
    error::ManifestError,
    model::{Manifest, SampleId, SampleInputs},
};

/// Parse a manifest file into a model::Manifest Struct
/// Manifest must be a csv with a header describing:
///
/// - sample: sample identifier
/// - snv: &PathBuf to vcf describing small variants
/// - cnv: (optional) &PathBuf to segment file describing genome-wide copy number changes  
/// - sv: (optional) &PathBuf to vcf describing structural variant breakends
pub fn parse_manifest(path: &PathBuf) -> Result<Manifest, ManifestError> {
    // Create a CSV reader (result)
    let reader_result = csv::ReaderBuilder::new()
        .flexible(false)
        .has_headers(true)
        .from_path(path);

    // Map the IO errors to our error type
    let mut reader = reader_result.map_err(|e| match e.kind() {
        csv::ErrorKind::Io(_) => ManifestError::ReadFailed {
            path: path.clone(),
            source: std::io::Error::other(e.to_string()),
        },
        _ => ManifestError::ParseFailed {
            path: path.clone(),
            source: e,
        },
    })?;

    let headers = reader
        .headers()
        .map_err(|e| ManifestError::ParseFailed {
            path: path.clone(),
            source: e,
        })?
        .clone();

    // Required columns
    let idx_sample_id =
        col_index(&headers, "sample").ok_or(ManifestError::MissingColumn { col: "sample" })?;
    let idx_snv_vcf =
        col_index(&headers, "snv").ok_or(ManifestError::MissingColumn { col: "snv" })?;

    // Optional columns
    let idx_sv = col_index(&headers, "sv");
    let idx_cnv = col_index(&headers, "cnv");

    // Setup HashSet to make sure same file is never read twice
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut samples: Vec<SampleInputs> = Vec::new();
    // Row numbers: +1 for header line; +1 because csv crate record indices start at 0.
    // We'll report human-friendly row numbers (1-based including header).
    for (i, result) in reader.records().enumerate() {
        let record = result.map_err(|e| ManifestError::ParseFailed {
            path: path.clone(),
            source: e,
        })?;
        let row_num = i + 2;

        let raw_id = get_required(&record, idx_sample_id, row_num, "sample")?;
        let sample = SampleId::new(raw_id).map_err(|reason| ManifestError::InvalidSampleId {
            row: row_num,
            reason,
        })?;

        if !seen_ids.insert(sample.as_str().to_owned()) {
            return Err(ManifestError::DuplicateSampleId {
                sample_id: sample.as_str().to_owned(),
            });
        }

        let snv = PathBuf::from(get_required(&record, idx_snv_vcf, row_num, "snv")?);

        let sv = idx_sv.and_then(|idx| get_optional_path(&record, idx));
        let cnv = idx_cnv.and_then(|idx| get_optional_path(&record, idx));

        // Filesystem checks: existence + is_file
        check_file(&sample, &snv)?;

        if let Some(ref p) = sv {
            check_file(&sample, p)?;
        }
        if let Some(ref p) = cnv {
            check_file(&sample, p)?;
        }

        samples.push(SampleInputs {
            sample,
            snv,
            sv,
            cnv,
        });
    }

    Ok(Manifest { samples })
}

/// Get column index of column with value matching col &str
fn col_index(headers: &StringRecord, col: &str) -> Option<usize> {
    // Pull index of first element where predicate == TRUE
    headers.iter().position(|s| s.trim() == col)
}

/// Check File Exists
fn check_file(sample: &SampleId, path: &Path) -> Result<(), ManifestError> {
    if !path.exists() {
        return Err(ManifestError::MissingFile {
            sample_id: sample.as_str().to_owned(),
            path: path.to_path_buf(),
        });
    }
    if !path.is_file() {
        return Err(ManifestError::NotAFile {
            sample_id: sample.as_str().to_owned(),
            path: path.to_path_buf(),
        });
    }
    Ok(())
}

fn get_required<'a>(
    record: &'a StringRecord,
    idx: usize,
    row: usize,
    col: &'static str, // Colunn is valid for the lifetime of StringRecord
) -> Result<&'a str, ManifestError> {
    let v = record.get(idx).map(str::trim).unwrap_or("");
    if v.is_empty() {
        return Err(ManifestError::MissingValue { row, col });
    }
    Ok(v)
}

fn get_optional_path(record: &StringRecord, idx: usize) -> Option<PathBuf> {
    let v = record.get(idx).map(str::trim).unwrap_or("");
    if v.is_empty() {
        None
    } else {
        Some(PathBuf::from(v))
    }
}
