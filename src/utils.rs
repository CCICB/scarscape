use anyhow::{bail, Context};
use csv::StringRecord;
use log::info;
use path_absolutize::Absolutize;
use serde::Deserialize;
use std::{collections::HashSet, fmt, hash::Hash, path::PathBuf};

pub fn log_error_and_panic(err: anyhow::Error) {
    let causes = err
        .chain()
        .skip(1)
        .map(|cause| format!("Caused by: {}", cause))
        .collect::<Vec<_>>()
        .join(" <= ");

    log::error!("Fatal Error: {}. {}", err, causes);
    std::process::exit(1)
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum Genome {
    /// the hg38/GRCh38 genome
    #[clap(name = "hg38")]
    GRCh38,
}

impl fmt::Display for Genome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Genome::GRCh38 => write!(f, "GRCh38 (hg38)"),
        }
    }
}

// struct GenomeFiles {
//     Genome: Genome,
//     // GenomeChromSizes: indexmap::map::IndexMap,
//     PathIGTCR: PathBuf,
//     PathGenome: PathBuf,
// }

#[derive(Debug, Deserialize)]
pub struct ManifestEntry {
    pub sample: String,
    pub sv: Option<PathBuf>,
    pub snv: Option<PathBuf>,
    pub cnv: Option<PathBuf>,
}

pub fn read_manifest(path: &PathBuf) -> Result<Vec<ManifestEntry>, anyhow::Error> {
    match path.exists() {
        true => log::info!("Reading manifest file: {:#?}", path),
        false => bail!("Could not find manifest file: {:#?}", path),
    }

    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .context("Failed to read manifest file")?;

    // Check Header Columns are as expected
    let headers = reader
        .headers()
        .context("Failed to read manifest header line")?;

    info!(
        "Found Manifest with {} columns: {}",
        headers.len(),
        headers
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    let missing_headers = missing_expected_headers(headers, &["sample", "sv", "snv", "cnv"]);
    if !missing_headers.is_empty() {
        bail!(
            "Manifest file is missing {} expected column: {}",
            missing_headers.len(),
            missing_headers.join(",")
        )
    }

    let mut n_entries = 0;
    let mut manifest: Vec<ManifestEntry> = vec![];
    for result in reader.deserialize() {
        let record: ManifestEntry = result.context("Failed to read manifest record")?;

        // Ensure sample name is more than 0 characters
        if record.sample.is_empty() {
            bail!(
                "Manifest contains empty sample identifier at line {} of [{:#?}]",
                n_entries + 2, // add 2: 1 for header, 1 because we iterate at end of loop
                path
            )
        }
        // Ensure record describes at least one datatype per sample
        if record.cnv.is_none() & record.snv.is_none() & record.sv.is_none() {
            bail!(
            "Manifest has no path to cnv, snv, or sv data for sample [{}]. Must contain at least one datatype",
            &record.sample.as_str()
          )
        }

        // Ensure all file paths point to real files (and they're valid file types) -- Continue this
        if let Some(p) = &record.sv {
            crate::vcf::check_vcf_exists_bgzipped_with_index(p)
                .context("Structural variant file in manifest")?
        }

        if let Some(p) = &record.snv {
            crate::vcf::check_vcf_exists_bgzipped_with_index(p)?
        }

        n_entries += 1;

        // Add record to manifest
        manifest.push(record);
    }

    let samples: Vec<String> = manifest.iter().map(|m| m.sample.clone()).collect();
    let dupsamples = find_duplicates(&samples);

    // Ensure there are no duplicate sample IDs (Throw error if there are)
    if !dupsamples.is_empty() {
        bail!(
            "Found {} duplicate sample IDs in manifest: [{}]",
            dupsamples.len(),
            dupsamples
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    }

    info!("Manifest describes {} unique samples", n_entries);
    Ok(manifest)
}

pub fn find_duplicates<T: Eq + Hash + fmt::Debug>(vec: &[T]) -> HashSet<&T> {
    let mut seen = HashSet::new();
    let mut duplicates = HashSet::new();
    for item in vec {
        if !seen.insert(item) {
            duplicates.insert(item);
        }
    }
    duplicates
}

fn missing_expected_headers(headers: &StringRecord, expected: &[&str]) -> Vec<String> {
    // Convert actual headers into a HashSet for efficient lookup.
    let headers_set: HashSet<&str> = headers.iter().collect();

    // Filter out expected headers that aren't in the actual headers.
    expected
        .iter()
        .filter(|&&header| !headers_set.contains(header))
        .map(|&header| header.to_string())
        .collect()
}

/// Converts the given output directory path to an absolute path and ensures the directory exists.
///
/// This function takes a reference to a [`PathBuf`], converts it into an absolute path using the
/// [`absolutize`] method, and creates the directory if it does not already exist.
///
/// # Arguments
///
/// * `path` - A reference to a [`PathBuf`] representing the directory path to convert and ensure exists.
///
/// # Returns
///
/// * `Ok(PathBuf)` containing the absolute path of the directory if successful.
/// * `Err(anyhow::Error)` if converting to an absolute path or creating the directory fails.
///
/// # Errors
///
/// This function will return an error in either of the following cases:
///
/// - The conversion of the given path to an absolute path fails.
/// - The directory creation fails (if the directory does not exist).
///
/// # Examples
///
/// ```rust
/// use std::path::PathBuf;
///
/// // Assume `make_output_directory` is imported from your module.
/// let path = PathBuf::from("relative/path");
/// let absolute_path = make_output_directory(&path)
///     .expect("Failed to create or retrieve the output directory");
/// println!("Absolute path: {:?}", absolute_path);
/// ```
///
/// [`absolutize`]: trait/YourTraitHere.html
pub fn make_output_directory(path: &PathBuf) -> Result<PathBuf, anyhow::Error> {
    let outdir: PathBuf = path
        .absolutize()
        .context("Failed to convert outdir to an absolute path")?
        .into();

    if !outdir.exists() {
        info!("Creating directory to store results {}", outdir.display());
        std::fs::create_dir(&outdir)?;
    }

    Ok(outdir)
}
