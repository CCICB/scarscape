use anyhow::bail;
use anyhow::Context;
use bedrs::Bed4;
use bedrs::Coordinates;
use bedrs::IntervalContainer;
use log::info;
use log::warn;
use noodles::core::position;
use noodles::core::Region;
use noodles::csi::BinningIndex;
use noodles::tabix;
use noodles::vcf;
use noodles_vcf::header::record::value::map::alternative_allele;
use noodles_vcf::variant::record::AlternateBases;
use noodles_vcf::variant::record::Filters;
use regex::Regex;
use serde::Deserialize;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::io::Read;
use std::path::Path;
use std::{fs::File, path::PathBuf};

use crate::utils::Genome;

/// Struct representing IGTCR (Immunoglobulin & T-cell receptor) regions.
/// It contains vectors for different receptor regions.
#[derive(Debug)]
pub struct RegionsIGTCR {
    pub igh: Vec<Region>,
    pub igk: Vec<Region>,
    pub igl: Vec<Region>,
    pub tra: Vec<Region>,
    pub trb: Vec<Region>,
    pub trd: Vec<Region>,
    pub trg: Vec<Region>,
}

impl RegionsIGTCR {
    /// Returns a single vector containing all B-cell receptor (BCR) regions (IGH, IGK, IGL).
    fn _get_ig_regions(&self) -> Vec<Region> {
        [self.igh.clone(), self.igk.clone(), self.igl.clone()].concat()
    }

    /// Returns a single vector containing all T-cell receptor (TCR) regions (TRA, TRB, TRD, TRG).
    fn _get_tcr_regions(&self) -> Vec<Region> {
        [
            self.tra.clone(),
            self.trb.clone(),
            self.trd.clone(),
            self.trg.clone(),
        ]
        .concat()
    }

    /// Returns the total number of IGTCR regions across all receptor types.
    fn len(&self) -> usize {
        self.igh.len()
            + self.igk.len()
            + self.igl.len()
            + self.tra.len()
            + self.trb.len()
            + self.trd.len()
            + self.trg.len()
    }

    /// Returns the total count of B-cell receptor (BCR) regions.
    fn count_bcr_regions(&self) -> usize {
        self.igh.len() + self.igk.len() + self.igl.len()
    }

    /// Returns the total count of T-cell receptor (TCR) regions.
    fn count_tcr_regions(&self) -> usize {
        self.tra.len() + self.trb.len() + self.trd.len() + self.trg.len()
    }
}

impl fmt::Display for RegionsIGTCR {
    /// Formats the IGTCR regions into a human-readable summary string.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "IGTCR Regions [total: {} | BCR: {} | TCR {} | igh: {} | igk: {} | igl: {} | tra: {} | trb: {} | trd: {} | trg: {}]",
            self.len(),
            self.count_bcr_regions(),
            self.count_tcr_regions(),
            self.igh.len(),
            self.igk.len(),
            self.igl.len(),
            self.tra.len(),
            self.trb.len(),
            self.trd.len(),
            self.trg.len()
        )
    }
}

/// Alias for an IntervalContainer holding Bed4 intervals.
type IntervalContainerBed4 = IntervalContainer<bedrs::Bed4<String, u64, String>, String, u64>;

/// Struct summarizing structural variant counts.
#[derive(Debug)]
pub struct StructuralVariantSummary {
    pub pass: u64,
    pub fail: u64,
    pub igh: u64, // Count of PASS variants in IGH region.
    pub igk: u64, // Count of PASS variants in IGK region.
    pub igl: u64, // Count of PASS variants in IGL region.
    pub tra: u64, // Count of PASS variants in TRA region.
    pub trb: u64, // Count of PASS variants in TRB region.
    pub trd: u64, // Count of PASS variants in TRD region.
    pub trg: u64, // Count of PASS variants in TRG region.
    pub total: u64,
}

impl StructuralVariantSummary {
    /// Returns the total number of PASS variants in all IGTCR regions.
    fn igtcr(&self) -> u64 {
        self.igh + self.igk + self.igl + self.tra + self.trb + self.trg + self.trd
    }

    /// Formats the structural variant summary as a CSV string.
    ///
    /// # Arguments
    ///
    /// * `sample` - A string slice that holds the sample identifier.
    ///
    /// # Returns
    ///
    /// A `String` containing the CSV-formatted summary.
    pub fn to_csv(&self, sample: &str) -> String {
        format!(
        "sample,total,pass,fail,igtcr,igh,igk,igl,tra,trb,trd,trg\n{},{},{},{},{},{},{},{},{},{},{},{}",
        sample,
        self.total,
        self.pass,
        self.fail,
        self.igtcr(),
        self.igh,
        self.igk,
        self.igl,
        self.tra,
        self.trb,
        self.trd,
        self.trg,
    )
    }

    /// Writes the CSV string representation of the summary to the specified file.
    ///
    /// # Arguments
    ///
    /// * `sample` - A string slice that holds the sample identifier.
    /// * `path` - A reference to the `PathBuf` representing the destination file path.
    ///
    /// # Errors
    ///
    /// Returns an error if the file already exists or if writing fails.
    pub fn write_to_csv(&self, sample: &str, path: &PathBuf) -> Result<(), anyhow::Error> {
        if path.exists() {
            bail!(
                "File [{}] already exists. Can NOT overwrite. Please manually delete and retry.",
                path.display()
            )
        }
        let csv_data = self.to_csv(sample);
        std::fs::write(path, csv_data)?;
        Ok(())
    }
}

impl fmt::Display for StructuralVariantSummary {
    /// Formats the structural variant summary into a human-readable string.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "SV stats [total: {} (pass: {}, fail: {}) | IGTCR: {} (igh: {}, igk: {}, igl: {}, tra: {}, trb: {}, trd: {}, trg: {})]",
            self.total,
            self.pass,
            self.fail,
            self.igtcr(),
            self.igh,
            self.igk,
            self.igl,
            self.tra,
            self.trb,
            self.trd,
            self.trg,
        )
    }
}

/// Counts structural variants in a VCF file and summarizes the results.
///
/// This function reads a bgzipped VCF file (with index), counts the total, PASS, and FAIL variants,
/// and further counts PASS variants in specific IGTCR regions.
///
/// # Arguments
///
/// * `path` - A reference to the `PathBuf` representing the VCF file path.
/// * `genome` - A reference to the `Genome` enum representing the reference genome used.
///
/// # Returns
///
/// Returns a `StructuralVariantSummary` containing the counts on success, or an error on failure.
pub fn sv_vcf_counts(
    path: &PathBuf,
    genome: &Genome,
) -> Result<StructuralVariantSummary, anyhow::Error> {
    // Read VCF file
    info!("Reading VCF file: {:#?}", path);

    // Check file exists and is an indexed BCF/VCF file so we can query it
    check_vcf_exists_bgzipped_with_index(path)?;

    // Check Index Exists
    let _index_path = assert_index_file_exists(path)?;
    let chroms = list_vcf_chromosomes_from_vcf_header(path)?;
    let vcf_uses_chr_prefix = chroms_include_chr(&chroms);

    // Read BED file into regions object describing TCR and BCR gene regions
    let igtcr_regions = get_igtcr_regions_bedrs(genome, vcf_uses_chr_prefix)
        .context("Failed to prepare IGTCR regions")?;

    let mut vcf_reader = vcf::io::indexed_reader::Builder::default().build_from_path(path)?;
    info!("Successfully to read VCF [{:#?}]", path);

    let mut ntotal: u64 = 0;
    let mut npass: u64 = 0;
    let mut nfail: u64 = 0;

    let header = vcf_reader.read_header()?;

    // Count SVs (General stats, PASS / FAIL / Total)
    for result in vcf_reader.records() {
        let record = result?;
        let pass = record
            .filters()
            .is_pass(&header)
            .expect("Failed to get sv filter");

        if pass {
            npass += 1
        } else {
            nfail += 1
        }
        ntotal += 1
    }

    // Count PASS SVs in IGTCR regions using vcftools region query
    let mut counts = count_pass_variants_in_regions_bedrs(igtcr_regions, path)?;

    // Ensure there is a value for all samples of interest
    counts.entry(String::from("IGH")).or_insert(0);
    counts.entry(String::from("IGK")).or_insert(0);
    counts.entry(String::from("IGL")).or_insert(0);
    counts.entry(String::from("TRA")).or_insert(0);
    counts.entry(String::from("TRB")).or_insert(0);
    counts.entry(String::from("TRD")).or_insert(0);
    counts.entry(String::from("TRG")).or_insert(0);

    let sv_summary = StructuralVariantSummary {
        pass: npass,
        fail: nfail,
        igh: *counts.get("IGH").unwrap(),
        igk: *counts.get("IGK").unwrap(),
        igl: *counts.get("IGL").unwrap(),
        tra: *counts.get("TRA").unwrap(),
        trb: *counts.get("TRB").unwrap(),
        trd: *counts.get("TRD").unwrap(),
        trg: *counts.get("TRG").unwrap(),
        total: ntotal,
    };

    Ok(sv_summary)
}

/// Enum representing different file types for VCF/BCF files.
#[derive(Debug)]
enum FileType {
    Vcf,
    VcfCompressedGzip,
    VcfCompressedBgzip,
    Bcf,
    BcfCompressedBgzip,
    BcfCompressedGzip,
    Unsure,
}

/// Determines the file type of the given path based on its filename and compression markers.
///
/// # Arguments
///
/// * `path` - A reference to a `Path` representing the file.
///
/// # Returns
///
/// A `FileType` enum variant indicating the file type.
fn path_to_filetype(path: &Path) -> FileType {
    let path_str = path.to_str().unwrap_or("");

    let vcf_gz_re =
        Regex::new(r"\.vcf\.b?gz$").expect("Failed to create regex for vcf (compressed)");
    let vcf_re = Regex::new(r"\.vcf$").expect("Failed to create regex for vcf (uncompressed)");
    let bcf_gz_re =
        Regex::new(r"\.bcf\.b?gz$").expect("Failed to create regex for bcf (compressed)");
    let bcf_re = Regex::new(r"\.bcf$").expect("Failed to create regex for bcf (uncompressed)");

    if vcf_gz_re.is_match(path_str) {
        match is_bgzf(path_str) {
            true => FileType::VcfCompressedBgzip,
            false => FileType::VcfCompressedGzip,
        }
    } else if vcf_re.is_match(path_str) {
        FileType::Vcf
    } else if bcf_re.is_match(path_str) {
        FileType::Bcf
    } else if bcf_gz_re.is_match(path_str) {
        match is_bgzf(path_str) {
            true => FileType::BcfCompressedBgzip,
            false => FileType::BcfCompressedGzip,
        }
    } else {
        warn!("Filetype uncertain {}", path.display());
        FileType::Unsure
    }
}

/// Checks whether the file at the given path is compressed with BGZF by reading its header.
///
/// # Arguments
///
/// * `file_path` - A string slice representing the file path.
///
/// # Returns
///
/// `true` if the file is BGZF-compressed, `false` otherwise.
fn is_bgzf(file_path: &str) -> bool {
    let mut file =
        File::open(file_path).expect("Failed to open file when checking bgzf compression");
    // Read the first 18 bytes: 10 bytes of fixed gzip header, 2 bytes of XLEN, and 6 bytes of extra field.
    let mut header = [0u8; 18];
    file.read_exact(&mut header)
        .expect("Failed to read header when checking bgzf compression");

    // Check for the standard gzip magic numbers: 0x1F, 0x8B, 0x08.
    if header[0] != 0x1F || header[1] != 0x8B || header[2] != 0x08 {
        // Not even a valid gzip file.
        return false;
    }

    // Bytes 10-11 store the XLEN (the length of the extra field).
    let xlen = u16::from_le_bytes([header[10], header[11]]);
    // BGZF extra field should be at least 6 bytes long and start with "BC".
    if xlen >= 6 && header[12] == b'B' && header[13] == b'C' {
        return true;
    }

    false
}

/// Asserts that an index file exists for the given VCF/BCF file and returns its path.
///
/// # Arguments
///
/// * `path_to_vcf` - A reference to a `Path` representing the VCF/BCF file.
///
/// # Returns
///
/// Returns a `PathBuf` pointing to the index file if it exists, or an error if it does not.
fn assert_index_file_exists(path_to_vcf: &Path) -> Result<PathBuf, anyhow::Error> {
    let path_tbi = PathBuf::from(format!("{}.tbi", path_to_vcf.to_string_lossy()));
    if !path_tbi.exists() {
        anyhow::bail!(
            "Could not find index file: [{}]. Please index the vcf by running `tabix -p vcf {}`",
            path_tbi.display(),
            path_to_vcf.display()
        );
    }
    Ok(path_tbi)
}

/// Lists the unique chromosome names from the contigs in a VCF header.
///
/// # Arguments
///
/// * `path` - A reference to a `PathBuf` representing the VCF file.
///
/// # Returns
///
/// A `HashSet<String>` containing the chromosome names extracted from the VCF header.
fn list_vcf_chromosomes_from_vcf_header(path: &PathBuf) -> Result<HashSet<String>, anyhow::Error> {
    let mut reader = vcf::io::reader::Builder::default().build_from_path(path)?;

    // Read the header
    let header = reader.read_header()?;

    // Extract unique chromosome names from the contigs in the header
    let chroms: HashSet<String> = header
        .contigs()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    Ok(chroms)
}

/// Lists the unique chromosome names by reading the tabix index file for a VCF.
///
/// # Arguments
///
/// * `path` - A reference to a `Path` representing the VCF file.
///
/// # Returns
///
/// A `HashSet<String>` containing the chromosome names from the tabix index header.
fn list_vcf_chromosomes_from_tbi(path: &Path) -> Result<HashSet<String>, anyhow::Error> {
    let path_tbi = assert_index_file_exists(path)?;

    // Open the tabix index file for your VCF (e.g., "sample.vcf.gz.tbi").
    let index = tabix::fs::read(path_tbi).context("Failed to open tabix file")?;

    let header = match index.header() {
        Some(x) => x,
        None => bail!("Could not find header in tabix file: {}", path.display()),
    };

    let refnames = header.reference_sequence_names();
    let chroms: HashSet<String> = refnames.iter().map(|n| n.to_string()).collect();

    Ok(chroms)
}

/// Checks if any of the chromosomes in the provided set have a "chr" prefix.
///
/// # Arguments
///
/// * `chroms` - A reference to a `HashSet<String>` containing chromosome names.
///
/// # Returns
///
/// `true` if one or more chromosome names start with "chr", otherwise `false`.
fn chroms_include_chr(chroms: &HashSet<String>) -> bool {
    let chroms_starting_with_chr: Vec<&str> = chroms
        .iter()
        .filter(|name| name.starts_with("chr"))
        .map(String::as_str)
        .collect();

    let n_chroms_starting_with_chr = chroms_starting_with_chr.len();

    if n_chroms_starting_with_chr > 0 {
        log::info!(
            "Found {} chromosomes starting with 'chr': {}",
            n_chroms_starting_with_chr,
            chroms_starting_with_chr.join(", ")
        );
        true
    } else {
        log::info!("Found no chromosomes starting with 'chr' in VCF",);
        false
    }
}

/// Filters the provided bedrs intervals to only include those on chromosomes present in the provided set.
///
/// # Arguments
///
/// * `intervals` - An `IntervalContainer` holding Bed4 intervals.
/// * `chroms` - A reference to a `HashSet<String>` of valid chromosome names.
///
/// # Returns
///
/// A new `IntervalContainer` with intervals restricted to the specified chromosomes.
fn filter_intervals_by_chroms(
    intervals: IntervalContainerBed4,
    chroms: &HashSet<String>,
) -> IntervalContainerBed4 {
    let regions_subset: IntervalContainer<bedrs::Bed4<String, u64, String>, String, u64> =
        intervals
            .iter()
            .filter(|bed| chroms.contains(bed.chr()))
            .cloned()
            .collect();

    regions_subset
}

/// Counts PASS variants within specified regions from a bgzipped VCF file.
///
/// This function queries the VCF file for each interval in the provided `IntervalContainer`
/// and counts variants that pass the filter criteria.
///
/// # Arguments
///
/// * `regions` - An `IntervalContainer` of Bed4 intervals defining the regions to query.
/// * `path_vcf` - A reference to the `PathBuf` of the VCF file.
///
/// # Returns
///
/// A `HashMap<String, u64>` where keys are region names (from the bed file) and values are the counts of PASS variants.
fn count_pass_variants_in_regions_bedrs(
    regions: IntervalContainer<bedrs::Bed4<String, u64, String>, String, u64>,
    path_vcf: &PathBuf,
) -> Result<HashMap<String, u64>, anyhow::Error> {
    let mut vcf_reader = vcf::io::indexed_reader::Builder::default().build_from_path(path_vcf)?;
    let header = vcf_reader.read_header()?;
    let chroms_in_vcf = list_vcf_chromosomes_from_tbi(path_vcf)?;

    // Filter regions for chromosomes present in the VCF to avoid query panics.
    let interval_subset = filter_intervals_by_chroms(regions, &chroms_in_vcf);

    // Set up HashMap to store results.
    let mut counts: HashMap<String, u64> = HashMap::new();

    // Run the query for each interval.
    for bed in interval_subset.iter() {
        let start = noodles::core::position::Position::try_from(bed.start() as usize + 1)?;
        let end = noodles::core::position::Position::try_from(bed.end() as usize)?;
        let region = noodles::core::region::Region::new(bed.chr().to_string(), start..=end);
        let query = vcf_reader.query(&header, &region)?;
        let name = bed.name();

        for result in query {
            let record = result.expect("failed to read variant in query subset");
            let is_pass = record
                .filters()
                .is_pass(&header)
                .expect("Failed to pull filter column");
            if is_pass {
                *counts.entry(name.to_owned()).or_insert(0) += 1;
            }
        }
    }

    Ok(counts)
}

/// Validates that the VCF/BCF file exists, is BGZF-compressed, and has an associated index file.
///
/// This function checks several conditions:
/// - The file exists at the given path.
/// - The file is a bgzipped VCF/BCF file (and not just gzipped).
/// - An index file exists for the given file.
///
/// # Arguments
///
/// * `path` - A reference to a `PathBuf` representing the file path.
///
/// # Returns
///
/// Returns `Ok(())` if all conditions are met, or an error otherwise.
pub fn check_vcf_exists_bgzipped_with_index(path: &PathBuf) -> Result<(), anyhow::Error> {
    // Check VCF exists.
    if !path.exists() {
        bail!("Could not find file: [{:#?}]", path);
    }
    // Check that VCF file is bgzipped and indexed so we can query it.
    let filetype = path_to_filetype(path);

    // Check filetype is indexed BCF/VCF file.
    match filetype {
        FileType::Bcf => bail!("BCF file must be compressed with bgzip and indexed with tabix. Try running `bgzip {}; tabix -p vcf {}.gz`", path.display(), path.display()),
        FileType::Vcf => bail!("VCF file must be compressed with bgzip and indexed with tabix. Try running `bgzip {}; tabix -p vcf {}.gz`", path.display(), path.display()),
        FileType::BcfCompressedBgzip => (),
        FileType::VcfCompressedBgzip => (),
        FileType::BcfCompressedGzip => bail!("BCF file must be compressed with bgzip not gzip. Problematic file: [{}]", path.display()),
        FileType::VcfCompressedGzip => bail!("VCF file must be compressed with bgzip not gzip. Problematic file: [{}]", path.display()),
        FileType::Unsure => bail!("Failed to recognise filetype from filename. Are you sure {} is a bgzipped VCF/BCF file?", path.display())
    };

    // Check if index is available.
    assert_index_file_exists(path)?;

    Ok(())
}

/// Returns the IGTCR regions bed file content as a string for the specified genome.
///
/// # Arguments
///
/// * `genome` - A reference to the `Genome` enum specifying the reference genome.
///
/// # Returns
///
/// A string slice containing the bed data.
pub fn genome_to_igtcr_regions_str(genome: &Genome) -> &str {
    match genome {
        Genome::GRCh38 => include_str!("data/igtcr_gene.38.bed4.bed"),
    }
}

/// Loads IGTCR regions from the embedded bed file into a bedrs `IntervalContainer`.
///
/// This function reads the bed file, parses each record into a Bed4 interval,
/// optionally strips the "chr" prefix, and then sorts the intervals.
///
/// # Arguments
///
/// * `genome` - A reference to the `Genome` enum specifying the reference genome.
/// * `include_chr_prefix` - A boolean indicating whether to keep the "chr" prefix in chromosome names.
///
/// # Returns
///
/// An `IntervalContainer` holding Bed4 intervals representing IGTCR regions, or an error if parsing fails.
pub fn get_igtcr_regions_bedrs(
    genome: &Genome,
    include_chr_prefix: bool,
) -> Result<IntervalContainerBed4, anyhow::Error> {
    let bed_str = genome_to_igtcr_regions_str(genome);
    // Create a CSV reader with tab as the delimiter.
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(bed_str.as_bytes());

    // Read bed file into a vector of Bed4 intervals.
    let mut regions: Vec<bedrs::Bed4<String, u64, String>> = vec![];

    // Deserialize each record into the ScarScapeBed4 struct.
    for result in reader.deserialize() {
        let record: ScarScapeBed4 =
            result.context("Failed to read IGTCR bed file. Please create a new github issue")?;

        let chrom_orig = record.chrom;
        let chrom = match include_chr_prefix {
            true => &chrom_orig,
            false => &chrom_orig
                .strip_prefix("chr")
                .context("failed to strip chr prefix from igtcr region")?
                .to_string(),
        };

        let interval: Bed4<String, u64, String> =
            bedrs::Bed4::new(chrom.clone(), record.start, record.end, record.name);
        regions.push(interval);
    }
    let mut interval_container = IntervalContainer::new(regions);

    // Sort the interval container.
    interval_container.sort();

    Ok(interval_container)
}

/// Struct representing a record from the IGTCR bed file used for deserialization.
#[derive(Deserialize, Debug)]
struct ScarScapeBed4 {
    /// Chromosome name.
    chrom: String,
    /// Start position.
    start: u64,
    /// End position.
    end: u64,
    /// Name field containing the receptor type.
    name: String,
}

pub struct SmallVariantSummary {
    small_variants: u64,
    small_variant_fail: u64,
    snvs: u64,
    doublets: u64,
    mnvs: u64,
    deletions: u64,
    insertions: u64,
    igtcr_small_varaints: u64,
    igh: u64,
    igk: u64,
    igl: u64,
    tra: u64,
    trb: u64,
    trd: u64,
    trg: u64,
    sbs96: SBS96,
}
impl SmallVariantSummary {
    // Compute total number of pass indels (insertions + deletions)
    fn indels(&self) -> u64 {
        self.insertions + self.deletions
    }

    /// Formats the small variant summary as a CSV string.
    ///
    /// # Arguments
    ///
    /// * `sample` - A string slice that holds the sample identifier.
    ///
    /// # Returns
    ///
    /// A `String` containing the CSV-formatted summary.
    pub fn to_csv(&self, sample: &str) -> String {
        format!(
            "sample,total,pass,fail,snv,doublet,mnv,indels,insertions,deletions,igtcr,igh,igk,igl,tra,trb,trd,trg\n\
             {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            sample,
            self.small_variants,
            self.small_variants,
            self.small_variant_fail,
            self.snvs,
            self.doublets,
            self.mnvs,
            self.indels(),
            self.insertions,
            self.deletions,
            self.igtcr_small_varaints,
            self.igh,
            self.igk,
            self.igl,
            self.tra,
            self.trb,
            self.trd,
            self.trg,
        )
    }

    /// Writes the CSV string representation of the summary to the specified file.
    ///
    /// # Arguments
    ///
    /// * `sample` - A string slice that holds the sample identifier.
    /// * `path` - A reference to the destination file path.
    ///
    /// # Errors
    ///
    /// Returns an error if the file already exists or if writing fails.
    pub fn write_to_csv(&self, sample: &str, path: &PathBuf) -> Result<(), anyhow::Error> {
        if path.exists() {
            bail!(
                "File [{}] already exists. Can NOT overwrite. Please manually delete and retry.",
                path.display()
            );
        }

        let csv_data = self.to_csv(sample);
        std::fs::write(path, csv_data)?;
        Ok(())
    }
}

impl fmt::Display for SmallVariantSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Small variant stats [total: {} (pass: {}, fail: {}) | SNV: {}, doublet: {}, MNV: {} | indels: {} (ins: {}, del: {}) | IGTCR: {} (IGH: {}, IGK: {}, IGL: {}, TRA: {}, TRB: {}, TRD: {}, TRG: {})]",
            self.small_variants,
            self.small_variants,
            self.small_variant_fail,
            self.snvs,
            self.doublets,
            self.mnvs,
            self.indels(),
            self.insertions,
            self.deletions,
            self.igtcr_small_varaints,
            self.igh,
            self.igk,
            self.igl,
            self.tra,
            self.trb,
            self.trd,
            self.trg,
        )
    }
}

pub struct SBS96 {
    counts: HashMap<String, u64>,
}

impl SBS96 {
    /// All 96 trinucleotide mutation contexts
    const CONTEXTS: [&'static str; 96] = [
        "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
        "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
        "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
        "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
        "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
        "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
        "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
        "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
        "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
        "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
        "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
        "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T",
    ];

    /// Create a new SBS96 with all context‐counts set to zero.
    pub fn new() -> Self {
        let mut counts = HashMap::with_capacity(Self::CONTEXTS.len());
        for &ctx in Self::CONTEXTS.iter() {
            counts.insert(ctx.to_string(), 0);
        }
        SBS96 { counts }
    }

    /// Increment the count for a given context (panics if invalid key).
    pub fn bump(&mut self, context: &str) {
        if let Some(count) = self.counts.get_mut(context) {
            *count += 1;
        } else {
            panic!("Invalid mutation context: {}", context);
        }
    }

    /// Get the count for a context (or zero if never initialized).
    pub fn get(&self, context: &str) -> u64 {
        *self.counts.get(context).unwrap_or(&0)
    }

    /// Iterate through all (context, count) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &u64)> {
        self.counts.iter()
    }
}

#[derive(Debug)]
enum SmallVariantTypes {
    SNV,
    DOUBLET,
    MNV,
    INSERTION,
    DELETION,
}

impl fmt::Display for SmallVariantTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SmallVariantTypes::SNV => "SNV",
            SmallVariantTypes::DOUBLET => "DOUBLET",
            SmallVariantTypes::MNV => "MNV",
            SmallVariantTypes::INSERTION => "INSERTION",
            SmallVariantTypes::DELETION => "DELETION",
        };
        write!(f, "{s}")
    }
}

pub fn small_variant_counts(
    path: &PathBuf, // Path to SNV VCF file
    genome: &Genome,
    ref_fasta: &PathBuf,
) -> Result<SmallVariantSummary, anyhow::Error> {
    // Read VCF file
    info!("Reading VCF file: {:#?}", path);

    // Check file exists and is an indexed BCF/VCF file so we can query it
    check_vcf_exists_bgzipped_with_index(path)?;

    // Check Index Exists
    let _index_path = assert_index_file_exists(path)?;
    let chroms = list_vcf_chromosomes_from_vcf_header(path)?;
    let vcf_uses_chr_prefix = chroms_include_chr(&chroms);

    // Read BED file into regions object describing TCR and BCR gene regions
    let igtcr_regions = get_igtcr_regions_bedrs(genome, vcf_uses_chr_prefix)
        .context("Failed to prepare IGTCR regions")?;

    // Read Small Variant VCF
    let mut vcf_reader = vcf::io::indexed_reader::Builder::default().build_from_path(path)?;
    info!("Successfully to read VCF [{path:#?}]");

    let mut npass: u64 = 0;
    let mut nfail: u64 = 0;
    let mut snvs: u64 = 0;
    let mut doublets: u64 = 0;
    let mut mnvs: u64 = 0;
    let mut insertions: u64 = 0;
    let mut deletions: u64 = 0;

    let header = vcf_reader.read_header()?;

    // Count Small Variants (General stats, PASS / FAIL / Total)
    for result in vcf_reader.records() {
        let record = result?;
        let pass = record
            .filters()
            .is_pass(&header)
            .expect("Failed to get sv filter");

        let variant = parse_vcf_record(&record)?;

        if variant.multiallelic {
            log::warn!("Small variant  VCF contains mutations with multiple alternative bases [{variant}]. Scarscape will only considering the first mutation. Run VCFs through VCF normalise to split multiallelics to have both alternative alleles counted")
        }

        if pass {
            npass += 1;
            match variant.class {
                SmallVariantTypes::SNV => snvs += 1,
                SmallVariantTypes::DOUBLET => doublets += 1,
                SmallVariantTypes::MNV => mnvs += 1,
                SmallVariantTypes::INSERTION => insertions += 1,
                SmallVariantTypes::DELETION => deletions += 1,
            }
        } else {
            nfail += 1
        }
    }

    // Count PASS small variants in IGTCR regions using vcftools region query
    let mut counts = count_pass_variants_in_regions_bedrs(igtcr_regions, path)?;

    // Ensure there is a value for all samples of interest
    counts.entry(String::from("IGH")).or_insert(0);
    counts.entry(String::from("IGK")).or_insert(0);
    counts.entry(String::from("IGL")).or_insert(0);
    counts.entry(String::from("TRA")).or_insert(0);
    counts.entry(String::from("TRB")).or_insert(0);
    counts.entry(String::from("TRD")).or_insert(0);
    counts.entry(String::from("TRG")).or_insert(0);

    // Sum igtcr small variants
    let pass_igtcr = counts.iter().map(|x| x.1).sum();

    let small_variant_summary = SmallVariantSummary {
        small_variants: npass,
        small_variant_fail: nfail,
        igtcr_small_varaints: pass_igtcr,
        snvs,
        doublets,
        mnvs,
        insertions,
        deletions,
        igh: *counts.get("IGH").unwrap(),
        igk: *counts.get("IGK").unwrap(),
        igl: *counts.get("IGL").unwrap(),
        tra: *counts.get("TRA").unwrap(),
        trb: *counts.get("TRB").unwrap(),
        trd: *counts.get("TRD").unwrap(),
        trg: *counts.get("TRG").unwrap(),
        sbs96: SBS96::new(),
    };

    Ok(small_variant_summary)
}

#[derive(Debug)]
pub struct Variant {
    chromosome: String,
    position: noodles::core::Position, // 1-based position (start)
    reference: String,
    alternative: String,
    delta: i64, // how big is the change to genome size (1 = inserting 1 base, -1 = deleting 1 base, etc)
    class: SmallVariantTypes,
    multiallelic: bool,
}
// Implement the `fmt::Display` trait for `Point`.
impl fmt::Display for Variant {
    // This trait requires the `fmt` function with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the write! macro to format the output.
        write!(
            f,
            "{}:{} {}>{} (delta: {}; class: {}; multiallelic:{})",
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
            self.delta,
            self.class,
            self.multiallelic
        )
    }
}

// Parses Variants from noodles VCF records - classifies VariantType. If multiallelic - drops all but the first allele.
pub fn parse_vcf_record(record: &noodles_vcf::Record) -> Result<Variant, anyhow::Error> {
    // Chromosome Name
    let chromosome = record.reference_sequence_name().to_string();

    // 1 based position
    let position = match record.variant_start() {
        Some(pos) => pos?,
        None => bail!("Couldn't get VCF entry position"),
    };

    // Reference base
    let reference = record.reference_bases().to_string();

    // Alternative base (pluck first alt for multiallelics)
    let altbases = record.alternate_bases();
    let multiallelic = altbases.len() > 1;
    let mut alternative: String = "".to_string();
    if !altbases.is_empty() {
        alternative = match altbases.iter().next() {
            Some(seq) => seq.unwrap_or("").to_string(),
            None => "".to_string(),
        };
    };

    // Classifying mutation type
    let alt_length = alternative.len();
    let ref_length = reference.len();
    let class = classify_small_variant(&reference, &alternative);
    let delta = alt_length as i64 - ref_length as i64;

    // Grab Allele Frequency
    Ok(Variant {
        chromosome,
        position,
        reference,
        alternative,
        delta,
        class,
        multiallelic,
    })
}

// Minimal starter classifier so the function compiles.
// Replace with your real rules.
fn classify_small_variant(reference: &str, alternative: &str) -> SmallVariantTypes {
    // If ALT was missing and we defaulted to "", don’t crash downstream.
    if reference.is_empty() || alternative.is_empty() {
        return SmallVariantTypes::SNV;
    }

    let r = reference.len();
    let a = alternative.len();

    if r == 1 && a == 1 {
        SmallVariantTypes::SNV
    } else if r == 2 && a == 2 {
        SmallVariantTypes::DOUBLET
    } else if r == a && r > 2 {
        SmallVariantTypes::MNV
    } else if a > r {
        SmallVariantTypes::INSERTION
    } else if r > a {
        SmallVariantTypes::DELETION
    } else {
        panic!("Unsure how to classify the small variant [{reference}>{alternative}]:");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Region;
    use std::collections::HashSet;
    use std::fs;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // --- Tests for RegionsIGTCR ---

    #[test]
    /// Tests the `get_igtcr_regions_bedrs` function to ensure the intervals are sorted,
    /// non-empty, and contain no overlapping regions after merging.
    fn test_get_igtcr_regions_bedrs() {
        let regions = get_igtcr_regions_bedrs(&Genome::GRCh38, true).unwrap();

        // Assert that the regions are sorted.
        assert!(regions.is_sorted(), "igtcr region bedfile should be sorted");

        // Assert that the regions are non-empty.
        assert!(
            !regions.is_empty(),
            "igtcr region bedfile should not be empty"
        );

        let merged_regions = regions.merge().unwrap();

        // Ensure no overlaps are present.
        assert_eq!(
            merged_regions.len(),
            regions.len(),
            "igtcr bed should not contain overlapping regions"
        )
    }

    #[test]
    fn test_regions_igtcr_methods_and_display() {
        // Create a dummy Region to use in all receptor vectors.
        let start = noodles::core::Position::try_from(100).unwrap();
        let end = noodles::core::Position::try_from(200).unwrap();
        let region = Region::new("chr1".to_string(), start..=end);
        let regions = RegionsIGTCR {
            igh: vec![region.clone()],
            igk: vec![region.clone()],
            igl: vec![region.clone()],
            tra: vec![region.clone()],
            trb: vec![region.clone()],
            trd: vec![region.clone()],
            trg: vec![region.clone()],
        };

        // Test lengths.
        assert_eq!(regions.len(), 7);
        assert_eq!(regions.count_bcr_regions(), 3);
        assert_eq!(regions.count_tcr_regions(), 4);

        // Test Display trait.
        let display_str = format!("{}", regions);
        assert!(display_str.contains("total: 7"));
        assert!(display_str.contains("BCR: 3"));
        assert!(display_str.contains("TCR 4"));
    }

    // --- Tests for StructuralVariantSummary ---

    #[test]
    fn test_structural_variant_summary_methods() {
        let summary = StructuralVariantSummary {
            pass: 10,
            fail: 5,
            igh: 2,
            igk: 1,
            igl: 1,
            tra: 2,
            trb: 1,
            trd: 1,
            trg: 2,
            total: 15,
        };

        // igtcr() should sum the counts for all receptor regions.
        assert_eq!(summary.igtcr(), 2 + 1 + 1 + 2 + 1 + 1 + 2);

        let csv_str = summary.to_csv("sample1");
        // Check that the CSV string contains the header and expected values.
        assert!(csv_str.contains("sample,total,pass,fail,igtcr,igh,igk,igl,tra,trb,trd,trg"));
        assert!(csv_str.contains("sample1,15,10,5,10,2,1,1,2,1,1,2"));

        // Test writing CSV to a temporary file.
        let tmp_file = NamedTempFile::new().unwrap();
        let path = tmp_file.path().to_path_buf();
        // Remove the file so that write_to_csv does not error on an existing file.
        fs::remove_file(&path).unwrap();
        assert!(summary.write_to_csv("sample1", &path).is_ok());
        let contents = fs::read_to_string(&path).unwrap();
        assert_eq!(contents, csv_str);
    }

    // --- Tests for is_bgzf and path_to_filetype ---

    #[test]
    fn test_is_bgzf_valid() {
        // Create a temporary file with a valid BGZF header.
        let mut tmp_file = NamedTempFile::new().unwrap();
        // Build an 18-byte header:
        // Bytes 0-2: standard gzip magic numbers.
        // Bytes 10-11: XLEN (set to 6).
        // Bytes 12-13: "BC" to indicate BGZF.
        let mut header = [0u8; 18];
        header[0] = 0x1F;
        header[1] = 0x8B;
        header[2] = 0x08;
        header[10] = 6;
        header[11] = 0;
        header[12] = b'B';
        header[13] = b'C';
        tmp_file.write_all(&header).unwrap();
        tmp_file.flush().unwrap();

        let file_path = tmp_file.path().to_str().unwrap();
        assert!(is_bgzf(file_path));
    }

    #[test]
    fn test_is_bgzf_invalid() {
        // Create a temporary file with a gzip header that is NOT BGZF.
        let mut tmp_file = NamedTempFile::new().unwrap();
        let mut header = [0u8; 18];
        header[0] = 0x1F;
        header[1] = 0x8B;
        header[2] = 0x08;
        header[10] = 6;
        header[11] = 0;
        // Instead of "BC", put different characters.
        header[12] = b'X';
        header[13] = b'Y';
        tmp_file.write_all(&header).unwrap();
        tmp_file.flush().unwrap();

        let file_path = tmp_file.path().to_str().unwrap();
        assert!(!is_bgzf(file_path));
    }

    #[test]
    fn test_path_to_filetype() {
        // For a file with extension ".vcf.bgz", we need to simulate a BGZF file.
        let mut tmp_bgzf = NamedTempFile::new().unwrap();
        let mut header = [0u8; 18];
        header[0] = 0x1F;
        header[1] = 0x8B;
        header[2] = 0x08;
        header[10] = 6;
        header[11] = 0;
        header[12] = b'B';
        header[13] = b'C';
        tmp_bgzf.write_all(&header).unwrap();
        tmp_bgzf.flush().unwrap();

        // Rename/copy the file to have a .vcf.bgz extension.
        let bgzf_path = tmp_bgzf.path().with_extension("vcf.bgz");
        fs::copy(tmp_bgzf.path(), &bgzf_path).unwrap();
        let filetype = path_to_filetype(&bgzf_path);
        match filetype {
            FileType::VcfCompressedBgzip => {}
            _ => panic!("Expected VcfCompressedBgzip, got {filetype:?} from file: {bgzf_path:#?}"),
        }

        // For a file with extension ".vcf.gz", simulate a file that is not BGZF.
        let mut tmp_gzip = NamedTempFile::new().unwrap();
        let mut header_gzip = [0u8; 18];
        header_gzip[0] = 0x1F;
        header_gzip[1] = 0x8B;
        header_gzip[2] = 0x08;
        header_gzip[10] = 6;
        header_gzip[11] = 0;
        // Set extra field to non-BGZF ("XY").
        header_gzip[12] = b'X';
        header_gzip[13] = b'Y';
        tmp_gzip.write_all(&header_gzip).unwrap();
        tmp_gzip.flush().unwrap();

        let gzip_path = tmp_gzip.path().with_extension("vcf.gz");
        fs::copy(tmp_gzip.path(), &gzip_path).unwrap();
        let filetype_gzip = path_to_filetype(&gzip_path);
        match filetype_gzip {
            FileType::VcfCompressedGzip => {}
            _ => panic!(
                "Expected VcfCompressedGzip, got {:?} from file: {:#?}",
                filetype_gzip, &gzip_path
            ),
        }
    }

    // --- Test for assert_index_file_exists ---

    #[test]
    fn test_assert_index_file_exists() {
        // Create a temporary file to act as the VCF.
        let tmp_vcf = NamedTempFile::new().unwrap();
        let vcf_path = tmp_vcf.path().to_path_buf();

        // Create an index file "<vcf_path>.tbi".
        let index_path = PathBuf::from(format!("{}.tbi", vcf_path.to_string_lossy()));
        fs::write(&index_path, "dummy index").unwrap();

        // Should return Ok with the index path.
        let result = assert_index_file_exists(&vcf_path);
        assert!(result.is_ok());
        let returned_path = result.unwrap();
        assert_eq!(returned_path, index_path);

        // Remove the index file and assert that an error is returned.
        fs::remove_file(&index_path).unwrap();
        let result_err = assert_index_file_exists(&vcf_path);
        assert!(result_err.is_err());
    }

    // --- Tests for chromosome and interval utilities ---

    #[test]
    fn test_chroms_include_chr() {
        let mut chroms = HashSet::new();
        chroms.insert("chr1".to_string());
        chroms.insert("chr2".to_string());
        assert!(chroms_include_chr(&chroms));

        let mut chroms_no_chr = HashSet::new();
        chroms_no_chr.insert("1".to_string());
        chroms_no_chr.insert("2".to_string());
        assert!(!chroms_include_chr(&chroms_no_chr));
    }

    #[test]
    fn test_filter_intervals_by_chroms() {
        // Create dummy intervals using bedrs::Bed4.
        use bedrs::Bed4;
        let interval1 = Bed4::new("chr1".to_string(), 100, 200, "region1".to_string());
        let interval2 = Bed4::new("chr2".to_string(), 150, 250, "region2".to_string());
        let interval3 = Bed4::new("chr3".to_string(), 300, 400, "region3".to_string());
        let intervals: IntervalContainerBed4 =
            vec![interval1.clone(), interval2.clone(), interval3.clone()]
                .into_iter()
                .collect();

        // Only allow "chr1" and "chr3".
        let mut allowed = HashSet::new();
        allowed.insert("chr1".to_string());
        allowed.insert("chr3".to_string());

        let filtered = filter_intervals_by_chroms(intervals, &allowed);
        assert_eq!(filtered.len(), 2);
        for interval in filtered.iter() {
            assert!(allowed.contains(interval.chr()));
        }
    }

    // --- Tests for IGTCR region loading functions ---

    #[test]
    fn test_genome_to_igtcr_regions_str() {
        let bed_str = genome_to_igtcr_regions_str(&Genome::GRCh38);
        assert!(!bed_str.is_empty(), "The bed string should not be empty");
    }

    #[test]
    fn test_get_igtcr_regions_bedrs_sort_and_merge() {
        // We already have a unit test for get_igtcr_regions_bedrs in the source.
        // This test will simply call that function and re-run some basic checks.
        let interval_container = get_igtcr_regions_bedrs(&Genome::GRCh38, true).unwrap();
        // Ensure the container is sorted.
        assert!(interval_container.is_sorted(), "Intervals should be sorted");
        // Ensure merging does not reduce the count (since the bed file is assumed to have no overlaps).
        let merged = interval_container.merge().unwrap();
        assert_eq!(
            merged.len(),
            interval_container.len(),
            "Merged intervals should have the same count as original intervals"
        );
    }

    // --- Test for check_vcf_exists_bgzipped_with_index ---
    #[test]
    fn test_check_vcf_exists_bgzipped_with_index() {
        // Create a temporary file with a valid BGZF header and a .vcf.bgz extension.
        let mut tmp_vcf = NamedTempFile::new().unwrap();
        let mut header = [0u8; 18];
        header[0] = 0x1F;
        header[1] = 0x8B;
        header[2] = 0x08;
        header[10] = 6;
        header[11] = 0;
        header[12] = b'B';
        header[13] = b'C';
        tmp_vcf.write_all(&header).unwrap();
        tmp_vcf.flush().unwrap();
        let vcf_path = tmp_vcf.path().with_extension("vcf.bgz");
        fs::copy(tmp_vcf.path(), &vcf_path).unwrap();

        // Create the corresponding index file "<vcf_path>.tbi".
        let index_path = PathBuf::from(format!("{}.tbi", vcf_path.to_string_lossy()));
        fs::write(&index_path, "dummy index").unwrap();

        // Should pass the check.
        let result = check_vcf_exists_bgzipped_with_index(&vcf_path);
        assert!(result.is_ok());
    }

    /// This test uses the real VCF file `tumor_sample.purple.sv.bgzipped.vcf.gz`
    /// (which must be accompanied by its tabix index file at
    /// `tumor_sample.purple.sv.bgzipped.vcf.gz.tbi` in the `testfiles` directory)
    /// and verifies that `sv_vcf_counts` correctly returns one PASS variant in the IGH region.
    #[test]
    fn test_sv_vcf_counts_tumor_sample() {
        // Construct the path to the test VCF file.
        let vcf_path = PathBuf::from("testfiles/tumor_sample.purple.sv.withIGHD.vcf.gz");

        // Ensure that the file exists. If it doesn't, the test should fail immediately.
        assert!(
            vcf_path.exists(),
            "The test VCF file {vcf_path:?} does not exist."
        );

        // Run the variant counting function.
        let summary = sv_vcf_counts(&vcf_path, &Genome::GRCh38)
            .expect("sv_vcf_counts failed on the tumor sample VCF");

        // We expect the file to contain one PASS variant overall.
        assert_eq!(summary.total, 325, "Expected one variant in total");
        assert_eq!(summary.pass, 163, "Expected one PASS variant");
        assert_eq!(summary.fail, 162, "Expected no FAIL variants");

        // Verify that the PASS variant is located in the IGH region.
        assert_eq!(summary.igh, 1, "Expected one PASS variant in IGH region");
        // All other regions should have a count of zero.
        assert_eq!(summary.igk, 0, "Expected zero PASS variants in IGK region");
        assert_eq!(summary.igl, 0, "Expected zero PASS variants in IGL region");
        assert_eq!(summary.tra, 0, "Expected zero PASS variants in TRA region");
        assert_eq!(summary.trb, 0, "Expected zero PASS variants in TRB region");
        assert_eq!(summary.trd, 0, "Expected zero PASS variants in TRD region");
        assert_eq!(summary.trg, 0, "Expected zero PASS variants in TRG region");
    }
}
