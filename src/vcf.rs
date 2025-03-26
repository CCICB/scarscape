use anyhow::bail;
use anyhow::Context;
use log::info;
use noodles::core;
use noodles::core::Region;
use noodles::vcf;
use noodles_vcf::variant::record::Filters;
use regex::Regex;
use serde::Deserialize;
use std::collections::HashSet;
use std::fmt;
use std::io::Read;
use std::path::Path;
use std::{fs::File, path::PathBuf};

use crate::utils::Genome;

pub struct RegionsIGTCR {
    igh: Vec<Region>,
    igk: Vec<Region>,
    igl: Vec<Region>,
    tra: Vec<Region>,
    trb: Vec<Region>,
    trd: Vec<Region>,
    trg: Vec<Region>,
}
impl RegionsIGTCR {
    // Get a single vector representing all BCR (IG) regions
    fn _get_ig_regions(&self) -> Vec<Region> {
        [self.igh.clone(), self.igk.clone(), self.igl.clone()].concat()
    }

    // Get a single vector representing all TCR (TR) regions
    fn _get_tcr_regions(&self) -> Vec<Region> {
        [
            self.tra.clone(),
            self.trb.clone(),
            self.trd.clone(),
            self.trg.clone(),
        ]
        .concat()
    }

    // Len should get total number of regions across all types
    fn len(&self) -> usize {
        self.igh.len()
            + self.igk.len()
            + self.igl.len()
            + self.tra.len()
            + self.trb.len()
            + self.trd.len()
            + self.trg.len()
    }

    fn count_bcr_regions(&self) -> usize {
        self.igh.len() + self.igk.len() + self.igl.len()
    }
    fn count_tcr_regions(&self) -> usize {
        self.tra.len() + self.trb.len() + self.trd.len() + self.trg.len()
    }
}

impl fmt::Display for RegionsIGTCR {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Write strictly the first element into the supplied output
        // stream: `f`. Returns `fmt::Result` which indicates whether the
        // operation succeeded or failed. Note that `write!` uses syntax which
        // is very similar to `println!`.
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

#[derive(Debug)]
pub struct StructuralVariantSummary {
    pass: u64,
    fail: u64,
    igh: u64, // Pass variants in igtcr regions
    igk: u64, // Pass variants in igtcr regions
    igl: u64, // Pass variants in igtcr regions
    tra: u64, // Pass variants in igtcr regions
    trb: u64, // Pass variants in igtcr regions
    trd: u64, // Pass variants in igtcr regions
    trg: u64, // Pass variants in igtcr regions
    total: u64,
}
impl StructuralVariantSummary {
    fn igtcr(&self) -> u64 {
        self.igh + self.igk + self.igl + self.tra + self.trb + self.trg + self.trd
    }

    /// Formats the summary as a CSV string.
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

    /// Writes the CSV string to the specified file path.
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
    // This trait requires `fmt` with this exact signature.
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
    let chroms = list_vcf_chromosomes(path)?;
    let vcf_uses_chr_prefix = chroms_include_chr(&chroms);

    // Read BED file into regions object describing TCR and BCR gene regions
    // let path_igtcr_region_bed = genome_to_igtcr_regions();
    let igtcr_regions = get_igtcr_regions(genome, vcf_uses_chr_prefix)
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
        // let chrom = record.reference_sequence_name();
        // let pos = record // 1 based inclusive start
        //     .variant_start()
        //     .context("Failed to get start coord")?;

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

    // Count PASS SVs in igtcr regions using vcftools region query
    let n_igh = count_pass_variants_in_regions(&igtcr_regions.igh, path)?;
    let n_igk = count_pass_variants_in_regions(&igtcr_regions.igk, path)?;
    let n_igl = count_pass_variants_in_regions(&igtcr_regions.igl, path)?;
    let n_tra = count_pass_variants_in_regions(&igtcr_regions.tra, path)?;
    let n_trb = count_pass_variants_in_regions(&igtcr_regions.trb, path)?;
    let n_trd = count_pass_variants_in_regions(&igtcr_regions.trd, path)?;
    let n_trg = count_pass_variants_in_regions(&igtcr_regions.trg, path)?;

    Ok(StructuralVariantSummary {
        pass: npass,
        fail: nfail,
        igh: n_igh,
        igk: n_igk,
        igl: n_igl,
        tra: n_tra,
        trb: n_trb,
        trd: n_trd,
        trg: n_trg,
        total: ntotal,
    })
}

enum FileType {
    Vcf,
    VcfCompressedGzip,
    VcfCompressedBgzip,
    Bcf,
    BcfCompressedBgzip,
    BcfCompressedGzip,
    Unsure,
}
fn path_to_filetype(path: &Path) -> FileType {
    let path_str = path.to_str().unwrap_or("");

    let vcf_gz_re =
        Regex::new(r"\.vcf\.[bg]z$").expect("Failed to create regex for vcf (compressed)");
    let vcf_re = Regex::new(r"\.vcf$").expect("Failed to create regex for vcf (uncompressed)");
    let bcf_gz_re =
        Regex::new(r"\.bcf\.[bg]z$").expect("Failed to create regex for bcf (compressed)");
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
        eprintln!("Ends with unsure");
        FileType::Unsure
    }
}

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

fn list_vcf_chromosomes(path: &PathBuf) -> Result<HashSet<String>, anyhow::Error> {
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

/// Filter a region vector to include only chromosomes in chroms
fn filter_regions(regions: &[Region], chroms: &HashSet<String>) -> Vec<Region> {
    // Filter for regions in VCF (otherwise query breaks)
    // let chroms_in_igtcr_bed = regions
    //     .iter()
    //     .map(|region| region.name().to_string())
    //     .collect::<HashSet<String>>();

    let regions_subset: Vec<Region> = regions
        .iter()
        .filter(|region| chroms.contains(region.name()))
        .cloned()
        .collect();

    // let n_regions_dropped = regions.len() - regions_subset.len();

    regions_subset
}

/// Count Pass Variants in Regions. Note that if regions are overlapping some variants could be double counted.
fn count_pass_variants_in_regions(
    regions: &[Region],
    path_vcf: &PathBuf,
) -> Result<u64, anyhow::Error> {
    let mut vcf_reader = vcf::io::indexed_reader::Builder::default().build_from_path(path_vcf)?;
    let header = vcf_reader.read_header()?;
    let chroms_in_vcf = list_vcf_chromosomes(path_vcf)?;
    let mut n_pass_variants = 0;

    // filter regions for chroms in vcf because otherwise query panics.
    let region_subset = filter_regions(regions, &chroms_in_vcf);
    // Run the actual query
    for region in region_subset {
        let query = vcf_reader.query(&header, &region)?;

        for result in query {
            let record = result.expect("failed to read variant in query subset");
            let is_pass = record
                .filters()
                .is_pass(&header)
                .expect("Failed to pull filter column");
            if is_pass {
                n_pass_variants += 1
            }
        }
    }
    Ok(n_pass_variants)
}

/// Checks that the VCF/BCF file at the given path exists, is bgzipped, and has an associated index.
///
/// This function performs several validations:
/// - It verifies that the file exists at the specified path.
/// - It determines the file type and ensures that the file is either a bgzipped VCF/BCF file or,
///   if not, returns an error with guidance on compressing and indexing the file correctly.
/// - It checks that an index file exists for the given file.
///
/// # Arguments
///
/// * `path` - A reference to a [`PathBuf`] representing the location of the VCF/BCF file.
///
/// # Returns
///
/// * `Ok(())` if the file exists, is bgzipped, and the index file is available.
/// * `Err(anyhow::Error)` if any of the validations fail. Errors are returned in cases such as:
///   - The file does not exist.
///   - The file is not bgzipped (e.g., it is gzipped instead).
///   - The file type is unrecognized or does not have an associated index.
///
pub fn check_vcf_exists_bgzipped_with_index(path: &PathBuf) -> Result<(), anyhow::Error> {
    // Check VCF exists
    if !path.exists() {
        bail!("Could not find file: [{:#?}]", path);
    }
    // Check that VCF file is bgzipped and indexed so we can query it
    let filetype = path_to_filetype(path);

    // Check filetype is indexed BCF/VCF file
    match filetype {
    FileType::Bcf => bail!("BCF file must be compressed with bgzip and indexed with tabix. Try running `bgzip {}; tabix -p vcf {}.gz`", path.display(), path.display()),
    FileType::Vcf => bail!("VCF file must be compressed with bgzip and indexed with tabix. Try running `bgzip {}; tabix -p vcf {}.gz`", path.display(), path.display()),
    FileType::BcfCompressedBgzip => (),
    FileType::VcfCompressedBgzip => (),
    FileType::BcfCompressedGzip => bail!("BCF file must be compressed with bgzip not gzip. Problematic file: [{}]", path.display()),
    FileType::VcfCompressedGzip => bail!("VCF file must be compressed with bgzip not gzip. Problematic file: [{}]", path.display()),
    FileType::Unsure => bail!("Failed to recognise filetype from filename. Are you sure {} is a bgzipped VCF/BCF file?", path.display())
  };

    // Check if index is available
    assert_index_file_exists(path)?;

    Ok(())
}

/// Return bed string
pub fn genome_to_igtcr_regions_str(genome: &Genome) -> &str {
    match genome {
        Genome::GRCh38 => include_str!("data/igtcr_gene.38.bed4.bed"),
    }
}

#[derive(Deserialize, Debug)]
struct Bed4 {
    // Define fields matching your JSON structure
    chrom: String,
    start: usize,
    end: usize,
    name: String, // ...
}

/// Load an IGTCR regions bed into a RegionsIGTCR struct, for use downstream in querying of SV VCFs
pub fn get_igtcr_regions(
    genome: &Genome,
    include_chr_prefix: bool,
) -> Result<RegionsIGTCR, anyhow::Error> {
    let bed_str = genome_to_igtcr_regions_str(genome);
    // Create a CSV reader with tab as the delimiter.
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(bed_str.as_bytes());

    // Read Bed into regions
    let mut igh_regions: Vec<Region> = vec![];
    let mut igk_regions: Vec<Region> = vec![];
    let mut igl_regions: Vec<Region> = vec![];
    let mut tra_regions: Vec<Region> = vec![];
    let mut trb_regions: Vec<Region> = vec![];
    let mut trd_regions: Vec<Region> = vec![];
    let mut trg_regions: Vec<Region> = vec![];

    // Deserialize each record into the Record struct.
    for result in reader.deserialize() {
        let record: Bed4 =
            result.context("Failed to read IGTCR bed file. Please create a new github issue")?;
        let start = core::Position::try_from(record.start)?;
        let end = core::Position::try_from(record.end)?;
        let name = record.name;
        let chrom_orig = record.chrom;

        let chrom = match include_chr_prefix {
            true => &chrom_orig,
            false => &chrom_orig
                .strip_prefix("chr")
                .context("failed to strip chr prefix from igtcr region")?
                .to_string(),
        };

        let region = Region::new(chrom.clone(), start..=end);

        if name.contains("|IGH|") {
            igh_regions.push(region.clone());
        } else if name.contains("|IGK|") {
            igk_regions.push(region.clone());
        } else if name.contains("|IGL|") {
            igl_regions.push(region.clone());
        } else if name.contains("|TRA|") {
            tra_regions.push(region.clone());
        } else if name.contains("|TRB|") {
            trb_regions.push(region.clone());
        } else if name.contains("|TRD|") {
            trd_regions.push(region.clone());
        } else if name.contains("|TRG|") {
            trg_regions.push(region.clone());
        } else {
            bail!("Could not identify receptor region  from name field of bed file. Problematic name: [{}]. Problematic region: [{}:{}-{}]", name, &chrom_orig, start, end);
        }
    }
    let regions = RegionsIGTCR {
        igh: igh_regions,
        igk: igk_regions,
        igl: igl_regions,
        tra: tra_regions,
        trb: trb_regions,
        trd: trd_regions,
        trg: trg_regions,
    };

    info!("Successfully loaded IG & TCR regions: {}", regions);

    Ok(regions)
}
