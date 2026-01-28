use crate::error::*;
use std::{fs::File, io::Read, path::Path};

/// Case insensitive extension matching
pub fn has_extension<P, I, S>(path: P, extensions: I) -> Result<bool>
where
    P: AsRef<Path>,
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let p = path.as_ref();
    let filename = p
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| Error::NonUtf8Path {
            path: p.to_path_buf(),
        })?;

    let filename = filename.to_ascii_lowercase();

    Ok(extensions.into_iter().any(|ext| {
        let ext = ext.as_ref().trim_start_matches('.').to_ascii_lowercase();
        filename == ext || filename.ends_with(&format!(".{ext}"))
    }))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfFileType {
    Vcf,
    VcfGzip,
    VcfBgzip,
    Bcf,
    BcfGzip,
    BcfBgzip,
}

impl std::fmt::Display for VcfFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            VcfFileType::Vcf => "Vcf",
            VcfFileType::VcfGzip => "VcfGzip",
            VcfFileType::VcfBgzip => "VcfBgzip",
            VcfFileType::Bcf => "Bcf",
            VcfFileType::BcfGzip => "BcfGzip",
            VcfFileType::BcfBgzip => "BcfBgzip",
        };
        write!(f, "{s}")
    }
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
pub fn path_to_filetype(path: &Path) -> Result<VcfFileType> {
    // Expected patterns
    // - vcf: .vcf, .vcf.gz, .vcf.bgz
    // - bcf: .bcf, .bcf.gz, .bcf.bgz
    let is_vcf_bgz = has_extension(path, ["vcf.bgz"])?;
    let is_vcf_gz = has_extension(path, ["vcf.gz"])?;
    let is_vcf = has_extension(path, ["vcf"])?;

    let is_bcf_bgz = has_extension(path, ["bcf.bgz"])?;
    let is_bcf_gz = has_extension(path, ["bcf.gz"])?;
    let is_bcf = has_extension(path, ["bcf"])?;

    // Compressed: use header sniff to decide bgzf vs gzip (for .gz/.bgz)
    if is_vcf_bgz || is_vcf_gz {
        return Ok(if is_bgzf(path)? {
            VcfFileType::VcfBgzip
        } else {
            VcfFileType::VcfGzip
        });
    }

    if is_bcf_bgz || is_bcf_gz {
        return Ok(if is_bgzf(path)? {
            VcfFileType::BcfBgzip
        } else {
            VcfFileType::BcfGzip
        });
    }

    if is_vcf {
        return Ok(VcfFileType::Vcf);
    }

    if is_bcf {
        return Ok(VcfFileType::Bcf);
    }

    Err(Error::UnsupportedExtension {
        path: path.to_path_buf(),
        expected: ".vcf, .vcf.gz, .vcf.bgz, .bcf, .bcf.gz, .bcf.bgz".to_string(),
    })
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
fn is_bgzf(path: &Path) -> Result<bool> {
    let mut file = File::open(path).map_err(|_| Error::OpenFailed {
        path: path.to_path_buf(),
    })?;

    let mut header = [0u8; 18];
    file.read_exact(&mut header)
        .map_err(|_| Error::ReadFailed {
            path: path.to_path_buf(),
        })?;

    // gzip magic
    if header[0] != 0x1F || header[1] != 0x8B || header[2] != 0x08 {
        return Ok(false);
    }

    let xlen = u16::from_le_bytes([header[10], header[11]]);

    Ok(xlen >= 6 && header[12] == b'B' && header[13] == b'C')
}
