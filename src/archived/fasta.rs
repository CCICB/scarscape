use std::path::PathBuf;

use anyhow::Error;
use noodles::{
    core::{Position, Region},
    fasta::record::Sequence,
};

/// Fetch Helper
pub fn fetch_seq(path: &PathBuf, region: &Region) -> Result<Sequence, Error> {
    // Parse file
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(path)?;

    let record = reader.query(region)?;
    let seq = record.sequence().to_owned();
    // log::info!("Fetched Sequence: {:#?}", seq);
    Ok(seq)
}

#[test]
fn seq_extract() {
    let genomepath = PathBuf::from("../testfiles/sbs96_contexts.fasta");
    assert!(genomepath.exists());
    let start = Position::new(3).unwrap();
    let end = Position::new(6).unwrap();
    let region = Region::new("chr1", start..=end);

    let seq_result = fetch_seq(&genomepath, &region);
    match seq_result {
        Err(err) => eprintln!("{:#?}", err),
        Ok(_) => eprintln!("No error"),
    }
    // eprintln!("{:#?}", seq);
}
