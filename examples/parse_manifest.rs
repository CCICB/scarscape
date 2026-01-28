// cargo run --example parse_manifest examples/data/manifest.csv
use scarscape::io::manifest::ManifestReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: read_manifest <manifest.tsv>");

    // Production entry: open from path and iterate domain records
    let reader = ManifestReader::from_path(path)?;

    for row in reader {
        let manifest_entry = row?;
        println!("{manifest_entry}");
    }

    Ok(())
}
