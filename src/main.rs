use std::{path::PathBuf, time::SystemTime};

use anyhow::{bail, Context};
use fern::colors::ColoredLevelConfig;
use log::{info, warn};
use scarscape::utils::Genome;

use clap::Parser;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// A TSV mapping samples to their corresponding mutation data files. Should contain 4 columns (with the following headers)
    /// [sample]: sample identifier
    /// [snv]: path to VCF file describing SNVs and INDELs
    /// [sv]: path to VCF file describing structural variants (1 entry per breakend)
    /// [cnv]: path to TSV describing segment copy number. Must contain the columns: chromosome, start (1-based), end (inclusive), copyNumber, minorAlleleCopyNumber
    #[arg(short, long, value_name = "manifest.tsv", verbatim_doc_comment)]
    manifest: PathBuf,

    /// Reference genome used for mutation calling.
    #[arg(value_enum, short, long, value_name = "hg38/hg19")]
    genome: Genome,

    /// Output folder. Will create if it does not exist
    #[arg(short, long, value_name = "scars", default_value = "scars")]
    outdir: PathBuf,
}

fn setup_logger() -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::new().info(fern::colors::Color::Green);

    // // Create name of logfile
    // let current_time = humantime::format_rfc3339_seconds(SystemTime::now())
    //     .to_string()
    //     .replace(":", "_");
    // let logfile = PathBuf::from(format!("scarscape.{current_time}.log"));

    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{} {} {}] {}",
                humantime::format_rfc3339_seconds(SystemTime::now()),
                colors.color(record.level()),
                record.target(),
                message
            ))
        })
        .level(log::LevelFilter::Info)
        .chain(std::io::stderr())
        // .chain(fern::log_file(logfile)?)
        .apply()?;
    Ok(())
}

fn main() {
    // Set up logging

    setup_logger().expect("Error setting up logging");
    // Log any bubbled up errors and panic
    if let Err(err) = run() {
        scarscape::utils::log_error_and_panic(err);
        std::process::exit(1)
    }
}

fn run() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    // Settings
    let genome: Genome = cli.genome;
    let path_manifest = cli.manifest;

    // Read Manifest
    let manifest = scarscape::utils::read_manifest(&path_manifest)?;

    // Set Up output Directory
    let outdir_raw = cli.outdir;
    let outdir = scarscape::utils::make_output_directory(&outdir_raw)?; // Returns absolutized path

    // Log Configuration
    info!(
        "User Supplied Setttings:  Genome: {}\nSV VCF: {}\nOutput Directory: {}",
        genome,
        path_manifest.display(),
        outdir.display()
    );

    // Iterate through manifest
    for manifest_entry in manifest {
        // Structural Variant Counts
        let _ = match manifest_entry.sv {
            Some(svpath) => {
                info!(
                    "Extracting SV counts for sample [{}]",
                    &manifest_entry.sample
                );
                let sv_summary = scarscape::vcf::sv_vcf_counts(&svpath, &genome)?;
                info!("Successfully computed SV statistics. {}", sv_summary);
                let outfile = &outdir.join(format!("{}.svcounts.csv", &manifest_entry.sample));
                info!("Writing SV counts to file: {}", outfile.display());

                if outfile.exists() {
                    bail!("Failed to write structural variant counts to [{}] as file already exists. Please manually remove and try again", outfile.display());
                };

                sv_summary
                    .write_to_csv(&manifest_entry.sample, outfile)
                    .context(format!(
                        "Failed to write structural variant counts to [{}]",
                        &outfile.display()
                    ))
            }
            None => {
                warn!(
                    "Skipping computation of SV counts for sample [{}] as no SV VCF is supplied",
                    &manifest_entry.sample
                );
                Ok(())
            }
        };
    }
    Ok(())
}
