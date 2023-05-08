use clap::Parser;
use std::path::Path;

/// Count k-mers in a set of fastq.gz files, with an option to assess cumulative subsets
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// k-mer length
    #[arg(short, default_value_t = 21)]
    k: u32,

    /// Maximum value for histogram
    #[arg(long, default_value_t = 10000)]
    histo_max: u64,

    /// Number of chunks to divide the data into
    #[arg(short, default_value_t = 10)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short, long, default_value_t = 0)]
    max_reads: u64,

    /// Directory and filename prefix for analysis output, for example out_dir/Nanomia-bijuga
    #[arg(short, long, default_value_t = String::from("sample") )]
    output: String,

    /// Input files, fastq
    #[arg(required = true)]
    input: Vec<String>,
}
fn main() {
    // Ingest command line arguments
    let args = Args::parse();

    // Print the program name and version
    println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    // Print the arguments
    println!("{:?}", args);

    // Parse the output path and create directories if necessary
    let path = Path::new(&args.output);
    let out_name = path.file_name().unwrap().to_str().unwrap(); // This is the prefix of the output files
    let directory = path.parent().unwrap().to_str().unwrap();
    let _ = std::fs::create_dir_all(directory);

    // Check that the arguments are valid
    assert!(
        args.k < 32,
        "k must be less than 32 due to use of 64 bit integers to encode kmers"
    );
    assert!(args.k > 0, "k must be greater than 0");
    assert!(args.k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tests() {
        assert_eq!(2 + 2, 4);
    }

}
