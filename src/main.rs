use clap::Parser;
use std::path::Path;
use std::io::BufRead;

// For new, just return everything before an N. But in the future may return
// a vector of integer encoded sequences that were separated by N.
fn seq_to_ints(seq:&String) -> Vec<Vec<u8>> {
    let mut ints: Vec<u8> = Vec::new();
    let mut frame: u8 = 0;
    for (i,c) in seq.chars().enumerate() {
        let base = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            'N' => 4,
            _ => 5,
        };
        if base > 3{
            break;
        }
        frame = (frame << 2) | base;
        if (i % 4 == 0) & (i>0) {
            ints.push(frame);
        }
    }
    vec![ints]
}


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

    
    // Ingest the fastq files
    let mut reads: Vec<Vec<u8>> = Vec::new();
    let mut n_reads_read = 0;
    let mut n_bases_read = 0;
    'processing_files: for file_name in args.input {
        let mut line_n = 0;
        // Open the file for buffered reading
        let file = std::fs::File::open(file_name).unwrap();
        let reader = std::io::BufReader::new(file);
        // Iterate over the lines of the file
        for line in reader.lines() {
            line_n += 1;
            if line_n % 4 == 2 {
                // This is a sequence line
                let line = line.unwrap();
                n_bases_read += line.len();
                let ints = seq_to_ints(&line);
                reads.extend(ints);
                n_reads_read += 1;
                
            }
            if args.max_reads > 0 && n_reads_read >= args.max_reads as usize {
                break 'processing_files;
            }
        }
    }

    let n_bases_ingested = reads.iter().map(|x| x.len()).sum::<usize>() * 4;

    // Print some stats
    println!("Read {} reads", n_reads_read);
    println!("Read {} bases", n_bases_read);
    println!("Ingested {} subreads", reads.len());
    println!("Ingested {} bases", n_bases_ingested);
    println!("Yield {}",  (n_bases_ingested as f64)/(n_bases_read as f64));


}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tests() {
        assert_eq!(2 + 2, 4);
    }

}
