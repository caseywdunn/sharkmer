use clap::Parser;
use rand::prelude::SliceRandom;
use std::collections::HashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;

fn xtest_ints_to_kmers() {
    let ints = vec![0b01101100, 0b00111001, 0b10100110];
    // original	                // reverse complement
    // 0b01101100_00111001_10   0b01_10010011_11000110  >
    // 0b101100_00111001_1010   0b0101_10010011_110001  >
    // 0b1100_00111001_101001   0b100101_10010011_1100  >
    // 0b00_00111001_10100110   0b01100101_10010011_11  <
    //
    let expected = vec![
        0b01_10010011_11000110,
        0b0101_10010011_110001,
        0b100101_10010011_1100,
        0b00_00111001_10100110,
    ];
    let actual = ints_to_kmers(ints, 9);
    assert_eq!(actual, expected);
}

fn revcomp_kmer(kmer: u64, k: u8) -> u64 {
    let mut revcomp = 0;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

// For new, just return everything before an N. But in the future may return
// a vector of integer encoded sequences that were separated by N.
fn seq_to_ints(seq: &str) -> Vec<Vec<u8>> {
    let mut ints: Vec<u8> = Vec::new();
    let mut frame: u8 = 0;
    for (i, c) in seq.chars().enumerate() {
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            'N' => 4,
            _ => 5,
        };
        if base > 3 {
            break;
        }
        frame = (frame << 2) | base;
        if ((i + 1) % 4 == 0) & (i > 0) {
            ints.push(frame);
        }
    }
    vec![ints]
}

fn ints_to_kmers(ints: Vec<u8>, k: u8) -> Vec<u64> {
    let mut kmers: Vec<u64> = Vec::new();
    let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
    let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
    let mut n_valid = 0; // number of valid bases in the frame
    let mask: u64 = (1 << (2 * k)) - 1;

    // Iterate over the bases
    for (i, &int) in ints.iter().enumerate() {
        // Iterate over the bases in the integer
        for j in 0..4 {
            let base = ((int >> (j * 2)) & 3) as u64;
            frame = (frame << 2) | base;
            revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));
            n_valid += 1;
            if n_valid >= k {
                let forward = frame & mask;
                let reverse = revframe & mask;
                if forward < reverse {
                    kmers.push(forward);
                } else {
                    kmers.push(reverse);
                }
            }
        }
    }
    kmers
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
    xtest_ints_to_kmers();
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
    println!(
        "Yield {}",
        (n_bases_ingested as f64) / (n_bases_read as f64)
    );

    // Randomize the order of the reads in place
    let mut rng = rand::thread_rng();
    reads.shuffle(&mut rng);

    // Create the hash table
    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();
    for read in reads {
        let kmers = ints_to_kmers(read, args.k as u8);
        for kmer in kmers {
            let count = kmer_counts.entry(kmer).or_insert(0);
            *count += 1;
        }
    }

    // Create the histogram
    let mut histo: Vec<u64> = vec![0; args.histo_max as usize];
    for (_, count) in kmer_counts.iter() {
        if *count < args.histo_max {
            histo[*count as usize] += 1;
        }
    }

    // Write the histogram to a file
    let mut file = std::fs::File::create(format!("{}.histo", args.output)).unwrap();
    for (i, count) in histo.iter().enumerate() {
        writeln!(file, "{}\t{}", i, count).unwrap();
    }
}

#[cfg(test)]
mod test;
