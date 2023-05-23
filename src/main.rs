use clap::Parser;
use rand::prelude::SliceRandom;
use std::collections::HashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use polars::prelude::*;


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
            // Get the base from the left side of the integer,
            // move it to the least two significant bits and mask it
            let base = ((int >> ((3 - j) * 2)) & 3) as u64;

            frame = (frame << 2) | base;
            revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));
            n_valid += 1;
            if n_valid >= k {
                let forward = frame & mask;
                let reverse = revframe & mask;

                // assert_eq!(forward, revcomp_kmer(reverse, k));

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

fn count_histogram(kmer_counts: &HashMap<u64, u64>, histo_max: u64) -> Vec<u64> {
    // Create a histogram of counts
    let mut histo: Vec<u64> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for count in kmer_counts.values() {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        } else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    histo
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
    print!("Ingested reads...");
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

    println!(" done");
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
    print!("Randomizing read order...");
    let mut rng = rand::thread_rng();
    reads.shuffle(&mut rng);
    println!(" done");

    // Create a hash table for each of n chunks of reads
    if reads.len() < args.n {
        // Throw an error
        panic!("Number of reads is less than number of chunks");
    }
    let mut chunk_size = reads.len() / args.n;

    // Create a vector of hash tables
    let mut chunk_kmer_counts: Vec<HashMap<u64, u64>> = Vec::new();
    for _ in 0..args.n {
        chunk_kmer_counts.push(HashMap::new());
    }

    print!("Hashing each chunk of reads...");
    // Iterate over the chunks
    for (i, chunk) in reads.chunks(chunk_size).enumerate() {
        println!("Processing chunk {}", i);
        // Create a hash table for this chunk
        let mut kmer_counts: HashMap<u64, u64> = HashMap::new();
        for read in chunk {
            let kmers = ints_to_kmers(read.to_vec(), args.k as u8);
            for kmer in kmers {
                let count = kmer_counts.entry(kmer).or_insert(0);
                *count += 1;
            }
        }
        chunk_kmer_counts[i] = kmer_counts;
    }
    println!(" done");

    // Create the histograms
    print!("Creating histograms...");
    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();

    // Create a polars dataframe with max_reads+2 rows and n columns, int32 type and fill it with zeros
    let mut histo_df = DataFrame::new(Vec::new()).unwrap();

    // Iterate over the chunks
    for chunk_kmer_count in chunk_kmer_counts {
        // Add the counts from chunk_kmer_count to the corresponding entries of kmer_counts, creating new entries as needed
        for (kmer, kmer_count) in chunk_kmer_count {
            let count = kmer_counts.entry(kmer).or_insert(0);
            *count += kmer_count;
        }

        // Create a histogram of counts
        let histo = count_histogram(&kmer_counts, args.histo_max);

        // Append the histogram to the dataframe
        let mut histo_series = Series::new("", histo);
        histo_df.with_column(histo_series).unwrap();

    }
    println!(" done");

    // Write the histograms to a tab delimited file, with the first column being the count
    // Skip the first row, which is the count of 0. Do not include a header
    print!("Writing histograms to file...");
    let mut file = std::fs::File::create(format!("{}.histo", out_name)).unwrap();

    CsvWriter::new(&file)
        .has_header(false)
        .with_delimiter(b'\t')
        .finish(histo_df);

    println!(" done");
}

#[cfg(test)]
mod test;
