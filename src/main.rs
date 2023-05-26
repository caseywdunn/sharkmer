use bloom::{BloomFilter, ASMS};
use clap::Parser;
use rand::prelude::SliceRandom;
use rayon::prelude::*;
//use std::collections::HashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use rustc_hash::FxHashMap;

// Create a structure with a hashmap for kmer counts and a u64 for the number of singleton kmers
struct KmerSummary {
    kmer_counts: FxHashMap<u64, u64>,
    n_singletons: u64,
}

// For new, just return everything before an N. But in the future may return
// a vector of integer encoded sequences that were separated by N.
fn seq_to_ints(seq: &str) -> Vec<Vec<u8>> {
    let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4);
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
    let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k as usize) + 1);
    let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
    let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
    let mut n_valid = 0; // number of valid bases in the frame
    let mask: u64 = (1 << (2 * k)) - 1;

    // Iterate over the bases
    for (_i, &int) in ints.iter().enumerate() {
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

fn count_histogram(kmer_counts: &FxHashMap<u64, u64>, histo_max: u64) -> Vec<u64> {
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
    let start = std::time::Instant::now();
    print!("Ingested reads...");
    std::io::stdout().flush().unwrap();
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
    println!("  Read {} reads", n_reads_read);
    println!("  Read {} bases", n_bases_read);
    println!("  Ingested {} subreads", reads.len());
    println!("  Ingested {} bases", n_bases_ingested);
    println!(
        "  Yield {}",
        (n_bases_ingested as f64) / (n_bases_read as f64)
    );
    println!("  Time to ingest reads: {:?}", start.elapsed());

    // Randomize the order of the reads in place
    print!("Randomizing read order...");
    std::io::stdout().flush().unwrap();
    let mut rng = rand::thread_rng();
    reads.shuffle(&mut rng);
    println!(" done");

    

    // Find kmers that occur multiple times with bloom filter
    let start = std::time::Instant::now();
    println!("Identifying kmers that occur more than once...");
    let mut pre_bloom = BloomFilter::with_rate(0.01, 4_294_967_295);
    let mut multi_bloom = BloomFilter::with_rate(0.01, 4_294_967_295);
    let mut n_kmers: u64 = 0;
    let mut n_multi_kmers: u64 = 0;

    // Print a progress bar, 2% per character
    println!("--------------------------------------------------- 100%");
    let reads_per_2_percent = (reads.len() / 50) as u64;

    for (n_reads, read) in (0_u64..).zip(reads.iter()) {
        if n_reads % reads_per_2_percent == 0 {
            print!(".");
            std::io::stdout().flush().unwrap();
        }

        let kmers = ints_to_kmers(read.to_vec(), args.k as u8);
        for kmer in kmers {
            n_kmers += 1;
            if pre_bloom.contains(&kmer) {
                multi_bloom.insert(&kmer);
                n_multi_kmers += 1;
            } else {
                pre_bloom.insert(&kmer);
            }
        }
    }
    println!(" done");
    println!("  Number of kmers processed: {}", n_kmers);
    println!("  Number of multi kmers: {}", n_multi_kmers);
    println!("  Number of once-off kmers: {}", n_kmers - n_multi_kmers);

    // Print the time taken to construct the bloom filter
    let duration = start.elapsed();
    println!("  Time to construct bloom filter: {:?}", duration);

    // Create a hash table for each of n chunks of reads
    if reads.len() < args.n {
        panic!("Number of reads is less than number of chunks");
    }
    let chunk_size = reads.len() / args.n;

    let start = std::time::Instant::now();
    print!("Hashing each chunk of reads...");
    std::io::stdout().flush().unwrap();
    // Iterate over the chunks
    let n: usize = args.n;
    let chunk_kmer_counts: Vec<KmerSummary> = (0..n)
        .into_par_iter()
        .map(|i| {
            let start = i * chunk_size;
            let end = (i + 1) * chunk_size;
            let mut kmer_counts: FxHashMap<u64, u64> = FxHashMap::default();
            let mut singles: u64 = 0;
            for read in reads[start..end].iter() {
                let kmers = ints_to_kmers(read.to_vec(), args.k as u8);
                for kmer in kmers {
                    if multi_bloom.contains(&kmer) {
                        let count = kmer_counts.entry(kmer).or_insert(0);
                        *count += 1;
                    } else {
                        singles += 1;
                    }
                }
            }

            KmerSummary {
                kmer_counts,
                n_singletons: singles,
            }
        })
        .collect();
    println!(" done, time: {:?}", start.elapsed());

    // Create the histograms
    print!("Creating histograms...");
    let start = std::time::Instant::now();
    std::io::stdout().flush().unwrap();
    let mut kmer_counts: FxHashMap<u64, u64> = FxHashMap::default();

    let mut histos: Vec<Vec<u64>> = Vec::with_capacity(args.n);

    // Iterate over the chunks
    for chunk_kmer_count in chunk_kmer_counts {
        for (kmer, kmer_count) in chunk_kmer_count.kmer_counts {
            let count = kmer_counts.entry(kmer).or_insert(0);
            *count += kmer_count;
        }

        let mut histo = count_histogram(&kmer_counts, args.histo_max);

        // Add the number of singletons to the histogram at index 1
        histo[1] += chunk_kmer_count.n_singletons;

        histos.push(histo.clone());
    }
    println!(" done, time: {:?}", start.elapsed());

    // Write the histograms to a tab delimited file, with the first column being the count
    // Skip the first row, which is the count of 0. Do not include a header
    print!("Writing histograms to file...");
    std::io::stdout().flush().unwrap();
    let mut file = std::fs::File::create(format!("{}.histo", out_name)).unwrap();
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);
        for histo in histos.iter() {
            line = format!("{}\t{}", line, histo[i]);
        }
        line = format!("{}\n", line);
        file.write_all(line.as_bytes()).unwrap();
    }

    // Write the final histogram to a file, ready for GenomeScope2 etc...
    std::io::stdout().flush().unwrap();
    let mut file = std::fs::File::create(format!("{}.final.histo", out_name)).unwrap();
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);

        line = format!("{}\t{}", line, histos[histos.len()-1][i]);

        line = format!("{}\n", line);
        file.write_all(line.as_bytes()).unwrap();
    }



    println!(" done");
}

#[cfg(test)]
mod test;
