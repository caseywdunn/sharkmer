use bloom::{BloomFilter, ASMS};
use clap::Parser;
use rand::prelude::SliceRandom;
use rayon::prelude::*;
use dashmap::DashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;

// Create a structure with a hashmap for kmer counts and a u64 for the number of singleton kmers
struct KmerSummary {
    kmer_counts: DashMap<u64, u64>,
    n_singletons: u64,
}

// Convert read to integer encoded subreads, split on N in original sequence
fn seq_to_ints(seq: &str) -> Vec<Vec<u8>> {
    let mut result: Vec<Vec<u8>> = Vec::new();
    let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4);
    let mut frame: u8 = 0;
    let mut position: usize = 0; // position in the sequence. not including Ns
    for c in seq.chars() {
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            'N' => {
                if !ints.is_empty() {
                    result.push(ints);
                    ints = Vec::with_capacity(seq.len() / 4);
                    frame = 0; // Reset frame before starting a new subread
                }
                continue;
            }
            _ => 5,
        };
        if base > 3 {
            break;
        }
        frame = (frame << 2) | base;
        if ((position + 1) % 4 == 0) & (position > 0) {
            ints.push(frame);
            frame = 0; // Reset frame after pushing to the vector
        }
        position += 1;
    }
    if !ints.is_empty() || result.is_empty() {
        // Don't miss the last part
        result.push(ints);
    }
    result
}

fn ints_to_kmers(ints: &Vec<u8>, k: u8) -> Vec<u64> {
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

fn count_histogram(kmer_counts: &DashMap<u64, u64>, histo_max: u64) -> Vec<u64> {
    // Create a histogram of counts
    let mut histo: Vec<u64> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for entry in kmer_counts.iter() {
        let count = *entry.value();
        if count <= histo_max {
            histo[count as usize] += 1;
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

    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Disable bloom filter, which runs faster but uses more memory
    #[arg(short, long, default_value_t = false)]
    disable_bloom: bool,

    /// Directory and filename prefix for analysis output, for example out_dir/Nanomia-bijuga
    #[arg(short, long, default_value_t = String::from("sample") )]
    output: String,

    /// Input files, fastq. If no files are specified, data will be read
    /// from stdin. This can be used to uncompress a gz file and send them
    /// to sharkmer.
    #[arg()]
    input: Option<Vec<String>>,
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

    // Set the number of threads for Rayon to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // Ingest the fastq files
    let start = std::time::Instant::now();
    print!("Ingested reads...");
    std::io::stdout().flush().unwrap();
    let mut reads: Vec<Vec<u8>> = Vec::new();
    let mut n_reads_read = 0;
    let mut n_bases_read = 0;

    match &args.input {
        Some(input_files) => {
            // read from one or more files
            'processing_files: for file_name in input_files.iter() {
                let mut line_n = 0;
                // Open the file for buffered reading
                let file_path = Path::new(&file_name);
                let file = std::fs::File::open(file_path).unwrap();
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
        }
        None => {
            // read from stdin
            let stdin = std::io::stdin();
            // Lock stdin for exclusive access
            let handle = stdin.lock();

            let mut line_n = 0;

            // Create a buffer for reading lines
            for line in handle.lines() {
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
                    break;
                }
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

    // Preallcoate the bloom filter
    let mut multi_bloom: Option<BloomFilter> = None;

    if !args.disable_bloom {
        // Find kmers that occur multiple times with bloom filter
        let start = std::time::Instant::now();
        println!("Building bloom filter to identify kmers that occur more than once...");
        // Resize the bloom filter
        
        multi_bloom = Some(BloomFilter::with_rate(0.01, 4_294_967_295));
        let mut n_kmers: u64 = 0;
        let mut n_multi_kmers: u64 = 0;
        if let Some(multi_bloom_ref) = multi_bloom.as_mut(){
            let mut pre_bloom = BloomFilter::with_rate(0.01, 4_294_967_295);
            // Print a progress bar, 2% per character
            println!("--------------------------------------------------- 100%");
            let reads_per_2_percent = (reads.len() / 50) as u64;

            for (n_reads, read) in (0_u64..).zip(reads.iter()) {
                if n_reads % reads_per_2_percent == 0 {
                    print!("#");
                    std::io::stdout().flush().unwrap();
                }

                let kmers = ints_to_kmers(read, args.k as u8);
                for kmer in kmers {
                    n_kmers += 1;
                    if pre_bloom.contains(&kmer) {
                        multi_bloom_ref.insert(&kmer);
                        n_multi_kmers += 1;
                    } else {
                        pre_bloom.insert(&kmer);
                    }
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
    } else {
        println!("Skipping bloom filter");
    }

    // Randomize the order of the reads in place
    print!("Randomizing read order...");
    std::io::stdout().flush().unwrap();
    let mut rng = rand::thread_rng();
    reads.shuffle(&mut rng);
    println!(" done");

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
            let kmer_counts: DashMap<u64, u64> = DashMap::new();
            let mut singles: u64 = 0;

            if args.disable_bloom {
                for read in reads[start..end].iter() {
                    let kmers = ints_to_kmers(read, args.k as u8);
                    for kmer in kmers {
                        kmer_counts.entry(kmer).and_modify(|count| *count += 1).or_insert(1);
                    }
                }
            } else {
                for read in reads[start..end].iter() {
                    let kmers = ints_to_kmers(read, args.k as u8);
                    for kmer in kmers {
                        if multi_bloom.as_ref().unwrap().contains(&kmer) {
                            kmer_counts.entry(kmer).and_modify(|count| *count += 1).or_insert(1);
                        } else {
                            singles += 1;
                        }
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
    let kmer_counts: DashMap<u64, u64> = DashMap::new();

    let mut histos: Vec<Vec<u64>> = Vec::with_capacity(args.n);

    // Iterate over the chunks
    for chunk_kmer_count in chunk_kmer_counts {
    
        let keys: Vec<u64> = chunk_kmer_count.kmer_counts.iter().map(|x| *x.key()).collect();

        (0..args.threads).into_par_iter().for_each(|i| {
            let start = i * keys.len() / args.threads;
            let end = (i + 1) * keys.len() / args.threads;
            for kmer in keys[start..end].iter() {
                if let Some(temp_ref) = chunk_kmer_count.kmer_counts.get(&kmer) {
                    let value = *temp_ref;
                    kmer_counts.entry(*kmer).and_modify(|count| *count += value).or_insert(value);
                } else {
                    panic!("Key not found");
                }
            }
        });

        // Take care of any remaining keys
        for kmer in keys[(args.threads * (keys.len() / args.threads))..].iter() {
            if let Some(temp_ref) = chunk_kmer_count.kmer_counts.get(&kmer) {
                let value = *temp_ref;
                kmer_counts.entry(*kmer).and_modify(|count| *count += value).or_insert(value);
            } else {
                panic!("Key not found");
            }
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

    println!(" done");

    // Write the final histogram to a file, ready for GenomeScope2 etc...
    print!("Writing final histogram to file...");
    let mut n_kmers: u64 = 0;
    let n_singleton_kmers: u64 = histos[histos.len() - 1][1];
    std::io::stdout().flush().unwrap();
    let mut file = std::fs::File::create(format!("{}.final.histo", out_name)).unwrap();
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);

        line = format!("{}\t{}", line, histos[histos.len() - 1][i]);
        n_kmers += histos[histos.len() - 1][i];

        line = format!("{}\n", line);
        file.write_all(line.as_bytes()).unwrap();
    }

    println!(" done");

    print!("Writing stats to file...");
    std::io::stdout().flush().unwrap();
    let mut file_stats = std::fs::File::create(format!("{}.stats", out_name)).unwrap();
    let mut line = format!("arguments\t{:?}\n", args);
    line = format!("{}kmer_length\t{}\n", line, args.k);
    line = format!("{}n_reads_read\t{}\n", line, n_reads_read);
    line = format!("{}n_bases_read\t{}\n", line, n_bases_read);
    line = format!("{}n_subreads_ingested\t{}\n", line, reads.len());
    line = format!("{}n_bases_ingested\t{}\n", line, n_bases_ingested);
    line = format!("{}n_kmers\t{}\n", line, n_kmers);
    line = format!("{}n_multi_kmers\t{}\n", line, n_kmers - n_singleton_kmers);
    line = format!("{}n_singleton_kmers\t{}\n", line, n_singleton_kmers);

    file_stats.write_all(line.as_bytes()).unwrap();
    println!(" done");
}

#[cfg(test)]
mod test;
