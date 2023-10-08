use bio::io::fasta;
use clap::Parser;
use rand::prelude::SliceRandom;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

mod kmer;
mod pcr;

pub enum ParameterValue {
    Int(u32),
    Str(String),
    Float(f64),
}

fn is_valid_nucleotide(c: char) -> bool {
    match c {
        'A' => true,
        'C' => true,
        'G' => true,
        'T' => true,
        'R' => true, // A or G
        'Y' => true, // C or T
        'S' => true, // C or G
        'W' => true, // A or T
        'K' => true, // G or T
        'M' => true, // A or C
        'B' => true, // C or G or T
        'D' => true, // A or G or T
        'H' => true, // A or C or T
        'V' => true, // A or C or G
        'N' => true, // A or C or G or T
        _ => false,
    }
}

pub fn parse_pcr_string(pcr_string: &str) -> Result<pcr::PCRParams, String> {
    // Split the string on underscores
    let split: Vec<&str> = pcr_string.split('_').collect();

    // Check that there are at least 4 elements
    if split.len() < 4 {
        return Err(format!(
            "Invalid pcr string, there are less than 4 elements separated by underscores: {}",
            pcr_string
        ));
    }

    let forward_seq = split[0].to_uppercase();
    let reverse_seq = split[1].to_uppercase();

    // Check that the forward and reverse primers contain only valid nucleotides
    for c in forward_seq.chars() {
        if !is_valid_nucleotide(c) {
            return Err(format!(
                "Invalid nucleotide {} in forward primer {}",
                c, split[0]
            ));
        }
    }
    for c in reverse_seq.chars() {
        if !is_valid_nucleotide(c) {
            return Err(format!(
                "Invalid nucleotide {} in reverse primer {}",
                c, split[1]
            ));
        }
    }

    // Check that the max-length is an integer
    let max_length: usize = match split[2].parse() {
        Ok(n) => n,
        Err(_) => return Err(format!("Invalid max-length: {}", split[2])),
    };

    let gene_name = split[3].to_string();

    let mut coverage = 3;
    let mut mismatches = 2;
    let mut trim = 15;

    // Loop over additional parameters, which are of the form key=value and are separated by underscores
    for item in split.iter().skip(4) {
        let key_value: Vec<&str> = item.split('=').collect();
        if key_value.len() != 2 {
            return Err(format!("Invalid parameter: {}", item));
        }

        let key = key_value[0];
        let value = key_value[1];

        match key {
            "coverage" => {
                coverage = value
                    .parse()
                    .map_err(|_| format!("Invalid value for {}: {}", key, value))?;
            }
            "mismatches" => {
                mismatches = value
                    .parse()
                    .map_err(|_| format!("Invalid value for {}: {}", key, value))?;
            }
            "trim" => {
                trim = value
                    .parse()
                    .map_err(|_| format!("Invalid value for {}: {}", key, value))?;
            }
            _ => {
                return Err(format!("Unexpected parameter: {}", key));
            }
        }
    }

    Ok(pcr::PCRParams {
        forward_seq,
        reverse_seq,
        max_length,
        gene_name,
        coverage,
        mismatches,
        trim,
    })
}

/// A collection of kmer counting and analysis tools
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// k-mer length
    #[arg(short, default_value_t = 21)]
    k: usize,

    /// Maximum value for histogram
    #[arg(long, default_value_t = 10000)]
    histo_max: u64,

    /// Number of chunks to divide the data into
    #[arg(short, default_value_t = 100)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short, long, default_value_t = 0)]
    max_reads: u64,

    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Name of the sample, could be species and sample ID eg Nanomia-bijuga-YPMIZ035039
    /// Will be used as the prefix for output files and some outout products
    #[arg(short, long, default_value_t = String::from("sample") )]
    sample: String,

    /// Directory for output files, defaults to current directory
    #[arg(short, long, default_value_t = String::from("./") )]
    outdir: String,

    /// Input files, fastq. If no files are specified, data will be read
    /// from stdin. This can be used to uncompress a gz file and send them
    /// to sharkmer.
    #[arg()]
    input: Option<Vec<String>>,

    /// Optional primer pairs for in silico PCR (sPCR). The format is:
    /// --pcr "forward_reverse_max-length_name"
    /// Where:
    ///   forward is the forward primer sequence in 5' to 3' orientation
    ///   reverse is the reverse primer sequence in 5' to 3' orientation
    ///     along the opposite strand as the forward primer, so that the
    ///     primers are in the same orientation that you would use in an
    ///     actual in vitro PCR reaction.
    ///   max-length is the maximum length of the PCR product, including
    ///     the primers.
    ///   name is a unique name for the primer pair or amplified gene
    ///     region. This will be used to specify amplified regions in
    ///     the output fasta file.
    /// More than one primer pair can be specified, for example:
    /// --pcr "forward1_reverse1_1000_name1" --pcr "forward2_reverse2_2000_name2"
    #[arg(short = 'p', long)]
    pcr: Vec<String>,

    /// Minimum coverage for kmer to be included in sPCR
    #[arg(short, long, default_value_t = 3)]
    coverage: u64,

    /// Verbosity
    #[arg(long, default_value_t = 0)]
    verbosity: usize,
}
fn main() {
    let start_run = std::time::Instant::now();

    // Ingest command line arguments
    let args = Args::parse();

    // Print the program name and version
    println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    // Print the arguments
    println!("{:?}", args);

    // Parse the outdir path and sample, create directories if necessary
    // output directory is outdir/sample
    let mut path = PathBuf::from(&args.outdir);
    path.push(&args.sample);

    // Create the output directory if it does not exist
    let directory = format!("{}/", path.to_str().unwrap());
    println!("Full output directory: {}", directory);
    std::fs::create_dir_all(&directory).unwrap();

    let k = args.k;

    // Check that the arguments are valid
    assert!(
        k < 32,
        "k must be less than 32 due to use of 64 bit integers to encode kmers"
    );
    assert!(k > 0, "k must be greater than 0");
    assert!(k % 2 == 1, "k must be odd");
    assert!(args.histo_max > 0, "histo_max must be greater than 0");
    assert!(args.n > 0, "n must be greater than 0");

    // Create an empty data frame for pcr runs
    let mut pcr_runs: Vec<pcr::PCRParams> = Vec::new();

    // Loop over the pcr strings, check that they are valid, and add each to the pcr_runs vector
    for pcr_string in args.pcr.iter() {
        let parsed_pcr = parse_pcr_string(pcr_string);
        match parsed_pcr {
            Ok(pcr_params) => {
                pcr_runs.push(pcr_params);
            }
            Err(err) => {
                panic!("Error parsing pcr string: {}", err);
            }
        }
    }

    // Check that there are no duplicate gene names
    let mut gene_names: Vec<String> = Vec::new();
    for pcr_params in pcr_runs.iter() {
        if gene_names.contains(&pcr_params.gene_name) {
            panic!("Duplicate gene name: {}", pcr_params.gene_name);
        } else {
            gene_names.push(pcr_params.gene_name.clone());
        }
    }

    // Set the number of threads for Rayon to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // Ingest the fastq files
    let start = std::time::Instant::now();
    print!("Ingesting reads...");
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
                        let ints = kmer::seq_to_ints(&line);
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
                    let ints = kmer::seq_to_ints(&line);
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

    let chunk_kmer_counts: Vec<kmer::KmerSummary> = (0..n)
        .into_par_iter()
        .map(|i| {
            let start = i * chunk_size;
            let end = (i + 1) * chunk_size;
            let mut kmer_counts: FxHashMap<u64, u64> = FxHashMap::default();

            for read in reads[start..end].iter() {
                let kmers = kmer::ints_to_kmers(read, &k);
                for kmer in kmers {
                    let count = kmer_counts.entry(kmer).or_insert(0);
                    *count += 1;
                }
            }

            kmer::KmerSummary {
                kmer_counts,
                n_singletons: 0,
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

        let mut histo = kmer::count_histogram(&kmer_counts, &args.histo_max);

        // Add the number of singletons to the histogram at index 1
        histo[1] += chunk_kmer_count.n_singletons;

        histos.push(histo.clone());
    }
    println!(" done, time: {:?}", start.elapsed());

    // Write the histograms to a tab delimited file, with the first column being the count
    // Skip the first row, which is the count of 0. Do not include a header
    print!("Writing histograms to file...");
    std::io::stdout().flush().unwrap();
    let mut file = std::fs::File::create(format!("{}{}.histo", directory, args.sample)).unwrap();
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
    let mut file =
        std::fs::File::create(format!("{}{}.final.histo", directory, args.sample)).unwrap();
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
    let mut file_stats =
        std::fs::File::create(format!("{}{}.stats", directory, args.sample)).unwrap();
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

    if !pcr_runs.is_empty() {
        println!("Running in silico PCR...");

        // Remove kmer_counts entries with less than coverage
        print!(
            "Removing kmers with coverage less than {}...",
            args.coverage
        );
        std::io::stdout().flush().unwrap();
        let mut kmer_counts_filtered: FxHashMap<u64, u64> = FxHashMap::default();
        let mut count_filtered_total: u64 = 0;
        let mut count_raw_total: u64 = 0;

        for (&kmer, &count) in &kmer_counts {
            if count >= args.coverage {
                kmer_counts_filtered.insert(kmer, count);
                count_filtered_total += count;
            }
            count_raw_total += count;
        }

        println!(
            "The total kmer count went from {} to {}",
            count_raw_total, count_filtered_total
        );
        println!(
            "The number of unique kmers went from {} to {}",
            kmer_counts.len(),
            kmer_counts_filtered.len()
        );

        for pcr_params in pcr_runs.iter() {
            let fasta = pcr::do_pcr(
                &kmer_counts_filtered,
                &{ args.k },
                &args.sample,
                args.verbosity,
                pcr_params,
            );

            println!("There are {} subassemblies", fasta.len());
            if !fasta.is_empty() {
                let fasta_path = format!(
                    "{}{}_{}.fasta",
                    directory, args.sample, pcr_params.gene_name
                );
                let mut fasta_writer =
                    fasta::Writer::new(std::fs::File::create(fasta_path).unwrap());
                for record in fasta {
                    fasta_writer.write_record(&record).unwrap();
                }
            }
        }

        println!("Done running in silico PCR");
    }

    println!("Total run time: {:?}", start_run.elapsed());
}

#[cfg(test)]
mod test;
