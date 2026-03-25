use anyhow::{Context, Result, bail, ensure};
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use std::io::BufRead;
use std::io::IsTerminal;
use std::io::Write;

use crate::cli::Args;
use crate::format::{format_bytes, format_count, format_duration};
use crate::kmer;
use crate::kmer::{Chunk, KmerCounts};

pub(crate) const FASTA_LINE_WIDTH: usize = 80;
pub(crate) const N_READS_PER_BATCH: u64 = 1000;

/// Result of querying ENA for an accession.
pub(crate) struct EnaResult {
    pub(crate) urls: Vec<String>,
    pub(crate) scientific_name: Option<String>,
}

/// Query ENA filereport API to get FASTQ URLs and metadata for an accession.
/// Returns URLs ordered by mate number (R1 first, then R2 if paired-end),
/// plus the scientific_name if available.
pub(crate) fn get_ena_fastq_urls(accession: &str) -> Result<EnaResult> {
    let url = format!(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=run_accession,fastq_ftp,scientific_name",
        accession
    );

    info!("Querying ENA for accession {}...", accession);
    let response = ureq::get(&url)
        .call()
        .with_context(|| format!("Failed to query ENA API for accession {}", accession))?;

    let body = response
        .into_string()
        .context("Failed to read ENA API response")?;

    // Response is TSV: header row, then data rows
    // Fields: run_accession\tfastq_ftp\tscientific_name
    let lines: Vec<&str> = body.lines().collect();
    ensure!(
        lines.len() >= 2,
        "ENA returned no results for accession '{}'. Check that the accession is valid.",
        accession
    );

    // Parse header to find column indices
    let headers: Vec<&str> = lines[0].split('\t').collect();
    let ftp_idx = headers
        .iter()
        .position(|h| *h == "fastq_ftp")
        .context("ENA response missing fastq_ftp column")?;
    let sci_name_idx = headers.iter().position(|h| *h == "scientific_name");

    let data_line = lines[1];
    let fields: Vec<&str> = data_line.split('\t').collect();
    ensure!(
        fields.len() > ftp_idx && !fields[ftp_idx].is_empty(),
        "ENA returned no FASTQ URLs for accession '{}'. The run may not have public FASTQ files.",
        accession
    );

    let urls: Vec<String> = fields[ftp_idx]
        .split(';')
        .map(|u| {
            if u.starts_with("ftp://") || u.starts_with("http") {
                u.to_string()
            } else {
                format!("http://{}", u)
            }
        })
        .collect();

    let scientific_name = sci_name_idx
        .and_then(|idx| fields.get(idx))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string());

    info!(
        "Found {} FASTQ file(s) for {}: {}",
        urls.len(),
        accession,
        urls.join(", ")
    );
    if let Some(ref name) = scientific_name {
        info!("Scientific name: {}", name);
    }

    Ok(EnaResult {
        urls,
        scientific_name,
    })
}

/// Warn if an output file already exists (it will be overwritten).
pub(crate) fn warn_if_exists(path: &str) {
    if std::path::Path::new(path).exists() {
        warn!("Overwriting existing file {}", path);
    }
}

/// Write a FASTA record with sequence lines wrapped at FASTA_LINE_WIDTH characters.
pub(crate) fn write_fasta_record(
    writer: &mut impl Write,
    record: &bio::io::fasta::Record,
) -> std::io::Result<()> {
    write!(writer, ">{}", record.id())?;
    if let Some(desc) = record.desc() {
        write!(writer, " {}", desc)?;
    }
    writeln!(writer)?;
    for chunk in record.seq().chunks(FASTA_LINE_WIDTH) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Validate a FASTQ record (4 lines). Returns an error if the record is malformed.
fn validate_fastq_record(
    header: &str,
    _sequence: &str,
    separator: &str,
    quality: &str,
    sequence_len: usize,
    record_num: u64,
) -> Result<()> {
    if header.starts_with('>') {
        bail!(
            "Input appears to be FASTA format, not FASTQ (record {} starts with '>'). \
             sharkmer requires FASTQ input with quality scores.",
            record_num + 1,
        );
    }
    ensure!(
        header.starts_with('@'),
        "FASTQ record {} has invalid header (expected '@', got '{}'): {}",
        record_num + 1,
        header.chars().next().unwrap_or(' '),
        header,
    );
    ensure!(
        separator.starts_with('+'),
        "FASTQ record {} has invalid separator line (expected '+', got '{}'): {}",
        record_num + 1,
        separator.chars().next().unwrap_or(' '),
        separator,
    );
    ensure!(
        quality.len() == sequence_len,
        "FASTQ record {} has mismatched sequence ({}) and quality ({}) lengths",
        record_num + 1,
        sequence_len,
        quality.len(),
    );
    Ok(())
}

/// Mutable state for FASTQ reading, shared across multiple input sources.
pub(crate) struct FastqReadState {
    pub(crate) chunks: Vec<kmer::Chunk>,
    pub(crate) chunk_index: usize,
    pub(crate) reads: Vec<kmer::Read>,
    pub(crate) n_reads_read: u64,
    pub(crate) n_bases_read: u64,
}

/// Read FASTQ records from a buffered reader, ingesting sequences into chunks.
///
/// Reads 4 lines at a time (header, sequence, separator, quality) and validates
/// the format. Returns true if max_reads was reached.
fn read_fastq<R: BufRead>(
    reader: R,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source_name: &str,
    progress: &ProgressBar,
) -> Result<bool> {
    let mut lines = reader.lines();
    let n_chunks = state.chunks.len();

    while let Some(header_result) = lines.next() {
        // Read 4 lines for one FASTQ record
        let header =
            header_result.with_context(|| format!("Failed to read header from {}", source_name))?;

        let sequence = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read sequence from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing sequence line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let separator = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read separator from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing separator line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let quality = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read quality from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing quality line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        // Validate the record
        let should_validate = state.n_reads_read == 0
            || (validate_every > 0 && state.n_reads_read % validate_every == 0);
        if should_validate {
            validate_fastq_record(
                &header,
                &sequence,
                &separator,
                &quality,
                sequence.len(),
                state.n_reads_read,
            )?;
        }

        // Process the sequence
        state.n_bases_read += sequence.len() as u64;
        let new_reads = kmer::seq_to_reads(&sequence)?;
        state.reads.extend(new_reads);
        state.n_reads_read += 1;

        // If we have read enough reads, ingest them into current chunk
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            state.chunks[state.chunk_index].ingest_reads(&state.reads)?;
            state.chunk_index += 1;
            if state.chunk_index == n_chunks {
                state.chunk_index = 0;
            }
            state.reads.clear();
            progress.set_position(state.n_reads_read);
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

/// Ingest FASTQ reads from all input sources (ENA, files, or stdin).
/// Returns the read state with populated chunks, and summary statistics.
pub(crate) fn ingest_reads(
    args: &Args,
    k: usize,
    mut cached_ena_result: Option<EnaResult>,
    show_progress: bool,
) -> Result<(FastqReadState, u64, u64, u64)> {
    let start = std::time::Instant::now();
    info!("Ingesting reads...");

    // When chunks == 0, allocate 1 internal chunk for kmer storage (needed for sPCR)
    let n_chunks = if args.chunks == 0 { 1 } else { args.chunks };
    let mut chunks: Vec<kmer::Chunk> = Vec::new();
    for _ in 0..n_chunks {
        chunks.push(Chunk::new(&k));
    }

    let mut state = FastqReadState {
        chunks,
        chunk_index: 0,
        reads: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
    };

    let max_reads = args.max_reads.unwrap_or(0);

    // Create progress indicator: bar with ETA if max_reads is known, spinner otherwise.
    let progress = if show_progress && max_reads > 0 {
        let pb = ProgressBar::new(max_reads);
        pb.set_style(
            ProgressStyle::with_template(
                "Ingesting reads {bar:30} {human_pos}/{human_len} [{per_sec}] [ETA {eta}]",
            )
            .expect("valid progress template"),
        );
        pb
    } else if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("Ingesting reads {spinner} {human_pos} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    if let Some(accession) = &args.ena {
        // Stream reads from ENA (use cached result if available)
        let ena_result = match cached_ena_result.take() {
            Some(r) => r,
            None => get_ena_fastq_urls(accession)?,
        };
        for url in &ena_result.urls {
            info!("Streaming from {}...", url);
            let response = ureq::get(url)
                .call()
                .with_context(|| format!("Failed to download {}", url))?;
            let reader =
                std::io::BufReader::new(flate2::read::GzDecoder::new(response.into_reader()));
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                url,
                &progress,
            )?;
            if reached_max {
                break;
            }
        }
    } else if let Some(input_files) = &args.input {
        // Read from one or more files
        for file_path in input_files.iter() {
            let file = std::fs::File::open(file_path)
                .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

            let file_name = file_path.to_string_lossy();
            let has_gz_ext = file_name.ends_with(".gz") || file_name.ends_with(".gzip");

            // Detect gzip magic bytes for files without .gz extension
            let use_gzip = if has_gz_ext {
                true
            } else {
                let mut buf_reader = std::io::BufReader::new(file);
                let magic = buf_reader.fill_buf().context("Failed to peek at file")?;
                let is_gzip = magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
                if is_gzip {
                    warn!(
                        "File '{}' appears to be gzipped (magic bytes detected) but lacks a .gz extension. Reading as gzipped.",
                        file_path.display()
                    );
                }
                // We need to re-open the file since BufReader consumed ownership
                drop(buf_reader);
                is_gzip
            };

            let file = std::fs::File::open(file_path)
                .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
            let reader: Box<dyn BufRead> = if use_gzip {
                Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
            } else {
                Box::new(std::io::BufReader::new(file))
            };
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                &file_name,
                &progress,
            )?;
            if reached_max {
                break;
            }
        }
    } else {
        // Read from stdin
        let stdin = std::io::stdin();
        ensure!(
            !stdin.is_terminal(),
            "No input files specified and stdin is a terminal.\n\
             Provide FASTQ files as arguments, use --ena, or pipe data via stdin.\n\
             Example: sharkmer -s sample -k 21 reads.fastq\n\
             Example: sharkmer -s sample --ena SRR5324768\n\
             Example: zcat reads.fastq.gz | sharkmer -s sample -k 21"
        );
        let handle = stdin.lock();

        read_fastq(
            handle,
            &mut state,
            max_reads,
            args.validate_every,
            "stdin",
            &progress,
        )?;
    }

    progress.finish_and_clear();

    // Ingest any remaining reads
    state.chunks[state.chunk_index].ingest_reads(&state.reads)?;
    state.reads.clear();

    let mut n_reads_ingested: u64 = 0;
    let mut n_bases_ingested: u64 = 0;
    let mut n_kmers_ingested: u64 = 0;
    for chunk in state.chunks.iter() {
        n_reads_ingested += chunk.get_n_reads();
        n_bases_ingested += chunk.get_n_bases();
        n_kmers_ingested += chunk.get_n_kmers();
    }

    info!(
        "Read {} reads, {} bases",
        format_count(state.n_reads_read),
        format_bytes(state.n_bases_read)
    );
    info!(
        "Ingested {} subreads, {} bases, {} kmers",
        format_count(n_reads_ingested),
        format_bytes(n_bases_ingested),
        format_count(n_kmers_ingested)
    );

    // Warn if very few reads for sPCR
    let has_pcr = !args.pcr_panel.is_empty()
        || !args.pcr_panel_file.is_empty()
        || !args.pcr_primers.is_empty();
    if has_pcr && state.n_reads_read < 10_000 {
        warn!(
            "Only {} reads ingested. sPCR typically needs many more reads to produce results.",
            state.n_reads_read
        );
    }
    info!("Time to ingest reads: {}", format_duration(start.elapsed()));

    if n_reads_ingested == 0 {
        bail!("No reads were ingested. Check that input files contain valid FASTQ records.");
    }

    Ok((state, n_reads_ingested, n_bases_ingested, n_kmers_ingested))
}

/// Consolidate chunks into a single KmerCounts table, optionally writing histograms.
/// Returns (kmer_counts, n_singleton_kmers).
pub(crate) fn consolidate_and_histogram(
    state: &mut FastqReadState,
    args: &Args,
    k: usize,
    sample: &str,
    directory: &str,
    n_kmers_ingested: u64,
    show_progress: bool,
) -> Result<(KmerCounts, Option<u64>)> {
    let spinner_msg = if args.chunks > 0 {
        "Consolidating kmer counts..."
    } else {
        "Merging kmer counts..."
    };
    info!("{}", spinner_msg);
    let spinner_style =
        ProgressStyle::with_template("{spinner:.cyan} {msg}").expect("valid spinner template");
    let consolidate_spinner = if show_progress {
        let sp = ProgressBar::new_spinner();
        sp.set_style(spinner_style);
        sp.set_message(spinner_msg.to_string());
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        sp
    } else {
        ProgressBar::hidden()
    };
    let start = std::time::Instant::now();

    let mut kmer_counts: KmerCounts = KmerCounts::new(&k);
    let mut n_singleton_kmers: Option<u64> = None;

    let histo_comment = format!(
        "# sharkmer {} k={} chunks={}",
        env!("CARGO_PKG_VERSION"),
        args.k,
        args.chunks
    );

    if args.chunks > 0 {
        // Incremental histogram mode: update histogram as chunks are merged
        let mut histos: Vec<kmer::Histogram> = Vec::with_capacity(args.chunks);
        let mut running_histo = kmer::Histogram::new(&args.histo_max);

        for chunk in state.chunks.drain(..) {
            kmer_counts.extend_with_histogram(chunk.get_kmer_counts(), &mut running_histo)?;
            drop(chunk);

            histos.push(running_histo.clone());
        }
        consolidate_spinner.finish_and_clear();
        info!(
            "Chunks consolidated, time: {}",
            format_duration(start.elapsed())
        );

        let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
        info!(
            "{} unique kmers with a total count of {} were found",
            kmer_counts.get_n_unique_kmers(),
            n_hashed_kmers
        );

        ensure!(
            n_hashed_kmers == n_kmers_ingested,
            "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
            n_hashed_kmers,
            n_kmers_ingested,
        );

        // Write incremental histograms with header
        info!("Writing histograms to file...");
        let histo_path = format!("{}{}.histo", directory, sample);
        warn_if_exists(&histo_path);
        let mut file =
            std::fs::File::create(&histo_path).context("Failed to create histogram file")?;

        // Comment line with version and parameters
        writeln!(file, "{}", histo_comment).context("Failed to write histogram comment")?;

        // Header row
        let mut header = "count".to_string();
        for i in 1..=histos.len() {
            header = format!("{}\tchunk_{}", header, i);
        }
        writeln!(file, "{}", header).context("Failed to write histogram header")?;

        // Data rows
        for i in 1..args.histo_max as usize + 2 {
            let mut line = format!("{}", i);
            for histo in histos.iter() {
                let histo_vec = kmer::Histogram::get_vector(histo)?;
                line = format!("{}\t{}", line, histo_vec[i]);
            }
            writeln!(file, "{}", line).context("Failed to write histogram data")?;
        }

        // Write the final histogram with header
        info!("Writing final histogram to file...");
        let last_histo = histos.last().context("No histograms were produced")?;
        let last_histo_vec = kmer::Histogram::get_vector(last_histo)?;

        let final_histo_path = format!("{}{}.final.histo", directory, sample);
        warn_if_exists(&final_histo_path);
        let mut file = std::fs::File::create(&final_histo_path)
            .context("Failed to create final histogram file")?;

        writeln!(file, "{}", histo_comment).context("Failed to write final histogram comment")?;
        writeln!(file, "count\tfrequency").context("Failed to write final histogram header")?;

        for (i, value) in last_histo_vec
            .iter()
            .enumerate()
            .skip(1)
            .take(args.histo_max as usize + 1)
        {
            writeln!(file, "{}\t{}", i, value).context("Failed to write final histogram data")?;
        }

        let singleton_count = *last_histo_vec
            .get(1)
            .context("Histogram vector too short to contain singleton count")?;
        n_singleton_kmers = Some(singleton_count);

        // Warn if singleton rate is very high (>95% of unique kmers)
        let n_unique = last_histo.get_n_unique_kmers();
        if n_unique > 0 {
            let singleton_rate = singleton_count as f64 / n_unique as f64;
            if singleton_rate > 0.95 {
                warn!(
                    "Very high singleton rate ({:.1}%). This may indicate very low coverage \
                     or contamination. sPCR results may be unreliable.",
                    singleton_rate * 100.0
                );
            }
        }

        let n_unique_kmers_histo: u64 = last_histo.get_n_unique_kmers();
        let n_kmers_histo: u64 = last_histo.get_n_kmers();

        debug!("{} unique kmers in histogram", n_unique_kmers_histo);
        debug!("{} kmers in histogram", n_kmers_histo);

        ensure!(
            n_kmers_histo == n_kmers_ingested,
            "The total count of kmers in the histogram ({}) does not equal the total expected count of kmers ({})",
            n_kmers_histo,
            n_kmers_ingested,
        );

        ensure!(
            n_unique_kmers_histo == kmer_counts.get_n_unique_kmers(),
            "The total count of unique kmers in the histogram ({}) does not equal the total count of hashed kmers ({})",
            n_unique_kmers_histo,
            kmer_counts.get_n_unique_kmers(),
        );
    } else {
        // No histogram mode: merge all chunks into a single kmer count table
        for chunk in state.chunks.drain(..) {
            kmer_counts.extend(chunk.get_kmer_counts())?;
            drop(chunk);
        }
        consolidate_spinner.finish_and_clear();
        info!(
            "Kmer counts merged, time: {}",
            format_duration(start.elapsed())
        );

        let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
        info!(
            "{} unique kmers with a total count of {} were found",
            kmer_counts.get_n_unique_kmers(),
            n_hashed_kmers
        );

        ensure!(
            n_hashed_kmers == n_kmers_ingested,
            "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
            n_hashed_kmers,
            n_kmers_ingested,
        );
    }

    Ok((kmer_counts, n_singleton_kmers))
}
