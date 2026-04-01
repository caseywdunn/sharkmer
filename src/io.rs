use anyhow::{Context, Result, bail, ensure};
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use std::io::BufRead;
use std::io::IsTerminal;
use std::io::Write;
use std::path::PathBuf;

use crate::cli::Args;
use crate::format::{format_bytes, format_count, format_duration};
use crate::kmer;
use crate::kmer::{Chunk, KmerCounts};

pub(crate) const FASTA_LINE_WIDTH: usize = 80;
pub(crate) const N_READS_PER_BATCH: u64 = 1000;

/// Describes how to re-acquire read data for Pass 2 (read threading).
pub(crate) enum ReadSourcePlan {
    /// Local files that can be reopened
    LocalFiles(Vec<PathBuf>),
    /// Remote files cached locally (paths to cached files)
    CachedRemote(Vec<PathBuf>),
    /// Remote files that must be re-downloaded (URLs)
    UncachedRemote(Vec<String>),
    /// Cannot re-read (stdin or --no-read-threading)
    Unavailable,
}

/// Plan for Pass 2 re-reading of FASTQ data.
pub(crate) struct ReadPlan {
    pub(crate) source: ReadSourcePlan,
    pub(crate) paired: bool,
    pub(crate) max_reads: u64,
}

/// A retained read for Pass 2 threading.
pub(crate) struct ReadRecord {
    /// Raw sequence (ASCII ACGT)
    pub(crate) sequence: String,
    /// Global read index (0-based, in order of ingestion)
    pub(crate) index: u64,
    /// Mate designation for paired-end reads
    pub(crate) mate: Mate,
}

/// Mate designation for paired-end reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mate {
    Unpaired,
    R1,
    R2,
}

/// A read retained during Pass 1 because it matched a primer Oligo.
#[allow(dead_code)]
pub(crate) struct RetainedRead {
    /// Raw sequence (ASCII ACGT)
    pub(crate) sequence: String,
    /// Index of the gene whose primer matched (into pcr_runs array)
    pub(crate) gene_index: usize,
}

/// Collects reads retained during Pass 1 primer Oligo matching.
pub(crate) struct RetainedReads {
    pub(crate) reads: Vec<RetainedRead>,
}

/// Matches kmers against pre-encoded primer Oligos during Pass 1.
/// For each kmer, checks whether high-order bits match any forward
/// primer Oligo or low-order bits match any reverse primer Oligo.
pub(crate) struct OligoMatcher {
    /// Per-gene: (forward Oligos shifted to high bits, mask, rc Oligos, rc mask)
    gene_matchers: Vec<GeneMatcher>,
}

struct GeneMatcher {
    /// Forward oligos shifted to kmer high-order bits
    forward_set: std::collections::HashSet<u64>,
    /// Mask for forward matching (covers high-order oligo bits)
    forward_mask: u64,
    /// Reverse complement oligos (unshifted, match low-order bits)
    rc_set: std::collections::HashSet<u64>,
    /// Mask for RC matching (covers low-order oligo bits)
    rc_mask: u64,
}

impl OligoMatcher {
    /// Build a matcher from pre-encoded primer Oligo sets.
    ///
    /// For matching against canonical kmers (as returned by `kmers_from_ascii`):
    /// - Forward primer Oligos are shifted to high-order bits and matched there
    /// - Both forward and reverse Oligos are also checked in RC form at
    ///   low-order bits, since the canonical kmer may be the RC of the
    ///   original read kmer
    pub(crate) fn new(oligo_sets: &[crate::pcr::PrimerOligoSet], k: usize) -> Self {
        let mut gene_matchers = Vec::new();

        for set in oligo_sets {
            let mut forward_set = std::collections::HashSet::new();
            let mut rc_set = std::collections::HashSet::new();

            // Forward oligos: shift to high-order bits of k-mer
            let fwd_len = set.forward_oligo_length;
            let fwd_shift = if k > fwd_len { 2 * (k - fwd_len) } else { 0 };
            let mut fwd_mask: u64 = 0;
            for _ in 0..(2 * fwd_len) {
                fwd_mask = (fwd_mask << 1) | 1;
            }
            fwd_mask <<= fwd_shift;

            for oligo in &set.forward_oligos {
                forward_set.insert(oligo.kmer << fwd_shift);
            }

            // Reverse oligos also shifted to high bits (primer at start of kmer)
            let rev_len = set.reverse_oligo_length;
            let rev_shift = if k > rev_len { 2 * (k - rev_len) } else { 0 };
            for oligo in &set.reverse_oligos {
                forward_set.insert(oligo.kmer << rev_shift);
            }

            // RC forms of both primer directions at low-order bits:
            // when the canonical kmer is the RC, the primer appears at the end
            let rc_mask_fwd = (1u64 << (2 * fwd_len)) - 1;
            let rc_mask_rev = (1u64 << (2 * rev_len)) - 1;
            // Use the wider mask (both should be same length after trimming,
            // but be safe)
            let rc_mask_val = rc_mask_fwd.max(rc_mask_rev);

            for oligo in &set.forward_oligos {
                rc_set.insert(crate::kmer::revcomp_kmer(&oligo.kmer, &fwd_len));
            }
            for oligo in &set.reverse_oligos {
                rc_set.insert(crate::kmer::revcomp_kmer(&oligo.kmer, &rev_len));
            }

            gene_matchers.push(GeneMatcher {
                forward_set,
                forward_mask: fwd_mask,
                rc_set,
                rc_mask: rc_mask_val,
            });
        }

        OligoMatcher { gene_matchers }
    }

    /// Check if a read contains any primer Oligo by examining its kmers.
    /// Returns the gene index if matched, or None.
    pub(crate) fn check_read(&self, sequence: &str, k: usize) -> Option<usize> {
        let kmers = match crate::kmer::encoding::kmers_from_ascii(sequence, k) {
            Ok(k) => k,
            Err(_) => return None,
        };

        for kmer in &kmers {
            for (gene_idx, matcher) in self.gene_matchers.iter().enumerate() {
                // Check forward orientation (oligo at high bits)
                if matcher.forward_set.contains(&(kmer & matcher.forward_mask)) {
                    return Some(gene_idx);
                }
                // Check RC orientation (oligo at low bits)
                if matcher.rc_set.contains(&(kmer & matcher.rc_mask)) {
                    return Some(gene_idx);
                }
            }
        }

        None
    }
}

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
    pub(crate) seqs: Vec<String>,
    pub(crate) n_reads_read: u64,
    pub(crate) n_bases_read: u64,
    /// Reads retained during Pass 1 because they matched a primer Oligo
    pub(crate) retained_reads: RetainedReads,
}

/// Read FASTQ records from a buffered reader, ingesting sequences into chunks.
///
/// Reads 4 lines at a time (header, sequence, separator, quality) and validates
/// the format. Returns true if max_reads was reached.
#[allow(clippy::too_many_arguments)]
fn read_fastq<R: BufRead>(
    reader: R,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source_name: &str,
    progress: &ProgressBar,
    oligo_matcher: Option<&OligoMatcher>,
    k: usize,
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
        state.seqs.push(sequence);
        state.n_reads_read += 1;

        // If we have read enough reads, ingest them into current chunk
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            // Check for primer Oligo matches before consuming sequences
            if let Some(matcher) = oligo_matcher {
                retain_primer_reads(state, matcher, k);
            }
            for seq in state.seqs.drain(..) {
                state.chunks[state.chunk_index].ingest_seq(&seq)?;
            }
            state.chunk_index += 1;
            if state.chunk_index == n_chunks {
                state.chunk_index = 0;
            }
            progress.set_position(state.n_reads_read);
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

/// Check sequences for primer Oligo matches and retain matching reads.
fn retain_primer_reads(state: &mut FastqReadState, oligo_matcher: &OligoMatcher, k: usize) {
    for seq in &state.seqs {
        if let Some(gene_index) = oligo_matcher.check_read(seq, k) {
            state.retained_reads.reads.push(RetainedRead {
                sequence: seq.clone(),
                gene_index,
            });
        }
    }
}

/// Ingest FASTQ reads from all input sources (ENA, files, or stdin).
/// Returns the read state with populated chunks, summary statistics, and a
/// `ReadPlan` describing how to re-read the same data in Pass 2.
pub(crate) fn ingest_reads(
    args: &Args,
    k: usize,
    mut cached_ena_result: Option<EnaResult>,
    cache_config: Option<&crate::cache::CacheConfig>,
    show_progress: bool,
    oligo_matcher: Option<&OligoMatcher>,
) -> Result<(FastqReadState, u64, u64, u64, ReadPlan)> {
    let start = std::time::Instant::now();
    info!("Ingesting reads...");

    // When chunks == 0, allocate 1 internal chunk for kmer storage (needed for sPCR)
    let n_chunks = if args.chunks == 0 { 1 } else { args.chunks };
    let chunks: Vec<kmer::Chunk> = (0..n_chunks).map(|_| Chunk::new(&k)).collect();

    let mut state = FastqReadState {
        chunks,
        chunk_index: 0,
        seqs: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
        retained_reads: RetainedReads { reads: Vec::new() },
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

    // Track the input source for Pass 2 re-reading
    let mut read_source_plan = ReadSourcePlan::Unavailable;

    if let Some(accession) = &args.ena {
        // Stream reads from ENA (use cached result if available)
        let ena_result = match cached_ena_result.take() {
            Some(r) => r,
            None => get_ena_fastq_urls(accession)?,
        };
        let mut cached_paths: Vec<PathBuf> = Vec::new();
        let mut is_cached = false;
        for url in &ena_result.urls {
            let reader: Box<dyn BufRead> = if let Some(cache) = cache_config {
                // Use cache: lookup or download
                let local_path = match cache.lookup(url, max_reads)? {
                    Some(path) => {
                        info!("Cache hit for {}", url);
                        path
                    }
                    None => {
                        info!("Cache miss for {}, downloading...", url);
                        cache.download_to_cache(url, max_reads)?
                    }
                };
                cached_paths.push(local_path.clone());
                is_cached = true;
                let file = std::fs::File::open(&local_path).with_context(|| {
                    format!("Failed to open cached file: {}", local_path.display())
                })?;
                Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
            } else {
                // No cache: stream directly
                info!("Streaming from {} (no cache)...", url);
                let response = ureq::get(url)
                    .call()
                    .with_context(|| format!("Failed to download {}", url))?;
                Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(
                    response.into_reader(),
                )))
            };
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                url,
                &progress,
                oligo_matcher,
                k,
            )?;
            if reached_max {
                break;
            }
        }
        read_source_plan = if is_cached {
            ReadSourcePlan::CachedRemote(cached_paths)
        } else {
            warn!("Read threading will require re-downloading reads from ENA (no cache)");
            ReadSourcePlan::UncachedRemote(ena_result.urls)
        };
    } else if let Some(input_files) = &args.input {
        if args.paired {
            // Paired-end: alternate reads from R1 and R2
            let reader1 = open_fastq_reader(&input_files[0])?;
            let reader2 = open_fastq_reader(&input_files[1])?;
            let name1 = input_files[0].to_string_lossy().to_string();
            let name2 = input_files[1].to_string_lossy().to_string();
            // Round max_reads up to even for balanced pairs
            let paired_max = if max_reads > 0 && max_reads % 2 != 0 {
                max_reads + 1
            } else {
                max_reads
            };
            read_fastq_paired(
                reader1,
                reader2,
                &mut state,
                paired_max,
                args.validate_every,
                &name1,
                &name2,
                &progress,
                oligo_matcher,
                k,
            )?;
        } else {
            // Read from one or more files sequentially
            for file_path in input_files.iter() {
                let reader = open_fastq_reader(file_path)?;
                let file_name = file_path.to_string_lossy();
                let reached_max = read_fastq(
                    reader,
                    &mut state,
                    max_reads,
                    args.validate_every,
                    &file_name,
                    &progress,
                    oligo_matcher,
                    k,
                )?;
                if reached_max {
                    break;
                }
            }
        }
        read_source_plan = ReadSourcePlan::LocalFiles(input_files.clone());
    } else {
        // Read from stdin (cannot re-read for Pass 2)
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
            oligo_matcher,
            k,
        )?;
        // read_source_plan remains Unavailable for stdin
    }

    progress.finish_and_clear();

    // Check remaining sequences for primer Oligo matches
    if let Some(matcher) = oligo_matcher {
        retain_primer_reads(&mut state, matcher, k);
    }

    // Ingest any remaining sequences
    for seq in state.seqs.drain(..) {
        state.chunks[state.chunk_index].ingest_seq(&seq)?;
    }

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
    if !state.retained_reads.reads.is_empty() {
        info!(
            "Retained {} reads matching primer Oligos for seed evaluation",
            state.retained_reads.reads.len()
        );
    }
    info!("Time to ingest reads: {}", format_duration(start.elapsed()));

    if n_reads_ingested == 0 {
        bail!("No reads were ingested. Check that input files contain valid FASTQ records.");
    }

    let read_plan = ReadPlan {
        source: read_source_plan,
        paired: args.paired,
        max_reads,
    };

    Ok((
        state,
        n_reads_ingested,
        n_bases_ingested,
        n_kmers_ingested,
        read_plan,
    ))
}

/// Open a FASTQ file as a buffered reader, with automatic gzip detection.
fn open_fastq_reader(file_path: &std::path::Path) -> Result<Box<dyn BufRead>> {
    let file_name = file_path.to_string_lossy();
    let has_gz_ext = file_name.ends_with(".gz") || file_name.ends_with(".gzip");

    let use_gzip = if has_gz_ext {
        true
    } else {
        let file = std::fs::File::open(file_path)
            .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
        let mut buf_reader = std::io::BufReader::new(file);
        let magic = buf_reader.fill_buf().context("Failed to peek at file")?;
        let is_gzip = magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
        drop(buf_reader);
        is_gzip
    };

    let file = std::fs::File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
    if use_gzip {
        Ok(Box::new(std::io::BufReader::new(
            flate2::read::GzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(std::io::BufReader::new(file)))
    }
}

/// Read FASTQ sequences from two files in alternating order (R1, R2, R1, R2, ...).
/// Returns true if max_reads was reached.
#[allow(clippy::too_many_arguments)]
fn read_fastq_paired<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source1_name: &str,
    source2_name: &str,
    progress: &ProgressBar,
    oligo_matcher: Option<&OligoMatcher>,
    k: usize,
) -> Result<bool> {
    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();
    let n_chunks = state.chunks.len();

    loop {
        // Read one record from R1
        let r1_done = read_one_fastq_record(&mut lines1, state, validate_every, source1_name)?;
        if r1_done {
            break;
        }
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            if let Some(matcher) = oligo_matcher {
                retain_primer_reads(state, matcher, k);
            }
            for seq in state.seqs.drain(..) {
                state.chunks[state.chunk_index].ingest_seq(&seq)?;
            }
            state.chunk_index = (state.chunk_index + 1) % n_chunks;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }

        // Read one record from R2
        let r2_done = read_one_fastq_record(&mut lines2, state, validate_every, source2_name)?;
        if r2_done {
            break;
        }
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            if let Some(matcher) = oligo_matcher {
                retain_primer_reads(state, matcher, k);
            }
            for seq in state.seqs.drain(..) {
                state.chunks[state.chunk_index].ingest_seq(&seq)?;
            }
            state.chunk_index = (state.chunk_index + 1) % n_chunks;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }
    }

    Ok(false)
}

/// Read a single FASTQ record (4 lines) from a lines iterator.
/// Returns true if EOF was reached (no more records).
fn read_one_fastq_record<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    state: &mut FastqReadState,
    validate_every: u64,
    source_name: &str,
) -> Result<bool> {
    let header = match lines.next() {
        Some(result) => {
            result.with_context(|| format!("Failed to read header from {}", source_name))?
        }
        None => return Ok(true), // EOF
    };

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

    let should_validate =
        state.n_reads_read == 0 || (validate_every > 0 && state.n_reads_read % validate_every == 0);
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

    state.n_bases_read += sequence.len() as u64;
    state.seqs.push(sequence);
    state.n_reads_read += 1;

    Ok(false) // not EOF
}

/// Re-read FASTQ sequences for Pass 2 (read threading).
/// Opens the same sources used in Pass 1 and collects sequences into `ReadRecord`s.
pub(crate) fn reread_sequences(plan: &ReadPlan, show_progress: bool) -> Result<Vec<ReadRecord>> {
    let start = std::time::Instant::now();
    info!("Pass 2: re-reading sequences for read threading...");

    let files: Vec<PathBuf> = match &plan.source {
        ReadSourcePlan::LocalFiles(paths) => paths.clone(),
        ReadSourcePlan::CachedRemote(paths) => paths.clone(),
        ReadSourcePlan::UncachedRemote(urls) => {
            // Re-download to temporary cache paths
            warn!("Re-downloading reads for Pass 2 (use --cache-dir to avoid this)");
            let mut paths = Vec::new();
            for url in urls {
                info!("Downloading {} for Pass 2...", url);
                let response = ureq::get(url)
                    .call()
                    .with_context(|| format!("Failed to download {} for Pass 2", url))?;
                // Write to a temporary file
                let tmp_path = std::env::temp_dir().join(format!(
                    "sharkmer_pass2_{}.fastq.gz",
                    url.len() // simple disambiguator
                ));
                let mut tmp_file = std::fs::File::create(&tmp_path).with_context(|| {
                    format!(
                        "Failed to create temp file for Pass 2: {}",
                        tmp_path.display()
                    )
                })?;
                std::io::copy(&mut response.into_reader(), &mut tmp_file)
                    .context("Failed to write Pass 2 temp file")?;
                paths.push(tmp_path);
            }
            paths
        }
        ReadSourcePlan::Unavailable => {
            return Ok(Vec::new());
        }
    };

    let progress = if show_progress && plan.max_reads > 0 {
        let pb = ProgressBar::new(plan.max_reads);
        pb.set_style(
            ProgressStyle::with_template("Pass 2 {bar:30} {human_pos}/{human_len} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("Pass 2 {spinner} {human_pos} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    let mut records: Vec<ReadRecord> = Vec::new();
    let mut index: u64 = 0;

    if plan.paired && files.len() == 2 {
        // Paired-end: alternate R1 and R2
        let reader1 = open_fastq_reader(&files[0])?;
        let reader2 = open_fastq_reader(&files[1])?;
        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            // Read R1
            match read_sequence_only(&mut lines1)? {
                Some(seq) => {
                    records.push(ReadRecord {
                        sequence: seq,
                        index,
                        mate: Mate::R1,
                    });
                    index += 1;
                }
                None => break,
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                break;
            }

            // Read R2
            match read_sequence_only(&mut lines2)? {
                Some(seq) => {
                    records.push(ReadRecord {
                        sequence: seq,
                        index,
                        mate: Mate::R2,
                    });
                    index += 1;
                }
                None => break,
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                break;
            }
        }
    } else {
        // Unpaired: read files sequentially
        for file_path in &files {
            let reader = open_fastq_reader(file_path)?;
            let mut lines = reader.lines();

            while let Some(seq) = read_sequence_only(&mut lines)? {
                records.push(ReadRecord {
                    sequence: seq,
                    index,
                    mate: Mate::Unpaired,
                });
                index += 1;
                if index % 10_000 == 0 {
                    progress.set_position(index);
                }
                if plan.max_reads > 0 && index >= plan.max_reads {
                    break;
                }
            }
            if plan.max_reads > 0 && index >= plan.max_reads {
                break;
            }
        }
    }

    progress.finish_and_clear();
    info!(
        "Pass 2: collected {} reads for threading in {}",
        format_count(index),
        format_duration(start.elapsed())
    );

    Ok(records)
}

/// Read a single FASTQ record and return only the sequence line.
/// Returns None at EOF.
fn read_sequence_only<R: BufRead>(lines: &mut std::io::Lines<R>) -> Result<Option<String>> {
    // Header
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ header")?;
        }
        None => return Ok(None),
    }
    // Sequence
    let sequence = match lines.next() {
        Some(result) => result.context("Failed to read FASTQ sequence")?,
        None => bail!("Truncated FASTQ record: missing sequence line"),
    };
    // Separator
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ separator")?;
        }
        None => bail!("Truncated FASTQ record: missing separator line"),
    }
    // Quality
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ quality")?;
        }
        None => bail!("Truncated FASTQ record: missing quality line"),
    }
    Ok(Some(sequence))
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
        let header: String = std::iter::once("count".to_string())
            .chain((1..=histos.len()).map(|i| format!("chunk_{}", i)))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(file, "{}", header).context("Failed to write histogram header")?;

        // Data rows — precompute histogram vectors to avoid re-computing per row
        let histo_vecs: Vec<Vec<u64>> = histos
            .iter()
            .map(kmer::Histogram::get_vector)
            .collect::<Result<_>>()?;
        for i in 1..args.histo_max as usize + 2 {
            let line: String = std::iter::once(i.to_string())
                .chain(histo_vecs.iter().map(|v| v[i].to_string()))
                .collect::<Vec<_>>()
                .join("\t");
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
