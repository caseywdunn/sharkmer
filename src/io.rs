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
pub(crate) struct RetainedRead {
    /// Raw sequence (ASCII ACGT)
    pub(crate) sequence: String,
}

/// Collects reads retained during Pass 1 primer Oligo matching.
pub(crate) struct RetainedReads {
    pub(crate) reads: Vec<RetainedRead>,
}

/// Two-stage primer Oligo filter for Pass 1 read retention.
///
/// Stage 1 (bloom filter): fast approximate check during kmer ingestion.
/// Uses a bloom filter with two hash functions for O(1) lookup with no
/// hash table overhead. May produce false positives, never false negatives.
///
/// Stage 2 (exact set): after ingestion, bloom-positive reads are verified
/// against an exact AHashSet to eliminate false positives.
///
/// When a panel mixes genes with different trim lengths (e.g. 15 for most
/// genes but 12 for a trimmed primer), the filter stores a separate
/// (mask, shift) pair per distinct trim length. `check_kmer` tries every
/// pair against the incoming kmer — previously only the last gene's pair
/// was retained, silently dropping reads matching earlier genes.
pub(crate) struct OligoFilter {
    /// Bloom filter for high-bit matching (oligo at kmer start).
    /// Entries are keyed by raw oligo kmer values (shift removed at
    /// insertion time), so a single bloom serves all trim lengths.
    high_bits: Vec<u64>,
    /// Bloom filter for low-bit matching (RC oligo at kmer end).
    /// Entries are keyed by raw rc oligo values.
    low_bits: Vec<u64>,
    /// One (mask, shift) pair per distinct trim length seen across genes.
    /// Typical panels have 1-3 entries here.
    high_masks: Vec<(u64, u32)>,
    /// One low-bit mask per distinct trim length seen across genes.
    /// Shift is always 0 for the low-bit side.
    low_masks: Vec<u64>,
    /// Exact set for high-bit verification, keyed by raw oligo values.
    high_exact: ahash::AHashSet<u64>,
    /// Exact set for low-bit verification, keyed by raw rc oligo values.
    low_exact: ahash::AHashSet<u64>,
}

/// Number of bits used for bloom filter indexing. 2^24 = 16M entries = 2MB per bitset.
/// Two hash functions reduce false positive rate. With ~50K Oligo variants,
/// occupancy is ~0.3%, giving FP rate ~(0.003)^2 ≈ 0.001%.
const BLOOM_BITS: u32 = 24;
const BLOOM_SIZE: usize = 1 << BLOOM_BITS;
#[allow(dead_code)]
const BLOOM_WORDS: usize = BLOOM_SIZE / 64;

impl OligoFilter {
    /// Build a filter from pre-encoded primer Oligo sets for all genes.
    ///
    /// Each oligo length seen across the input contributes one entry to
    /// `high_masks` and `low_masks`; duplicate lengths are deduplicated.
    /// The bloom filters and exact sets store raw oligo values (shift
    /// removed), so a single shared set of each covers all trim lengths.
    #[allow(dead_code)]
    pub(crate) fn new(oligo_sets: &[crate::pcr::PrimerOligoSet], k: usize) -> Self {
        let mut high_bits = vec![0u64; BLOOM_WORDS];
        let mut low_bits = vec![0u64; BLOOM_WORDS];
        let mut high_exact = ahash::AHashSet::new();
        let mut low_exact = ahash::AHashSet::new();
        let mut seen_lengths: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut high_masks: Vec<(u64, u32)> = Vec::new();
        let mut low_masks: Vec<u64> = Vec::new();

        let register_length = |len: usize,
                               seen: &mut std::collections::HashSet<usize>,
                               high: &mut Vec<(u64, u32)>,
                               low: &mut Vec<u64>| {
            if len == 0 || !seen.insert(len) {
                return;
            }
            let shift = 2 * (k - len);
            let size_mask = (1u64 << (2 * len)) - 1;
            high.push((size_mask << shift, shift as u32));
            low.push(size_mask);
        };

        for set in oligo_sets {
            let fwd_len = set.forward_oligo_length;
            if fwd_len > 0 {
                register_length(fwd_len, &mut seen_lengths, &mut high_masks, &mut low_masks);
                for oligo in &set.forward_oligos {
                    // Bloom and exact sets are keyed by the raw oligo value
                    // (no shift baked in), so the same entry is checked at
                    // any trim length via the corresponding (mask, shift).
                    bloom_insert(&mut high_bits, oligo.kmer);
                    high_exact.insert(oligo.kmer);
                    let rc = crate::kmer::revcomp_kmer(&oligo.kmer, &fwd_len);
                    bloom_insert(&mut low_bits, rc);
                    low_exact.insert(rc);
                }
            }

            let rev_len = set.reverse_oligo_length;
            if rev_len > 0 {
                register_length(rev_len, &mut seen_lengths, &mut high_masks, &mut low_masks);
                for oligo in &set.reverse_oligos {
                    bloom_insert(&mut high_bits, oligo.kmer);
                    high_exact.insert(oligo.kmer);
                    let rc = crate::kmer::revcomp_kmer(&oligo.kmer, &rev_len);
                    bloom_insert(&mut low_bits, rc);
                    low_exact.insert(rc);
                }
            }
        }

        OligoFilter {
            high_bits,
            low_bits,
            high_masks,
            low_masks,
            high_exact,
            low_exact,
        }
    }

    /// Fast approximate check during kmer ingestion (bloom filter).
    /// May return false positives, never false negatives. Tries every
    /// registered trim-length mask; early-exits on the first hit.
    #[inline]
    pub(crate) fn check_kmer(&self, kmer: u64) -> bool {
        for &(mask, shift) in &self.high_masks {
            let high_val = (kmer & mask) >> shift;
            if bloom_contains(&self.high_bits, high_val) {
                return true;
            }
        }
        for &mask in &self.low_masks {
            if bloom_contains(&self.low_bits, kmer & mask) {
                return true;
            }
        }
        false
    }

    /// Exact verification of a read against the Oligo set (AHashSet).
    /// Called on bloom-positive reads to eliminate false positives.
    pub(crate) fn verify_read(&self, sequence: &str, k: usize) -> bool {
        let kmers = match crate::kmer::encoding::kmers_from_ascii(sequence, k) {
            Ok(k) => k,
            Err(_) => return false,
        };
        for kmer in &kmers {
            for &(mask, shift) in &self.high_masks {
                if self.high_exact.contains(&((kmer & mask) >> shift)) {
                    return true;
                }
            }
            for &mask in &self.low_masks {
                if self.low_exact.contains(&(kmer & mask)) {
                    return true;
                }
            }
        }
        false
    }
}

/// Two hash functions for the bloom filter. Both derived from the same
/// u64 value using different bit mixing to produce independent indices.
#[inline]
fn bloom_h1(val: u64) -> usize {
    (val as usize) & (BLOOM_SIZE - 1)
}

#[inline]
fn bloom_h2(val: u64) -> usize {
    // Rotate bits and XOR to get a second independent hash
    let mixed = (val >> 17) ^ (val.wrapping_mul(0x9E3779B97F4A7C15));
    (mixed as usize) & (BLOOM_SIZE - 1)
}

/// Insert a value into the bloom filter using two hash functions.
#[inline]
#[allow(dead_code)]
fn bloom_insert(bits: &mut [u64], val: u64) {
    let i1 = bloom_h1(val);
    bits[i1 / 64] |= 1u64 << (i1 % 64);
    let i2 = bloom_h2(val);
    bits[i2 / 64] |= 1u64 << (i2 % 64);
}

/// Check if a value may be in the bloom filter (both bits must be set).
#[inline]
fn bloom_contains(bits: &[u64], val: u64) -> bool {
    let i1 = bloom_h1(val);
    if (bits[i1 / 64] & (1u64 << (i1 % 64))) == 0 {
        return false;
    }
    let i2 = bloom_h2(val);
    (bits[i2 / 64] & (1u64 << (i2 % 64))) != 0
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

/// Build a context-rich error message for an I/O error that occurred mid-stream
/// during FASTQ reading. Distinguishes abrupt connection drops (the common ENA
/// failure mode) from actual FASTQ parse issues, and suggests the cache as a
/// workaround when streaming fails.
fn stream_io_error(
    line_role: &str,
    err: &std::io::Error,
    source_name: &str,
    n_reads_read: u64,
) -> anyhow::Error {
    use std::io::ErrorKind;
    let is_stream_drop = matches!(
        err.kind(),
        ErrorKind::UnexpectedEof
            | ErrorKind::ConnectionReset
            | ErrorKind::ConnectionAborted
            | ErrorKind::BrokenPipe
            | ErrorKind::TimedOut
            | ErrorKind::Interrupted
    );
    let is_remote = source_name.starts_with("http://") || source_name.starts_with("https://");

    if is_stream_drop {
        if is_remote {
            anyhow::anyhow!(
                "Stream from {} dropped while reading {} line of record {} (I/O error: {} — kind {:?}).\n\
                 This is usually a transient network interruption, not a bad FASTQ file.\n\
                 Retry the run. If it repeats, use the read cache (the default for --ena); \
                 cached downloads are verified by SHA-256 and do not suffer mid-stream drops.",
                source_name,
                line_role,
                n_reads_read + 1,
                err,
                err.kind(),
            )
        } else {
            anyhow::anyhow!(
                "Local read stream ended unexpectedly while reading {} line of record {} in {} \
                 (I/O error: {} — kind {:?}). The file may be truncated or corrupted.",
                line_role,
                n_reads_read + 1,
                source_name,
                err,
                err.kind(),
            )
        }
    } else {
        anyhow::anyhow!(
            "Failed to read {} line of record {} in {}: {} (kind {:?})",
            line_role,
            n_reads_read + 1,
            source_name,
            err,
            err.kind(),
        )
    }
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
    oligo_filter: Option<&OligoFilter>,
) -> Result<bool> {
    let mut lines = reader.lines();
    let n_chunks = state.chunks.len();

    while let Some(header_result) = lines.next() {
        // Read 4 lines for one FASTQ record
        let header = header_result
            .map_err(|e| stream_io_error("header", &e, source_name, state.n_reads_read))?;

        let sequence = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("sequence", &e, source_name, state.n_reads_read))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing sequence line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let separator = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("separator", &e, source_name, state.n_reads_read))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing separator line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let quality = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("quality", &e, source_name, state.n_reads_read))?
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
            drain_batch(state, oligo_filter, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

/// Drain accumulated sequences into the current chunk, optionally checking
/// each kmer against the Oligo filter during ingestion. Matching reads
/// are retained for seed evaluation.
fn drain_batch(
    state: &mut FastqReadState,
    oligo_filter: Option<&OligoFilter>,
    n_chunks: usize,
) -> Result<()> {
    if let Some(filter) = oligo_filter {
        for seq in state.seqs.drain(..) {
            let matched = state.chunks[state.chunk_index].ingest_seq_with_filter(&seq, filter)?;
            if matched {
                state
                    .retained_reads
                    .reads
                    .push(RetainedRead { sequence: seq });
            }
        }
    } else {
        for seq in state.seqs.drain(..) {
            state.chunks[state.chunk_index].ingest_seq(&seq)?;
        }
    }
    state.chunk_index = (state.chunk_index + 1) % n_chunks;
    Ok(())
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
    oligo_filter: Option<&OligoFilter>,
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

    let mut max_reads = args.max_reads.unwrap_or(0);

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
                oligo_filter,
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
            // Round max_reads up to even so Pass 1 consumes balanced
            // pairs. The rounded value is also recorded in the ReadPlan
            // below so Pass 2 (re-read for threading) consumes the same
            // number of reads; leaving the unrounded value here would
            // cause Pass 2 to stop one R2 short of Pass 1 and drop the
            // final read pair silently.
            if max_reads > 0 && max_reads % 2 != 0 {
                max_reads += 1;
            }
            read_fastq_paired(
                reader1,
                reader2,
                &mut state,
                max_reads,
                args.validate_every,
                &name1,
                &name2,
                &progress,
                oligo_filter,
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
                    oligo_filter,
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
            oligo_filter,
        )?;
        // read_source_plan remains Unavailable for stdin
    }

    progress.finish_and_clear();

    // Ingest any remaining sequences (with Oligo filter check if active)
    let n_chunks = state.chunks.len();
    drain_batch(&mut state, oligo_filter, n_chunks)?;

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
    // Verify bloom-positive reads against exact AHashSet and report counts
    if !state.retained_reads.reads.is_empty() {
        let bloom_count = state.retained_reads.reads.len();
        if let Some(filter) = oligo_filter {
            state
                .retained_reads
                .reads
                .retain(|r| filter.verify_read(&r.sequence, k));
            let verified_count = state.retained_reads.reads.len();
            info!(
                "Retained reads for seed evaluation: {} bloom hits, {} verified",
                format_count(bloom_count as u64),
                format_count(verified_count as u64)
            );
        } else {
            info!(
                "Retained {} reads for seed evaluation",
                format_count(bloom_count as u64)
            );
        }
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

    // Open the file once. If we need to peek at the gzip magic bytes, do so
    // through the BufReader: fill_buf() loads bytes into the internal buffer
    // without consuming them, so subsequent reads (including through a
    // GzDecoder wrapper) see the stream from position 0. This avoids a
    // redundant second open that showed up as two syscalls per FASTQ file.
    let file = std::fs::File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
    let mut buf_reader = std::io::BufReader::new(file);

    let use_gzip = if has_gz_ext {
        true
    } else {
        let magic = buf_reader.fill_buf().context("Failed to peek at file")?;
        magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b
    };

    if use_gzip {
        Ok(Box::new(std::io::BufReader::new(
            flate2::read::GzDecoder::new(buf_reader),
        )))
    } else {
        Ok(Box::new(buf_reader))
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
    oligo_filter: Option<&OligoFilter>,
) -> Result<bool> {
    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();
    let n_chunks = state.chunks.len();
    let mut r1_records: u64 = 0;
    let mut r2_records: u64 = 0;
    let mismatch_exit: bool;

    loop {
        // Read one record from R1
        let r1_done = read_one_fastq_record(&mut lines1, state, validate_every, source1_name)?;
        if r1_done {
            // R1 exhausted. Verify R2 is also at EOF so we don't silently drop
            // unmatched R2 records (truncated R1, mismatched inputs).
            let r2_extra = read_one_fastq_record(&mut lines2, state, validate_every, source2_name)?;
            mismatch_exit = !r2_extra;
            if !r2_extra {
                r2_records += 1;
            }
            break;
        }
        r1_records += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            drain_batch(state, oligo_filter, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }

        // Read one record from R2
        let r2_done = read_one_fastq_record(&mut lines2, state, validate_every, source2_name)?;
        if r2_done {
            // R2 ended mid-pair: R1 already consumed one more record than R2.
            mismatch_exit = true;
            break;
        }
        r2_records += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            drain_batch(state, oligo_filter, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }
    }

    if mismatch_exit && r1_records != r2_records {
        warn!(
            "Paired-end input length mismatch: {} has {} reads, {} has {} reads. \
             Extra reads in the longer file were not processed.",
            source1_name, r1_records, source2_name, r2_records
        );
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
            result.map_err(|e| stream_io_error("header", &e, source_name, state.n_reads_read))?
        }
        None => return Ok(true), // EOF
    };

    let sequence = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("sequence", &e, source_name, state.n_reads_read))?
        }
        None => bail!(
            "Truncated FASTQ record at record {} in {}: missing sequence line",
            state.n_reads_read + 1,
            source_name,
        ),
    };

    let separator = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("separator", &e, source_name, state.n_reads_read))?
        }
        None => bail!(
            "Truncated FASTQ record at record {} in {}: missing separator line",
            state.n_reads_read + 1,
            source_name,
        ),
    };

    let quality = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("quality", &e, source_name, state.n_reads_read))?
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

    // Temp files held here are auto-deleted when this vec drops at end of
    // function. Each NamedTempFile has a unique path (random suffix), so
    // concurrent sharkmer invocations and same-length URLs cannot collide.
    let mut _pass2_tempfiles: Vec<tempfile::NamedTempFile> = Vec::new();
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
                // Create a unique temp file (random suffix) in the system temp
                // directory. The NamedTempFile guard is kept alive in
                // `_pass2_tempfiles` below so the file persists through Pass 2
                // reading and is deleted on function return.
                let mut tmp = tempfile::Builder::new()
                    .prefix("sharkmer_pass2_")
                    .suffix(".fastq.gz")
                    .tempfile()
                    .context("Failed to create temp file for Pass 2")?;
                std::io::copy(&mut response.into_reader(), tmp.as_file_mut())
                    .context("Failed to write Pass 2 temp file")?;
                paths.push(tmp.path().to_path_buf());
                _pass2_tempfiles.push(tmp);
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
        // Paired-end: alternate R1 and R2. Track per-side record counts
        // so a truncated mate file can be surfaced to the user rather
        // than silently producing an unbalanced read set for Pass 2
        // threading. This mirrors the check in read_fastq_paired for
        // Pass 1.
        let reader1 = open_fastq_reader(&files[0])?;
        let reader2 = open_fastq_reader(&files[1])?;
        let name1 = files[0].to_string_lossy().to_string();
        let name2 = files[1].to_string_lossy().to_string();
        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();
        let mut r1_records: u64 = 0;
        let mut r2_records: u64 = 0;
        let mismatch_exit: bool;

        'outer: loop {
            // Read R1
            match read_sequence_only(&mut lines1)? {
                Some(seq) => {
                    records.push(ReadRecord {
                        sequence: seq,
                        index,
                        mate: Mate::R1,
                    });
                    index += 1;
                    r1_records += 1;
                }
                None => {
                    // R1 exhausted. Check whether R2 also reached EOF —
                    // if not, R2 has extra unmatched records.
                    let r2_extra = read_sequence_only(&mut lines2)?.is_some();
                    mismatch_exit = r2_extra;
                    if r2_extra {
                        r2_records += 1; // counted so the warning numbers are meaningful
                    }
                    break 'outer;
                }
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                mismatch_exit = false;
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
                    r2_records += 1;
                }
                None => {
                    // R2 ended mid-pair: R1 already consumed one more.
                    mismatch_exit = true;
                    break 'outer;
                }
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                mismatch_exit = false;
                break;
            }
        }

        if mismatch_exit && r1_records != r2_records {
            warn!(
                "Pass 2 paired-end input length mismatch: {} has {} reads, {} has {} reads. \
                 Extra reads in the longer file were not re-read for threading.",
                name1, r1_records, name2, r2_records
            );
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

    let estimated_capacity: usize = state.chunks.iter().map(|c| c.get_kmer_counts().len()).sum();
    let mut kmer_counts: KmerCounts = KmerCounts::new_with_capacity(&k, estimated_capacity);
    let mut n_singleton_kmers: Option<u64> = None;

    let histo_comment = format!(
        "# sharkmer {} k={} chunks={}",
        env!("CARGO_PKG_VERSION"),
        args.k,
        args.chunks
    );

    if args.chunks > 0 {
        // Incremental histogram mode: update histogram as chunks are merged.
        // Snapshot only the Vec<u64> at each step (not the full Histogram struct
        // with its FxHashMap).
        let mut histo_vecs: Vec<Vec<u64>> = Vec::with_capacity(args.chunks);
        let mut running_histo = kmer::Histogram::new(&args.histo_max);

        for chunk in state.chunks.drain(..) {
            kmer_counts.extend_with_histogram(chunk.get_kmer_counts(), &mut running_histo)?;
            drop(chunk);

            histo_vecs.push(running_histo.get_vector()?);
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
            .chain((1..=histo_vecs.len()).map(|i| format!("chunk_{}", i)))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(file, "{}", header).context("Failed to write histogram header")?;

        // Data rows — histo_vecs already contains the precomputed vectors
        for i in 1..args.histo_max as usize + 2 {
            let line: String = std::iter::once(i.to_string())
                .chain(histo_vecs.iter().map(|v| v[i].to_string()))
                .collect::<Vec<_>>()
                .join("\t");
            writeln!(file, "{}", line).context("Failed to write histogram data")?;
        }

        // Write the final histogram with header
        info!("Writing final histogram to file...");
        let last_histo_vec = histo_vecs.last().context("No histograms were produced")?;

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
        let n_unique = running_histo.get_n_unique_kmers();
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

        let n_unique_kmers_histo: u64 = running_histo.get_n_unique_kmers();
        let n_kmers_histo: u64 = running_histo.get_n_kmers();

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcr::{Oligo, PrimerOligoSet};

    #[test]
    fn test_bloom_insert_and_contains() {
        let mut bits = vec![0u64; BLOOM_WORDS];
        let val = 0xDEAD_BEEF_u64;

        assert!(!bloom_contains(&bits, val));
        bloom_insert(&mut bits, val);
        assert!(bloom_contains(&bits, val));
    }

    #[test]
    fn test_bloom_no_false_negatives() {
        let mut bits = vec![0u64; BLOOM_WORDS];
        let values: Vec<u64> = (0..1000).collect();
        for &v in &values {
            bloom_insert(&mut bits, v);
        }
        for &v in &values {
            assert!(
                bloom_contains(&bits, v),
                "bloom must not have false negatives for {}",
                v
            );
        }
    }

    #[test]
    fn test_oligo_filter_empty() {
        let filter = OligoFilter::new(&[], 19);
        // Empty filter should never match
        assert!(!filter.check_kmer(0));
        assert!(!filter.check_kmer(u64::MAX));
    }

    #[test]
    fn test_oligo_filter_verify_read_eliminates_bloom_false_positive() {
        // Build a filter with a known oligo
        let oligo_val: u64 = 0b1001_1000; // "GCGA" in 2-bit encoding, length 4
        let set = PrimerOligoSet {
            gene_name: "test".to_string(),
            forward_oligos: vec![Oligo {
                length: 4,
                kmer: oligo_val,
            }],
            reverse_oligos: Vec::new(),
            forward_oligo_length: 4,
            reverse_oligo_length: 0,
        };
        let k = 5;
        let filter = OligoFilter::new(&[set], k);

        // A read that does NOT contain "GCGA" should fail exact verification
        // even if bloom is a false positive
        assert!(!filter.verify_read("AAAAA", k));
    }

    /// When a panel mixes genes with different trim lengths, the filter must
    /// detect oligos from BOTH lengths, not just whichever gene was iterated
    /// last. Prior to the fix, check_kmer only matched the last gene's
    /// trim length; kmers placed at any other position silently missed.
    #[test]
    fn test_oligo_filter_mixed_trim_lengths() {
        let k: usize = 19;

        // Two distinct non-trivial oligos at different trim lengths. The
        // bit patterns are chosen so that a truncation of either cannot
        // match the other: oligo_a's top 12 bases differ from oligo_b,
        // and oligo_b as a suffix of oligo_a's bit pattern would require
        // specific alignment that doesn't occur here.
        let oligo_a_value: u64 = 0x2AF3_B1C9; // 30 bits, length-15
        let oligo_b_value: u64 = 0x00AB_CDEF; // 24 bits, length-12

        // Set A declared FIRST, set B declared LAST. Under the old code,
        // only set_b's (mask, shift) would be retained and oligo_a at its
        // length-15 position would be unreachable from check_kmer.
        let set_a = PrimerOligoSet {
            gene_name: "A".to_string(),
            forward_oligos: vec![Oligo {
                length: 15,
                kmer: oligo_a_value,
            }],
            reverse_oligos: Vec::new(),
            forward_oligo_length: 15,
            reverse_oligo_length: 0,
        };
        let set_b = PrimerOligoSet {
            gene_name: "B".to_string(),
            forward_oligos: vec![Oligo {
                length: 12,
                kmer: oligo_b_value,
            }],
            reverse_oligos: Vec::new(),
            forward_oligo_length: 12,
            reverse_oligo_length: 0,
        };

        let filter = OligoFilter::new(&[set_a, set_b], k);

        // Two distinct trim lengths must be registered, not one.
        assert_eq!(
            filter.high_masks.len(),
            2,
            "filter must retain a mask per distinct trim length"
        );

        // A kmer with oligo_a at the length-15 position (high 30 bits),
        // arbitrary bits elsewhere. check_kmer must find it via the
        // length-15 mask; the old code silently missed this.
        let kmer_with_a = (oligo_a_value << (2 * (k - 15))) | 0b01;
        assert!(
            filter.check_kmer(kmer_with_a),
            "oligo_a at its length-15 position must match (old code missed \
             this when set_b's trim length overwrote set_a's mask)"
        );

        // And oligo_b at its length-12 position still works.
        let kmer_with_b = (oligo_b_value << (2 * (k - 12))) | 0b10;
        assert!(
            filter.check_kmer(kmer_with_b),
            "oligo_b at its length-12 position must match"
        );
    }
}
