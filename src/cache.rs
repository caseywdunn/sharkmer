use crate::kmer::budget::ResourceBudget;
use crate::output::MetadataTemporary as CacheTemporary;
use anyhow::{Context, Result};
use fs2::FileExt;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Sidecar metadata stored alongside each cached file.
#[derive(Serialize, Deserialize, Debug)]
pub(crate) struct CacheMeta {
    /// The original URL this file was downloaded from.
    url: String,
    /// SHA-256 hex digest of the cached data file.
    sha256: String,
    /// Whether the download completed fully (true) or was truncated at max_reads.
    complete: bool,
    /// Number of FASTQ reads stored in the cached file.
    #[serde(default)]
    n_reads: u64,
}

mod metadata;
mod transaction;

/// Configuration and operations for the read cache.
pub(crate) struct CacheConfig {
    pub(crate) cache_dir: PathBuf,
    _lock_file: File,
}

impl CacheConfig {
    /// Create a new cache configuration.
    /// Uses `cache_dir_override` if provided, otherwise the platform default
    /// (`~/.cache/sharkmer/reads/` on Linux, `~/Library/Caches/sharkmer/reads/` on macOS).
    pub(crate) fn new(cache_dir_override: Option<&Path>) -> Result<Self> {
        let requested_cache_dir = if let Some(dir) = cache_dir_override {
            dir.to_path_buf()
        } else {
            let base = dirs::cache_dir().context("Could not determine platform cache directory")?;
            base.join("sharkmer").join("reads")
        };
        fs::create_dir_all(&requested_cache_dir).with_context(|| {
            format!(
                "Failed to create cache directory: {}",
                requested_cache_dir.display()
            )
        })?;
        let cache_dir = fs::canonicalize(&requested_cache_dir).with_context(|| {
            format!(
                "Failed to canonicalize cache directory: {}",
                requested_cache_dir.display()
            )
        })?;
        let lock_file = acquire_cache_lock(&cache_dir)?;
        Ok(CacheConfig {
            cache_dir,
            _lock_file: lock_file,
        })
    }

    /// Look up a URL in the cache.
    /// Returns `Some(path)` if a valid cached file exists with enough reads,
    /// `None` on cache miss. A complete download (EOF reached) satisfies any
    /// `max_reads` value. A partial download satisfies only if it contains
    /// at least `max_reads` reads.
    ///
    /// **Cache trust model:** the SHA-256 checksum is verified on every
    /// lookup. Ambiguous entries are preserved rather than removed, because
    /// their contents cannot safely be attributed to sharkmer.
    pub(crate) fn lookup(&self, url: &str, max_reads: u64) -> Result<Option<PathBuf>> {
        self.lookup_with_budget(url, max_reads, None)
    }

    /// Validate cached source metadata under the exact counter's shared ledger.
    pub(crate) fn lookup_bounded(
        &self,
        url: &str,
        max_reads: u64,
        budget: &Arc<ResourceBudget>,
    ) -> Result<Option<PathBuf>> {
        self.lookup_with_budget(url, max_reads, Some(budget))
    }

    fn lookup_with_budget(
        &self,
        url: &str,
        max_reads: u64,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Option<PathBuf>> {
        transaction::recover(&self.cache_dir, budget)?;
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);
        let meta = match validate_entry_with_budget(&key, &data_path, &meta_path, budget)? {
            Some(meta) => meta,
            None => return Ok(None),
        };

        // Check if cached reads are sufficient for the requested max_reads.
        // A complete download (EOF reached) satisfies any request.
        if !meta.complete {
            let insufficient = if max_reads == 0 {
                true // wants all reads, but cache is partial
            } else {
                meta.n_reads < max_reads
            };
            if insufficient {
                info!(
                    "Cache has {} reads but {} requested; retaining it until a replacement download succeeds for {}",
                    meta.n_reads,
                    if max_reads == 0 {
                        "all".to_string()
                    } else {
                        max_reads.to_string()
                    },
                    url
                );
                return Ok(None);
            }
        }

        Ok(Some(data_path))
    }

    /// Download a URL to the cache directory, storing at most `max_reads` FASTQ
    /// reads (0 = unlimited). The remote gzipped FASTQ is decompressed on the
    /// fly, and only the retained reads are re-compressed into the cache file.
    /// Returns the path to the cached file.
    pub(crate) fn download_to_cache(&self, url: &str, max_reads: u64) -> Result<PathBuf> {
        self.download(url, max_reads, None)
    }

    /// Bound remote cache-miss records under the exact counter's shared input
    /// allowance before allocating a line or replacing a verified cache entry.
    pub(crate) fn download_to_cache_bounded(
        &self,
        url: &str,
        max_reads: u64,
        record_bytes: usize,
        budget: &Arc<ResourceBudget>,
    ) -> Result<PathBuf> {
        self.download(url, max_reads, Some((record_bytes, budget)))
    }

    fn download(
        &self,
        url: &str,
        max_reads: u64,
        resources: Option<(usize, &Arc<ResourceBudget>)>,
    ) -> Result<PathBuf> {
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);
        let budget = resources.map(|(_, budget)| budget);
        transaction::recover(&self.cache_dir, budget)?;
        let previous = validate_entry_with_budget(&key, &data_path, &meta_path, budget)?.is_some();

        info!(
            "Downloading {} to cache (max_reads: {})...",
            url,
            if max_reads == 0 {
                "unlimited".to_string()
            } else {
                max_reads.to_string()
            }
        );

        // Input pipeline: HTTP → MultiGzDecoder → BufReader
        let response =
            crate::http::get(url).with_context(|| format!("Failed to download {}", url))?;
        let transaction = transaction::Transaction::begin(&self.cache_dir, &key, previous, budget)?;
        let mut temporary = CacheTemporary::create(transaction.directory(), ".new-data-", budget)?;
        let buf_reader = crate::codec::Decoder::new(response, budget)?;

        // Finalize the encoder and its buffer before consulting the owned
        // writer. Its retained cause preserves a typed quota refusal even
        // when compression or buffered I/O returns a generic io::Error.
        let written = (|| {
            let mut gz_writer = crate::codec::Encoder::new(temporary.writer(), budget)?;
            let selection =
                write_cache_records(buf_reader, &mut gz_writer, url, max_reads, resources)?;
            gz_writer
                .finish()
                .context("Failed to finish and flush gzip cache file")?;
            Ok(selection)
        })();
        let (n_reads, complete) = temporary.finish(written)?;

        info!(
            "Cached {} reads for {} (complete: {})",
            n_reads, url, complete
        );

        // Compute SHA-256 of the cached file
        let sha256 = compute_sha256(temporary.path())?;
        let meta = metadata::Owned::generated(url, sha256, complete, n_reads, budget)?;
        let meta_temporary = metadata::temporary(transaction.directory(), &meta, budget)?;
        transaction.publish(temporary, meta_temporary)?;

        Ok(data_path)
    }

    /// Clear verified cache entries while preserving the cache directory and lock.
    pub(crate) fn clear(cache_dir_override: Option<&Path>) -> Result<()> {
        let requested_cache_dir = if let Some(dir) = cache_dir_override {
            dir.to_path_buf()
        } else {
            let base = dirs::cache_dir().context("Could not determine platform cache directory")?;
            base.join("sharkmer").join("reads")
        };
        fs::create_dir_all(&requested_cache_dir).with_context(|| {
            format!(
                "Failed to create cache directory: {}",
                requested_cache_dir.display()
            )
        })?;
        let cache_dir = fs::canonicalize(&requested_cache_dir).with_context(|| {
            format!(
                "Failed to canonicalize cache directory: {}",
                requested_cache_dir.display()
            )
        })?;
        let _lock_file = acquire_cache_lock(&cache_dir)?;
        transaction::recover(&cache_dir, None)?;
        let cleared = clear_owned_entries(&cache_dir)?;
        info!(
            "Cleared {} verified cache entr{}; preserving unowned, malformed, modified, and symlinked paths.",
            cleared,
            if cleared == 1 { "y" } else { "ies" }
        );
        Ok(())
    }

    fn data_path(&self, key: &str) -> PathBuf {
        self.cache_dir.join(format!("{}.fastq.gz", key))
    }

    fn meta_path(&self, key: &str) -> PathBuf {
        self.cache_dir.join(format!("{}.meta.yaml", key))
    }
}

/// Preserve cache selection/normalization while bounding exact-run records.
fn write_cache_records<R: BufRead, W: Write>(
    reader: R,
    output: &mut W,
    url: &str,
    max_reads: u64,
    resources: Option<(usize, &Arc<ResourceBudget>)>,
) -> Result<(u64, bool)> {
    if let Some((limit, budget)) = resources {
        let mut records = crate::fastq::FastqReader::new(reader, limit, Arc::clone(budget))?;
        let mut n_reads = 0;
        let complete = loop {
            if !records.next(url, n_reads)? {
                break true;
            }
            for field in records.fields() {
                output.write_all(field)?;
                output.write_all(b"\n")?;
            }
            n_reads = n_reads
                .checked_add(1)
                .context("Cache record count overflow")?;
            if max_reads > 0 && n_reads >= max_reads {
                break false;
            }
        };
        return Ok((n_reads, complete));
    }
    let mut lines = reader.lines();
    let mut n_reads: u64 = 0;
    let mut complete = false;

    loop {
        // Read header line (first line of a FASTQ record)
        let header = match lines.next() {
            Some(Ok(line)) => line,
            Some(Err(e)) => {
                return Err(e).with_context(|| format!("Failed reading FASTQ from {}", url));
            }
            None => {
                complete = true;
                break;
            }
        };

        // Read remaining 3 lines of the FASTQ record
        let sequence = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url))?
            .with_context(|| format!("Failed reading FASTQ from {}", url))?;
        let separator = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url))?
            .with_context(|| format!("Failed reading FASTQ from {}", url))?;
        let quality = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url))?
            .with_context(|| format!("Failed reading FASTQ from {}", url))?;

        // Write the 4 lines to the cache file
        writeln!(output, "{}", header)?;
        writeln!(output, "{}", sequence)?;
        writeln!(output, "{}", separator)?;
        writeln!(output, "{}", quality)?;

        n_reads += 1;

        if max_reads > 0 && n_reads >= max_reads {
            break;
        }
    }
    Ok((n_reads, complete))
}

fn cache_key(url: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(url.as_bytes());
    format!("{:x}", hasher.finalize())
}

/// Compute the SHA-256 hex digest of a file.
fn compute_sha256(path: &Path) -> Result<String> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open {} for checksum", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 8192];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

#[cfg(test)]
fn read_meta(path: &Path) -> Result<Option<metadata::Owned>> {
    metadata::read(path, None)
}

#[cfg(test)]
fn write_meta(path: &Path, meta: &CacheMeta) -> Result<()> {
    let contents = serde_yaml_ng::to_string(meta)
        .map_err(|e| anyhow::anyhow!("Failed to serialize cache sidecar: {}", e))?;
    fs::write(path, contents)
        .with_context(|| format!("Failed to write sidecar {}", path.display()))?;
    Ok(())
}

fn acquire_cache_lock(cache_dir: &Path) -> Result<File> {
    let lock_path = cache_dir.join(".sharkmer-cache.lock");
    reject_non_regular_existing_path(&lock_path, "cache lock")?;
    let lock_file = OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)
        .with_context(|| format!("Failed to open cache lock {}", lock_path.display()))?;
    lock_file
        .lock_exclusive()
        .with_context(|| format!("Failed to lock cache directory {}", cache_dir.display()))?;
    Ok(lock_file)
}

fn reject_non_regular_existing_path(path: &Path, label: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            anyhow::ensure!(
                metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
                "Refusing non-regular {} {}. Use a fresh cache directory or --no-cache.",
                label,
                path.display()
            );
            Ok(())
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect {}", path.display())),
    }
}

fn entry_file_state(path: &Path, label: &str) -> Result<bool> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            anyhow::ensure!(
                metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
                "Refusing non-regular cache {} {}. Use a fresh cache directory or --no-cache.",
                label,
                path.display()
            );
            Ok(true)
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(error) => {
            Err(error).with_context(|| format!("Failed to inspect cache {}", path.display()))
        }
    }
}

fn validate_entry(
    key: &str,
    data_path: &Path,
    meta_path: &Path,
) -> Result<Option<metadata::Owned>> {
    validate_entry_with_budget(key, data_path, meta_path, None)
}

fn validate_entry_with_budget(
    key: &str,
    data_path: &Path,
    meta_path: &Path,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Option<metadata::Owned>> {
    let data_exists = entry_file_state(data_path, "data")?;
    let meta_exists = entry_file_state(meta_path, "metadata")?;
    if !data_exists && !meta_exists {
        return Ok(None);
    }
    anyhow::ensure!(
        data_exists && meta_exists,
        "Refusing ambiguous cache entry for {}. Use a fresh cache directory or --no-cache.",
        key
    );
    let meta = metadata::read(meta_path, budget)
        .with_context(|| {
            format!(
                "Unable to validate cache metadata for {}; the cache entry is preserved",
                key
            )
        })?
        .context("Cache metadata disappeared during validation")?;
    anyhow::ensure!(
        cache_key(&meta.url) == key,
        "Refusing cache entry with metadata that does not own key {}. Use a fresh cache directory or --no-cache.",
        key
    );
    anyhow::ensure!(
        meta.sha256.len() == 64
            && meta
                .sha256
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()),
        "Refusing cache entry with an invalid checksum receipt for {}. Use a fresh cache directory or --no-cache.",
        key
    );
    let actual_sha256 = compute_sha256(data_path)?;
    anyhow::ensure!(
        actual_sha256 == meta.sha256,
        "Refusing modified cache entry for {}. Use a fresh cache directory or --no-cache.",
        key
    );
    Ok(Some(meta))
}

fn clear_owned_entries(cache_dir: &Path) -> Result<usize> {
    let mut cleared = 0usize;
    for entry in fs::read_dir(cache_dir)
        .with_context(|| format!("Failed to read cache directory {}", cache_dir.display()))?
    {
        let entry = entry
            .with_context(|| format!("Failed to read cache entry in {}", cache_dir.display()))?;
        let path = entry.path();
        let file_name = match path.file_name().and_then(|name| name.to_str()) {
            Some(name) => name,
            None => continue,
        };
        let key = match file_name.strip_suffix(".fastq.gz") {
            Some(key) if key.len() == 64 && key.bytes().all(|byte| byte.is_ascii_hexdigit()) => key,
            _ => match file_name.strip_suffix(".meta.yaml") {
                Some(key)
                    if key.len() == 64 && key.bytes().all(|byte| byte.is_ascii_hexdigit()) =>
                {
                    match entry_file_state(&path, "metadata") {
                        Ok(true) => {}
                        Ok(false) => continue,
                        Err(error) => {
                            warn!(
                                "Preserving unverified cache-shaped entry {}: {}",
                                path.display(),
                                error
                            );
                            continue;
                        }
                    }
                    let data_path = cache_dir.join(format!("{key}.fastq.gz"));
                    match entry_file_state(&data_path, "data") {
                        Ok(true) => {}
                        Ok(false) => warn!(
                            "Preserving unverified cache-shaped entry {}: matching data is missing",
                            path.display()
                        ),
                        Err(error) => warn!(
                            "Preserving unverified cache-shaped entry {}: {}",
                            path.display(),
                            error
                        ),
                    }
                    continue;
                }
                _ => continue,
            },
        };
        let meta_path = cache_dir.join(format!("{key}.meta.yaml"));
        match validate_entry(key, &path, &meta_path) {
            Ok(Some(_)) => {}
            Ok(None) => continue,
            Err(error) => {
                warn!(
                    "Preserving unverified cache-shaped entry {}: {}",
                    path.display(),
                    error
                );
                continue;
            }
        }
        fs::remove_file(&path)
            .with_context(|| format!("Failed to remove verified cache data {}", path.display()))?;
        fs::remove_file(&meta_path).with_context(|| {
            format!(
                "Failed to remove verified cache metadata {}",
                meta_path.display()
            )
        })?;
        cleared += 1;
    }
    Ok(cleared)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::TcpListener;

    const RECORD_ONE: &str = "@read-one\nACGT\n+\n!!!!\n";
    const RECORD_TWO: &str = "@read-two\nTGCA\n+\n####\n";

    fn cache_read_budget(bytes: u64) -> Arc<ResourceBudget> {
        use crate::kmer::budget::{MEMORY_CLASSES, MemoryClass, ResourceLimits};
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::InputBuffer as usize] = bytes;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    fn cache_metadata_budget(input_bytes: u64, disk: u64) -> Arc<ResourceBudget> {
        cache_budget(input_bytes, disk)
    }

    fn cache_http_budget(input_bytes: u64, disk: u64) -> Arc<ResourceBudget> {
        cache_metadata_budget(input_bytes, disk)
    }

    fn cache_budget(input_bytes: u64, disk: u64) -> Arc<ResourceBudget> {
        use crate::kmer::budget::{MEMORY_CLASSES, MemoryClass, ResourceLimits};
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::InputBuffer as usize] = input_bytes;
        classes[MemoryClass::RetainedResults as usize] = 8 << 20;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: disk,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn bounded_cache_records_preserve_normalization_selection_and_completion() {
        let input = b"@read-one\r\nACGT\r\n+\r\n!!!!\r\n@read-two\nTGCA\n+\n####";
        for maximum in [0, 1, 2, 3] {
            let budget = cache_read_budget(512);
            let mut legacy = Vec::new();
            let expected =
                write_cache_records(&input[..], &mut legacy, "source", maximum, None).unwrap();
            let mut bounded = Vec::new();
            let actual = write_cache_records(
                &input[..],
                &mut bounded,
                "source",
                maximum,
                Some((128, &budget)),
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(bounded, legacy);
            assert_eq!(
                actual,
                (
                    if maximum == 1 { 1 } else { 2 },
                    maximum == 0 || maximum > 2
                )
            );
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().peak_reserved_memory_bytes > 0);
        }
        for input in [&b"@bad\nACGT\n+\n"[..], &b"@bad\nACGT\n+\n!\xff!!\n"[..]] {
            let budget = cache_read_budget(512);
            assert!(write_cache_records(input, &mut Vec::new(), "source", 0, None).is_err());
            assert!(
                write_cache_records(input, &mut Vec::new(), "source", 0, Some((128, &budget)))
                    .is_err()
            );
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn bounded_http_cache_preserves_gzip_members_and_publishes_verified_records() {
        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let mut body = gzip_member(&RECORD_ONE.replace('\n', "\r\n"));
        body.extend(gzip_member(RECORD_TWO));
        let (url, server) = serve_once(body);
        let budget = cache_http_budget(512, 1 << 20);
        let path = config
            .download_to_cache_bounded(&url, 0, 128, &budget)
            .unwrap();
        server.join().unwrap().unwrap();
        let mut actual = String::new();
        flate2::read::MultiGzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut actual)
            .unwrap();
        assert_eq!(actual, format!("{RECORD_ONE}{RECORD_TWO}"));
        let meta = read_meta(&config.meta_path(&cache_key(&url)))
            .unwrap()
            .unwrap();
        assert_eq!(meta.n_reads, 2);
        assert!(meta.complete);
        assert_eq!(meta.sha256, compute_sha256(&path).unwrap());
        assert_eq!(config.lookup(&url, 0).unwrap(), Some(path));
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            fs::metadata(config.data_path(&cache_key(&url)))
                .unwrap()
                .len()
                + fs::metadata(config.meta_path(&cache_key(&url)))
                    .unwrap()
                    .len()
                + transaction::CONTROL_PEAK_BYTES
        );
    }

    #[test]
    fn compressed_cache_disk_refusal_preserves_prior_or_absent_entry() {
        use crate::kmer::budget::{ResourceExhausted, ResourceKind};
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        encoder.write_all(RECORD_ONE.as_bytes()).unwrap();
        let compressed = encoder.finish().unwrap();
        // Zero bytes and a missing gzip trailer refuse the compressed data;
        // allowing exactly the data refuses metadata while it is still owned.
        for prior in [false, true] {
            for disk in [
                0,
                transaction::CONTEXT_DISK_BYTES,
                transaction::CONTEXT_DISK_BYTES + compressed.len() as u64 - 1,
                transaction::CONTEXT_DISK_BYTES + compressed.len() as u64,
            ] {
                let directory = tempfile::tempdir().unwrap();
                let config = test_config(directory.path());
                let (url, server) = serve_once(gzip_member(RECORD_ONE));
                let data_path = config.data_path(&cache_key(&url));
                let meta_path = config.meta_path(&cache_key(&url));
                let previous = prior.then(|| {
                    create_cache_entry(&config, &url, false, 1);
                    (fs::read(&data_path).unwrap(), fs::read(&meta_path).unwrap())
                });
                let foreign = directory.path().join("foreign.txt");
                fs::write(&foreign, b"unowned data\n").unwrap();
                let budget = cache_http_budget(512, disk);
                let error = config
                    .download_to_cache_bounded(&url, 0, 128, &budget)
                    .unwrap_err();
                server.join().unwrap().unwrap();
                assert_eq!(
                    error.downcast_ref::<ResourceExhausted>().unwrap().resource,
                    ResourceKind::TemporaryDisk,
                    "{error:#}"
                );
                match previous {
                    Some((data, meta)) => {
                        assert_eq!(fs::read(&data_path).unwrap(), data);
                        assert_eq!(fs::read(&meta_path).unwrap(), meta);
                    }
                    None => {
                        assert!(!data_path.exists());
                        assert!(!meta_path.exists());
                    }
                }
                assert_eq!(fs::read(&foreign).unwrap(), b"unowned data\n");
                assert_eq!(
                    fs::read_dir(directory.path()).unwrap().count(),
                    if prior { 4 } else { 2 }
                );
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
                assert!(budget.snapshot().peak_reserved_temporary_disk_bytes <= disk);
                assert!(budget.snapshot().valid);
            }
        }
    }

    #[test]
    fn cache_pair_acceptance_releases_both_closed_temporary_charges() {
        let directory = tempfile::tempdir().unwrap();
        let budget = cache_metadata_budget(512, 1 << 20);
        let url = "http://example.com/reads.fastq.gz";
        let key = cache_key(url);
        let transaction =
            transaction::Transaction::begin(directory.path(), &key, false, Some(&budget)).unwrap();
        let mut data =
            CacheTemporary::create(transaction.directory(), ".new-data-", Some(&budget)).unwrap();
        data.writer().write_all(b"compressed cache bytes").unwrap();
        let selection = data.finish(Ok((1_u64, true))).unwrap();
        assert_eq!(selection, (1, true));
        let meta = CacheMeta {
            url: "http://example.com/reads.fastq.gz".to_owned(),
            sha256: compute_sha256(data.path()).unwrap(),
            complete: true,
            n_reads: 1,
        };
        let metadata = metadata::temporary(transaction.directory(), &meta, Some(&budget)).unwrap();
        let bytes =
            fs::metadata(data.path()).unwrap().len() + fs::metadata(metadata.path()).unwrap().len();
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            bytes + transaction::CONTEXT_DISK_BYTES
        );
        #[cfg(target_os = "linux")]
        for entry in fs::read_dir("/proc/self/fd").unwrap() {
            if let Ok(open) = fs::read_link(entry.unwrap().path()) {
                assert_ne!(open, data.path());
                assert_ne!(open, metadata.path());
            }
        }
        let data_path = directory.path().join(format!("{key}.fastq.gz"));
        let meta_path = directory.path().join(format!("{key}.meta.yaml"));
        transaction.publish(data, metadata).unwrap();
        assert_eq!(fs::read(&data_path).unwrap(), b"compressed cache bytes");
        assert_eq!(read_meta(&meta_path).unwrap().unwrap().sha256, meta.sha256);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            bytes + transaction::CONTROL_PEAK_BYTES
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn bounded_http_cache_refusal_preserves_previous_generation_and_releases_input() {
        use crate::kmer::budget::{MemoryClass, ResourceExhausted, ResourceKind};
        for record_limit in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let config = test_config(directory.path());
            let input = if record_limit {
                format!("@{}\nACGT\n+\n!!!!\n", "x".repeat(4096))
            } else {
                RECORD_ONE.to_owned()
            };
            let (url, server) = serve_once(gzip_member(&input));
            let path = create_cache_entry(&config, &url, false, 1);
            let meta_path = config.meta_path(&cache_key(&url));
            let data = fs::read(&path).unwrap();
            let meta = fs::read(&meta_path).unwrap();
            let foreign = directory.path().join("foreign.txt");
            fs::write(&foreign, "unowned data\n").unwrap();
            let budget = cache_http_budget(if record_limit { 512 } else { 8 }, 1 << 20);
            let error = config
                .download_to_cache_bounded(&url, 0, 128, &budget)
                .unwrap_err();
            server.join().unwrap().unwrap();
            if record_limit {
                assert!(
                    error.to_string().contains("128 byte buffer limit"),
                    "{error:#}"
                );
            } else {
                assert_eq!(
                    error.downcast_ref::<ResourceExhausted>().unwrap().resource,
                    ResourceKind::StageMemory(MemoryClass::InputBuffer)
                );
            }
            assert_eq!(fs::read(&path).unwrap(), data);
            assert_eq!(fs::read(&meta_path).unwrap(), meta);
            assert_eq!(fs::read_to_string(&foreign).unwrap(), "unowned data\n");
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 4);
        }
    }

    fn test_config(cache_directory: &Path) -> CacheConfig {
        CacheConfig::new(Some(cache_directory)).unwrap()
    }

    fn gzip_member(contents: &str) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(contents.as_bytes()).unwrap();
        encoder.finish().unwrap()
    }

    fn serve_once(body: Vec<u8>) -> (String, std::thread::JoinHandle<std::io::Result<()>>) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let server = std::thread::spawn(move || -> std::io::Result<()> {
            let (mut stream, _) = listener.accept()?;
            stream.set_read_timeout(Some(std::time::Duration::from_secs(1)))?;
            let mut reader = std::io::BufReader::new(stream.try_clone()?);
            loop {
                let mut header_line = String::new();
                let bytes_read = BufRead::read_line(&mut reader, &mut header_line)?;
                if bytes_read == 0 || header_line == "\r\n" {
                    break;
                }
            }
            write!(
                stream,
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            )?;
            stream.write_all(&body)?;
            Ok(())
        });
        (format!("http://{address}/reads.fastq.gz"), server)
    }

    #[test]
    fn test_cache_key_deterministic() {
        let url = "http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR571/ERR571460/ERR571460.fastq.gz";
        assert_eq!(cache_key(url), cache_key(url));
    }

    #[test]
    fn test_cache_key_different_urls() {
        let url1 = "http://example.com/reads_1.fastq.gz";
        let url2 = "http://example.com/reads_2.fastq.gz";
        assert_ne!(cache_key(url1), cache_key(url2));
    }

    #[test]
    fn test_cache_key_is_hex_sha256() {
        let key = cache_key("http://example.com/test.fastq.gz");
        assert_eq!(key.len(), 64); // SHA-256 hex = 64 chars
        assert!(key.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn test_meta_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let meta_path = dir.path().join("test.meta.yaml");

        let meta = CacheMeta {
            url: "http://example.com/test.fastq.gz".to_string(),
            sha256: "abc123".to_string(),
            complete: true,
            n_reads: 500,
        };

        write_meta(&meta_path, &meta).unwrap();
        let loaded = read_meta(&meta_path).unwrap().unwrap();

        assert_eq!(loaded.url, meta.url);
        assert_eq!(loaded.sha256, meta.sha256);
        assert_eq!(loaded.complete, meta.complete);
        assert_eq!(loaded.n_reads, meta.n_reads);
    }

    #[test]
    fn test_read_meta_missing_file() {
        let result = read_meta(Path::new("/nonexistent/path.meta.yaml")).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_lookup_miss_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());
        let result = config
            .lookup("http://example.com/test.fastq.gz", 0)
            .unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_compute_sha256() {
        let dir = tempfile::tempdir().unwrap();
        let file_path = dir.path().join("test.txt");
        std::fs::write(&file_path, b"hello world").unwrap();

        let hash = compute_sha256(&file_path).unwrap();
        // Known SHA-256 of "hello world"
        assert_eq!(
            hash,
            "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"
        );
    }

    #[test]
    fn test_lookup_miss_bad_checksum() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());

        let url = "http://example.com/test.fastq.gz";
        let key = cache_key(url);

        // Create a data file and sidecar with wrong checksum
        let data_path = config.data_path(&key);
        let meta_path = config.meta_path(&key);

        std::fs::write(&data_path, b"some data").unwrap();
        let meta = CacheMeta {
            url: url.to_string(),
            sha256: "wrong_checksum".to_string(),
            complete: true,
            n_reads: 0,
        };
        write_meta(&meta_path, &meta).unwrap();

        let error = config.lookup(url, 0).unwrap_err();
        assert!(error.to_string().contains("invalid checksum receipt"));
        assert!(data_path.exists());
        assert!(meta_path.exists());
    }

    #[test]
    fn test_lookup_hit_valid_checksum() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());

        let url = "http://example.com/test.fastq.gz";
        let key = cache_key(url);
        let data_path = config.data_path(&key);
        let meta_path = config.meta_path(&key);

        let data = b"some cached data";
        std::fs::write(&data_path, data).unwrap();
        let sha256 = compute_sha256(&data_path).unwrap();

        let meta = CacheMeta {
            url: url.to_string(),
            sha256,
            complete: true,
            n_reads: 0,
        };
        write_meta(&meta_path, &meta).unwrap();

        let result = config.lookup(url, 0).unwrap();
        assert_eq!(result, Some(data_path));
    }

    /// Helper: create a cache entry with given complete/n_reads values.
    fn create_cache_entry(
        config: &CacheConfig,
        url: &str,
        complete: bool,
        n_reads: u64,
    ) -> PathBuf {
        let key = cache_key(url);
        let data_path = config.data_path(&key);
        let meta_path = config.meta_path(&key);

        let data = b"cached data";
        std::fs::write(&data_path, data).unwrap();
        let sha256 = compute_sha256(&data_path).unwrap();

        let meta = CacheMeta {
            url: url.to_string(),
            sha256,
            complete,
            n_reads,
        };
        write_meta(&meta_path, &meta).unwrap();
        data_path
    }

    #[test]
    fn test_lookup_hit_complete_any_max_reads() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, true, 100);

        // Complete download satisfies any max_reads
        assert_eq!(config.lookup(url, 500).unwrap(), Some(data_path.clone()));
        assert_eq!(config.lookup(url, 0).unwrap(), Some(data_path));
    }

    #[test]
    fn test_lookup_hit_sufficient_reads() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, false, 1000);

        // Partial with 1000 reads satisfies request for 500
        assert_eq!(config.lookup(url, 500).unwrap(), Some(data_path));
    }

    #[test]
    fn test_lookup_miss_insufficient_reads() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, false, 100);
        let key = cache_key(url);
        let meta_path = config.meta_path(&key);

        assert_eq!(config.lookup(url, 500).unwrap(), None);
        assert!(data_path.exists());
        assert!(meta_path.exists());
    }

    #[test]
    fn test_lookup_miss_partial_unlimited() {
        let dir = tempfile::tempdir().unwrap();
        let config = test_config(dir.path());
        let url = "http://example.com/test.fastq.gz";
        create_cache_entry(&config, url, false, 1000);

        // Partial cache never satisfies unlimited (max_reads=0)
        assert_eq!(config.lookup(url, 0).unwrap(), None);
    }

    #[test]
    fn test_meta_backward_compat() {
        // Old sidecar YAML without n_reads should deserialize with n_reads=0
        let yaml = "url: http://example.com/test.fastq.gz\nsha256: abc123\ncomplete: true\n";
        let meta: CacheMeta = serde_yaml_ng::from_str(yaml).unwrap();
        assert_eq!(meta.n_reads, 0);
        assert!(meta.complete);
    }

    #[test]
    fn test_clear_nonexistent_dir() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        CacheConfig::clear(Some(&cache_directory)).unwrap();
        assert!(cache_directory.exists());
        assert!(cache_directory.join(".sharkmer-cache.lock").exists());
    }

    #[test]
    fn clear_removes_only_verified_entries_and_preserves_foreign_paths() {
        let dir = tempfile::tempdir().unwrap();
        let cache_dir = dir.path().join("cache");
        let config = test_config(&cache_dir);
        let url = "http://example.com/clear.fastq.gz";
        let data_path = create_cache_entry(&config, url, true, 1);
        let meta_path = config.meta_path(&cache_key(url));
        let sentinel = cache_dir.join("foreign.txt");
        let nested = cache_dir.join("foreign-dir");
        fs::write(&sentinel, b"keep").unwrap();
        fs::create_dir(&nested).unwrap();
        drop(config);
        CacheConfig::clear(Some(&cache_dir)).unwrap();
        assert!(!data_path.exists());
        assert!(!meta_path.exists());
        assert_eq!(fs::read(&sentinel).unwrap(), b"keep");
        assert!(nested.exists());
        assert!(cache_dir.join(".sharkmer-cache.lock").is_file());
    }

    #[test]
    fn clear_accepts_legacy_verified_metadata_without_read_count() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        let config = test_config(&cache_directory);
        let url = "http://example.com/legacy.fastq.gz";
        let key = cache_key(url);
        let data_path = config.data_path(&key);
        fs::write(&data_path, b"legacy cache data").unwrap();
        let checksum = compute_sha256(&data_path).unwrap();
        let meta_path = config.meta_path(&key);
        fs::write(
            &meta_path,
            format!("url: {url}\nsha256: {checksum}\ncomplete: true\n"),
        )
        .unwrap();
        drop(config);

        CacheConfig::clear(Some(&cache_directory)).unwrap();

        assert!(!data_path.exists());
        assert!(!meta_path.exists());
    }

    #[test]
    fn clear_preserves_modified_and_orphaned_entries() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        let config = test_config(&cache_directory);
        let modified_url = "http://example.com/modified.fastq.gz";
        let modified_path = create_cache_entry(&config, modified_url, true, 1);
        fs::write(&modified_path, b"modified").unwrap();
        let orphan_url = "http://example.com/orphan.fastq.gz";
        let orphan_path = config.data_path(&cache_key(orphan_url));
        fs::write(&orphan_path, b"orphan").unwrap();
        let orphan_meta_url = "http://example.com/orphan-meta.fastq.gz";
        let orphan_meta_path = config.meta_path(&cache_key(orphan_meta_url));
        fs::write(&orphan_meta_path, b"orphan metadata").unwrap();
        drop(config);

        CacheConfig::clear(Some(&cache_directory)).unwrap();

        assert!(modified_path.exists());
        assert!(orphan_path.exists());
        assert!(orphan_meta_path.exists());
    }

    #[test]
    fn mismatched_metadata_url_is_preserved_and_refused() {
        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let url = "http://example.com/requested.fastq.gz";
        let key = cache_key(url);
        let data_path = config.data_path(&key);
        fs::write(&data_path, b"foreign cache data").unwrap();
        let checksum = compute_sha256(&data_path).unwrap();
        let meta_path = config.meta_path(&key);
        write_meta(
            &meta_path,
            &CacheMeta {
                url: "http://example.com/other.fastq.gz".to_string(),
                sha256: checksum,
                complete: true,
                n_reads: 1,
            },
        )
        .unwrap();

        let error = config.lookup(url, 1).unwrap_err();

        assert!(error.to_string().contains("does not own key"));
        assert!(data_path.exists());
        assert!(meta_path.exists());
    }

    #[test]
    fn cache_lease_stays_live_for_replay() {
        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let url = "http://example.com/replay.fastq.gz";
        let data_path = create_cache_entry(&config, url, false, 1);
        let first_pass = config.lookup(url, 1).unwrap().unwrap();
        let first_bytes = fs::read(&first_pass).unwrap();
        let second_pass = config.lookup(url, 1).unwrap().unwrap();

        assert_eq!(data_path, first_pass);
        assert_eq!(first_pass, second_pass);
        assert_eq!(first_bytes, fs::read(&second_pass).unwrap());
    }

    #[test]
    fn cache_lease_blocks_a_second_configuration_until_drop() {
        use std::sync::mpsc;

        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let cache_directory = config.cache_dir.clone();
        let (sender, receiver) = mpsc::channel();
        let handle = std::thread::spawn(move || {
            sender.send("waiting").unwrap();
            let second = CacheConfig::new(Some(&cache_directory)).unwrap();
            sender.send("locked").unwrap();
            drop(second);
        });
        assert_eq!(receiver.recv().unwrap(), "waiting");
        assert!(receiver.try_recv().is_err());
        drop(config);
        assert_eq!(
            receiver
                .recv_timeout(std::time::Duration::from_secs(1))
                .unwrap(),
            "locked"
        );
        handle.join().unwrap();
    }

    #[test]
    fn clear_waits_for_the_active_cache_lease() {
        use std::sync::mpsc;

        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let cache_directory = config.cache_dir.clone();
        let (sender, receiver) = mpsc::channel();
        let handle = std::thread::spawn(move || {
            sender.send("waiting").unwrap();
            CacheConfig::clear(Some(&cache_directory)).unwrap();
            sender.send("cleared").unwrap();
        });
        assert_eq!(receiver.recv().unwrap(), "waiting");
        assert!(receiver.try_recv().is_err());
        drop(config);
        assert_eq!(
            receiver
                .recv_timeout(std::time::Duration::from_secs(1))
                .unwrap(),
            "cleared"
        );
        handle.join().unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn cache_aliases_share_the_same_lock() {
        use std::os::unix::fs::symlink;

        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        let alias_directory = directory.path().join("alias");
        let config = test_config(&cache_directory);
        symlink(&cache_directory, &alias_directory).unwrap();
        let alias_lock = OpenOptions::new()
            .read(true)
            .write(true)
            .open(alias_directory.join(".sharkmer-cache.lock"))
            .unwrap();

        assert!(alias_lock.try_lock_exclusive().is_err());
        drop(config);
        alias_lock.try_lock_exclusive().unwrap();
        alias_lock.unlock().unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn symlinked_cache_entry_is_preserved_and_refused() {
        use std::os::unix::fs::symlink;

        let directory = tempfile::tempdir().unwrap();
        let config = test_config(directory.path());
        let url = "http://example.com/symlink.fastq.gz";
        let key = cache_key(url);
        let target = directory.path().join("target");
        fs::write(&target, b"keep").unwrap();
        let data_path = config.data_path(&key);
        symlink(&target, &data_path).unwrap();

        let error = config.lookup(url, 1).unwrap_err();
        assert!(error.to_string().contains("non-regular cache data"));
        assert_eq!(fs::read(&target).unwrap(), b"keep");
    }

    #[test]
    fn corrupt_second_gzip_member_is_not_published_to_cache() {
        let dir = tempfile::tempdir().unwrap();
        let cache_dir = dir.path().join("cache");
        let config = test_config(&cache_dir);
        let mut corrupt_second_member = gzip_member(RECORD_TWO);
        *corrupt_second_member.last_mut().unwrap() ^= 0xff;
        let mut body = gzip_member(RECORD_ONE);
        body.extend(corrupt_second_member);
        let (url, server) = serve_once(body);

        let result = config.download_to_cache(&url, 0);
        server.join().unwrap().unwrap();
        assert!(result.is_err());
        let key = cache_key(&url);
        assert!(!config.data_path(&key).exists());
        assert!(!config.meta_path(&key).exists());
    }

    #[test]
    fn failed_replacement_keeps_the_previous_verified_generation() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        let config = test_config(&cache_directory);
        let mut corrupt_second_member = gzip_member(RECORD_TWO);
        *corrupt_second_member.last_mut().unwrap() ^= 0xff;
        let mut body = gzip_member(RECORD_ONE);
        body.extend(corrupt_second_member);
        let (url, server) = serve_once(body);
        let old_path = create_cache_entry(&config, &url, false, 1);
        let old_data = fs::read(&old_path).unwrap();
        let old_meta = fs::read(config.meta_path(&cache_key(&url))).unwrap();

        assert!(config.download_to_cache(&url, 0).is_err());
        server.join().unwrap().unwrap();

        assert_eq!(fs::read(&old_path).unwrap(), old_data);
        assert_eq!(
            fs::read(config.meta_path(&cache_key(&url))).unwrap(),
            old_meta
        );
    }

    #[test]
    fn foreign_metadata_destination_is_preserved_before_publication() {
        let directory = tempfile::tempdir().unwrap();
        let url = "http://example.com/new.fastq.gz";
        let key = cache_key(url);
        let transaction =
            transaction::Transaction::begin(directory.path(), &key, false, None).unwrap();
        let data_path = directory.path().join(format!("{key}.fastq.gz"));
        let meta_path = directory.path().join(format!("{key}.meta.yaml"));
        let mut data = CacheTemporary::create(transaction.directory(), ".new-data-", None).unwrap();
        data.writer().write_all(b"new cache data").unwrap();
        data.finish(Ok(())).unwrap();
        let meta = metadata::temporary(
            transaction.directory(),
            &CacheMeta {
                url: url.to_owned(),
                sha256: compute_sha256(data.path()).unwrap(),
                complete: true,
                n_reads: 1,
            },
            None,
        )
        .unwrap();
        fs::create_dir(&meta_path).unwrap();
        let error = transaction.publish(data, meta).unwrap_err();
        assert!(format!("{error:#}").contains("Non-regular"));
        assert!(!data_path.exists());
        assert!(meta_path.is_dir());
    }

    #[test]
    fn successful_replacement_updates_a_valid_insufficient_entry() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        let config = test_config(&cache_directory);
        let mut body = gzip_member(RECORD_ONE);
        body.extend(gzip_member(RECORD_TWO));
        let (url, server) = serve_once(body);
        let old_path = create_cache_entry(&config, &url, false, 1);
        let old_data = fs::read(&old_path).unwrap();

        let new_path = config.download_to_cache(&url, 2).unwrap();
        server.join().unwrap().unwrap();
        let meta = read_meta(&config.meta_path(&cache_key(&url)))
            .unwrap()
            .unwrap();

        assert_eq!(old_path, new_path);
        assert_ne!(fs::read(&new_path).unwrap(), old_data);
        assert_eq!(meta.n_reads, 2);
        assert!(config.lookup(&url, 2).unwrap().is_some());
    }

    #[test]
    fn cache_download_limit_does_not_read_an_unneeded_corrupt_member() {
        let dir = tempfile::tempdir().unwrap();
        let cache_dir = dir.path().join("cache");
        let config = test_config(&cache_dir);
        let mut corrupt_second_member = gzip_member(RECORD_TWO);
        *corrupt_second_member.last_mut().unwrap() ^= 0xff;
        let mut body = gzip_member(RECORD_ONE);
        body.extend(corrupt_second_member);
        let (url, server) = serve_once(body);

        let data_path = config.download_to_cache(&url, 1).unwrap();
        server.join().unwrap().unwrap();
        let meta_path = config.meta_path(&cache_key(&url));
        let meta = read_meta(&meta_path).unwrap().unwrap();

        assert_eq!(config.lookup(&url, 1).unwrap(), Some(data_path));
        assert_eq!(meta.n_reads, 1);
        assert!(!meta.complete);
    }
}
