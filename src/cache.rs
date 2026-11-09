use anyhow::{Context, Result};
use fs2::FileExt;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};

/// Sidecar metadata stored alongside each cached file.
#[derive(Serialize, Deserialize, Debug)]
struct CacheMeta {
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

struct CacheBackup {
    data_path: tempfile::TempPath,
    meta_path: tempfile::TempPath,
}

impl CacheBackup {
    fn restore(self, data_path: &Path, meta_path: &Path) -> Result<()> {
        self.meta_path
            .persist(meta_path)
            .map_err(|error| error.error)
            .with_context(|| format!("Failed to restore cache metadata {}", meta_path.display()))?;
        self.data_path
            .persist(data_path)
            .map_err(|error| error.error)
            .with_context(|| format!("Failed to restore cache data {}", data_path.display()))?;
        Ok(())
    }
}

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
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);
        let meta = match validate_entry(&key, &data_path, &meta_path)? {
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
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);
        validate_entry(&key, &data_path, &meta_path)?;

        // Create a unique temp file in the cache directory. Using tempfile
        // gives us a random-suffixed name (avoiding collisions between
        // concurrent sharkmer invocations downloading the same URL) and
        // auto-deletes on error paths. Keeping it in the cache dir ensures
        // the final rename stays on the same filesystem and is atomic.
        let temporary = tempfile::Builder::new()
            .prefix(&format!("{}.", key))
            .suffix(".gz.tmp")
            .tempfile_in(&self.cache_dir)
            .with_context(|| {
                format!(
                    "Failed to create temp file in cache directory: {}",
                    self.cache_dir.display()
                )
            })?;
        let (temporary_file, temporary_path) = temporary.into_parts();

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
        let response = ureq::get(url)
            .call()
            .with_context(|| format!("Failed to download {}", url))?;
        let gz_reader = flate2::read::MultiGzDecoder::new(response.into_reader());
        let buf_reader = std::io::BufReader::new(gz_reader);
        let mut lines = buf_reader.lines();

        // Output pipeline: File → BufWriter → GzEncoder
        let mut gz_writer = flate2::write::GzEncoder::new(
            std::io::BufWriter::new(temporary_file),
            flate2::Compression::fast(),
        );

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
                .ok_or_else(|| {
                    anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url)
                })?
                .with_context(|| format!("Failed reading FASTQ from {}", url))?;
            let separator = lines
                .next()
                .ok_or_else(|| {
                    anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url)
                })?
                .with_context(|| format!("Failed reading FASTQ from {}", url))?;
            let quality = lines
                .next()
                .ok_or_else(|| {
                    anyhow::anyhow!("Truncated FASTQ at record {} in {}", n_reads + 1, url)
                })?
                .with_context(|| format!("Failed reading FASTQ from {}", url))?;

            // Write the 4 lines to the cache file
            writeln!(gz_writer, "{}", header)?;
            writeln!(gz_writer, "{}", sequence)?;
            writeln!(gz_writer, "{}", separator)?;
            writeln!(gz_writer, "{}", quality)?;

            n_reads += 1;

            if max_reads > 0 && n_reads >= max_reads {
                break;
            }
        }

        // Finish gzip stream (writes trailer) and flush to disk
        let buf_writer = gz_writer
            .finish()
            .context("Failed to finish gzip stream for cache file")?;
        buf_writer
            .into_inner()
            .map_err(|e| anyhow::anyhow!("Failed to flush cache file: {}", e))?;

        info!(
            "Cached {} reads for {} (complete: {})",
            n_reads, url, complete
        );

        // Compute SHA-256 of the cached file
        let sha256 = compute_sha256(&temporary_path)?;
        let meta = CacheMeta {
            url: url.to_string(),
            sha256,
            complete,
            n_reads,
        };
        let meta_temporary = write_meta_temporary(&self.cache_dir, &meta)?;
        let existing_meta = validate_entry(&key, &data_path, &meta_path)?;
        let backup = if existing_meta.is_some() {
            Some(backup_entry(&self.cache_dir, &data_path, &meta_path)?)
        } else {
            None
        };
        publish_entry(
            temporary_path,
            meta_temporary,
            &data_path,
            &meta_path,
            backup,
        )?;

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

fn read_meta(path: &Path) -> Result<Option<CacheMeta>> {
    let contents = match fs::read_to_string(path) {
        Ok(contents) => contents,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("Failed to read sidecar {}", path.display()));
        }
    };
    let meta: CacheMeta = serde_yaml_ng::from_str(&contents)
        .map_err(|e| anyhow::anyhow!("Failed to parse sidecar {}: {}", path.display(), e))?;
    Ok(Some(meta))
}

#[cfg(test)]
fn write_meta(path: &Path, meta: &CacheMeta) -> Result<()> {
    let contents = serde_yaml_ng::to_string(meta)
        .map_err(|e| anyhow::anyhow!("Failed to serialize cache sidecar: {}", e))?;
    fs::write(path, contents)
        .with_context(|| format!("Failed to write sidecar {}", path.display()))?;
    Ok(())
}

fn write_meta_temporary(cache_dir: &Path, meta: &CacheMeta) -> Result<tempfile::NamedTempFile> {
    let mut temporary = tempfile::Builder::new()
        .prefix(".sharkmer-cache-meta-")
        .suffix(".tmp")
        .tempfile_in(cache_dir)
        .with_context(|| {
            format!(
                "Failed to create temporary cache metadata in {}",
                cache_dir.display()
            )
        })?;
    serde_yaml_ng::to_writer(temporary.as_file_mut(), meta)
        .context("Failed to write temporary cache metadata")?;
    temporary
        .as_file_mut()
        .flush()
        .context("Failed to flush temporary cache metadata")?;
    temporary
        .as_file_mut()
        .sync_all()
        .context("Failed to sync temporary cache metadata")?;
    Ok(temporary)
}

fn backup_entry(cache_dir: &Path, data_path: &Path, meta_path: &Path) -> Result<CacheBackup> {
    Ok(CacheBackup {
        data_path: hard_link_temporary(cache_dir, data_path, ".sharkmer-cache-data-backup-")?,
        meta_path: hard_link_temporary(cache_dir, meta_path, ".sharkmer-cache-meta-backup-")?,
    })
}

fn hard_link_temporary(
    cache_dir: &Path,
    source_path: &Path,
    prefix: &str,
) -> Result<tempfile::TempPath> {
    let temporary = tempfile::Builder::new()
        .prefix(prefix)
        .tempfile_in(cache_dir)
        .with_context(|| format!("Failed to create cache backup in {}", cache_dir.display()))?;
    let temporary_path = temporary.into_temp_path();
    fs::remove_file(&temporary_path).with_context(|| {
        format!(
            "Failed to prepare cache backup {}",
            temporary_path.display()
        )
    })?;
    fs::hard_link(source_path, &temporary_path).with_context(|| {
        format!(
            "Failed to link cache backup from {} to {}",
            source_path.display(),
            temporary_path.display()
        )
    })?;
    Ok(temporary_path)
}

fn publish_entry(
    data_temporary: tempfile::TempPath,
    meta_temporary: tempfile::NamedTempFile,
    data_path: &Path,
    meta_path: &Path,
    backup: Option<CacheBackup>,
) -> Result<()> {
    data_temporary
        .persist(data_path)
        .with_context(|| format!("Failed to persist cache file to {}", data_path.display()))?;
    if let Err(error) = meta_temporary
        .persist(meta_path)
        .map_err(|error| error.error)
        .with_context(|| format!("Failed to persist cache sidecar to {}", meta_path.display()))
    {
        if let Some(backup) = backup {
            backup.restore(data_path, meta_path)?;
        } else {
            fs::remove_file(data_path).with_context(|| {
                format!(
                    "Failed to remove newly published cache data after metadata failure {}",
                    data_path.display()
                )
            })?;
        }
        return Err(error);
    }
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

fn validate_entry(key: &str, data_path: &Path, meta_path: &Path) -> Result<Option<CacheMeta>> {
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
    let meta = read_meta(meta_path)
        .map_err(|error| {
            anyhow::anyhow!(
                "Refusing malformed cache entry for {}. Use a fresh cache directory or --no-cache: {error:#}",
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
    fn fresh_entry_metadata_publish_failure_removes_new_data() {
        let directory = tempfile::tempdir().unwrap();
        let cache_directory = directory.path().join("cache");
        fs::create_dir(&cache_directory).unwrap();
        let data_path = cache_directory.join("new.fastq.gz");
        let meta_path = cache_directory.join("blocked.meta.yaml");
        let mut data_temporary_file = tempfile::Builder::new()
            .tempfile_in(&cache_directory)
            .unwrap();
        data_temporary_file.write_all(b"new cache data").unwrap();
        data_temporary_file.flush().unwrap();
        let data_temporary = data_temporary_file.into_temp_path();
        let meta_temporary = write_meta_temporary(
            &cache_directory,
            &CacheMeta {
                url: "http://example.com/new.fastq.gz".to_string(),
                sha256: "0".repeat(64),
                complete: true,
                n_reads: 1,
            },
        )
        .unwrap();
        fs::create_dir(&meta_path).unwrap();

        let error = publish_entry(data_temporary, meta_temporary, &data_path, &meta_path, None)
            .unwrap_err();

        assert!(error.to_string().contains("persist cache sidecar"));
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
