use anyhow::{Context, Result};
use log::{info, warn};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
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

/// Configuration and operations for the read cache.
pub(crate) struct CacheConfig {
    pub(crate) cache_dir: PathBuf,
}

impl CacheConfig {
    /// Create a new cache configuration.
    /// Uses `cache_dir_override` if provided, otherwise the platform default
    /// (`~/.cache/sharkmer/reads/` on Linux, `~/Library/Caches/sharkmer/reads/` on macOS).
    pub(crate) fn new(cache_dir_override: Option<&Path>) -> Result<Self> {
        let cache_dir = if let Some(dir) = cache_dir_override {
            dir.to_path_buf()
        } else {
            let base = dirs::cache_dir().context("Could not determine platform cache directory")?;
            base.join("sharkmer").join("reads")
        };
        Ok(CacheConfig { cache_dir })
    }

    /// Look up a URL in the cache.
    /// Returns `Some(path)` if a valid cached file exists with enough reads,
    /// `None` on cache miss. A complete download (EOF reached) satisfies any
    /// `max_reads` value. A partial download satisfies only if it contains
    /// at least `max_reads` reads.
    pub(crate) fn lookup(&self, url: &str, max_reads: u64) -> Result<Option<PathBuf>> {
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);

        if !data_path.exists() {
            return Ok(None);
        }

        // Read sidecar
        let meta = match read_meta(&meta_path) {
            Ok(Some(m)) => m,
            Ok(None) => {
                warn!(
                    "Cache data file exists without sidecar, treating as miss: {}",
                    data_path.display()
                );
                return Ok(None);
            }
            Err(e) => {
                warn!("Failed to read cache sidecar, treating as miss: {}", e);
                return Ok(None);
            }
        };

        // Verify checksum
        let actual_sha256 = compute_sha256(&data_path)
            .with_context(|| format!("Failed to compute SHA-256 of {}", data_path.display()))?;

        if actual_sha256 != meta.sha256 {
            warn!(
                "Cache checksum mismatch for {} (expected {}, got {}), re-downloading",
                url, meta.sha256, actual_sha256
            );
            evict_stale(&data_path, &meta_path);
            return Ok(None);
        }

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
                    "Cache has {} reads but {} requested, re-downloading {}",
                    meta.n_reads,
                    if max_reads == 0 {
                        "all".to_string()
                    } else {
                        max_reads.to_string()
                    },
                    url
                );
                evict_stale(&data_path, &meta_path);
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

        // Ensure cache directory exists
        std::fs::create_dir_all(&self.cache_dir).with_context(|| {
            format!(
                "Failed to create cache directory: {}",
                self.cache_dir.display()
            )
        })?;

        // Create a unique temp file in the cache directory. Using tempfile
        // gives us a random-suffixed name (avoiding collisions between
        // concurrent sharkmer invocations downloading the same URL) and
        // auto-deletes on error paths. Keeping it in the cache dir ensures
        // the final rename stays on the same filesystem and is atomic.
        let tmp = tempfile::Builder::new()
            .prefix(&format!("{}.", key))
            .suffix(".gz.tmp")
            .tempfile_in(&self.cache_dir)
            .with_context(|| {
                format!(
                    "Failed to create temp file in cache directory: {}",
                    self.cache_dir.display()
                )
            })?;
        let (tmp_file, tmp_path) = tmp.into_parts();

        info!(
            "Downloading {} to cache (max_reads: {})...",
            url,
            if max_reads == 0 {
                "unlimited".to_string()
            } else {
                max_reads.to_string()
            }
        );

        // Input pipeline: HTTP → GzDecoder → BufReader
        let response = ureq::get(url)
            .call()
            .with_context(|| format!("Failed to download {}", url))?;
        let gz_reader = flate2::read::GzDecoder::new(response.into_reader());
        let buf_reader = std::io::BufReader::new(gz_reader);
        let mut lines = buf_reader.lines();

        // Output pipeline: File → BufWriter → GzEncoder
        let mut gz_writer = flate2::write::GzEncoder::new(
            std::io::BufWriter::new(tmp_file),
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
        let sha256 = compute_sha256(&tmp_path)?;

        // Atomic rename (same filesystem — tmp_path was created in cache_dir).
        // persist() consumes the TempPath guard so its auto-delete drop no
        // longer fires on the now-renamed file.
        tmp_path
            .persist(&data_path)
            .with_context(|| format!("Failed to persist cache file to {}", data_path.display()))?;

        // Write sidecar
        let meta = CacheMeta {
            url: url.to_string(),
            sha256,
            complete,
            n_reads,
        };
        write_meta(&meta_path, &meta)?;

        Ok(data_path)
    }

    /// Delete the entire cache directory.
    pub(crate) fn clear(cache_dir_override: Option<&Path>) -> Result<()> {
        let cache_dir = if let Some(dir) = cache_dir_override {
            dir.to_path_buf()
        } else {
            let base = dirs::cache_dir().context("Could not determine platform cache directory")?;
            base.join("sharkmer").join("reads")
        };

        if cache_dir.exists() {
            info!("Removing cache directory: {}", cache_dir.display());
            std::fs::remove_dir_all(&cache_dir).with_context(|| {
                format!("Failed to remove cache directory: {}", cache_dir.display())
            })?;
            info!("Cache cleared.");
        } else {
            info!(
                "Cache directory does not exist, nothing to clear: {}",
                cache_dir.display()
            );
        }
        Ok(())
    }

    fn data_path(&self, key: &str) -> PathBuf {
        self.cache_dir.join(format!("{}.fastq.gz", key))
    }

    fn meta_path(&self, key: &str) -> PathBuf {
        self.cache_dir.join(format!("{}.meta.yaml", key))
    }
}

/// Compute a deterministic cache key from a URL (hex-encoded SHA-256).
/// Remove a stale cache entry (data file and sidecar). Logs a warning if
/// removal fails, but does not error — the next download attempt will
/// overwrite the stale file, and failing silently here would leave the
/// cache in an inconsistent state without any user-visible signal.
fn evict_stale(data_path: &Path, meta_path: &Path) {
    for path in [data_path, meta_path] {
        if let Err(e) = std::fs::remove_file(path) {
            if e.kind() != std::io::ErrorKind::NotFound {
                warn!(
                    "Failed to remove stale cache entry {}: {}",
                    path.display(),
                    e
                );
            }
        }
    }
}

fn cache_key(url: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(url.as_bytes());
    format!("{:x}", hasher.finalize())
}

/// Compute the SHA-256 hex digest of a file.
fn compute_sha256(path: &Path) -> Result<String> {
    let mut file = std::fs::File::open(path)
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

/// Read sidecar metadata from a JSON file.
fn read_meta(path: &Path) -> Result<Option<CacheMeta>> {
    if !path.exists() {
        return Ok(None);
    }
    let contents = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read sidecar {}", path.display()))?;
    let meta: CacheMeta = serde_yaml_ng::from_str(&contents)
        .map_err(|e| anyhow::anyhow!("Failed to parse sidecar {}: {}", path.display(), e))?;
    Ok(Some(meta))
}

/// Write sidecar metadata to a JSON file.
fn write_meta(path: &Path, meta: &CacheMeta) -> Result<()> {
    let contents = serde_yaml_ng::to_string(meta)
        .map_err(|e| anyhow::anyhow!("Failed to serialize cache sidecar: {}", e))?;
    std::fs::write(path, contents)
        .with_context(|| format!("Failed to write sidecar {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };
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
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };

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

        // Lookup should return None (checksum mismatch)
        let result = config.lookup(url, 0).unwrap();
        assert!(result.is_none());

        // Files should be cleaned up
        assert!(!data_path.exists());
        assert!(!meta_path.exists());
    }

    #[test]
    fn test_lookup_hit_valid_checksum() {
        let dir = tempfile::tempdir().unwrap();
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };

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
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, true, 100);

        // Complete download satisfies any max_reads
        assert_eq!(config.lookup(url, 500).unwrap(), Some(data_path.clone()));
        assert_eq!(config.lookup(url, 0).unwrap(), Some(data_path));
    }

    #[test]
    fn test_lookup_hit_sufficient_reads() {
        let dir = tempfile::tempdir().unwrap();
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, false, 1000);

        // Partial with 1000 reads satisfies request for 500
        assert_eq!(config.lookup(url, 500).unwrap(), Some(data_path));
    }

    #[test]
    fn test_lookup_miss_insufficient_reads() {
        let dir = tempfile::tempdir().unwrap();
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };
        let url = "http://example.com/test.fastq.gz";
        let data_path = create_cache_entry(&config, url, false, 100);
        let key = cache_key(url);
        let meta_path = config.meta_path(&key);

        // Partial with 100 reads does not satisfy request for 500
        assert_eq!(config.lookup(url, 500).unwrap(), None);
        // Stale files should be cleaned up
        assert!(!data_path.exists());
        assert!(!meta_path.exists());
    }

    #[test]
    fn test_lookup_miss_partial_unlimited() {
        let dir = tempfile::tempdir().unwrap();
        let config = CacheConfig {
            cache_dir: dir.path().to_path_buf(),
        };
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
        // Should succeed without error
        CacheConfig::clear(Some(Path::new("/tmp/sharkmer_test_nonexistent_cache_dir"))).unwrap();
    }

    #[test]
    fn test_clear_existing_dir() {
        let dir = tempfile::tempdir().unwrap();
        let cache_dir = dir.path().join("cache");
        std::fs::create_dir_all(&cache_dir).unwrap();
        std::fs::write(cache_dir.join("test.txt"), b"data").unwrap();

        CacheConfig::clear(Some(&cache_dir)).unwrap();
        assert!(!cache_dir.exists());
    }
}
