use anyhow::{Context, Result};
use log::{info, warn};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::{Path, PathBuf};

/// Sidecar metadata stored alongside each cached file.
#[derive(Serialize, Deserialize, Debug)]
struct CacheMeta {
    /// The original URL this file was downloaded from.
    url: String,
    /// SHA-256 hex digest of the cached data file.
    sha256: String,
    /// Whether the download completed fully (true) or was interrupted (false).
    complete: bool,
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
    /// Returns `Some(path)` if a valid cached file exists, `None` on cache miss.
    pub(crate) fn lookup(&self, url: &str) -> Result<Option<PathBuf>> {
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
            // Remove stale cache entry
            let _ = std::fs::remove_file(&data_path);
            let _ = std::fs::remove_file(&meta_path);
            return Ok(None);
        }

        Ok(Some(data_path))
    }

    /// Download a URL to the cache directory. Returns the path to the cached file.
    pub(crate) fn download_to_cache(&self, url: &str) -> Result<PathBuf> {
        let key = cache_key(url);
        let data_path = self.data_path(&key);
        let meta_path = self.meta_path(&key);
        let tmp_path = data_path.with_extension("gz.tmp");

        // Ensure cache directory exists
        std::fs::create_dir_all(&self.cache_dir).with_context(|| {
            format!(
                "Failed to create cache directory: {}",
                self.cache_dir.display()
            )
        })?;

        info!("Downloading {} to cache...", url);

        // Download to temp file
        let response = ureq::get(url)
            .call()
            .with_context(|| format!("Failed to download {}", url))?;

        let mut reader = response.into_reader();
        let mut tmp_file = std::fs::File::create(&tmp_path)
            .with_context(|| format!("Failed to create temp file: {}", tmp_path.display()))?;

        let bytes_written = std::io::copy(&mut reader, &mut tmp_file)
            .with_context(|| format!("Failed to write cached data for {}", url))?;

        drop(tmp_file);

        info!(
            "Downloaded {} bytes for {}",
            crate::format::format_bytes(bytes_written),
            url
        );

        // Compute SHA-256 of the downloaded file
        let sha256 = compute_sha256(&tmp_path)?;

        // Atomic rename
        std::fs::rename(&tmp_path, &data_path).with_context(|| {
            format!(
                "Failed to rename {} to {}",
                tmp_path.display(),
                data_path.display()
            )
        })?;

        // Write sidecar
        let meta = CacheMeta {
            url: url.to_string(),
            sha256,
            complete: true,
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
    let meta: CacheMeta = serde_yml::from_str(&contents)
        .map_err(|e| anyhow::anyhow!("Failed to parse sidecar {}: {}", path.display(), e))?;
    Ok(Some(meta))
}

/// Write sidecar metadata to a JSON file.
fn write_meta(path: &Path, meta: &CacheMeta) -> Result<()> {
    let contents = serde_yml::to_string(meta)
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
        };

        write_meta(&meta_path, &meta).unwrap();
        let loaded = read_meta(&meta_path).unwrap().unwrap();

        assert_eq!(loaded.url, meta.url);
        assert_eq!(loaded.sha256, meta.sha256);
        assert_eq!(loaded.complete, meta.complete);
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
        let result = config.lookup("http://example.com/test.fastq.gz").unwrap();
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
        };
        write_meta(&meta_path, &meta).unwrap();

        // Lookup should return None (checksum mismatch)
        let result = config.lookup(url).unwrap();
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
        };
        write_meta(&meta_path, &meta).unwrap();

        let result = config.lookup(url).unwrap();
        assert_eq!(result, Some(data_path));
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
