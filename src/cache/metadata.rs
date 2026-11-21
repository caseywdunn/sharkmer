//! Bounded immutable cache sidecars and their parser/emitter ownership.

use super::*;
use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::{MemoryClass, MemoryReservation};
use anyhow::ensure;

// A cache record contains one URL and one digest. Keep the same generous byte
// ceiling as the ENA client's bounded response, and enforce it on both ends.
const MAX_BYTES: usize = 10 * 1024 * 1024;

/// Typed strings drop before the reservation for their surviving capacities.
pub(super) struct Owned {
    data: CacheMeta,
    _owner: Option<MemoryReservation>,
}
impl std::ops::Deref for Owned {
    type Target = CacheMeta;
    fn deref(&self) -> &Self::Target {
        &self.data
    }
}
impl Owned {
    pub(super) fn generated(
        url: &str,
        sha256: String,
        complete: bool,
        n_reads: u64,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Self> {
        ensure!(
            url.len() <= MAX_BYTES,
            "Cache metadata URL exceeds the sidecar byte limit"
        );
        let length = url
            .len()
            .checked_add(sha256.capacity())
            .context("Cache metadata size overflow")?;
        let owner = budget
            .map(|budget| budget.reserve(MemoryClass::RetainedResults, u64::try_from(length)?))
            .transpose()?;
        let mut text = String::new();
        text.try_reserve_exact(url.len())
            .context("Unable to allocate cache metadata URL")?;
        text.push_str(url);
        ensure!(
            text.capacity() <= url.len(),
            "Cache URL exceeded its admitted capacity"
        );
        Ok(Self {
            data: CacheMeta {
                url: text,
                sha256,
                complete,
                n_reads,
            },
            _owner: owner,
        })
    }
}

pub(super) fn read(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<Option<Owned>> {
    let file = match File::open(path) {
        Ok(file) => file,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("Failed to read sidecar {}", path.display()));
        }
    };
    let length = usize::try_from(file.metadata()?.len())?;
    ensure!(
        length <= MAX_BYTES,
        "Cache sidecar exceeds the {MAX_BYTES} byte limit: {}",
        path.display()
    );
    let mut input =
        BudgetVec::with_capacity_in_class(length, budget, MemoryClass::RetainedResults)?;
    let mut reader = file.take(MAX_BYTES as u64 + 1);
    let mut buffer = [0; 8192];
    loop {
        let n = reader.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        ensure!(
            input
                .len()
                .checked_add(n)
                .is_some_and(|length| length <= MAX_BYTES),
            "Cache sidecar grew beyond the {MAX_BYTES} byte limit"
        );
        for byte in &buffer[..n] {
            input.push(*byte)?;
        }
    }
    drop(reader);
    let _workspace = crate::output::parsing::workspace::Workspace::new(&input, budget)?;
    // The only dynamically owned model fields are two scalar Strings. Even if
    // both alias the same scalar, their total decoded bytes fit this bound.
    let bound = input
        .len()
        .checked_add(1)
        .and_then(|n| n.checked_mul(2))
        .context("Cache metadata scalar bound overflow")?;
    let mut owner = budget
        .map(|budget| budget.reserve(MemoryClass::RetainedResults, u64::try_from(bound)?))
        .transpose()?;
    let data: CacheMeta = serde_yaml_ng::from_slice(&input)
        .with_context(|| format!("Failed to parse sidecar {}", path.display()))?;
    let actual = data
        .url
        .capacity()
        .checked_add(data.sha256.capacity())
        .context("Cache metadata capacity overflow")?;
    ensure!(
        actual <= bound,
        "Decoded cache metadata exceeded its admitted bound"
    );
    if let Some(owner) = &mut owner {
        owner.resize(u64::try_from(actual)?)?;
    }
    Ok(Some(Owned {
        data,
        _owner: owner,
    }))
}

struct CountBytes(usize);
impl Write for CountBytes {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let length = self
            .0
            .checked_add(bytes.len())
            .ok_or_else(|| std::io::Error::other("Cache metadata byte count overflow"))?;
        if length > MAX_BYTES {
            return Err(std::io::Error::other(
                "Cache metadata exceeds the sidecar byte limit",
            ));
        }
        self.0 = length;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

pub(super) fn temporary(
    directory: &Path,
    meta: &CacheMeta,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<crate::output::MetadataTemporary> {
    let yaml = crate::stats::yaml::prepare_cache(meta, budget)?;
    // Both passes borrow the same immutable, derived model. Refuse oversized
    // generated metadata before even creating a temporary file.
    yaml.write(&mut CountBytes(0))
        .context("Unable to admit cache metadata bytes")?;
    let mut temporary =
        crate::output::MetadataTemporary::create(directory, ".sharkmer-cache-meta-", budget)?;
    let written = yaml.write(temporary.writer());
    temporary
        .finish(written)
        .context("Failed to write temporary cache metadata")?;
    Ok(temporary)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};

    fn budget(disk: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: disk,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn cached_model_preserves_legacy_defaults_escapes_aliases_and_string_ownership() {
        for input in [
            "url: http://example.com/source\nsha256: '0000'\ncomplete: false\n",
            "url: &value 'http://example.com/shared'\nsha256: *value\ncomplete: true\nn_reads: 123\nignored: [1, 2]\n",
            "url: \"http://example.com/\\u03b1\"\nsha256: 'abcd'\ncomplete: true\nn_reads: 0\n",
        ] {
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("metadata.yaml");
            fs::write(&path, input).unwrap();
            let ledger = budget(0);
            let expected: CacheMeta = serde_yaml_ng::from_str(input).unwrap();
            let actual = read(&path, Some(&ledger)).unwrap().unwrap();
            assert_eq!(
                serde_yaml_ng::to_string(&*actual).unwrap(),
                serde_yaml_ng::to_string(&expected).unwrap()
            );
            assert_eq!(
                ledger.snapshot().reserved_memory_bytes,
                (actual.url.capacity() + actual.sha256.capacity()) as u64
            );
            drop(actual);
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
            assert_eq!(Arc::strong_count(&ledger), 1);
            assert!(ledger.snapshot().valid);
        }
    }

    #[test]
    fn loader_refusal_releases_input_and_preserves_every_metadata_byte() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("metadata.yaml");
        let input = "url: http://example.com/source\nsha256: 'abcd'\ncomplete: true\n";
        fs::write(&path, input).unwrap();
        let ledger = budget(0);
        let mut refused = 0;
        let mut accepted = 0;
        for available in [0, 1, 64, 512, 4096, 65536, 1 << 20] {
            let competing = ledger
                .reserve(MemoryClass::RetainedResults, (8 << 20) - available)
                .unwrap();
            let before = ledger.snapshot().reserved_memory_bytes;
            match read(&path, Some(&ledger)) {
                Ok(meta) => {
                    accepted += 1;
                    drop(meta);
                }
                Err(error) => {
                    refused += 1;
                    assert!(
                        error.downcast_ref::<ResourceExhausted>().is_some(),
                        "{error:#}"
                    );
                }
            }
            assert_eq!(ledger.snapshot().reserved_memory_bytes, before);
            assert_eq!(fs::read(&path).unwrap(), input.as_bytes());
            drop(competing);
        }
        assert!(refused > 0 && accepted > 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }

    #[test]
    fn oversized_file_refuses_before_buffer_allocation_and_stays_untouched() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("metadata.yaml");
        File::create(&path)
            .unwrap()
            .set_len(MAX_BYTES as u64 + 1)
            .unwrap();
        let ledger = budget(0);
        let competing = ledger
            .reserve(MemoryClass::RetainedResults, 8 << 20)
            .unwrap();
        let error = read(&path, Some(&ledger)).err().unwrap();
        assert!(error.to_string().contains("byte limit"));
        assert_eq!(fs::metadata(&path).unwrap().len(), MAX_BYTES as u64 + 1);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 8 << 20);
        drop(competing);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }

    #[test]
    fn oversized_producer_refuses_before_creating_any_temporary() {
        let directory = tempfile::tempdir().unwrap();
        let meta = CacheMeta {
            url: "x".repeat(MAX_BYTES),
            sha256: "0".repeat(64),
            complete: true,
            n_reads: 1,
        };
        assert!(temporary(directory.path(), &meta, None).is_err());
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
    }

    #[test]
    fn metadata_disk_refusal_preserves_its_type_and_cleans_the_closed_temporary() {
        let directory = tempfile::tempdir().unwrap();
        let ledger = budget(0);
        let meta = Owned::generated(
            "http://example.com/source",
            "0".repeat(64),
            true,
            2,
            Some(&ledger),
        )
        .unwrap();
        let baseline = ledger.snapshot().reserved_memory_bytes;
        let error = temporary(directory.path(), &meta, Some(&ledger))
            .err()
            .unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, baseline);
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        drop(meta);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }

    #[test]
    fn published_metadata_keeps_exact_bytes_and_guard_until_directory_sync() {
        let directory = tempfile::tempdir().unwrap();
        let ledger = budget(1 << 20);
        let meta = Owned::generated(
            "http://example.com/source",
            "0".repeat(64),
            true,
            2,
            Some(&ledger),
        )
        .unwrap();
        let expected = serde_yaml_ng::to_string(&*meta).unwrap();
        let temporary = temporary(directory.path(), &meta, Some(&ledger)).unwrap();
        assert_eq!(fs::read(temporary.path()).unwrap(), expected.as_bytes());
        assert_eq!(
            ledger.snapshot().reserved_temporary_disk_bytes,
            expected.len() as u64
        );
        #[cfg(target_os = "linux")]
        for entry in fs::read_dir("/proc/self/fd").unwrap() {
            if let Ok(path) = fs::read_link(entry.unwrap().path()) {
                assert_ne!(path, temporary.path());
            }
        }
        let path = directory.path().join("accepted.yaml");
        let mut published = temporary.publish(&path, true).unwrap();
        assert_eq!(
            ledger.snapshot().reserved_temporary_disk_bytes,
            expected.len() as u64
        );
        crate::output::sync_directory(directory.path()).unwrap();
        published.release();
        drop(meta);
        let parsed = read(&path, Some(&ledger)).unwrap().unwrap();
        assert_eq!(serde_yaml_ng::to_string(&*parsed).unwrap(), expected);
        drop(parsed);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
        assert!(ledger.snapshot().valid);
    }

    #[test]
    fn bounded_lookup_preserves_underlying_resource_error_and_valid_entry_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let config = CacheConfig::new(Some(directory.path())).unwrap();
        let url = "http://example.com/source.fastq.gz";
        let key = cache_key(url);
        let data = config.data_path(&key);
        let meta = config.meta_path(&key);
        fs::write(&data, b"verified cache bytes").unwrap();
        write_meta(
            &meta,
            &CacheMeta {
                url: url.to_owned(),
                sha256: compute_sha256(&data).unwrap(),
                complete: true,
                n_reads: 1,
            },
        )
        .unwrap();
        let before = fs::read(&meta).unwrap();
        let ledger = budget(0);
        let competing = ledger
            .reserve(MemoryClass::RetainedResults, 8 << 20)
            .unwrap();
        let error = config.lookup_bounded(url, 0, &ledger).unwrap_err();
        assert!(
            error.downcast_ref::<ResourceExhausted>().is_some(),
            "{error:#}"
        );
        assert_eq!(fs::read(&meta).unwrap(), before);
        assert_eq!(fs::read(&data).unwrap(), b"verified cache bytes");
        drop(competing);
        assert_eq!(config.lookup_bounded(url, 0, &ledger).unwrap(), Some(data));
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }
}
