//! Sorted full-key disk counts with a single charged page cache.
//!
//! Native format 1 uses a 64-byte header, 4,096-byte pages and a 96-byte
//! completion footer. Each page holds 254 little-endian (u64 key, u64 count)
//! entries followed by SHA-256(header, page index, payload). A completed-stream
//! SHA-256 binds header, pages, totals and selected-input identity. Last-page
//! padding is zero. No hash/fingerprint substitutes for full-key membership.

use super::budget::{MemoryClass, MemoryReservation, ResourceBudget, ResourceExhausted};
use super::identity::{CanonicalKey, KmerSize};
use super::store::{CountMetadata, CountStore};
use super::temporary::OwnedDiskFile;
use anyhow::{Context, Result, ensure};
use rand::RngCore;
use sha2::{Digest, Sha256};
use std::num::NonZeroU64;
use std::path::Path;
use std::sync::{Arc, Mutex};

const HEADER_BYTES: usize = 64;
const FOOTER_BYTES: usize = 96;
pub(crate) const PAGE_BYTES: usize = 4096;
const PAYLOAD_BYTES: usize = PAGE_BYTES - 32;
pub(crate) const PAGE_ENTRIES: usize = PAYLOAD_BYTES / 16;
const MAGIC: &[u8; 8] = b"SHKCTP01";
const COMPLETE: &[u8; 8] = b"SHKEND01";

fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    let mut value = [0; 8];
    value.copy_from_slice(&bytes[offset..offset + 8]);
    u64::from_le_bytes(value)
}

fn header(k: KmerSize, identity: [u8; 16]) -> [u8; HEADER_BYTES] {
    let mut bytes = [0; HEADER_BYTES];
    bytes[..8].copy_from_slice(MAGIC);
    bytes[8..12].copy_from_slice(&1u32.to_le_bytes());
    bytes[12] = k.get() as u8;
    bytes[13] = 64;
    bytes[14] = 1; // min(forward, reverse complement)
    bytes[15] = 1; // uppercase ACGTN; N breaks windows
    bytes[16..20].copy_from_slice(&(PAGE_ENTRIES as u32).to_le_bytes());
    bytes[20..24].copy_from_slice(&(PAGE_BYTES as u32).to_le_bytes());
    bytes[24..40].copy_from_slice(&identity);
    bytes
}

fn page_checksum(header: &[u8; HEADER_BYTES], index: u64, payload: &[u8]) -> [u8; 32] {
    let mut hash = Sha256::new();
    hash.update(b"sharkmer-count-page-v1\0");
    hash.update(header);
    hash.update(index.to_le_bytes());
    hash.update(payload);
    hash.finalize().into()
}

struct PageBuffer {
    bytes: Vec<u8>,
    budget: Arc<ResourceBudget>,
    // The allocation is freed before its ownership charge.
    _reservation: MemoryReservation,
}

impl PageBuffer {
    fn new(budget: &Arc<ResourceBudget>) -> Result<Self> {
        let mut reservation = budget.reserve(MemoryClass::LookupCache, PAGE_BYTES as u64)?;
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(PAGE_BYTES)
            .context("Unable to reserve count-file page")?;
        reservation.resize(u64::try_from(bytes.capacity())?)?;
        bytes.resize(PAGE_BYTES, 0);
        Ok(Self {
            bytes,
            budget: Arc::clone(budget),
            _reservation: reservation,
        })
    }
}

/// Bounded streaming writer. Input must already be strictly sorted and reduced.
/// The same charged page allocation moves into the immutable reader.
pub(crate) struct SortedCountWriter {
    file: OwnedDiskFile,
    page: PageBuffer,
    header: [u8; HEADER_BYTES],
    k: KmerSize,
    used: usize,
    pages: u64,
    keys: u64,
    occurrences: u64,
    last_key: Option<u64>,
    stream_hash: Sha256,
    max_lookup_batch: usize,
    poisoned: bool,
}

impl SortedCountWriter {
    /// Reserve one page before allocation and stage an incomplete owned file.
    pub(crate) fn new(
        directory: &Path,
        k: KmerSize,
        budget: &Arc<ResourceBudget>,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        ensure!(
            max_lookup_batch > 0,
            "Count lookup batch reservation must be positive"
        );
        let page = PageBuffer::new(budget)?;
        let mut file = OwnedDiskFile::new(directory, budget)?;
        let mut identity = [0; 16];
        rand::rngs::OsRng
            .try_fill_bytes(&mut identity)
            .context("Unable to initialize count-file identity")?;
        let header = header(k, identity);
        file.append(&header)?;
        let mut stream_hash = Sha256::new();
        stream_hash.update(header);
        Ok(Self {
            file,
            page,
            header,
            k,
            used: 0,
            pages: 0,
            keys: 0,
            occurrences: 0,
            last_key: None,
            stream_hash,
            max_lookup_batch,
            poisoned: false,
        })
    }

    fn flush_page(&mut self) -> Result<()> {
        ensure!(!self.poisoned, "Count-file writer has failed");
        if self.used == 0 {
            return Ok(());
        }
        let next_page = self
            .pages
            .checked_add(1)
            .context("Count page total overflow")?;
        let checksum = page_checksum(&self.header, self.pages, &self.page.bytes[..PAYLOAD_BYTES]);
        self.page.bytes[PAYLOAD_BYTES..].copy_from_slice(&checksum);
        if let Err(error) = self.file.append(&self.page.bytes) {
            if error.downcast_ref::<ResourceExhausted>().is_none() {
                self.poisoned = true;
            }
            return Err(error);
        }
        self.stream_hash.update(&self.page.bytes);
        self.pages = next_page;
        self.used = 0;
        self.page.bytes.fill(0);
        Ok(())
    }

    /// Accept one exact positive count without collecting the complete table.
    /// Invalid identities/order/overflow leave the accepted prefix unchanged.
    pub(crate) fn push(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        ensure!(!self.poisoned, "Count-file writer has failed");
        ensure!(key.size() == self.k, "Count key k differs from file k");
        ensure!(
            self.last_key.is_none_or(|previous| previous < key.bits()),
            "Count-file keys must be strictly increasing"
        );
        let keys = self
            .keys
            .checked_add(1)
            .context("Count cardinality overflow")?;
        let occurrences = self
            .occurrences
            .checked_add(count.get())
            .context("Exact count occurrence total overflow")?;
        if self.used == PAGE_ENTRIES {
            self.flush_page()?;
        }
        let offset = self.used * 16;
        self.page.bytes[offset..offset + 8].copy_from_slice(&key.bits().to_le_bytes());
        self.page.bytes[offset + 8..offset + 16].copy_from_slice(&count.get().to_le_bytes());
        self.used += 1;
        self.keys = keys;
        self.occurrences = occurrences;
        self.last_key = Some(key.bits());
        Ok(())
    }

    /// Complete and sync before validating the entire stream for publication.
    pub(crate) fn finish(mut self, selected_input_sha256: [u8; 32]) -> Result<DiskCountStore> {
        self.flush_page()?;
        let mut footer = [0; FOOTER_BYTES];
        footer[..8].copy_from_slice(COMPLETE);
        footer[8..16].copy_from_slice(&self.keys.to_le_bytes());
        footer[16..24].copy_from_slice(&self.occurrences.to_le_bytes());
        footer[24..32].copy_from_slice(&self.pages.to_le_bytes());
        footer[32..64].copy_from_slice(&selected_input_sha256);
        self.stream_hash.update(&footer[..64]);
        footer[64..].copy_from_slice(&self.stream_hash.finalize());
        self.file.append(&footer)?;
        self.file.sync()?;
        DiskCountStore::open_owned(self.file, self.page, self.max_lookup_batch)
    }
}

struct DiskState {
    file: OwnedDiskFile,
    page: PageBuffer,
    cached_page: Option<u64>,
    failed: bool,
    #[cfg(test)]
    page_reads: u64,
}

/// Exact immutable full-key lookup with one 4 KiB cache, shared across workers.
/// A visit calls its visitor outside the file mutex, allowing reentrant lookup.
pub(crate) struct DiskCountStore {
    metadata: CountMetadata,
    header: [u8; HEADER_BYTES],
    file_bytes: u64,
    state: Mutex<DiskState>,
}

fn load_page(
    file: &mut OwnedDiskFile,
    page: &mut PageBuffer,
    header: &[u8; HEADER_BYTES],
    index: u64,
) -> Result<()> {
    file.verify_extent()?;
    let offset = index
        .checked_mul(PAGE_BYTES as u64)
        .and_then(|n| n.checked_add(HEADER_BYTES as u64))
        .context("Count page offset overflow")?;
    file.read_exact_at(offset, &mut page.bytes)?;
    let checksum = page_checksum(header, index, &page.bytes[..PAYLOAD_BYTES]);
    ensure!(
        page.bytes[PAYLOAD_BYTES..] == checksum,
        "Count page checksum mismatch"
    );
    Ok(())
}

fn entry(page: &PageBuffer, index: usize, k: KmerSize) -> Result<(CanonicalKey, NonZeroU64)> {
    ensure!(index < PAGE_ENTRIES, "Count page entry index exceeds page");
    let bits = read_u64(&page.bytes, index * 16);
    let key = CanonicalKey::new(bits, k)?;
    ensure!(key.bits() == bits, "Stored count key is not canonical");
    let count =
        NonZeroU64::new(read_u64(&page.bytes, index * 16 + 8)).context("Stored count is zero")?;
    Ok((key, count))
}

impl DiskCountStore {
    fn open_owned(
        mut file: OwnedDiskFile,
        mut page: PageBuffer,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        file.verify_extent()?;
        let file_bytes = file.len()?;
        ensure!(
            file_bytes >= (HEADER_BYTES + FOOTER_BYTES) as u64,
            "Incomplete count file"
        );
        let mut persisted_header = [0; HEADER_BYTES];
        file.read_exact_at(0, &mut persisted_header)?;
        let k = KmerSize::new(usize::from(persisted_header[12]))?;
        let mut identity = [0; 16];
        identity.copy_from_slice(&persisted_header[24..40]);
        ensure!(
            persisted_header == header(k, identity),
            "Unsupported count-file identity or format"
        );
        let mut footer = [0; FOOTER_BYTES];
        file.read_exact_at(file_bytes - FOOTER_BYTES as u64, &mut footer)?;
        ensure!(&footer[..8] == COMPLETE, "Count file is not complete");
        let keys = read_u64(&footer, 8);
        let occurrences = read_u64(&footer, 16);
        let pages = read_u64(&footer, 24);
        ensure!(
            pages == keys.div_ceil(PAGE_ENTRIES as u64),
            "Count-file page cardinality mismatch"
        );
        let expected_bytes = pages
            .checked_mul(PAGE_BYTES as u64)
            .and_then(|n| n.checked_add((HEADER_BYTES + FOOTER_BYTES) as u64))
            .context("Count-file size overflow")?;
        ensure!(
            file_bytes == expected_bytes,
            "Count-file extent/cardinality mismatch"
        );
        let mut stream_hash = Sha256::new();
        stream_hash.update(persisted_header);
        let mut previous = None;
        let mut seen = 0u64;
        let mut total = 0u64;
        for index in 0..pages {
            load_page(&mut file, &mut page, &persisted_header, index)?;
            stream_hash.update(&page.bytes);
            let used = (keys - seen).min(PAGE_ENTRIES as u64) as usize;
            for offset in 0..used {
                let (key, count) = entry(&page, offset, k)?;
                ensure!(
                    previous.is_none_or(|bits| bits < key.bits()),
                    "Stored count keys are not strictly sorted"
                );
                total = total
                    .checked_add(count.get())
                    .context("Stored occurrence total overflow")?;
                seen += 1;
                previous = Some(key.bits());
            }
            ensure!(
                page.bytes[used * 16..PAYLOAD_BYTES]
                    .iter()
                    .all(|byte| *byte == 0),
                "Count-page padding is nonzero"
            );
        }
        ensure!(
            seen == keys && total == occurrences,
            "Count-file occurrence/cardinality mismatch"
        );
        stream_hash.update(&footer[..64]);
        ensure!(
            stream_hash.finalize().as_slice() == &footer[64..],
            "Completed count-stream checksum mismatch"
        );
        let mut selected_input_sha256 = [0; 32];
        selected_input_sha256.copy_from_slice(&footer[32..64]);
        let metadata = CountMetadata {
            format_version: 1,
            k,
            count_width_bits: 64,
            backend: "checksummed_sorted_file",
            hash_algorithm: "none",
            canonicalization: "min_forward_reverse_complement",
            preprocessing: "uppercase_acgtn_n_breaks_windows",
            overflow_policy: "error_preserve_previous_counts",
            selected_input_sha256,
            distinct_keys: keys,
            occurrences,
            max_lookup_batch,
        };
        Ok(Self {
            metadata,
            header: persisted_header,
            file_bytes,
            state: Mutex::new(DiskState {
                file,
                page,
                cached_page: pages.checked_sub(1),
                failed: false,
                #[cfg(test)]
                page_reads: pages,
            }),
        })
    }

    /// Report the live file extent, distinct from RAM or filesystem allocation.
    pub(crate) fn file_bytes(&self) -> u64 {
        self.file_bytes
    }

    /// Confirm removal before a producer or output transaction completes.
    pub(crate) fn close(self) -> Result<()> {
        let state = self
            .state
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        state.file.close()
    }

    /// Report the single owned cache allocation (not a process-memory bound).
    pub(crate) fn cache_bytes(&self) -> usize {
        self.state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .page
            .bytes
            .capacity()
    }

    /// Peak charged count pages, including independent concurrent scan buffers.
    pub(crate) fn peak_page_bytes(&self) -> u64 {
        self.state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .page
            .budget
            .snapshot()
            .peak_reserved_by_class[MemoryClass::LookupCache as usize]
    }

    pub(crate) fn entry_at(&self, index: u64) -> Result<(CanonicalKey, NonZeroU64)> {
        ensure!(
            index < self.metadata.distinct_keys,
            "Count entry index exceeds store"
        );
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        ensure!(!state.failed, "Count-file access previously failed");
        let page_index = index / PAGE_ENTRIES as u64;
        if state.cached_page != Some(page_index) {
            state.cached_page = None;
            let DiskState { file, page, .. } = &mut *state;
            if let Err(error) = load_page(file, page, &self.header, page_index) {
                state.failed = true;
                return Err(error);
            }
            state.cached_page = Some(page_index);
            #[cfg(test)]
            {
                state.page_reads += 1;
            }
        }
        let result = entry(
            &state.page,
            (index % PAGE_ENTRIES as u64) as usize,
            self.metadata.k,
        );
        if result.is_err() {
            state.failed = true;
        }
        result
    }
}

#[cfg(test)]
mod tests;

impl CountStore for DiskCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }

    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        ensure!(
            key.size() == self.metadata.k,
            "Count lookup k differs from store k"
        );
        let (mut left, mut right) = (0, self.metadata.distinct_keys);
        while left < right {
            let middle = left + (right - left) / 2;
            let (stored, count) = self.entry_at(middle)?;
            match stored.bits().cmp(&key.bits()) {
                std::cmp::Ordering::Equal => return Ok(Some(count)),
                std::cmp::Ordering::Less => left = middle + 1,
                std::cmp::Ordering::Greater => right = middle,
            }
        }
        Ok(None)
    }

    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        if self.metadata.distinct_keys == 0 {
            return Ok(());
        }
        let budget = {
            let state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            ensure!(!state.failed, "Count-file access previously failed");
            Arc::clone(&state.page.budget)
        };
        // Each scan retains one charged page independently of other scans and
        // the shared lookup cache. Visitors always run outside the file lock.
        let mut page = PageBuffer::new(&budget)?;
        let pages = self.metadata.distinct_keys.div_ceil(PAGE_ENTRIES as u64);
        let mut seen = 0;
        for index in 0..pages {
            {
                let mut state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                ensure!(!state.failed, "Count-file access previously failed");
                if let Err(error) = load_page(&mut state.file, &mut page, &self.header, index) {
                    state.failed = true;
                    return Err(error);
                }
                #[cfg(test)]
                {
                    state.page_reads += 1;
                }
            }
            let used = (self.metadata.distinct_keys - seen).min(PAGE_ENTRIES as u64) as usize;
            for offset in 0..used {
                let (key, count) = match entry(&page, offset, self.metadata.k) {
                    Ok(value) => value,
                    Err(error) => {
                        self.state
                            .lock()
                            .unwrap_or_else(std::sync::PoisonError::into_inner)
                            .failed = true;
                        return Err(error);
                    }
                };
                visitor(key, count)?;
                seen += 1;
            }
        }
        Ok(())
    }
}
