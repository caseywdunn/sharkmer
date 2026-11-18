//! Byte-packed full-key columns with exact u64 counts and bounded growth.
//!
//! Count zero marks an empty bucket; every observed abundance remains positive.
//! Keys retain ceil(2*k/8) bytes. Seven trailing bytes permit safe eight-byte
//! loads at every bucket, followed by a width mask. Membership compares the
//! complete key; neither fingerprints nor compressed counts establish presence.

use super::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use super::identity::{CanonicalKey, KmerSize, OrientedKey};
use super::store::{CountMetadata, CountStore, TABLE_HASHER_NAME, TableHasher};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::hash::BuildHasher;
use std::num::NonZeroU64;
use std::sync::Arc;

struct Columns {
    keys: Vec<u8>,
    counts: Vec<u64>,
    width: usize,
    mask: u64,
    // Both allocations drop before their shared reservation.
    _reservation: MemoryReservation,
}

impl Columns {
    fn allocation_bytes(buckets: usize, width: usize) -> Result<usize> {
        ensure!(
            buckets >= 16 && buckets.is_power_of_two(),
            "Invalid packed bucket count"
        );
        ensure!((1..=8).contains(&width), "Invalid packed key width");
        buckets
            .checked_mul(width + size_of::<u64>())
            .and_then(|bytes| bytes.checked_add(7))
            .filter(|bytes| *bytes <= isize::MAX as usize)
            .context("Packed count allocation size overflow")
    }

    fn new(buckets: usize, width: usize, budget: &Arc<ResourceBudget>) -> Result<Self> {
        let bytes = Self::allocation_bytes(buckets, width)?;
        let reservation = budget.reserve(MemoryClass::CountSlots, u64::try_from(bytes)?)?;
        let key_bytes = buckets
            .checked_mul(width)
            .and_then(|n| n.checked_add(7))
            .context("Packed key allocation overflow")?;
        let mut keys = Vec::new();
        keys.try_reserve_exact(key_bytes)
            .context("Unable to reserve packed key column")?;
        let mut counts = Vec::new();
        counts
            .try_reserve_exact(buckets)
            .context("Unable to reserve exact count column")?;
        ensure!(
            keys.capacity() == key_bytes && counts.capacity() == buckets,
            "Packed allocation differs from reservation"
        );
        keys.resize(key_bytes, 0);
        counts.resize(buckets, 0);
        Ok(Self {
            keys,
            counts,
            width,
            mask: u64::MAX >> (8 * (8 - width)),
            _reservation: reservation,
        })
    }

    #[inline]
    fn bits(&self, bucket: usize) -> u64 {
        let offset = bucket * self.width;
        let mut bytes = [0; 8];
        bytes.copy_from_slice(&self.keys[offset..offset + 8]);
        u64::from_le_bytes(bytes) & self.mask
    }

    #[inline]
    fn find(&self, bits: u64, hash: u64) -> usize {
        let mask = self.counts.len() - 1;
        let mut bucket = hash as usize & mask;
        // Growth keeps at least one eighth of buckets empty. A full-key
        // comparison is required even when hashes/probe origins collide.
        while self.counts[bucket] != 0 && self.bits(bucket) != bits {
            bucket = (bucket + 1) & mask;
        }
        bucket
    }

    #[inline]
    fn put(&mut self, bucket: usize, bits: u64, count: u64) {
        let offset = bucket * self.width;
        self.keys[offset..offset + self.width].copy_from_slice(&bits.to_le_bytes()[..self.width]);
        self.counts[bucket] = count;
    }

    fn capacity(&self) -> usize {
        self.counts.len() / 8 * 7
    }
    fn bytes(&self) -> usize {
        self.keys.capacity() + self.counts.capacity() * size_of::<u64>()
    }
}

/// Measurements emitted only when packed storage was actually used.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct PackedMetrics {
    pub(crate) key_bytes_per_bucket: usize,
    pub(crate) count_bytes_per_bucket: usize,
    pub(crate) buckets: usize,
    pub(crate) padding_bytes: usize,
}

/// One owner preserves every accepted count across rejected updates/growth.
pub(crate) struct PackedCountBuilder {
    columns: Columns,
    budget: Arc<ResourceBudget>,
    hasher: TableHasher,
    k: KmerSize,
    distinct: usize,
    occurrences: u64,
    peak_table_bytes: usize,
    max_lookup_batch: usize,
}

impl PackedCountBuilder {
    /// Reserve both columns before allocation; resizing retains the old charge.
    pub(crate) fn with_budget(
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        ensure!(
            max_lookup_batch > 0,
            "Packed lookup batch reservation must be positive"
        );
        let columns = Columns::new(16, (2 * k.get()).div_ceil(8), &budget)?;
        let peak_table_bytes = columns.bytes();
        Ok(Self {
            columns,
            budget,
            hasher: TableHasher::default(),
            k,
            distinct: 0,
            occurrences: 0,
            peak_table_bytes,
            max_lookup_batch,
        })
    }

    fn grow(&mut self) -> Result<()> {
        let buckets = self
            .columns
            .counts
            .len()
            .checked_mul(2)
            .context("Packed bucket growth overflow")?;
        let mut replacement = Columns::new(buckets, self.columns.width, &self.budget)?;
        for (bucket, &count) in self.columns.counts.iter().enumerate() {
            if count != 0 {
                let bits = self.columns.bits(bucket);
                let destination = replacement.find(bits, BuildHasher::hash_one(&self.hasher, bits));
                replacement.put(destination, bits, count);
            }
        }
        let peak = self
            .columns
            .bytes()
            .checked_add(replacement.bytes())
            .context("Packed growth peak overflow")?;
        self.peak_table_bytes = self.peak_table_bytes.max(peak);
        self.columns = replacement;
        Ok(())
    }

    /// Checked exact addition. Key zero is valid; count zero is never observed.
    #[inline]
    pub(crate) fn add(&mut self, key: CanonicalKey, increment: NonZeroU64) -> Result<()> {
        ensure!(
            key.size() == self.k,
            "Count key k differs from packed store k"
        );
        let Some(total) = self.occurrences.checked_add(increment.get()) else {
            anyhow::bail!("Exact k-mer occurrence total overflow");
        };
        let bits = key.bits();
        let hash = BuildHasher::hash_one(&self.hasher, bits);
        let mut bucket = self.columns.find(bits, hash);
        let existing = self.columns.counts[bucket];
        if existing == 0 && self.distinct == self.columns.capacity() {
            self.grow()?;
            bucket = self.columns.find(bits, hash);
        }
        let Some(updated) = existing.checked_add(increment.get()) else {
            anyhow::bail!("Exact k-mer count overflow");
        };
        self.columns.put(bucket, bits, updated);
        self.distinct += usize::from(existing == 0);
        self.occurrences = total;
        Ok(())
    }

    /// Accepted occurrence total without a second global scan.
    pub(crate) fn occurrences(&self) -> u64 {
        self.occurrences
    }

    /// Move the columns into an immutable assembly store without consolidation.
    pub(crate) fn finalize(self, selected_input_sha256: [u8; 32]) -> PackedCountStore {
        let metadata = CountMetadata {
            format_version: 1,
            k: self.k,
            count_width_bits: 64,
            backend: "packed_columns_memory",
            hash_algorithm: TABLE_HASHER_NAME,
            canonicalization: "min_forward_reverse_complement",
            preprocessing: "uppercase_acgtn_n_breaks_windows",
            overflow_policy: "error_preserve_previous_counts",
            selected_input_sha256,
            distinct_keys: self.distinct as u64,
            occurrences: self.occurrences,
            max_lookup_batch: self.max_lookup_batch,
        };
        PackedCountStore {
            builder: self,
            metadata,
        }
    }
}

/// Full-key exact membership over immutable packed columns.
pub(crate) struct PackedCountStore {
    builder: PackedCountBuilder,
    metadata: CountMetadata,
}

impl PackedCountStore {
    pub(crate) fn capacity(&self) -> usize {
        self.builder.columns.capacity()
    }
    pub(crate) fn table_bytes(&self) -> usize {
        self.builder.columns.bytes()
    }
    pub(crate) fn peak_table_bytes(&self) -> usize {
        self.builder.peak_table_bytes
    }
    pub(crate) fn metrics(&self) -> PackedMetrics {
        PackedMetrics {
            key_bytes_per_bucket: self.builder.columns.width,
            count_bytes_per_bucket: size_of::<u64>(),
            buckets: self.builder.columns.counts.len(),
            padding_bytes: 7,
        }
    }

    fn lookup_bits(&self, bits: u64) -> Option<NonZeroU64> {
        let bucket = self
            .builder
            .columns
            .find(bits, BuildHasher::hash_one(&self.builder.hasher, bits));
        NonZeroU64::new(self.builder.columns.counts[bucket])
    }

    /// Concrete scans let assembly inline visitors without copying columns.
    pub(crate) fn visit_entries<F: FnMut(CanonicalKey, NonZeroU64) -> Result<()> + ?Sized>(
        &self,
        visitor: &mut F,
    ) -> Result<()> {
        for (bucket, &count) in self.builder.columns.counts.iter().enumerate() {
            if let Some(count) = NonZeroU64::new(count) {
                visitor(
                    CanonicalKey::new(self.builder.columns.bits(bucket), self.metadata.k)?,
                    count,
                )?;
            }
        }
        Ok(())
    }
}

impl CountStore for PackedCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }
    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        ensure!(
            key.size() == self.metadata.k,
            "Count lookup k differs from packed store k"
        );
        Ok(self.lookup_bits(key.bits()))
    }
    fn lookup_oriented(&self, key: OrientedKey) -> Result<Option<NonZeroU64>> {
        ensure!(
            key.size() == self.metadata.k,
            "Count lookup k differs from packed store k"
        );
        Ok(self
            .lookup_bits(key.bits())
            .or_else(|| self.lookup_bits(key.reverse_complement().bits())))
    }
    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        self.visit_entries(visitor)
    }
}

#[cfg(test)]
mod tests;
