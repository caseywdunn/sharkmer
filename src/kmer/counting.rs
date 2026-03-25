// kmer/counting.rs
//! Kmer counting with feature-gated hash map backends.

#[cfg(all(feature = "ahashmap", feature = "fxhashmap"))]
compile_error!("Features 'ahashmap' and 'fxhashmap' are mutually exclusive. Enable only one.");

use anyhow::{Result, bail};

#[cfg(feature = "fxhashmap")]
use rustc_hash::FxHashMap;

#[cfg(feature = "ahashmap")]
use std::collections::HashMap as StdHashMap;

use super::encoding::{Read, revcomp_kmer};
use super::histogram::Histogram;

/// Trait abstracting over hash map implementations for kmer counting.
#[allow(dead_code)]
trait KmerMap {
    fn new_map() -> Self;
    fn new_map_with_capacity(capacity: usize) -> Self;
    fn insert_or_add(&mut self, key: u64, count: u64);
    fn get_value(&self, key: u64) -> Option<&u64>;
    fn get_count(&self, key: u64) -> u64;
    fn contains(&self, key: u64) -> bool;
    fn retain_above(&mut self, min_count: u64);
}

#[cfg(feature = "fxhashmap")]
impl KmerMap for rustc_hash::FxHashMap<u64, u64> {
    fn new_map() -> Self {
        FxHashMap::default()
    }
    fn new_map_with_capacity(capacity: usize) -> Self {
        FxHashMap::with_capacity_and_hasher(capacity, Default::default())
    }
    fn insert_or_add(&mut self, key: u64, count: u64) {
        let c = self.entry(key).or_insert(0);
        *c += count;
    }
    fn get_value(&self, key: u64) -> Option<&u64> {
        self.get(&key)
    }
    fn get_count(&self, key: u64) -> u64 {
        *self.get(&key).unwrap_or(&0)
    }
    fn contains(&self, key: u64) -> bool {
        self.contains_key(&key)
    }
    fn retain_above(&mut self, min_count: u64) {
        self.retain(|_, count| *count >= min_count);
    }
}

#[cfg(feature = "ahashmap")]
type AHashMap<K, V> = StdHashMap<K, V, ahash::RandomState>;

#[cfg(feature = "ahashmap")]
impl KmerMap for AHashMap<u64, u64> {
    fn new_map() -> Self {
        AHashMap::with_hasher(ahash::RandomState::new())
    }
    fn new_map_with_capacity(capacity: usize) -> Self {
        AHashMap::with_capacity_and_hasher(capacity, ahash::RandomState::new())
    }
    fn insert_or_add(&mut self, key: u64, count: u64) {
        let c = self.entry(key).or_insert(0);
        *c += count;
    }
    fn get_value(&self, key: u64) -> Option<&u64> {
        self.get(&key)
    }
    fn get_count(&self, key: u64) -> u64 {
        *self.get(&key).unwrap_or(&0)
    }
    fn contains(&self, key: u64) -> bool {
        self.contains_key(&key)
    }
    fn retain_above(&mut self, min_count: u64) {
        self.retain(|_, count| *count >= min_count);
    }
}

#[cfg(feature = "fxhashmap")]
type MapType = rustc_hash::FxHashMap<u64, u64>;

#[cfg(feature = "ahashmap")]
type MapType = AHashMap<u64, u64>;

pub struct KmerCounts {
    kmers: MapType,
    k: usize,
}

impl KmerCounts {
    pub fn new(k: &usize) -> KmerCounts {
        KmerCounts {
            kmers: MapType::new_map(),
            k: *k,
        }
    }

    #[allow(dead_code)]
    pub fn new_with_capacity(k: &usize, capacity: usize) -> KmerCounts {
        KmerCounts {
            kmers: MapType::new_map_with_capacity(capacity),
            k: *k,
        }
    }

    pub fn ingest_reads(&mut self, reads: &[Read]) -> Result<()> {
        for read in reads {
            for kmer in read.get_kmers(&self.k)? {
                self.kmers.insert_or_add(kmer, 1);
            }
        }
        Ok(())
    }

    /// Adds a kmer and its count to the KmerCounts object.
    pub fn insert(&mut self, kmer: &u64, count: &u64) {
        self.kmers.insert_or_add(*kmer, *count);
    }

    /// Adds the counts from another KmerCounts object to this one.
    pub fn extend(&mut self, other: &KmerCounts) -> Result<()> {
        if self.k != other.k {
            bail!("Cannot extend KmerCounts with different k");
        }

        for (kmer, count) in other.iter() {
            self.kmers.insert_or_add(*kmer, *count);
        }
        Ok(())
    }

    /// Adds the counts from another KmerCounts object to this one,
    /// incrementally updating a histogram to reflect the new counts.
    /// This avoids a full table scan to rebuild the histogram after each merge.
    pub fn extend_with_histogram(
        &mut self,
        other: &KmerCounts,
        histo: &mut Histogram,
    ) -> Result<()> {
        if self.k != other.k {
            bail!("Cannot extend KmerCounts with different k");
        }

        for (kmer, new_count) in other.iter() {
            let old_count = self.kmers.get_count(*kmer);
            self.kmers.insert_or_add(*kmer, *new_count);
            histo.move_count(old_count, old_count + *new_count);
        }
        Ok(())
    }

    /// Get the count of the canonical form of a kmer (min of forward and revcomp).
    pub fn get_canonical_count(&self, kmer: &u64) -> u64 {
        let revcomp = revcomp_kmer(kmer, &self.k);
        let canonical = if *kmer < revcomp { *kmer } else { revcomp };
        self.kmers.get_count(canonical)
    }

    #[allow(dead_code)]
    pub fn get(&self, kmer: &u64) -> Option<&u64> {
        self.kmers.get_value(*kmer)
    }

    /// Look up a kmer by checking both orientations (forward and reverse complement).
    /// Returns the count if found in either orientation, None if absent.
    pub fn get_canonical(&self, kmer: &u64) -> Option<&u64> {
        self.kmers
            .get_value(*kmer)
            .or_else(|| self.kmers.get_value(revcomp_kmer(kmer, &self.k)))
    }

    #[allow(dead_code)]
    pub fn get_count(&self, kmer: &u64) -> u64 {
        self.kmers.get_count(*kmer)
    }

    #[allow(dead_code)]
    pub fn contains(&self, kmer: &u64) -> bool {
        self.kmers.contains(*kmer)
    }

    #[allow(dead_code)]
    pub fn remove_low_count_kmers(&mut self, min_count: &u64) {
        self.kmers.retain_above(*min_count);
    }

    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u64)> {
        self.kmers.iter()
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn get_k(&self) -> usize {
        self.k
    }

    pub fn get_n_kmers(&self) -> u64 {
        self.kmers.values().sum()
    }

    pub fn get_n_unique_kmers(&self) -> u64 {
        self.kmers.len() as u64
    }

    pub fn counts(&self) -> Vec<u64> {
        self.kmers.values().copied().collect()
    }

    pub fn kmers(&self) -> Vec<u64> {
        self.kmers.keys().copied().collect()
    }

    pub fn get_max_count(&self) -> u64 {
        *self.kmers.values().max().unwrap_or(&0)
    }

    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    /// Create a filtered view of this KmerCounts that lazily excludes entries
    /// below `min_count`. No data is copied — lookups check the threshold on the fly.
    pub fn filtered_view(&self, min_count: u64) -> FilteredKmerCounts<'_> {
        FilteredKmerCounts {
            inner: self,
            min_count,
        }
    }
}

/// A read-only view of a `KmerCounts` that filters out entries below a minimum count.
/// No data is copied — the threshold is checked at lookup time.
pub struct FilteredKmerCounts<'a> {
    inner: &'a KmerCounts,
    min_count: u64,
}

impl<'a> FilteredKmerCounts<'a> {
    pub fn get_k(&self) -> usize {
        self.inner.k
    }

    /// Look up a kmer by checking both orientations, returning the count
    /// only if it meets the minimum threshold.
    pub fn get_canonical(&self, kmer: &u64) -> Option<u64> {
        self.inner.get_canonical(kmer).and_then(|&count| {
            if count >= self.min_count {
                Some(count)
            } else {
                None
            }
        })
    }

    /// Get the count of the canonical form of a kmer, returning 0 if below threshold.
    pub fn get_canonical_count(&self, kmer: &u64) -> u64 {
        let count = self.inner.get_canonical_count(kmer);
        if count >= self.min_count { count } else { 0 }
    }

    /// Iterate over all entries, including those below threshold.
    /// Callers that need filtering should check counts themselves
    /// (as find_oligos_in_kmers already does).
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u64)> {
        self.inner.iter()
    }
}

impl Clone for KmerCounts {
    fn clone(&self) -> Self {
        KmerCounts {
            kmers: self.kmers.clone(),
            k: self.k,
        }
    }
}
