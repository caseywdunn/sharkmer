// kmer/counting.rs
//! Kmer counting with feature-gated hash map backends.

#[cfg(all(feature = "ahashmap", feature = "fxhashmap"))]
compile_error!("Features 'ahashmap' and 'fxhashmap' are mutually exclusive. Enable only one.");

use anyhow::{Result, bail};

#[cfg(feature = "fxhashmap")]
use rustc_hash::FxHashMap;

#[cfg(feature = "ahashmap")]
use std::collections::HashMap as StdHashMap;

use super::encoding::{Read, kmers_from_ascii, revcomp_kmer};
use super::histogram::Histogram;

/// Trait abstracting over hash map implementations for kmer counting.
#[allow(dead_code)]
trait KmerMap {
    fn new_map() -> Self;
    fn new_map_with_capacity(capacity: usize) -> Self;
    fn insert_or_add(&mut self, key: u64, count: u32);
    /// Insert a kmer or add to its existing count, returning
    /// `(old_count, new_count_after_saturation)`. The new count is
    /// `old.saturating_add(count)`, so at `u32::MAX` further additions are
    /// silently dropped. Callers that maintain a parallel histogram must
    /// use `new_count_after_saturation` as the new bin rather than
    /// `old + count`, otherwise the histogram drifts out of sync with the
    /// stored count once any kmer saturates.
    fn insert_or_add_get_counts(&mut self, key: u64, count: u32) -> (u32, u32);
    fn get_value(&self, key: u64) -> Option<&u32>;
    fn get_count(&self, key: u64) -> u32;
    fn contains(&self, key: u64) -> bool;
    fn retain_above(&mut self, min_count: u32);
}

#[cfg(feature = "fxhashmap")]
impl KmerMap for rustc_hash::FxHashMap<u64, u32> {
    fn new_map() -> Self {
        FxHashMap::default()
    }
    fn new_map_with_capacity(capacity: usize) -> Self {
        FxHashMap::with_capacity_and_hasher(capacity, Default::default())
    }
    fn insert_or_add(&mut self, key: u64, count: u32) {
        let c = self.entry(key).or_insert(0);
        *c = c.saturating_add(count);
    }
    fn insert_or_add_get_counts(&mut self, key: u64, count: u32) -> (u32, u32) {
        let c = self.entry(key).or_insert(0);
        let old = *c;
        let new = old.saturating_add(count);
        *c = new;
        (old, new)
    }
    fn get_value(&self, key: u64) -> Option<&u32> {
        self.get(&key)
    }
    fn get_count(&self, key: u64) -> u32 {
        *self.get(&key).unwrap_or(&0)
    }
    fn contains(&self, key: u64) -> bool {
        self.contains_key(&key)
    }
    fn retain_above(&mut self, min_count: u32) {
        self.retain(|_, count| *count >= min_count);
    }
}

#[cfg(feature = "ahashmap")]
type AHashMap<K, V> = StdHashMap<K, V, ahash::RandomState>;

#[cfg(feature = "ahashmap")]
impl KmerMap for AHashMap<u64, u32> {
    fn new_map() -> Self {
        AHashMap::with_hasher(ahash::RandomState::new())
    }
    fn new_map_with_capacity(capacity: usize) -> Self {
        AHashMap::with_capacity_and_hasher(capacity, ahash::RandomState::new())
    }
    fn insert_or_add(&mut self, key: u64, count: u32) {
        let c = self.entry(key).or_insert(0);
        *c = c.saturating_add(count);
    }
    fn insert_or_add_get_counts(&mut self, key: u64, count: u32) -> (u32, u32) {
        let c = self.entry(key).or_insert(0);
        let old = *c;
        let new = old.saturating_add(count);
        *c = new;
        (old, new)
    }
    fn get_value(&self, key: u64) -> Option<&u32> {
        self.get(&key)
    }
    fn get_count(&self, key: u64) -> u32 {
        *self.get(&key).unwrap_or(&0)
    }
    fn contains(&self, key: u64) -> bool {
        self.contains_key(&key)
    }
    fn retain_above(&mut self, min_count: u32) {
        self.retain(|_, count| *count >= min_count);
    }
}

#[cfg(feature = "fxhashmap")]
type MapType = rustc_hash::FxHashMap<u64, u32>;

#[cfg(feature = "ahashmap")]
type MapType = AHashMap<u64, u32>;

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

    pub fn new_with_capacity(k: &usize, capacity: usize) -> KmerCounts {
        KmerCounts {
            kmers: MapType::new_map_with_capacity(capacity),
            k: *k,
        }
    }

    #[allow(dead_code)]
    pub fn ingest_reads(&mut self, reads: &[Read]) -> Result<()> {
        for read in reads {
            for kmer in read.get_kmers(&self.k)? {
                self.kmers.insert_or_add(kmer, 1);
            }
        }
        Ok(())
    }

    /// Ingest kmers directly from an ASCII sequence, bypassing Read struct encoding.
    pub fn ingest_seq(&mut self, seq: &str) -> Result<()> {
        for kmer in kmers_from_ascii(seq, self.k)? {
            self.kmers.insert_or_add(kmer, 1);
        }
        Ok(())
    }

    /// Ingest kmers from an ASCII sequence, checking each kmer against an
    /// Oligo filter. Returns true if any kmer matched the filter.
    pub fn ingest_seq_with_filter(
        &mut self,
        seq: &str,
        filter: &crate::io::OligoFilter,
    ) -> Result<bool> {
        let mut matched = false;
        for kmer in kmers_from_ascii(seq, self.k)? {
            self.kmers.insert_or_add(kmer, 1);
            if !matched && filter.check_kmer(kmer) {
                matched = true;
            }
        }
        Ok(matched)
    }

    /// Adds a kmer and its count to the KmerCounts object.
    pub fn insert(&mut self, kmer: &u64, count: &u32) {
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

        // Track whether any kmer hit u32::MAX during this merge so the
        // operator sees a single warning rather than silent saturation.
        let mut any_saturated = false;
        for (kmer, new_count) in other.iter() {
            let (old_count, stored_count) = self.kmers.insert_or_add_get_counts(*kmer, *new_count);
            // Use the actual stored count (post-saturation) as the new
            // histogram bin, not old + new. Without this, a saturating
            // kmer's histogram entry drifts above u32::MAX while its real
            // count is capped.
            histo.move_count(old_count as u64, stored_count as u64);
            if stored_count == u32::MAX && old_count < u32::MAX {
                any_saturated = true;
            }
        }
        if any_saturated {
            log::warn!(
                "One or more kmer counts reached u32::MAX during ingestion. \
                 Counts are now capped; histogram reflects the capped values. \
                 This is only expected on extremely high-copy repeat kmers."
            );
        }
        Ok(())
    }

    /// Get the count of the canonical form of a kmer (min of forward and revcomp).
    pub fn get_canonical_count(&self, kmer: &u64) -> u32 {
        let revcomp = revcomp_kmer(kmer, &self.k);
        let canonical = if *kmer < revcomp { *kmer } else { revcomp };
        self.kmers.get_count(canonical)
    }

    #[allow(dead_code)]
    pub fn get(&self, kmer: &u64) -> Option<&u32> {
        self.kmers.get_value(*kmer)
    }

    /// Look up a kmer by checking both orientations (forward and reverse complement).
    /// Returns the count if found in either orientation, None if absent.
    pub fn get_canonical(&self, kmer: &u64) -> Option<&u32> {
        self.kmers
            .get_value(*kmer)
            .or_else(|| self.kmers.get_value(revcomp_kmer(kmer, &self.k)))
    }

    #[allow(dead_code)]
    pub fn get_count(&self, kmer: &u64) -> u32 {
        self.kmers.get_count(*kmer)
    }

    #[allow(dead_code)]
    pub fn contains(&self, kmer: &u64) -> bool {
        self.kmers.contains(*kmer)
    }

    #[allow(dead_code)]
    pub fn remove_low_count_kmers(&mut self, min_count: &u32) {
        self.kmers.retain_above(*min_count);
    }

    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u32)> {
        self.kmers.iter()
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    pub fn get_k(&self) -> usize {
        self.k
    }

    /// Total count of all kmers (sum of all counts). Returns u64 to avoid overflow.
    pub fn get_n_kmers(&self) -> u64 {
        self.kmers.values().map(|&v| v as u64).sum()
    }

    pub fn get_n_unique_kmers(&self) -> u64 {
        self.kmers.len() as u64
    }

    pub fn counts(&self) -> Vec<u32> {
        self.kmers.values().copied().collect()
    }

    pub fn kmers(&self) -> Vec<u64> {
        self.kmers.keys().copied().collect()
    }

    pub fn get_max_count(&self) -> u32 {
        *self.kmers.values().max().unwrap_or(&0)
    }

    pub fn get_median_count(&self) -> u32 {
        let mut counts: Vec<u32> = self.kmers.values().copied().collect();
        if counts.is_empty() {
            return 0;
        }
        // select_nth_unstable is expected-linear vs sort's O(n log n).
        // For odd len, the nth element is the median directly.
        // For even len, partition at `mid` (placing the upper half to the
        // right of index mid), then the max of the lower half is the
        // (mid-1)th order statistic — the other middle value. Averaged
        // with rounding-down via half-sums to avoid overflow.
        let mid = counts.len() / 2;
        if counts.len() % 2 == 0 {
            counts.select_nth_unstable(mid);
            let upper_min = counts[mid];
            let lower_max = *counts[..mid]
                .iter()
                .max()
                .expect("non-empty lower half when len >= 2");
            (lower_max / 2) + (upper_min / 2)
        } else {
            *counts.select_nth_unstable(mid).1
        }
    }

    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    /// Create a filtered view of this KmerCounts that lazily excludes entries
    /// below `min_count`. No data is copied — lookups check the threshold on the fly.
    pub fn filtered_view(&self, min_count: u32) -> FilteredKmerCounts<'_> {
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
    min_count: u32,
}

impl<'a> FilteredKmerCounts<'a> {
    pub fn get_k(&self) -> usize {
        self.inner.k
    }

    /// Look up a kmer by checking both orientations, returning the count
    /// only if it meets the minimum threshold.
    pub fn get_canonical(&self, kmer: &u64) -> Option<u32> {
        self.inner.get_canonical(kmer).and_then(|&count| {
            if count >= self.min_count {
                Some(count)
            } else {
                None
            }
        })
    }

    /// Get the count of the canonical form of a kmer, returning 0 if below threshold.
    pub fn get_canonical_count(&self, kmer: &u64) -> u32 {
        let count = self.inner.get_canonical_count(kmer);
        if count >= self.min_count { count } else { 0 }
    }

    /// Iterate over all entries, including those below threshold.
    /// Callers that need filtering should check counts themselves
    /// (as find_oligos_in_kmers already does).
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u32)> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_and_basic_ops() {
        let kc = KmerCounts::new(&5);
        assert_eq!(kc.get_k(), 5);
        assert!(kc.is_empty());
        assert_eq!(kc.len(), 0);
        assert_eq!(kc.get_n_kmers(), 0);
    }

    #[test]
    fn test_insert_and_get() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&42, &3);
        assert!(!kc.is_empty());
        assert_eq!(kc.get_count(&42), 3);
        assert!(kc.contains(&42));
        assert!(!kc.contains(&99));
    }

    #[test]
    fn test_insert_accumulates() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&42, &3);
        kc.insert(&42, &7);
        assert_eq!(kc.get_count(&42), 10);
        assert_eq!(kc.len(), 1);
    }

    #[test]
    fn test_saturating_add() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &u32::MAX);
        kc.insert(&1, &1);
        assert_eq!(kc.get_count(&1), u32::MAX);
    }

    #[test]
    fn test_extend_merges() {
        let mut a = KmerCounts::new(&5);
        a.insert(&1, &10);
        a.insert(&2, &20);

        let mut b = KmerCounts::new(&5);
        b.insert(&2, &5);
        b.insert(&3, &15);

        a.extend(&b).unwrap();
        assert_eq!(a.get_count(&1), 10);
        assert_eq!(a.get_count(&2), 25);
        assert_eq!(a.get_count(&3), 15);
    }

    #[test]
    fn test_extend_different_k_fails() {
        let mut a = KmerCounts::new(&5);
        let b = KmerCounts::new(&7);
        assert!(a.extend(&b).is_err());
    }

    #[test]
    fn test_median_odd() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &10);
        kc.insert(&2, &20);
        kc.insert(&3, &30);
        assert_eq!(kc.get_median_count(), 20);
    }

    #[test]
    fn test_median_even() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &10);
        kc.insert(&2, &20);
        // Median of [10, 20] = (10/2 + 20/2) = 15
        assert_eq!(kc.get_median_count(), 15);
    }

    #[test]
    fn test_median_empty() {
        let kc = KmerCounts::new(&5);
        assert_eq!(kc.get_median_count(), 0);
    }

    #[test]
    fn test_get_max_count() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &5);
        kc.insert(&2, &100);
        kc.insert(&3, &50);
        assert_eq!(kc.get_max_count(), 100);
    }

    #[test]
    fn test_remove_low_count_kmers() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &1);
        kc.insert(&2, &5);
        kc.insert(&3, &10);
        kc.remove_low_count_kmers(&5);
        assert!(!kc.contains(&1));
        assert!(kc.contains(&2));
        assert!(kc.contains(&3));
    }

    #[test]
    fn test_filtered_view_hides_low_counts() {
        let mut kc = KmerCounts::new(&5);
        kc.insert(&1, &2);
        kc.insert(&2, &10);

        let fv = kc.filtered_view(5);
        assert!(fv.get_canonical(&1).is_none());
        assert_eq!(fv.get_canonical(&2), Some(10));
        assert_eq!(fv.get_canonical_count(&1), 0);
        assert_eq!(fv.get_canonical_count(&2), 10);
    }

    #[test]
    fn test_ingest_seq() {
        let mut kc = KmerCounts::new(&3);
        // ACGT has 2 3-mers: ACG and CGT. These are reverse complements
        // of each other, so canonical storage keeps only 1 entry with count 2.
        kc.ingest_seq("ACGT").unwrap();
        assert_eq!(kc.get_n_unique_kmers(), 1);
        assert_eq!(kc.get_n_kmers(), 2);
    }

    #[test]
    fn test_extend_with_histogram() {
        let mut a = KmerCounts::new(&5);
        a.insert(&1, &3);

        let mut b = KmerCounts::new(&5);
        b.insert(&1, &2);
        b.insert(&2, &5);

        let histo_max: u64 = 100;
        let mut histo = Histogram::new(&histo_max);
        // Pre-populate histogram to match a's state
        histo.move_count(0, 3); // kmer 1 has count 3
        a.extend_with_histogram(&b, &mut histo).unwrap();

        assert_eq!(a.get_count(&1), 5);
        assert_eq!(a.get_count(&2), 5);
    }
}
