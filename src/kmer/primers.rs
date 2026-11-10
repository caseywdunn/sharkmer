//! Directional primer candidates, separate from canonical global counts.

use super::identity::{KmerSize, OrientedKey};
use ahash::AHashMap;
use anyhow::{Context, Result};

/// A local directional collection. Reverse complements are distinct keys and
/// counts retain the exact store's precision.
pub struct PrimerKmers {
    k: KmerSize,
    entries: AHashMap<u64, u64>,
}

impl PrimerKmers {
    /// Construct an empty collection with validated identity width.
    pub fn new(k: &usize) -> Result<Self> {
        Ok(Self {
            k: KmerSize::new(*k)?,
            entries: AHashMap::new(),
        })
    }

    /// Add directional abundance without saturation or mutation on overflow.
    pub fn insert(&mut self, bits: &u64, count: &u64) -> Result<()> {
        OrientedKey::new(*bits, self.k)?;
        let sum = self
            .entries
            .get(bits)
            .copied()
            .unwrap_or(0)
            .checked_add(*count)
            .context("Primer count overflow")?;
        self.entries.insert(*bits, sum);
        Ok(())
    }

    /// Read a directional entry, without canonicalizing it.
    pub fn get(&self, bits: &u64) -> Option<&u64> {
        self.entries.get(bits)
    }

    /// Visit the retained directional entries in unspecified order.
    pub fn iter(&self) -> impl Iterator<Item = (&u64, &u64)> {
        self.entries.iter()
    }

    /// Return the local candidate cardinality.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Report whether no primer candidates were retained.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return this collection's k-mer length.
    #[cfg(test)]
    pub fn get_k(&self) -> usize {
        self.k.get()
    }

    /// Collect local directional keys for deterministic graph seeding.
    pub fn kmers(&self) -> Vec<u64> {
        self.entries.keys().copied().collect()
    }

    /// Return the largest local count, or zero for an empty collection.
    pub fn get_max_count(&self) -> u64 {
        self.entries.values().copied().max().unwrap_or(0)
    }

    /// Return the integer median, averaging even cardinalities without overflow.
    pub fn get_median_count(&self) -> u64 {
        let mut counts: Vec<u64> = self.entries.values().copied().collect();
        if counts.is_empty() {
            return 0;
        }
        let mid = counts.len() / 2;
        counts.select_nth_unstable(mid);
        let upper = counts[mid];
        if mid * 2 != counts.len() {
            upper
        } else {
            let lower = counts[..mid].iter().copied().max().unwrap_or(upper);
            lower + (upper - lower) / 2
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directional_counts_keep_strands_and_u64_precision() {
        let mut kmers = PrimerKmers::new(&3).unwrap();
        kmers.insert(&0, &u64::MAX).unwrap();
        kmers.insert(&63, &(u64::MAX - 2)).unwrap();
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers.get(&0), Some(&u64::MAX));
        assert_eq!(kmers.get(&63), Some(&(u64::MAX - 2)));
        assert_eq!(kmers.get_median_count(), u64::MAX - 1);
        assert!(kmers.insert(&0, &1).is_err());
        assert!(kmers.insert(&64, &1).is_err());
        assert_eq!(kmers.get(&0), Some(&u64::MAX));
        assert_eq!(kmers.len(), 2);
    }
}
