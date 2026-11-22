//! Directional primer candidates, separate from canonical global counts.

use super::budget::ResourceBudget;
use super::budget::collections::{BudgetMap, BudgetVec};
use super::identity::{KmerSize, OrientedKey};
use anyhow::{Context, Result};
use std::sync::Arc;

/// A local directional collection. Reverse complements are distinct keys and
/// counts retain the exact store's precision.
pub struct PrimerKmers {
    k: KmerSize,
    entries: BudgetMap<u64, u64>,
}

impl PrimerKmers {
    /// Construct an empty collection with validated identity width.
    #[cfg(test)]
    pub fn new(k: &usize) -> Result<Self> {
        Self::with_budget(k, None)
    }

    /// Construct directional candidates sharing the active target ceiling.
    pub(crate) fn with_budget(k: &usize, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        Ok(Self {
            k: KmerSize::new(*k)?,
            entries: BudgetMap::new(budget),
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
        self.entries.insert(*bits, sum)?;
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

    /// Return the largest local count, or zero for an empty collection.
    pub fn get_max_count(&self) -> u64 {
        self.entries
            .iter()
            .map(|(_, count)| *count)
            .max()
            .unwrap_or(0)
    }

    /// Return the integer median, averaging even cardinalities without overflow.
    pub fn get_median_count(&self) -> Result<u64> {
        let mut counts = BudgetVec::collect(
            self.entries.iter().map(|(_, count)| *count),
            self.entries.budget(),
        )?;
        if counts.is_empty() {
            return Ok(0);
        }
        let mid = counts.len() / 2;
        counts.select_nth_unstable(mid);
        let upper = counts[mid];
        if mid * 2 != counts.len() {
            Ok(upper)
        } else {
            let lower = counts[..mid].iter().copied().max().unwrap_or(upper);
            Ok(lower + (upper - lower) / 2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primer_table_and_median_refusal_preserve_exact_candidates_and_retry() {
        use crate::kmer::budget::{MemoryClass, ResourceExhausted, ResourceLimits};
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1000,
            process_overhead_bytes: 0,
            class_bytes: [1000; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let mut candidates = PrimerKmers::with_budget(&5, Some(&budget)).unwrap();
        for key in 0..14 {
            candidates.insert(&key, &(key + 1)).unwrap();
        }
        let before = budget.snapshot().reserved_memory_bytes;
        let held = budget.reserve(MemoryClass::TargetGraphs, 300).unwrap();
        assert!(
            candidates
                .insert(&14, &15)
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(candidates.len(), 14);
        assert_eq!(budget.snapshot().reserved_memory_bytes, before + 300);
        assert!((0..14).all(|key| candidates.get(&key) == Some(&(key + 1))));
        drop(held);
        candidates.insert(&14, &15).unwrap();
        let held = budget.reserve(MemoryClass::TargetGraphs, 400).unwrap();
        assert!(
            candidates
                .get_median_count()
                .unwrap_err()
                .downcast_ref::<ResourceExhausted>()
                .is_some()
        );
        assert_eq!(candidates.len(), 15);
        assert_eq!(candidates.get_max_count(), 15);
        drop(held);
        assert_eq!(candidates.get_median_count().unwrap(), 8);
        drop(candidates);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn directional_counts_keep_strands_and_u64_precision() {
        let mut kmers = PrimerKmers::new(&3).unwrap();
        kmers.insert(&0, &u64::MAX).unwrap();
        kmers.insert(&63, &(u64::MAX - 2)).unwrap();
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers.get(&0), Some(&u64::MAX));
        assert_eq!(kmers.get(&63), Some(&(u64::MAX - 2)));
        assert_eq!(kmers.get_median_count().unwrap(), u64::MAX - 1);
        assert!(kmers.insert(&0, &1).is_err());
        assert!(kmers.insert(&64, &1).is_err());
        assert_eq!(kmers.get(&0), Some(&u64::MAX));
        assert_eq!(kmers.len(), 2);
    }
}
