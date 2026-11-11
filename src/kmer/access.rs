//! Fallible assembly access without changing legacy count precision.

use super::counting::KmerCounts;
use super::identity::{CanonicalKey, KmerSize, OrientedKey};
use super::store::CountStore;
use anyhow::Result;

enum Source<'a> {
    Legacy(&'a KmerCounts),
    Exact(&'a dyn CountStore),
}

/// A borrowed assembly view. Legacy counts retain their released u32 policy;
/// exact stores retain u64 values and propagate access failures.
pub struct FilteredKmerCounts<'a> {
    source: Source<'a>,
    minimum: u64,
}

impl<'a> FilteredKmerCounts<'a> {
    /// Adapt the legacy counter without claiming it is an exact count store.
    pub fn legacy(counts: &'a KmerCounts, minimum: u32) -> Self {
        Self {
            source: Source::Legacy(counts),
            minimum: u64::from(minimum),
        }
    }

    /// Borrow a finalized exact store without collecting another global table.
    pub fn from_store(store: &'a dyn CountStore, minimum: u64) -> Self {
        Self {
            source: Source::Exact(store),
            minimum,
        }
    }

    /// Tighten the existing view's threshold without copying count data.
    pub fn with_minimum(&self, minimum: u64) -> Self {
        Self {
            source: match self.source {
                Source::Legacy(counts) => Source::Legacy(counts),
                Source::Exact(store) => Source::Exact(store),
            },
            minimum: self.minimum.max(minimum),
        }
    }

    /// Return the stored k-mer length.
    pub fn get_k(&self) -> usize {
        match self.source {
            Source::Legacy(counts) => counts.get_k(),
            Source::Exact(store) => store.metadata().k.get(),
        }
    }

    /// Resolve either orientation, separating absence from access failure.
    pub fn get_canonical(&self, bits: &u64) -> Result<Option<u64>> {
        let count = match self.source {
            Source::Legacy(counts) => counts.get_canonical(bits).copied().map(u64::from),
            Source::Exact(store) => store
                .lookup_oriented(OrientedKey::new(*bits, KmerSize::new(self.get_k())?)?)?
                .map(|count| count.get()),
        };
        Ok(count.filter(|count| *count >= self.minimum))
    }

    /// Resolve canonical identity; absent or threshold-filtered keys return zero.
    /// Access and identity errors remain errors.
    pub fn get_canonical_count(&self, bits: &u64) -> Result<u64> {
        let count = match self.source {
            Source::Legacy(counts) => u64::from(counts.get_canonical_count(bits)),
            Source::Exact(store) => store
                .lookup(CanonicalKey::new(*bits, KmerSize::new(self.get_k())?)?)?
                .map_or(0, |count| count.get()),
        };
        Ok(if count >= self.minimum { count } else { 0 })
    }

    /// Visit threshold-admitted entries once, propagating visitor/store errors.
    pub fn visit(&self, visitor: &mut dyn FnMut(u64, u64) -> Result<()>) -> Result<()> {
        let mut admitted = |bits, count| {
            if count >= self.minimum {
                visitor(bits, count)?;
            }
            Ok(())
        };
        match self.source {
            Source::Legacy(counts) => {
                for (&bits, &count) in counts.iter() {
                    admitted(bits, u64::from(count))?;
                }
                Ok(())
            }
            Source::Exact(store) => {
                store.visit(&mut |key, count| admitted(key.bits(), count.get()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::store::MemoryCountBuilder;
    use std::num::NonZeroU64;

    #[test]
    fn exact_view_preserves_large_counts_strands_thresholds_and_errors() {
        let k = KmerSize::new(3).unwrap();
        let mut builder = MemoryCountBuilder::new(k, 4096, 8).unwrap();
        let large = u64::from(u32::MAX) + 100;
        builder
            .add(
                CanonicalKey::new(0, k).unwrap(),
                NonZeroU64::new(large).unwrap(),
            )
            .unwrap();
        builder
            .add(
                CanonicalKey::new(1, k).unwrap(),
                NonZeroU64::new(2).unwrap(),
            )
            .unwrap();
        let store = builder.finalize([0; 32]);
        let view = FilteredKmerCounts::from_store(&store, 3);
        assert_eq!(view.get_canonical(&0).unwrap(), Some(large));
        assert_eq!(view.get_canonical(&63).unwrap(), Some(large));
        assert_eq!(view.get_canonical_count(&63).unwrap(), large);
        assert_eq!(view.get_canonical(&1).unwrap(), None);
        assert_eq!(view.get_canonical_count(&2).unwrap(), 0);
        assert!(view.get_canonical(&64).is_err());
        let mut observed = Vec::new();
        view.visit(&mut |key, count| {
            observed.push((key, count));
            Ok(())
        })
        .unwrap();
        assert_eq!(observed, [(0, large)]);
        assert!(
            view.visit(&mut |_, _| anyhow::bail!("visitor failed"))
                .is_err()
        );
        assert_eq!(
            view.with_minimum(large + 1).get_canonical(&0).unwrap(),
            None
        );
        assert_eq!(view.with_minimum(1).get_canonical(&1).unwrap(), None);
    }
}
