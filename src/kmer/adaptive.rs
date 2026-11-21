//! One-pass exact counting that spills an accepted memory prefix on pressure.
//!
//! The transition streams weighted entries, then retries the unaccepted update.
//! It never reopens the original input or constructs a second global table.

use super::budget::{MemoryClass, ResourceBudget, ResourceExhausted, ResourceKind};
use super::external::{ExternalCountBuilder, ExternalMetrics};
use super::identity::{CanonicalKey, KmerSize};
use super::sorted_file::DiskCountStore;
use super::store::{CountStore, MemoryCountBuilder, MemoryCountStore};
use anyhow::{Result, ensure};
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};
use std::sync::Arc;

enum State {
    Memory(MemoryCountBuilder),
    External(ExternalCountBuilder),
}

/// Actual final storage; fitting inputs never allocate external sort buffers.
pub(crate) enum AdaptiveStore {
    Memory(MemoryCountStore),
    Disk(DiskCountStore, ExternalMetrics),
}

/// Serial development candidate for lossless automatic spilling.
pub(crate) struct AdaptiveCountBuilder {
    state: Option<State>,
    k: KmerSize,
    directory: PathBuf,
    budget: Arc<ResourceBudget>,
    sort_bytes: usize,
    max_lookup_batch: usize,
    occurrences: u64,
    peak_table_bytes: usize,
}

impl AdaptiveCountBuilder {
    /// Keep the supplied budget live across table transfer and disk production.
    pub(crate) fn new(
        directory: &Path,
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        sort_bytes: usize,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        ensure!(
            sort_bytes >= 16,
            "External sort buffer must hold a key/count pair"
        );
        let memory = MemoryCountBuilder::with_budget(k, Arc::clone(&budget), max_lookup_batch)?;
        Ok(Self {
            state: Some(State::Memory(memory)),
            k,
            directory: directory.to_path_buf(),
            budget,
            sort_bytes,
            max_lookup_batch,
            occurrences: 0,
            peak_table_bytes: 0,
        })
    }

    /// Checked total of accepted updates, including the transferred prefix.
    pub(crate) fn occurrences(&self) -> u64 {
        self.occurrences
    }

    /// Declared key size, retained when moving this producer to another owner.
    pub(crate) fn size(&self) -> KmerSize {
        self.k
    }

    /// The same ledger must follow this producer through owner changes.
    pub(crate) fn resource_budget(&self) -> &Arc<ResourceBudget> {
        &self.budget
    }

    /// Current in-memory allocation; a spilled/failed producer has no table.
    pub(crate) fn memory_table_bytes(&self) -> Option<usize> {
        match &self.state {
            Some(State::Memory(builder)) => Some(builder.table_bytes()),
            _ => None,
        }
    }

    /// Spill only typed count-growth pressure; all other errors propagate.
    #[inline]
    pub(crate) fn add(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        let result = match self.state.as_mut() {
            Some(State::Memory(builder)) => builder.add(key, count),
            Some(State::External(builder)) => builder.add(key, count),
            None => anyhow::bail!("Adaptive count producer previously failed"),
        };
        if let Err(error) = result {
            let count_pressure = matches!(self.state, Some(State::Memory(_)))
                && error
                    .downcast_ref::<ResourceExhausted>()
                    .is_some_and(|error| {
                        matches!(
                            error.resource,
                            ResourceKind::StageMemory(MemoryClass::CountSlots)
                                | ResourceKind::RunMemory
                        )
                    });
            if !count_pressure {
                return Err(error);
            }
            // Taking state first poisons publication if creating, transferring,
            // verifying or retrying fails. The previous add accepted no update.
            self.spill_and_retry(key, count)?;
        }
        self.occurrences = match self.state.as_ref() {
            Some(State::Memory(builder)) => builder.occurrences(),
            Some(State::External(builder)) => builder.occurrences(),
            None => anyhow::bail!("Adaptive count producer previously failed"),
        };
        Ok(())
    }

    fn spill_and_retry(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        let Some(State::Memory(builder)) = self.state.take() else {
            anyhow::bail!("Adaptive spill requires a live memory prefix");
        };
        let prefix = builder.finalize([0; 32]);
        self.peak_table_bytes = prefix.peak_table_bytes();
        ensure!(
            prefix.metadata().occurrences == self.occurrences,
            "Adaptive prefix total differs"
        );
        let mut external = ExternalCountBuilder::new(
            &self.directory,
            self.k,
            Arc::clone(&self.budget),
            self.sort_bytes,
            self.max_lookup_batch,
        )?;
        prefix.visit(&mut |key, count| external.add(key, count))?;
        ensure!(
            external.occurrences() == self.occurrences,
            "Adaptive transfer lost counts"
        );
        // Free the original table before admitting the next selected window.
        drop(prefix);
        external.add(key, count)?;
        self.state = Some(State::External(external));
        Ok(())
    }

    /// Move final storage and bind the selected input once, after complete ingest.
    pub(crate) fn finish(self, digest: [u8; 32]) -> Result<AdaptiveStore> {
        match self.state {
            Some(State::Memory(builder)) => Ok(AdaptiveStore::Memory(builder.finalize(digest))),
            Some(State::External(builder)) => {
                let (store, mut metrics) = builder.finish(digest)?;
                metrics.peak_table_bytes = self.peak_table_bytes;
                Ok(AdaptiveStore::Disk(store, metrics))
            }
            None => anyhow::bail!("Adaptive count producer previously failed"),
        }
    }
}

#[cfg(test)]
mod tests;

impl CountStore for AdaptiveStore {
    fn metadata(&self) -> &super::store::CountMetadata {
        match self {
            Self::Memory(store) => store.metadata(),
            Self::Disk(store, _) => store.metadata(),
        }
    }

    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        match self {
            Self::Memory(store) => store.lookup(key),
            Self::Disk(store, _) => store.lookup(key),
        }
    }

    fn lookup_oriented(&self, key: super::identity::OrientedKey) -> Result<Option<NonZeroU64>> {
        match self {
            Self::Memory(store) => store.lookup_oriented(key),
            Self::Disk(store, _) => store.lookup_oriented(key),
        }
    }

    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        match self {
            Self::Memory(store) => store.visit(visitor),
            Self::Disk(store, _) => store.visit(visitor),
        }
    }
}
