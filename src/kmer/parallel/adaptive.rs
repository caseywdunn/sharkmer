//! Independently spillable owners share count, sort, page, queue and disk limits.

use super::*;
use std::path::Path;

impl CountPool<AdaptiveCountBuilder> {
    /// Move a complete accepted prefix into one owner without copying its
    /// table. The pool census includes that prefix before accepting new keys.
    pub(crate) fn from_prefix(prefix: AdaptiveCountBuilder) -> Result<Self> {
        let k = prefix.size();
        let budget = Arc::clone(prefix.resource_budget());
        let occurrences = prefix.occurrences();
        let mut prefix = Some(prefix);
        let mut pool = Self::with_factory(k, budget, 1, || {
            prefix.take().context("Count prefix already transferred")
        })?;
        pool.occurrences = occurrences;
        Ok(pool)
    }

    /// Divide the shared sort allocation between owners; never multiply it.
    pub(crate) fn adaptive(
        directory: &Path,
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        workers: usize,
        sort_bytes: usize,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        ensure!(workers > 0, "Parallel counting requires at least one owner");
        let each_sort = sort_bytes / workers / 16 * 16;
        ensure!(
            each_sort >= 16,
            "Shared sort buffer cannot supply every counting owner"
        );
        Self::with_factory(k, Arc::clone(&budget), workers, || {
            AdaptiveCountBuilder::new(
                directory,
                k,
                Arc::clone(&budget),
                each_sort,
                max_lookup_batch,
            )
        })
    }

    /// Join every owner and retain full-key memory/disk shards without merging.
    pub(crate) fn finish_adaptive(
        self,
        digest: [u8; 32],
    ) -> Result<(AdaptiveShardedCountStore, ParallelMetrics)> {
        let budget = Arc::clone(&self.budget);
        let completed = self.finish_owned(digest)?;
        let store = AdaptiveShardedCountStore::new(
            completed.stores,
            budget,
            completed.metadata_reservation,
            digest,
        )?;
        Ok((store, completed.metrics))
    }
}

/// Actual final storage composition, distinct from biological sequence output.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct AdaptiveShardMetrics {
    pub(crate) memory_shards: usize,
    pub(crate) disk_shards: usize,
    pub(crate) count_file_bytes: u64,
    pub(crate) lookup_cache_bytes: u64,
    pub(crate) sorted_run_count: u64,
    pub(crate) merge_count: u64,
}

/// Disjoint owner storage; membership probes exactly one complete-key shard.
pub(crate) struct AdaptiveShardedCountStore {
    stores: Vec<AdaptiveStore>,
    metadata: CountMetadata,
    metrics: AdaptiveShardMetrics,
    peak_table_bytes: usize,
    _metadata_reservation: MemoryReservation,
}

impl AdaptiveShardedCountStore {
    fn new(
        stores: Vec<AdaptiveStore>,
        budget: Arc<ResourceBudget>,
        reservation: MemoryReservation,
        digest: [u8; 32],
    ) -> Result<Self> {
        let mut metadata = stores
            .first()
            .context("No completed adaptive owners")?
            .metadata()
            .clone();
        metadata.backend = "adaptive_sharded_exact";
        metadata.hash_algorithm = super::super::store::TABLE_HASHER_NAME;
        metadata.selected_input_sha256 = digest;
        metadata.distinct_keys = 0;
        metadata.occurrences = 0;
        let mut metrics = AdaptiveShardMetrics {
            memory_shards: 0,
            disk_shards: 0,
            count_file_bytes: 0,
            lookup_cache_bytes: 0,
            sorted_run_count: 0,
            merge_count: 0,
        };
        let add = |sum: &mut u64, value: u64| -> Result<()> {
            *sum = sum
                .checked_add(value)
                .context("Adaptive shard summary overflow")?;
            Ok(())
        };
        for store in &stores {
            ensure!(
                store.metadata().k == metadata.k
                    && store.metadata().selected_input_sha256 == digest,
                "Adaptive shard identity mismatch"
            );
            add(&mut metadata.distinct_keys, store.metadata().distinct_keys)?;
            add(&mut metadata.occurrences, store.metadata().occurrences)?;
            match store {
                AdaptiveStore::Memory(_) => metrics.memory_shards += 1,
                AdaptiveStore::Disk(store, external) => {
                    metrics.disk_shards += 1;
                    add(&mut metrics.count_file_bytes, store.file_bytes())?;
                    add(
                        &mut metrics.lookup_cache_bytes,
                        u64::try_from(store.cache_bytes())?,
                    )?;
                    add(&mut metrics.sorted_run_count, external.sorted_runs)?;
                    add(&mut metrics.merge_count, external.merges)?;
                }
            }
        }
        let peak_table_bytes = usize::try_from(
            budget.snapshot().peak_reserved_by_class[MemoryClass::CountSlots as usize],
        )?;
        Ok(Self {
            stores,
            metadata,
            metrics,
            peak_table_bytes,
            _metadata_reservation: reservation,
        })
    }

    /// Actual owned memory table capacity; disk shards add no memory capacity.
    pub(crate) fn capacity(&self) -> usize {
        self.stores
            .iter()
            .map(|store| match store {
                AdaptiveStore::Memory(store) => store.capacity(),
                AdaptiveStore::Disk(..) => 0,
            })
            .sum()
    }

    /// Current owned count-table buffers, excluding disk and lookup caches.
    pub(crate) fn table_bytes(&self) -> usize {
        self.stores
            .iter()
            .map(|store| match store {
                AdaptiveStore::Memory(store) => store.table_bytes(),
                AdaptiveStore::Disk(..) => 0,
            })
            .sum()
    }

    /// Combined live table peak from the shared ledger, including resize overlap.
    pub(crate) fn peak_table_bytes(&self) -> usize {
        self.peak_table_bytes
    }

    /// Actual count storage composition and external work.
    pub(crate) fn metrics(&self) -> AdaptiveShardMetrics {
        self.metrics
    }

    /// Try every disk close before returning an error; memory guards drop too.
    pub(crate) fn close(self) -> Result<()> {
        let mut failure = None;
        for store in self.stores {
            if let AdaptiveStore::Disk(store, _) = store
                && let Err(error) = store.close()
                && failure.is_none()
            {
                failure = Some(error);
            }
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    }
}

impl CountStore for AdaptiveShardedCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }

    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        ensure!(
            key.size() == self.metadata.k,
            "Count lookup k differs from adaptive shards"
        );
        self.stores[owner(key.bits(), self.stores.len())].lookup(key)
    }

    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        for store in &self.stores {
            store.visit(visitor)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
