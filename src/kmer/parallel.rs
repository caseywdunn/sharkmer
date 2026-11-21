//! Bounded coordinator/owner ingestion without per-read allocations or table merge.

use super::adaptive::{AdaptiveCountBuilder, AdaptiveStore};
use super::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use super::identity::{CanonicalKey, KmerSize, OrientedKey};
use super::store::{CountMetadata, CountStore, MemoryCountBuilder, MemoryCountStore};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::num::NonZeroU64;
use std::sync::{Arc, mpsc};
use std::thread::JoinHandle;

const BATCH_KEYS: usize = 1024;
const WORK_QUEUE_BATCHES: usize = 1;
const BUFFERS_PER_WORKER: usize = WORK_QUEUE_BATCHES + 2;
// Conservative fixed channel/descriptor overhead, separate from key allocations.
const CHANNEL_BYTES_PER_WORKER: usize = 4096;

fn owner(bits: u64, workers: usize) -> usize {
    if workers == 1 {
        return 0;
    }
    let mut value = bits;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    ((value ^ (value >> 31)) % workers as u64) as usize
}

/// Statically selected owner implementation; ordinary memory queues retain
/// concrete table calls without virtual dispatch in the counting loop.
pub(crate) trait OwnerCounter: Send + 'static {
    type Store: CountStore;
    fn accept(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()>;
    fn complete(self, digest: [u8; 32]) -> Result<Self::Store>;
}

impl OwnerCounter for MemoryCountBuilder {
    type Store = MemoryCountStore;
    #[inline]
    fn accept(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        self.add(key, count)
    }
    fn complete(self, digest: [u8; 32]) -> Result<Self::Store> {
        Ok(self.finalize(digest))
    }
}

impl OwnerCounter for AdaptiveCountBuilder {
    type Store = AdaptiveStore;
    #[inline]
    fn accept(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        self.add(key, count)
    }
    fn complete(self, digest: [u8; 32]) -> Result<Self::Store> {
        self.finish(digest)
    }
}

struct Worker<B: OwnerCounter> {
    work: Option<mpsc::SyncSender<Vec<u64>>>,
    recycled: Option<mpsc::Receiver<Vec<u64>>>,
    current: Vec<u64>,
    handle: Option<JoinHandle<Result<B>>>,
}

/// Owners that retain memory-only tables, used by the existing comparison.
pub(crate) type ParallelCountBuilder = CountPool<MemoryCountBuilder>;
/// Owners that transfer their own accepted prefix to disk on memory pressure.
pub(crate) type AdaptiveParallelCountBuilder = CountPool<AdaptiveCountBuilder>;

/// Emitted only when background counting actually ran.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct ParallelMetrics {
    pub(crate) counting_workers: usize,
    pub(crate) keys_per_batch: usize,
    pub(crate) queued_batches_per_worker: usize,
    pub(crate) reusable_key_buffers: usize,
    pub(crate) owned_worker_key_bytes: u64,
    pub(crate) reserved_worker_bytes: u64,
}

/// Worker-owned exact tables share one resize/accounting budget.
pub(crate) struct CountPool<B: OwnerCounter> {
    k: KmerSize,
    workers: Vec<Worker<B>>,
    budget: Arc<ResourceBudget>,
    occurrences: u64,
    failed: bool,
    metrics: ParallelMetrics,
    // Drop joins every worker and drops its channels before this guard releases.
    _pool_reservation: MemoryReservation,
}

impl<B: OwnerCounter> CountPool<B> {
    /// Reserve the fixed key pools, worker vector and channel overhead.
    pub(crate) fn reservation_bytes(workers: usize) -> Result<u64> {
        let each = BUFFERS_PER_WORKER
            .checked_mul(BATCH_KEYS)
            .and_then(|n| n.checked_mul(size_of::<u64>()))
            .and_then(|n| n.checked_add(size_of::<Worker<B>>() + CHANNEL_BYTES_PER_WORKER))
            .context("Counting worker byte size overflow")?;
        u64::try_from(
            workers
                .checked_mul(each)
                .context("Counting queue byte size overflow")?,
        )?
        .checked_add(Self::descriptor_bytes(workers)?)
        .context("Counting queue byte size overflow")
    }

    fn descriptor_bytes(workers: usize) -> Result<u64> {
        u64::try_from(
            workers
                .checked_mul(size_of::<B::Store>())
                .context("Count-shard descriptor size overflow")?,
        )
        .context("Count-shard descriptor size overflow")
    }

    fn with_factory(
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        workers: usize,
        mut create: impl FnMut() -> Result<B>,
    ) -> Result<Self> {
        ensure!(workers > 0, "Parallel counting requires at least one owner");
        let bytes = Self::reservation_bytes(workers)?
            .checked_sub(Self::descriptor_bytes(workers)?)
            .context("Counting queue reservation underflow")?;
        let reservation = budget.reserve(MemoryClass::WorkerQueues, bytes)?;
        let mut ports = Vec::new();
        ports
            .try_reserve_exact(workers)
            .context("Unable to reserve counting worker descriptors")?;
        ensure!(
            ports.capacity() == workers,
            "Counting descriptor allocation differs from reservation"
        );
        let mut result = Self {
            k,
            workers: ports,
            budget,
            occurrences: 0,
            failed: false,
            metrics: ParallelMetrics {
                counting_workers: workers,
                keys_per_batch: BATCH_KEYS,
                queued_batches_per_worker: WORK_QUEUE_BATCHES,
                reusable_key_buffers: workers
                    .checked_mul(BUFFERS_PER_WORKER)
                    .context("Counting buffer count overflow")?,
                owned_worker_key_bytes: 0,
                reserved_worker_bytes: bytes,
            },
            _pool_reservation: reservation,
        };
        let one = NonZeroU64::new(1).context("Positive count increment is zero")?;
        for index in 0..workers {
            let builder = create()?;
            let allocate = || -> Result<Vec<u64>> {
                let mut buffer = Vec::new();
                buffer
                    .try_reserve_exact(BATCH_KEYS)
                    .context("Unable to reserve counting key buffer")?;
                ensure!(
                    buffer.capacity() == BATCH_KEYS,
                    "Counting buffer allocation differs from reservation"
                );
                Ok(buffer)
            };
            let current = allocate()?;
            let (work_sender, work_receiver) = mpsc::sync_channel::<Vec<u64>>(WORK_QUEUE_BATCHES);
            let (recycle_sender, recycle_receiver) =
                mpsc::sync_channel::<Vec<u64>>(BUFFERS_PER_WORKER - 1);
            for _ in 1..BUFFERS_PER_WORKER {
                recycle_sender
                    .send(allocate()?)
                    .context("Unable to initialize counting buffer pool")?;
            }
            let handle = std::thread::Builder::new()
                .name(format!("count-owner-{index}"))
                .spawn(move || {
                    let mut builder = builder;
                    while let Ok(mut buffer) = work_receiver.recv() {
                        for &bits in &buffer {
                            // The coordinator checked k and accepted a typed
                            // canonical key before it entered this owned pool.
                            builder.accept(CanonicalKey::from_stored(bits, k), one)?;
                        }
                        buffer.clear();
                        // During finalization the coordinator closes recycling while
                        // owners still drain every already-queued batch.
                        let _ = recycle_sender.send(buffer);
                    }
                    Ok(builder)
                })
                .context("Unable to start counting owner")?;
            result.workers.push(Worker {
                work: Some(work_sender),
                recycled: Some(recycle_receiver),
                current,
                handle: Some(handle),
            });
        }
        result.metrics.owned_worker_key_bytes = u64::try_from(
            workers
                .checked_mul(BUFFERS_PER_WORKER * BATCH_KEYS * size_of::<u64>())
                .context("Counting key-pool byte size overflow")?,
        )?;
        Ok(result)
    }

    /// Checked total of coordinator-accepted keys, verified again after joins.
    pub(crate) fn occurrences(&self) -> u64 {
        self.occurrences
    }

    fn dispatch(worker: &mut Worker<B>, recycle: bool) -> Result<()> {
        if worker.current.is_empty() {
            return Ok(());
        }
        worker
            .work
            .as_ref()
            .context("Counting queue is closed")?
            .send(std::mem::take(&mut worker.current))
            .context("Counting owner stopped before accepting a batch")?;
        if recycle {
            worker.current = worker
                .recycled
                .as_ref()
                .context("Counting buffer pool is closed")?
                .recv()
                .context("Counting owner stopped before returning a buffer")?;
            ensure!(
                worker.current.is_empty() && worker.current.capacity() == BATCH_KEYS,
                "Invalid recycled counting buffer"
            );
        }
        Ok(())
    }

    fn owner_error(worker: &mut Worker<B>, channel_error: anyhow::Error) -> anyhow::Error {
        drop(worker.work.take());
        drop(worker.recycled.take());
        let Some(handle) = worker.handle.take() else {
            return channel_error;
        };
        match handle.join() {
            Ok(Err(error)) => error.context("Counting owner failed"),
            Err(_) => anyhow::anyhow!("Counting owner panicked"),
            Ok(Ok(builder)) => {
                drop(builder);
                channel_error.context("Counting owner ended unexpectedly")
            }
        }
    }

    /// Route one canonical window without allocating a new read or batch.
    #[inline]
    pub(crate) fn add(&mut self, key: CanonicalKey) -> Result<()> {
        ensure!(!self.failed, "Parallel count producer previously failed");
        ensure!(
            key.size() == self.k,
            "Count key k differs from counting owners"
        );
        let Some(total) = self.occurrences.checked_add(1) else {
            anyhow::bail!("Exact count occurrence total overflow");
        };
        let index = owner(key.bits(), self.workers.len());
        let worker = &mut self.workers[index];
        worker.current.push(key.bits());
        if worker.current.len() == BATCH_KEYS
            && let Err(error) = Self::dispatch(worker, true)
        {
            self.failed = true;
            return Err(Self::owner_error(worker, error));
        }
        self.occurrences = total;
        Ok(())
    }

    fn close_channels(&mut self) {
        for worker in &mut self.workers {
            drop(worker.work.take());
            drop(worker.recycled.take());
        }
    }

    fn finish_owned(mut self, digest: [u8; 32]) -> Result<CompletedOwners<B::Store>> {
        ensure!(!self.failed, "Parallel count producer previously failed");
        for worker in &mut self.workers {
            if let Err(error) = Self::dispatch(worker, false) {
                return Err(Self::owner_error(worker, error));
            }
        }
        self.close_channels();
        let mut metadata_reservation = self.budget.reserve(
            MemoryClass::WorkerQueues,
            Self::descriptor_bytes(self.workers.len())?,
        )?;
        let mut stores = Vec::new();
        stores
            .try_reserve_exact(self.workers.len())
            .context("Unable to reserve completed count shards")?;
        metadata_reservation.resize(u64::try_from(
            stores
                .capacity()
                .checked_mul(size_of::<B::Store>())
                .context("Count-shard descriptor size overflow")?,
        )?)?;
        for worker in &mut self.workers {
            let builder = worker
                .handle
                .take()
                .context("Counting owner handle is missing")?
                .join()
                .map_err(|_| anyhow::anyhow!("Counting owner panicked"))??;
            stores.push(builder.complete(digest)?);
        }
        let occurrences = stores.iter().try_fold(0u64, |sum, store| {
            sum.checked_add(store.metadata().occurrences)
                .context("Count-shard occurrence overflow")
        })?;
        ensure!(
            occurrences == self.occurrences,
            "Counting owners differ from accepted windows"
        );
        Ok(CompletedOwners {
            stores,
            metadata_reservation,
            metrics: self.metrics,
        })
    }
}

struct CompletedOwners<S> {
    stores: Vec<S>,
    metadata_reservation: MemoryReservation,
    metrics: ParallelMetrics,
}

impl CountPool<MemoryCountBuilder> {
    /// Start memory-only owners under their shared table and queue budgets.
    pub(crate) fn new(
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        workers: usize,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        Self::with_factory(k, Arc::clone(&budget), workers, || {
            MemoryCountBuilder::with_budget(k, Arc::clone(&budget), max_lookup_batch)
        })
    }

    /// Drain and join owners, moving their concrete tables without consolidation.
    pub(crate) fn finish(self, digest: [u8; 32]) -> Result<(ShardedCountStore, ParallelMetrics)> {
        let budget = Arc::clone(&self.budget);
        let complete = self.finish_owned(digest)?;
        let store = ShardedCountStore::new(
            complete.stores,
            budget,
            complete.metadata_reservation,
            digest,
        )?;
        Ok((store, complete.metrics))
    }
}

impl<B: OwnerCounter> Drop for CountPool<B> {
    fn drop(&mut self) {
        self.close_channels();
        for worker in &mut self.workers {
            if let Some(handle) = worker.handle.take() {
                // Even on another owner's error, join before returning the
                // shared key-pool reservation or letting the process continue.
                match handle.join() {
                    Ok(Ok(builder)) => drop(builder),
                    Ok(Err(error)) => log::debug!("Counting owner stopped: {error:#}"),
                    Err(_) => log::error!("Counting owner panicked during cleanup"),
                }
            }
        }
    }
}

/// Borrowable exact shards; lookup selects one full-key table, visits retain all.
pub(crate) struct ShardedCountStore {
    stores: Vec<MemoryCountStore>,
    metadata: CountMetadata,
    peak_table_bytes: usize,
    _metadata_reservation: MemoryReservation,
}

impl ShardedCountStore {
    fn new(
        stores: Vec<MemoryCountStore>,
        budget: Arc<ResourceBudget>,
        metadata_reservation: MemoryReservation,
        digest: [u8; 32],
    ) -> Result<Self> {
        let mut metadata = stores
            .first()
            .context("No completed counting owners")?
            .metadata()
            .clone();
        metadata.backend = "sharded_hashbrown_memory";
        metadata.selected_input_sha256 = digest;
        metadata.distinct_keys = 0;
        metadata.occurrences = 0;
        for store in &stores {
            ensure!(
                store.metadata().k == metadata.k,
                "Count-shard identity mismatch"
            );
            metadata.distinct_keys = metadata
                .distinct_keys
                .checked_add(store.metadata().distinct_keys)
                .context("Count-shard cardinality overflow")?;
            metadata.occurrences = metadata
                .occurrences
                .checked_add(store.metadata().occurrences)
                .context("Count-shard occurrence overflow")?;
        }
        let peak_table_bytes = usize::try_from(
            budget.snapshot().peak_reserved_by_class[MemoryClass::CountSlots as usize],
        )
        .context("Count-shard byte total exceeds address space")?;
        Ok(Self {
            stores,
            metadata,
            peak_table_bytes,
            _metadata_reservation: metadata_reservation,
        })
    }

    pub(crate) fn capacity(&self) -> usize {
        self.stores.iter().map(MemoryCountStore::capacity).sum()
    }
    pub(crate) fn table_bytes(&self) -> usize {
        self.stores.iter().map(MemoryCountStore::table_bytes).sum()
    }
    pub(crate) fn peak_table_bytes(&self) -> usize {
        self.peak_table_bytes
    }

    /// Concrete scan path allows assembly visitors to inline across shards.
    pub(crate) fn visit_entries<F: FnMut(CanonicalKey, NonZeroU64) -> Result<()> + ?Sized>(
        &self,
        visitor: &mut F,
    ) -> Result<()> {
        for store in &self.stores {
            store.visit_entries(visitor)?;
        }
        Ok(())
    }
}

impl CountStore for ShardedCountStore {
    fn metadata(&self) -> &CountMetadata {
        &self.metadata
    }
    fn lookup(&self, key: CanonicalKey) -> Result<Option<NonZeroU64>> {
        ensure!(
            key.size() == self.metadata.k,
            "Count lookup k differs from shard k"
        );
        self.stores[owner(key.bits(), self.stores.len())].lookup(key)
    }
    fn lookup_oriented(&self, key: OrientedKey) -> Result<Option<NonZeroU64>> {
        if self.stores.len() == 1 {
            return self.stores[0].lookup_oriented(key);
        }
        self.lookup(key.canonical())
    }
    fn visit(&self, visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>) -> Result<()> {
        self.visit_entries(visitor)
    }
}

mod adaptive;
pub(crate) use adaptive::{AdaptiveShardMetrics, AdaptiveShardedCountStore};

#[cfg(test)]
mod tests;
