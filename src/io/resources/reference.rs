//! Frozen numeric formulas from the verified dedup-workspace checkpoint.
//! Original src/io.rs SHA-256: 1c34c3bfb8649fd0955f1b1d4d4f77c721cd8815f82fdafba243e77be3e16a64.
//! Only returns/types/names are adapted; builder construction and OS census
//! are excluded so these tests compare the unchanged declared resource limits.
use super::*;

pub(super) fn total(
    mut classes: [u64; crate::kmer::budget::MEMORY_CLASSES],
    disk_bytes: u64,
    graph_bytes: u64,
) -> Result<ResourceLimits> {
    classes[MemoryClass::RetainedResults as usize] = 0;
    classes[MemoryClass::TargetGraphs as usize] = graph_bytes;
    let memory_bytes = classes.iter().try_fold(0u64, |sum, bytes| {
        sum.checked_add(*bytes)
            .context("Count/assembly resource budget overflow")
    })?;
    // Results share the existing combined ceiling; this class ceiling
    // does not add another allowance to the configured run memory.
    classes[MemoryClass::RetainedResults as usize] = memory_bytes;
    Ok(ResourceLimits {
        memory_bytes,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk_bytes,
        temporary_open_files: 128,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn table(
    _k: usize,
    table_bytes: u64,
    record_bytes: usize,
    disk_bytes: u64,
    replay_bytes: u64,
    _packed: bool,
    graph_buffer_bytes: u64,
) -> Result<ResourceLimits> {
    let input_bytes = u64::try_from(record_bytes)?
        .checked_mul(4)
        .context("FASTQ buffer budget overflow")?;
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    let key_bytes = u64::try_from(KMER_BATCH_CAPACITY * size_of::<CanonicalKey>())?;
    classes[MemoryClass::CountSlots as usize] = table_bytes;
    classes[MemoryClass::InputBuffer as usize] = input_bytes;
    classes[MemoryClass::WorkerQueues as usize] = key_bytes;
    classes[MemoryClass::ReplayBuffer as usize] = replay_bytes;
    total(classes, disk_bytes, graph_buffer_bytes)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn parallel(
    k: usize,
    table_bytes: u64,
    record_bytes: usize,
    threads: usize,
    disk_bytes: u64,
    replay_bytes: u64,
    graph_buffer_bytes: u64,
) -> Result<ResourceLimits> {
    if threads <= 1 {
        return table(
            k,
            table_bytes,
            record_bytes,
            disk_bytes,
            replay_bytes,
            false,
            graph_buffer_bytes,
        );
    }
    let workers = threads - 1;
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    let key_bytes = u64::try_from(KMER_BATCH_CAPACITY * size_of::<CanonicalKey>())?;
    classes[MemoryClass::CountSlots as usize] = table_bytes;
    classes[MemoryClass::ReplayBuffer as usize] = replay_bytes;
    classes[MemoryClass::InputBuffer as usize] = u64::try_from(record_bytes)?
        .checked_mul(4)
        .context("FASTQ buffer budget overflow")?;
    // A fixed development ceiling prevents a thread request from creating
    // an unbounded key pool. Whole-process/thread overhead is still separate.
    classes[MemoryClass::WorkerQueues as usize] = key_bytes
        .checked_add(ParallelCountBuilder::reservation_bytes(workers)?.min(64 << 20))
        .context("Counting worker reservation overflow")?;
    total(classes, disk_bytes, graph_buffer_bytes)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn spill(
    _k: usize,
    table_bytes: u64,
    sort_bytes: u64,
    record_bytes: usize,
    disk_bytes: u64,
    _directory: &std::path::Path,
    replay_bytes: u64,
    workers: usize,
    identity_worker: bool,
    owner_transfer_threshold: Option<usize>,
    graph_buffer_bytes: u64,
) -> Result<ResourceLimits> {
    ensure!(
        owner_transfer_threshold.is_none_or(|bytes| bytes > 0 && workers == 0 && identity_worker),
        "Owner transfer requires a positive threshold and one identity worker"
    );
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    let key_bytes = u64::try_from(KMER_BATCH_CAPACITY * size_of::<CanonicalKey>())?;
    classes[MemoryClass::CountSlots as usize] = table_bytes;
    classes[MemoryClass::InputBuffer as usize] = u64::try_from(record_bytes)?
        .checked_mul(4)
        .context("FASTQ buffer budget overflow")?;
    let owners = u64::try_from(workers.max(1))?;
    let future_workers = workers.max(usize::from(owner_transfer_threshold.is_some()));
    let pool_bytes = if future_workers == 0 {
        0
    } else {
        AdaptiveParallelCountBuilder::reservation_bytes(future_workers)?.min(64 << 20)
    };
    let identity_bytes = if identity_worker {
        IdentityWriter::reservation_bytes()?
    } else {
        0
    };
    classes[MemoryClass::WorkerQueues as usize] = ExternalCountBuilder::run_slot_bytes()?
        .checked_mul(owners)
        .and_then(|bytes| bytes.checked_add(key_bytes))
        .and_then(|bytes| bytes.checked_add(pool_bytes))
        .and_then(|bytes| bytes.checked_add(identity_bytes))
        .context("Count worker reservation overflow")?;
    classes[MemoryClass::SortBuffer as usize] = sort_bytes;
    classes[MemoryClass::LookupCache as usize] = ExternalCountBuilder::cache_limit_bytes()
        .checked_mul(owners)
        .context("Count page reservation overflow")?;
    classes[MemoryClass::ReplayBuffer as usize] = replay_bytes;
    total(classes, disk_bytes, graph_buffer_bytes)
}
