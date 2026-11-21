use super::*;
use crate::kmer::budget::{ResourceExhausted, ResourceLimits};
use std::collections::BTreeMap;

fn resources(sort: usize, disk: u64) -> Arc<ResourceBudget> {
    resources_with_files(sort, disk, 128)
}

fn resources_with_files(sort: usize, disk: u64, files: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::SortBuffer as usize] = sort as u64;
    classes[MemoryClass::WorkerQueues as usize] = ExternalCountBuilder::run_slot_bytes().unwrap();
    classes[MemoryClass::LookupCache as usize] = ExternalCountBuilder::cache_limit_bytes();
    ResourceBudget::new(ResourceLimits {
        memory_bytes: classes.iter().sum(),
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk,
        temporary_open_files: files,
    })
    .unwrap()
}

#[test]
fn multi_run_reduction_matches_weighted_oracle_without_growing_the_sort_buffer() {
    let k = KmerSize::new(8).unwrap();
    for capacity in [1, 7, 31] {
        let directory = tempfile::tempdir().unwrap();
        let budget = resources(capacity * 16, 1 << 24);
        let mut builder =
            ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), capacity * 16, 8)
                .unwrap();
        let pointer = builder.pairs.as_ptr();
        let mut expected = BTreeMap::<u64, u64>::new();
        for index in 0..600u64 {
            let key = CanonicalKey::new(((index % 300) * 557) % 65536, k).unwrap();
            let count = if index == 0 {
                u64::from(u32::MAX) + 17
            } else {
                1 + index % 7
            };
            *expected.entry(key.bits()).or_default() += count;
            builder.add(key, NonZeroU64::new(count).unwrap()).unwrap();
            assert_eq!(builder.pairs.as_ptr(), pointer);
            assert_eq!(builder.pairs.capacity(), capacity);
            assert_eq!(builder.runs.len(), RUN_SLOTS);
            assert!(std::fs::read_dir(directory.path()).unwrap().count() <= RUN_SLOTS);
        }
        assert_eq!(builder.occurrences(), expected.values().sum());
        let (store, metrics) = builder.finish([7; 32]).unwrap();
        let mut observed = BTreeMap::new();
        store
            .visit(&mut |key, count| {
                assert!(observed.insert(key.bits(), count.get()).is_none());
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected);
        for (&bits, &count) in &expected {
            assert_eq!(
                store
                    .lookup(CanonicalKey::new(bits, k).unwrap())
                    .unwrap()
                    .unwrap()
                    .get(),
                count
            );
        }
        assert_eq!(store.metadata().selected_input_sha256, [7; 32]);
        assert_eq!(store.metadata().occurrences, expected.values().sum());
        assert_eq!(metrics.sorted_runs, 600u64.div_ceil(capacity as u64));
        assert_eq!(metrics.merges, metrics.sorted_runs - 1);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::SortBuffer as usize],
            0
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::WorkerQueues as usize],
            0
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, PAGE_BYTES as u64);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
        drop(store);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    }
}

#[test]
fn overflow_and_wrong_k_leave_the_accepted_prefix_available() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(16, 1 << 20);
    let k = KmerSize::new(3).unwrap();
    let mut builder =
        ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 16, 8).unwrap();
    let key = CanonicalKey::new(0, k).unwrap();
    builder
        .add(key, NonZeroU64::new(u64::MAX).unwrap())
        .unwrap();
    for wrong in [
        key,
        CanonicalKey::new(0, KmerSize::new(5).unwrap()).unwrap(),
    ] {
        let error = builder.add(wrong, NonZeroU64::new(1).unwrap()).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
        assert_eq!(builder.occurrences(), u64::MAX);
        assert!(!builder.failed);
    }
    let (store, _) = builder.finish([4; 32]).unwrap();
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), u64::MAX);
}

#[test]
fn failed_spill_cannot_publish_an_accepted_subset_and_cleans_owned_files() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(16, 4257);
    let k = KmerSize::new(5).unwrap();
    let mut builder =
        ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 16, 8).unwrap();
    for bits in [0, 1] {
        builder
            .add(
                CanonicalKey::new(bits, k).unwrap(),
                NonZeroU64::new(1).unwrap(),
            )
            .unwrap();
    }
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
    let error = builder
        .add(
            CanonicalKey::new(2, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap_err();
    assert!(error.downcast_ref::<ResourceExhausted>().is_some());
    assert_eq!(builder.occurrences(), 2);
    assert!(builder.failed);
    assert_eq!(builder.runs[0].as_ref().unwrap().metadata().occurrences, 1);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
    assert!(builder.finish([0; 32]).is_err());
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn empty_streams_are_bound_to_input_and_invalid_buffers_allocate_nothing() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(32, 1 << 20);
    let k = KmerSize::new(31).unwrap();
    assert!(ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 15, 8).is_err());
    assert!(ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 32, 0).is_err());
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    let builder =
        ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 32, 8).unwrap();
    let (store, metrics) = builder.finish([11; 32]).unwrap();
    assert_eq!(store.metadata().distinct_keys, 0);
    assert_eq!(store.metadata().occurrences, 0);
    assert_eq!(store.metadata().selected_input_sha256, [11; 32]);
    assert_eq!(metrics.sorted_runs, 0);
    assert_eq!(metrics.merges, 0);
    drop(store);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
}

#[cfg(unix)]
#[test]
fn failed_merge_input_cleanup_invalidates_publication_and_retains_disk_charge() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(16, 1 << 20);
    let k = KmerSize::new(5).unwrap();
    let mut builder =
        ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 16, 8).unwrap();
    for bits in [0, 1] {
        builder
            .add(
                CanonicalKey::new(bits, k).unwrap(),
                NonZeroU64::new(1).unwrap(),
            )
            .unwrap();
    }
    let original = std::fs::read_dir(directory.path())
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let surviving = directory.path().join("surviving-run");
    std::fs::rename(&original, &surviving).unwrap();
    std::fs::create_dir(&original).unwrap();
    let error = builder
        .add(
            CanonicalKey::new(2, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap_err();
    assert!(error.downcast_ref::<ResourceExhausted>().is_none());
    assert!(format!("{error:#}").contains("remains charged"));
    assert!(builder.failed);
    assert!(builder.finish([0; 32]).is_err());
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(
        budget.snapshot().reserved_temporary_disk_bytes,
        std::fs::metadata(&surviving).unwrap().len()
    );
    assert!(budget.snapshot().reserved_temporary_disk_bytes > 0);
    std::fs::remove_dir(original).unwrap();
    std::fs::remove_file(surviving).unwrap();
}

#[test]
fn merge_descriptor_exhaustion_is_typed_and_cannot_publish_a_partial_prefix() {
    use crate::kmer::budget::TemporaryFilesExhausted;

    let directory = tempfile::tempdir().unwrap();
    let k = KmerSize::new(8).unwrap();
    let budget = resources_with_files(16, 1 << 20, 2);
    let mut builder =
        ExternalCountBuilder::new(directory.path(), k, Arc::clone(&budget), 16, 8).unwrap();
    builder
        .add(
            CanonicalKey::new(1, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap();
    builder
        .add(
            CanonicalKey::new(2, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap();
    let error = builder
        .add(
            CanonicalKey::new(3, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap_err();
    assert_eq!(
        error.downcast_ref::<TemporaryFilesExhausted>(),
        Some(&TemporaryFilesExhausted {
            required_files: 3,
            limit_files: 2,
        })
    );
    assert_eq!(builder.occurrences(), 2);
    assert!(builder.finish([0; 32]).is_err());
    assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
    assert_eq!(budget.snapshot().peak_reserved_temporary_open_files, 2);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    assert!(budget.snapshot().valid);
}
