use super::*;
use crate::kmer::budget::ResourceLimits;
use std::collections::BTreeMap;

fn resources(table: u64, sort: u64, disk: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::CountSlots as usize] = table;
    classes[MemoryClass::SortBuffer as usize] = sort;
    classes[MemoryClass::WorkerQueues as usize] = ExternalCountBuilder::run_slot_bytes().unwrap();
    classes[MemoryClass::LookupCache as usize] = ExternalCountBuilder::cache_limit_bytes();
    ResourceBudget::new(ResourceLimits {
        memory_bytes: classes.iter().sum(),
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk,
        temporary_open_files: 128,
    })
    .unwrap()
}

#[test]
fn fitting_input_retains_memory_without_allocating_disk_or_sort_buffers() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(600, 64, 0);
    let k = KmerSize::new(7).unwrap();
    let mut builder =
        AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 64, 8).unwrap();
    for index in 0..100 {
        builder
            .add(
                CanonicalKey::new(index % 5, k).unwrap(),
                NonZeroU64::new(3).unwrap(),
            )
            .unwrap();
    }
    let AdaptiveStore::Memory(store) = builder.finish([9; 32]).unwrap() else {
        panic!("unnecessary spill")
    };
    assert_eq!(store.metadata().occurrences, 300);
    assert_eq!(store.metadata().selected_input_sha256, [9; 32]);
    assert_eq!(
        budget.snapshot().peak_reserved_by_class[MemoryClass::SortBuffer as usize],
        0
    );
    assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 0);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    drop(store);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn weighted_prefix_and_unaccepted_update_transfer_once_and_release_table() {
    for width in [3, 8, 19, 31] {
        let directory = tempfile::tempdir().unwrap();
        let budget = resources(600, 64, 1 << 24);
        let k = KmerSize::new(width).unwrap();
        let mut builder =
            AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 64, 8).unwrap();
        let mut expected = BTreeMap::<u64, u64>::new();
        for index in 0..180u64 {
            let key = CanonicalKey::new(index % (1u64 << (2 * width).min(8)), k).unwrap();
            let count = if index == 0 {
                u64::from(u32::MAX) + 21
            } else {
                1 + index % 13
            };
            builder.add(key, NonZeroU64::new(count).unwrap()).unwrap();
            *expected.entry(key.bits()).or_default() += count;
            assert_eq!(builder.occurrences(), expected.values().sum());
        }
        assert!(matches!(builder.state, Some(State::External(_))));
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::CountSlots as usize],
            0
        );
        let AdaptiveStore::Disk(store, metrics) = builder.finish([17; 32]).unwrap() else {
            panic!("missing spill")
        };
        assert!(metrics.peak_table_bytes > 0);
        assert!(metrics.peak_table_bytes <= 600);
        assert_eq!(store.metadata().selected_input_sha256, [17; 32]);
        let mut observed = BTreeMap::new();
        store
            .visit(&mut |key, count| {
                assert!(observed.insert(key.bits(), count.get()).is_none());
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected);
        assert_eq!(store.metadata().occurrences, expected.values().sum());
        store.close().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    }
}

#[test]
fn wrong_k_and_overflow_do_not_trigger_spill_or_lose_valid_counts() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(600, 64, 0);
    let k = KmerSize::new(8).unwrap();
    let mut builder =
        AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 64, 8).unwrap();
    let key = CanonicalKey::new(0, k).unwrap();
    builder
        .add(key, NonZeroU64::new(u64::MAX).unwrap())
        .unwrap();
    assert!(builder.add(key, NonZeroU64::new(1).unwrap()).is_err());
    assert!(
        builder
            .add(
                CanonicalKey::new(0, KmerSize::new(3).unwrap()).unwrap(),
                NonZeroU64::new(1).unwrap()
            )
            .is_err()
    );
    assert_eq!(builder.occurrences(), u64::MAX);
    assert!(matches!(builder.state, Some(State::Memory(_))));
    let AdaptiveStore::Memory(store) = builder.finish([1; 32]).unwrap() else {
        panic!("invalid spill")
    };
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), u64::MAX);
    assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 0);
}

#[test]
fn failed_transfer_poisoning_prevents_partial_publication_and_cleans_files() {
    for disk in [0, 1 << 20] {
        let directory = tempfile::tempdir().unwrap();
        let budget = resources(600, 64, disk);
        let k = KmerSize::new(8).unwrap();
        let mut builder =
            AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 64, 8).unwrap();
        // Also exercise physical filesystem failure, distinct from disk limits.
        if disk != 0 {
            std::fs::remove_dir(directory.path()).unwrap();
        }
        let mut failed = false;
        for bits in 0..100 {
            match builder.add(
                CanonicalKey::new(bits, k).unwrap(),
                NonZeroU64::new(1).unwrap(),
            ) {
                Ok(()) => (),
                Err(error) => {
                    if disk == 0 {
                        assert_eq!(
                            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
                            ResourceKind::TemporaryDisk
                        );
                    } else {
                        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
                    }
                    failed = true;
                    break;
                }
            }
        }
        assert!(failed);
        assert!(builder.state.is_none());
        assert!(
            builder
                .add(
                    CanonicalKey::new(0, k).unwrap(),
                    NonZeroU64::new(1).unwrap()
                )
                .is_err()
        );
        assert!(builder.finish([0; 32]).is_err());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        if directory.path().exists() {
            assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn combined_memory_pressure_spills_while_other_stage_ownership_stays_live() {
    let directory = tempfile::tempdir().unwrap();
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::CountSlots as usize] = 500_000;
    classes[MemoryClass::SortBuffer as usize] = 4096;
    classes[MemoryClass::WorkerQueues as usize] = ExternalCountBuilder::run_slot_bytes().unwrap();
    classes[MemoryClass::LookupCache as usize] = ExternalCountBuilder::cache_limit_bytes();
    classes[MemoryClass::ReadEvidence as usize] = 400_000;
    let budget = ResourceBudget::new(ResourceLimits {
        memory_bytes: 500_000,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 1 << 24,
        temporary_open_files: 128,
    })
    .unwrap();
    let evidence = budget.reserve(MemoryClass::ReadEvidence, 400_000).unwrap();
    let k = KmerSize::new(19).unwrap();
    let mut builder =
        AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 4096, 8).unwrap();
    for bits in 0..2200 {
        // The 1,793rd unique key needs simultaneous old/new table allocations
        // of 104,480 bytes: below CountSlots, above the remaining run budget.
        if bits == 1792 {
            assert!(matches!(builder.state, Some(State::Memory(_))));
        }
        builder
            .add(
                CanonicalKey::new(bits, k).unwrap(),
                NonZeroU64::new(1).unwrap(),
            )
            .unwrap();
        if bits == 1792 {
            assert!(matches!(builder.state, Some(State::External(_))));
        }
    }
    assert_eq!(
        budget.snapshot().reserved_by_class[MemoryClass::CountSlots as usize],
        0
    );
    assert_eq!(
        budget.snapshot().reserved_by_class[MemoryClass::ReadEvidence as usize],
        400_000
    );
    assert!(budget.snapshot().peak_reserved_memory_bytes <= 500_000);
    let AdaptiveStore::Disk(store, _) = builder.finish([5; 32]).unwrap() else {
        panic!("missing spill")
    };
    assert_eq!(store.metadata().distinct_keys, 2200);
    assert_eq!(store.metadata().occurrences, 2200);
    for bits in 0..2200 {
        assert_eq!(
            store
                .lookup(CanonicalKey::new(bits, k).unwrap())
                .unwrap()
                .unwrap()
                .get(),
            1
        );
    }
    store.close().unwrap();
    assert_eq!(budget.snapshot().reserved_memory_bytes, 400_000);
    drop(evidence);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
}
