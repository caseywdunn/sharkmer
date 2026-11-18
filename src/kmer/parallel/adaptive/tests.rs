use super::*;
use crate::kmer::budget::{ResourceExhausted, ResourceKind, ResourceLimits};
use crate::kmer::external::ExternalCountBuilder;
use std::collections::BTreeMap;

fn resources(workers: usize, table: u64, sort: u64, disk: u64) -> Arc<ResourceBudget> {
    resources_with_files(workers, table, sort, disk, 128)
}

fn resources_with_files(
    workers: usize,
    table: u64,
    sort: u64,
    disk: u64,
    files: u64,
) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::CountSlots as usize] = table;
    classes[MemoryClass::SortBuffer as usize] = sort;
    classes[MemoryClass::WorkerQueues as usize] =
        AdaptiveParallelCountBuilder::reservation_bytes(workers).unwrap()
            + ExternalCountBuilder::run_slot_bytes().unwrap() * workers as u64;
    classes[MemoryClass::LookupCache as usize] =
        ExternalCountBuilder::cache_limit_bytes() * workers as u64;
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
fn skewed_owner_spills_without_losing_queued_tails_or_duplicating_other_shards() {
    let k = KmerSize::new(19).unwrap();
    for workers in [1, 3, 7] {
        let directory = tempfile::tempdir().unwrap();
        let budget = resources(workers, 4096, 4096, 1 << 25);
        let mut builder = AdaptiveParallelCountBuilder::adaptive(
            directory.path(),
            k,
            Arc::clone(&budget),
            workers,
            4096,
            3,
        )
        .unwrap();
        let mut expected = BTreeMap::new();
        // All but one owner retain a tiny table. Only owner zero needs to spill.
        for selected in 1..workers {
            let bits = (0..100000)
                .find(|&bits| owner(bits, workers) == selected)
                .unwrap();
            builder.add(CanonicalKey::new(bits, k).unwrap()).unwrap();
            expected.insert(bits, 1u64);
        }
        let mut unique = 0;
        for bits in 0..100000 {
            if owner(bits, workers) != 0 {
                continue;
            }
            let key = CanonicalKey::new(bits, k).unwrap();
            for _ in 0..(1 + bits % 3) {
                builder.add(key).unwrap();
                *expected.entry(bits).or_default() += 1;
            }
            unique += 1;
            if unique == 2107 {
                break;
            }
        }
        assert_eq!(builder.occurrences(), expected.values().sum());
        let (store, metrics) = builder.finish_adaptive([7; 32]).unwrap();
        assert_eq!(metrics.counting_workers, workers);
        assert_eq!(store.metrics().disk_shards, 1);
        assert_eq!(store.metrics().memory_shards, workers - 1);
        assert_eq!(store.metadata().selected_input_sha256, [7; 32]);
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
        let absent = CanonicalKey::new(200000, k).unwrap();
        assert_eq!(store.lookup(absent).unwrap(), None);
        let key = CanonicalKey::new(*expected.first_key_value().unwrap().0, k).unwrap();
        let mut output = [None; 3];
        store
            .lookup_batch(&[key, absent, key], &mut output)
            .unwrap();
        assert_eq!(output[0], output[2]);
        assert!(output[0].is_some() && output[1].is_none());
        assert!(
            store
                .lookup(CanonicalKey::new(0, KmerSize::new(3).unwrap()).unwrap())
                .is_err()
        );
        assert!(
            store
                .visit(&mut |_, _| anyhow::bail!("visitor stopped"))
                .is_err()
        );
        assert!(budget.snapshot().peak_reserved_by_class[MemoryClass::SortBuffer as usize] <= 4096);
        assert!(budget.snapshot().peak_reserved_by_class[MemoryClass::CountSlots as usize] <= 4096);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::CountSlots as usize],
            store.table_bytes() as u64
        );
        store.close().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        assert_eq!(Arc::strong_count(&budget), 1);
    }
}

#[test]
fn fitting_parallel_input_allocates_neither_sort_buffers_nor_disk_files() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(3, 1 << 20, 4096, 0);
    let k = KmerSize::new(19).unwrap();
    let mut builder = AdaptiveParallelCountBuilder::adaptive(
        directory.path(),
        k,
        Arc::clone(&budget),
        3,
        4096,
        8,
    )
    .unwrap();
    for bits in 0..33 {
        builder.add(CanonicalKey::new(bits, k).unwrap()).unwrap();
    }
    let (store, _) = builder.finish_adaptive([1; 32]).unwrap();
    assert_eq!(store.metrics().disk_shards, 0);
    assert_eq!(store.metrics().memory_shards, 3);
    assert_eq!(store.metadata().occurrences, 33);
    assert_eq!(
        budget.snapshot().peak_reserved_by_class[MemoryClass::SortBuffer as usize],
        0
    );
    assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 0);
    store.close().unwrap();
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn failed_spill_or_drop_joins_all_owners_and_returns_owned_resources() {
    for disk in [0, 1 << 24] {
        let directory = tempfile::tempdir().unwrap();
        let budget = resources(3, 4096, 4096, disk);
        let k = KmerSize::new(19).unwrap();
        let mut builder = AdaptiveParallelCountBuilder::adaptive(
            directory.path(),
            k,
            Arc::clone(&budget),
            3,
            4096,
            8,
        )
        .unwrap();
        let mut failure = None;
        for bits in 0..5000 {
            if let Err(error) = builder.add(CanonicalKey::new(bits, k).unwrap()) {
                failure = Some(error);
                break;
            }
        }
        if disk == 0 {
            let error = match failure {
                Some(error) => {
                    drop(builder);
                    error
                }
                None => builder.finish_adaptive([0; 32]).err().unwrap(),
            };
            assert_eq!(
                error.downcast_ref::<ResourceExhausted>().unwrap().resource,
                ResourceKind::TemporaryDisk
            );
        } else {
            assert!(failure.is_none());
            drop(builder);
        }
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        assert_eq!(Arc::strong_count(&budget), 1);
    }
}

#[test]
fn partial_constructor_failure_joins_started_adaptive_owners() {
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(3, 288, 4096, 1 << 20);
    assert!(
        AdaptiveParallelCountBuilder::adaptive(
            directory.path(),
            KmerSize::new(19).unwrap(),
            Arc::clone(&budget),
            3,
            4096,
            8
        )
        .is_err()
    );
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(Arc::strong_count(&budget), 1);
}

#[test]
fn descriptor_pressure_in_one_owner_joins_every_owner_and_releases_all_files() {
    use crate::kmer::budget::TemporaryFilesExhausted;

    let directory = tempfile::tempdir().unwrap();
    let k = KmerSize::new(19).unwrap();
    let budget = resources_with_files(3, 4096, 48, 1 << 22, 2);
    let mut builder =
        AdaptiveParallelCountBuilder::adaptive(directory.path(), k, Arc::clone(&budget), 3, 48, 3)
            .unwrap();
    let mut failure = None;
    for bits in 0..20000 {
        if let Err(error) = builder.add(CanonicalKey::new(bits, k).unwrap()) {
            failure = Some(error);
            break;
        }
    }
    let error = if let Some(error) = failure {
        drop(builder);
        error
    } else {
        builder.finish_adaptive([0; 32]).err().unwrap()
    };
    assert_eq!(
        error.downcast_ref::<TemporaryFilesExhausted>(),
        Some(&TemporaryFilesExhausted {
            required_files: 3,
            limit_files: 2,
        })
    );
    let snapshot = budget.snapshot();
    assert_eq!(snapshot.reserved_temporary_open_files, 0);
    assert_eq!(snapshot.peak_reserved_temporary_open_files, 2);
    assert_eq!(snapshot.reserved_temporary_disk_bytes, 0);
    assert_eq!(snapshot.reserved_memory_bytes, 0);
    assert!(snapshot.valid);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
}

#[test]
fn moved_prefix_retains_empty_weighted_and_external_counts_in_its_original_ledger() {
    let k = KmerSize::new(19).unwrap();
    for spill in [false, true] {
        for weighted in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let budget = resources(1, if spill { 4096 } else { 1 << 20 }, 4096, 1 << 24);
            let mut prefix =
                AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 4096, 3)
                    .unwrap();
            let mut expected = BTreeMap::new();
            if weighted {
                for bits in 0..1000 {
                    let count = if bits == 0 {
                        u64::from(u32::MAX) + 17
                    } else {
                        2 + bits % 3
                    };
                    prefix
                        .add(
                            CanonicalKey::new(bits, k).unwrap(),
                            NonZeroU64::new(count).unwrap(),
                        )
                        .unwrap();
                    expected.insert(bits, count);
                }
            }
            let total = expected.values().sum::<u64>();
            let mut owner = AdaptiveParallelCountBuilder::from_prefix(prefix).unwrap();
            assert_eq!(owner.occurrences(), total);
            assert!(
                owner
                    .add(CanonicalKey::new(1, KmerSize::new(3).unwrap()).unwrap())
                    .is_err()
            );
            for bits in 0..1307 {
                owner.add(CanonicalKey::new(bits, k).unwrap()).unwrap();
                *expected.entry(bits).or_default() += 1;
            }
            let (store, metrics) = owner.finish_adaptive([9; 32]).unwrap();
            assert_eq!(metrics.counting_workers, 1);
            assert_eq!(store.metadata().selected_input_sha256, [9; 32]);
            assert_eq!(store.metadata().occurrences, expected.values().sum());
            let mut actual = BTreeMap::new();
            store
                .visit(&mut |key, count| {
                    actual.insert(key.bits(), count.get());
                    Ok(())
                })
                .unwrap();
            assert_eq!(actual, expected);
            store.close().unwrap();
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
            assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
            assert!(budget.snapshot().valid);
            assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn owner_prefix_overflow_and_startup_reservation_failure_preserve_error_semantics() {
    let k = KmerSize::new(19).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let budget = resources(1, 1 << 20, 4096, 1 << 24);
    let mut prefix =
        AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 4096, 3).unwrap();
    let key = CanonicalKey::new(1, k).unwrap();
    prefix.add(key, NonZeroU64::new(u64::MAX).unwrap()).unwrap();
    let mut owner = AdaptiveParallelCountBuilder::from_prefix(prefix).unwrap();
    assert!(owner.add(key).is_err());
    assert_eq!(owner.occurrences(), u64::MAX);
    let (store, _) = owner.finish_adaptive([0; 32]).unwrap();
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), u64::MAX);
    store.close().unwrap();
    let prefix =
        AdaptiveCountBuilder::new(directory.path(), k, Arc::clone(&budget), 4096, 3).unwrap();
    let blocked = budget
        .reserve(
            MemoryClass::WorkerQueues,
            budget.class_limit(MemoryClass::WorkerQueues),
        )
        .unwrap();
    let error = AdaptiveParallelCountBuilder::from_prefix(prefix)
        .err()
        .unwrap();
    assert!(error.downcast_ref::<ResourceExhausted>().is_some());
    drop(blocked);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
}
