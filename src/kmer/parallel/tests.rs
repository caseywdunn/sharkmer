use super::*;
use crate::kmer::access::FilteredKmerCounts;
use crate::kmer::budget::{ResourceExhausted, ResourceLimits};
use std::collections::BTreeMap;

fn resources(workers: usize, table_bytes: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::CountSlots as usize] = table_bytes;
    classes[MemoryClass::WorkerQueues as usize] =
        ParallelCountBuilder::reservation_bytes(workers).unwrap();
    ResourceBudget::new(ResourceLimits {
        memory_bytes: classes.iter().sum(),
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 0,
        temporary_open_files: 128,
    })
    .unwrap()
}

#[test]
fn owner_shards_preserve_all_windows_tails_and_absent_keys_without_consolidation() {
    let k = KmerSize::new(8).unwrap();
    for workers in [1, 3, 7] {
        let budget = resources(workers, 8 << 20);
        let mut builder = ParallelCountBuilder::new(k, Arc::clone(&budget), workers, 3).unwrap();
        let mut expected = BTreeMap::<u64, u64>::new();
        for index in 0..10037u64 {
            let key = CanonicalKey::new(index.wrapping_mul(7919) % 65536, k).unwrap();
            *expected.entry(key.bits()).or_default() += 1;
            builder.add(key).unwrap();
            assert!(
                builder
                    .workers
                    .iter()
                    .all(|worker| worker.current.capacity() == BATCH_KEYS)
            );
        }
        assert_eq!(builder.occurrences(), 10037);
        let (store, metrics) = builder.finish([17; 32]).unwrap();
        let mut observed = BTreeMap::new();
        store
            .visit(&mut |key, count| {
                assert!(observed.insert(key.bits(), count.get()).is_none());
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected);
        assert_eq!(store.metadata().occurrences, 10037);
        assert_eq!(store.metadata().selected_input_sha256, [17; 32]);
        assert_eq!(metrics.counting_workers, workers);
        assert_eq!(metrics.reusable_key_buffers, BUFFERS_PER_WORKER * workers);
        assert_eq!(
            metrics.owned_worker_key_bytes,
            (BUFFERS_PER_WORKER * workers * BATCH_KEYS * size_of::<u64>()) as u64
        );
        for (&bits, &count) in &expected {
            let key = CanonicalKey::new(bits, k).unwrap();
            assert_eq!(store.lookup(key).unwrap().unwrap().get(), count);
            assert_eq!(
                store
                    .stores
                    .iter()
                    .filter(|shard| shard.lookup(key).unwrap().is_some())
                    .count(),
                1
            );
        }
        for bits in 0..65536 {
            let key = CanonicalKey::new(bits, k).unwrap();
            if !expected.contains_key(&key.bits()) {
                assert_eq!(store.lookup(key).unwrap(), None);
                break;
            }
        }
        let view = FilteredKmerCounts::from_sharded(&store, 2);
        for (&bits, &count) in expected.iter().take(31) {
            assert_eq!(
                view.get_canonical_count(&bits).unwrap(),
                if count >= 2 { count } else { 0 }
            );
        }
        assert!(
            store
                .lookup(CanonicalKey::new(0, KmerSize::new(7).unwrap()).unwrap())
                .is_err()
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::CountSlots as usize],
            store.table_bytes() as u64
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::WorkerQueues as usize],
            ParallelCountBuilder::descriptor_bytes(workers).unwrap()
        );
        assert!(store.peak_table_bytes() as u64 <= 8 << 20);
        drop(store);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
    }
}

#[test]
fn dropping_a_producer_drains_and_joins_before_returning_its_reservations() {
    let budget = resources(3, 8 << 20);
    let k = KmerSize::new(8).unwrap();
    let mut builder = ParallelCountBuilder::new(k, Arc::clone(&budget), 3, 3).unwrap();
    for index in 0..5007 {
        builder
            .add(CanonicalKey::new(index % 65536, k).unwrap())
            .unwrap();
    }
    drop(builder);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert!(budget.snapshot().valid);
    assert_eq!(Arc::strong_count(&budget), 1);
}

#[test]
fn worker_growth_errors_keep_their_type_and_cannot_publish_partial_counts() {
    let budget = resources(1, 4096);
    let k = KmerSize::new(8).unwrap();
    let mut builder = ParallelCountBuilder::new(k, Arc::clone(&budget), 1, 3).unwrap();
    let mut failure = None;
    for index in 0..10000 {
        if let Err(error) = builder.add(CanonicalKey::new(index % 65536, k).unwrap()) {
            failure = Some(error);
            break;
        }
    }
    let error = match failure {
        Some(error) => {
            drop(builder);
            error
        }
        None => builder.finish([0; 32]).err().unwrap(),
    };
    assert!(
        error.downcast_ref::<ResourceExhausted>().is_some(),
        "{error:#}"
    );
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert!(budget.snapshot().valid);
    assert_eq!(Arc::strong_count(&budget), 1);
}

#[test]
fn invalid_identity_preserves_valid_prefix_and_constructor_failure_joins_started_owners() {
    let k = KmerSize::new(8).unwrap();
    let budget = resources(3, 8 << 20);
    let mut builder = ParallelCountBuilder::new(k, Arc::clone(&budget), 3, 3).unwrap();
    let key = CanonicalKey::new(0, k).unwrap();
    for _ in 0..7 {
        builder.add(key).unwrap();
    }
    assert!(
        builder
            .add(CanonicalKey::new(0, KmerSize::new(7).unwrap()).unwrap())
            .is_err()
    );
    let (store, _) = builder.finish([9; 32]).unwrap();
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), 7);
    drop(store);
    let budget = resources(3, 288);
    assert!(ParallelCountBuilder::new(k, Arc::clone(&budget), 3, 3).is_err());
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(Arc::strong_count(&budget), 1);
}
