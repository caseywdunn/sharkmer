use super::*;
use crate::kmer::access::FilteredKmerCounts;
use crate::kmer::budget::{ResourceExhausted, ResourceLimits};
use std::collections::BTreeMap;

fn budget(bytes: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::CountSlots as usize] = bytes;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: bytes,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 0,
        temporary_open_files: 128,
    })
    .unwrap()
}

#[test]
fn packed_columns_match_exact_counts_across_all_key_widths_and_moving_finalization() {
    for size in 1..=31 {
        let k = KmerSize::new(size).unwrap();
        let pool = budget(8 << 20);
        let mut builder = PackedCountBuilder::with_budget(k, Arc::clone(&pool), 16).unwrap();
        let mut expected = BTreeMap::<u64, u64>::new();
        let mask = (1u64 << (2 * size)) - 1;
        let mut state = 19u64;
        for index in 0..10037 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let key = CanonicalKey::new(state & mask, k).unwrap();
            let count = if index % 97 == 0 {
                u64::from(u32::MAX) + 7
            } else {
                1
            };
            builder.add(key, NonZeroU64::new(count).unwrap()).unwrap();
            *expected.entry(key.bits()).or_default() += count;
        }
        let key_pointer = builder.columns.keys.as_ptr();
        let count_pointer = builder.columns.counts.as_ptr();
        let store = builder.finalize([29; 32]);
        assert_eq!(store.builder.columns.keys.as_ptr(), key_pointer);
        assert_eq!(store.builder.columns.counts.as_ptr(), count_pointer);
        assert_eq!(store.metadata().selected_input_sha256, [29; 32]);
        assert_eq!(store.metadata().distinct_keys, expected.len() as u64);
        assert_eq!(store.metadata().occurrences, expected.values().sum());
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
            let reverse = OrientedKey::new(bits, k).unwrap().reverse_complement();
            assert_eq!(
                store.lookup_oriented(reverse).unwrap().unwrap().get(),
                count
            );
        }
        let view = FilteredKmerCounts::from_packed(&store, u64::from(u32::MAX) + 1);
        for (&bits, &count) in expected.iter().take(19) {
            assert_eq!(
                view.get_canonical_count(&bits).unwrap(),
                if count > u64::from(u32::MAX) {
                    count
                } else {
                    0
                }
            );
        }
        let width = (2 * size).div_ceil(8);
        assert_eq!(store.metrics().key_bytes_per_bucket, width);
        assert_eq!(
            store.table_bytes(),
            store.metrics().buckets * (width + 8) + 7
        );
        assert_eq!(
            pool.snapshot().reserved_memory_bytes,
            store.table_bytes() as u64
        );
        assert_eq!(
            pool.snapshot().peak_reserved_memory_bytes,
            store.peak_table_bytes() as u64
        );
        drop(store);
        assert_eq!(pool.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&pool), 1);
    }
}

#[test]
fn colliding_full_keys_and_wrapped_probes_never_alias_absent_keys_or_key_zero() {
    let pool = budget(4096);
    for width in 1..=8 {
        let mut columns = Columns::new(16, width, &pool).unwrap();
        // All keys deliberately share the last probe origin, crossing index 0.
        for bits in 0..14 {
            let bucket = columns.find(bits, 15);
            assert_eq!(columns.counts[bucket], 0);
            columns.put(bucket, bits, u64::from(u32::MAX) + bits + 1);
        }
        for bits in 0..14 {
            let bucket = columns.find(bits, 15);
            assert_eq!(columns.bits(bucket), bits);
            assert_eq!(columns.counts[bucket], u64::from(u32::MAX) + bits + 1);
        }
        let absent = columns.find(255, 15);
        assert_eq!(columns.counts[absent], 0);
        assert!(columns.bits(15) == 0 && columns.counts[15] > 0);
    }
    assert_eq!(pool.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn rejected_resize_preserves_prefix_and_existing_updates_and_can_retry_after_release() {
    let k = KmerSize::new(19).unwrap();
    let old = Columns::allocation_bytes(16, 5).unwrap() as u64;
    let replacement = Columns::allocation_bytes(32, 5).unwrap() as u64;
    let pool = budget(old + replacement);
    let obstruction = pool.reserve(MemoryClass::CountSlots, 1).unwrap();
    let mut builder = PackedCountBuilder::with_budget(k, Arc::clone(&pool), 4).unwrap();
    for bits in 0..14 {
        builder
            .add(
                CanonicalKey::new(bits, k).unwrap(),
                NonZeroU64::new(1).unwrap(),
            )
            .unwrap();
    }
    let peak = pool.snapshot().peak_reserved_memory_bytes;
    let next = CanonicalKey::new(14, k).unwrap();
    assert!(
        builder
            .add(next, NonZeroU64::new(1).unwrap())
            .unwrap_err()
            .downcast_ref::<ResourceExhausted>()
            .is_some()
    );
    assert_eq!(builder.occurrences(), 14);
    assert_eq!(pool.snapshot().reserved_memory_bytes, old + 1);
    assert_eq!(pool.snapshot().peak_reserved_memory_bytes, peak);
    builder
        .add(
            CanonicalKey::new(0, k).unwrap(),
            NonZeroU64::new(2).unwrap(),
        )
        .unwrap();
    drop(obstruction);
    builder.add(next, NonZeroU64::new(1).unwrap()).unwrap();
    let store = builder.finalize([0; 32]);
    assert_eq!(store.metadata().occurrences, 17);
    assert_eq!(
        store
            .lookup(CanonicalKey::new(0, k).unwrap())
            .unwrap()
            .unwrap()
            .get(),
        3
    );
    for bits in 1..15 {
        assert_eq!(
            store
                .lookup(CanonicalKey::new(bits, k).unwrap())
                .unwrap()
                .unwrap()
                .get(),
            1
        );
    }
    assert_eq!(store.peak_table_bytes() as u64, old + replacement);
}

#[test]
fn overflow_wrong_k_and_batch_errors_leave_accepted_counts_and_outputs_intact() {
    let k = KmerSize::new(31).unwrap();
    let pool = budget(4096);
    let mut builder = PackedCountBuilder::with_budget(k, Arc::clone(&pool), 2).unwrap();
    let key = CanonicalKey::new(0, k).unwrap();
    builder
        .add(key, NonZeroU64::new(u64::MAX).unwrap())
        .unwrap();
    assert!(builder.add(key, NonZeroU64::new(1).unwrap()).is_err());
    assert!(
        builder
            .add(
                CanonicalKey::new(1, KmerSize::new(30).unwrap()).unwrap(),
                NonZeroU64::new(1).unwrap()
            )
            .is_err()
    );
    let store = builder.finalize([0; 32]);
    assert_eq!(store.metadata().occurrences, u64::MAX);
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), u64::MAX);
    assert_eq!(
        store.lookup(CanonicalKey::new(1, k).unwrap()).unwrap(),
        None
    );
    let wrong = CanonicalKey::new(0, KmerSize::new(30).unwrap()).unwrap();
    assert!(store.lookup(wrong).is_err());
    let mut output = [NonZeroU64::new(7); 3];
    assert!(store.lookup_batch(&[key, key, key], &mut output).is_err());
    assert_eq!(output, [NonZeroU64::new(7); 3]);
    assert!(store.lookup_batch(&[key, wrong], &mut output[..2]).is_err());
    assert_eq!(output, [NonZeroU64::new(7); 3]);
    store.lookup_batch(&[key, key], &mut output[..2]).unwrap();
    assert_eq!(output[..2], [NonZeroU64::new(u64::MAX); 2]);
    assert!(
        store
            .visit(&mut |_, _| anyhow::bail!("visitor failed"))
            .is_err()
    );
    assert!(Columns::allocation_bytes(usize::MAX, 8).is_err());
}
