use super::*;
use crate::kmer::budget::{ResourceKind, ResourceLimits};
use std::collections::BTreeSet;
use std::io::{Seek, SeekFrom, Write};

fn test_budget(disk: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::LookupCache as usize] = 4 * PAGE_BYTES as u64;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: 4 * PAGE_BYTES as u64,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk,
        temporary_open_files: 128,
    })
    .unwrap()
}

fn rows(number: usize) -> Vec<(CanonicalKey, NonZeroU64)> {
    let k = KmerSize::new(8).unwrap();
    let mut keys = BTreeSet::new();
    for bits in 0..65536 {
        keys.insert(CanonicalKey::new(bits, k).unwrap());
        if keys.len() == number {
            break;
        }
    }
    assert_eq!(keys.len(), number);
    keys.into_iter()
        .enumerate()
        .map(|(i, key)| (key, NonZeroU64::new(1 + (i % 7) as u64).unwrap()))
        .collect()
}

fn store(
    directory: &Path,
    budget: &Arc<ResourceBudget>,
    rows: &[(CanonicalKey, NonZeroU64)],
) -> DiskCountStore {
    let mut writer =
        SortedCountWriter::new(directory, KmerSize::new(8).unwrap(), budget, 3).unwrap();
    for &(key, count) in rows {
        writer.push(key, count).unwrap();
    }
    writer.finish([9; 32]).unwrap()
}

#[test]
fn disk_lookup_full_keys_batches_reentrant_visits_and_shared_workers_match() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    let all = rows(1200);
    let kept = all
        .iter()
        .enumerate()
        .filter_map(|(i, value)| (i % 3 != 1).then_some(*value))
        .collect::<Vec<_>>();
    let store = store(directory.path(), &budget, &kept);
    for (i, &(key, count)) in all.iter().enumerate() {
        assert_eq!(store.lookup(key).unwrap(), (i % 3 != 1).then_some(count));
    }
    assert!(
        store
            .lookup(CanonicalKey::new(0, KmerSize::new(7).unwrap()).unwrap())
            .is_err()
    );
    let mut batch = [None; 3];
    store
        .lookup_batch(&[kept[10].0, all[1].0, kept[10].0], &mut batch)
        .unwrap();
    assert_eq!(batch, [Some(kept[10].1), None, Some(kept[10].1)]);
    let before = batch;
    assert!(store.lookup_batch(&[kept[0].0; 4], &mut [None; 4]).is_err());
    assert!(store.lookup_batch(&[kept[0].0], &mut []).is_err());
    assert!(
        store
            .lookup_batch(
                &[CanonicalKey::new(0, KmerSize::new(7).unwrap()).unwrap(); 3],
                &mut batch
            )
            .is_err()
    );
    assert_eq!(batch, before);
    let mut observed = Vec::new();
    store
        .visit(&mut |key, count| {
            assert_eq!(store.lookup(key)?, Some(count));
            observed.push((key, count));
            Ok(())
        })
        .unwrap();
    assert_eq!(observed, kept);
    let mut called = 0;
    assert!(
        store
            .visit(&mut |_, _| {
                called += 1;
                anyhow::bail!("visitor failure")
            })
            .is_err()
    );
    assert_eq!(called, 1);
    std::thread::scope(|scope| {
        for offset in 0..4 {
            let store = &store;
            let kept = &kept;
            scope.spawn(move || {
                for &(key, count) in kept.iter().skip(offset).step_by(17) {
                    assert_eq!(store.lookup(key).unwrap(), Some(count));
                }
            });
        }
    });
    assert_eq!(
        store.metadata().occurrences,
        kept.iter().map(|(_, n)| n.get()).sum::<u64>()
    );
    assert_eq!(store.metadata().selected_input_sha256, [9; 32]);
    assert_eq!(store.cache_bytes(), PAGE_BYTES);
    assert_eq!(
        budget.snapshot().peak_reserved_memory_bytes,
        2 * PAGE_BYTES as u64
    );
    assert!(store.file_bytes() > 3 * PAGE_BYTES as u64);
    assert_eq!(
        budget.snapshot().reserved_temporary_disk_bytes,
        store.file_bytes()
    );
    drop(store);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
}

#[test]
fn empty_files_and_exact_u64_counts_preserve_identity_and_overflow_prefix() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    for k in [1, 31] {
        let k = KmerSize::new(k).unwrap();
        let empty = SortedCountWriter::new(directory.path(), k, &budget, 2)
            .unwrap()
            .finish([1; 32])
            .unwrap();
        assert_eq!(empty.metadata().distinct_keys, 0);
        assert_eq!(empty.metadata().occurrences, 0);
        assert_eq!(
            empty.lookup(CanonicalKey::new(0, k).unwrap()).unwrap(),
            None
        );
        assert_eq!(empty.file_bytes(), (HEADER_BYTES + FOOTER_BYTES) as u64);
        drop(empty);
    }
    let k = KmerSize::new(8).unwrap();
    let mut writer = SortedCountWriter::new(directory.path(), k, &budget, 2).unwrap();
    let key = CanonicalKey::new(0, k).unwrap();
    writer
        .push(key, NonZeroU64::new(u64::MAX).unwrap())
        .unwrap();
    for bad in [
        key,
        CanonicalKey::new(1, k).unwrap(),
        CanonicalKey::new(0, KmerSize::new(3).unwrap()).unwrap(),
    ] {
        let error = writer.push(bad, NonZeroU64::new(1).unwrap()).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
    }
    let store = writer.finish([3; 32]).unwrap();
    assert_eq!(store.lookup(key).unwrap().unwrap().get(), u64::MAX);
    assert_eq!(store.metadata().occurrences, u64::MAX);
    assert_eq!(store.metadata().distinct_keys, 1);
}

#[test]
fn page_growth_and_footer_limits_fail_without_leaking_live_files_or_memory() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(HEADER_BYTES as u64 + PAGE_BYTES as u64);
    let mut writer =
        SortedCountWriter::new(directory.path(), KmerSize::new(8).unwrap(), &budget, 3).unwrap();
    for &(key, count) in &rows(10) {
        writer.push(key, count).unwrap();
    }
    let error = writer.finish([0; 32]).err().unwrap();
    assert_eq!(
        error.downcast_ref::<ResourceExhausted>().unwrap().resource,
        ResourceKind::TemporaryDisk
    );
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
    let budget = test_budget(HEADER_BYTES as u64);
    let mut writer =
        SortedCountWriter::new(directory.path(), KmerSize::new(8).unwrap(), &budget, 3).unwrap();
    let items = rows(PAGE_ENTRIES + 1);
    for &(key, count) in &items[..PAGE_ENTRIES] {
        writer.push(key, count).unwrap();
    }
    assert!(
        writer
            .push(items[PAGE_ENTRIES].0, items[PAGE_ENTRIES].1)
            .is_err()
    );
    assert_eq!(writer.keys, PAGE_ENTRIES as u64);
    assert_eq!(writer.file.len().unwrap(), HEADER_BYTES as u64);
    drop(writer);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
}

fn reopen_after(store: DiskCountStore, mutate: impl FnOnce(&Path)) -> Result<DiskCountStore> {
    let state = store.state.into_inner().unwrap();
    mutate(state.file.path()?);
    DiskCountStore::open_owned(state.file, state.page, 3)
}

fn overwrite(path: &Path, offset: u64, bytes: &[u8]) {
    let mut file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
    file.seek(SeekFrom::Start(offset)).unwrap();
    file.write_all(bytes).unwrap();
    file.sync_all().unwrap();
}

#[test]
fn malformed_or_incomplete_streams_are_rejected_before_exposing_a_store() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    for offset in [
        0,
        13,
        HEADER_BYTES + 8,
        HEADER_BYTES + PAYLOAD_BYTES,
        HEADER_BYTES + PAGE_BYTES,
        HEADER_BYTES + PAGE_BYTES + 32,
    ] {
        let store = store(directory.path(), &budget, &rows(10));
        let error = reopen_after(store, |path| {
            let mut data = std::fs::read(path).unwrap();
            data[offset] ^= 1;
            overwrite(path, 0, &data);
        })
        .err()
        .unwrap();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
    }
    for extra in [false, true] {
        let store = store(directory.path(), &budget, &rows(10));
        assert!(
            reopen_after(store, |path| {
                let file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
                let n = file.metadata().unwrap().len();
                file.set_len(if extra { n + 1 } else { n - 1 }).unwrap();
            })
            .is_err()
        );
    }
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
}

#[test]
fn valid_page_checksum_does_not_waive_canonical_count_order_or_padding_checks() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    for (case, expected) in [(0, "zero"), (1, "sorted"), (2, "padding"), (3, "exceeds k")] {
        let store = store(directory.path(), &budget, &rows(10));
        let error = reopen_after(store, |path| {
            let mut bytes = std::fs::read(path).unwrap();
            let payload = &mut bytes[HEADER_BYTES..HEADER_BYTES + PAYLOAD_BYTES];
            match case {
                0 => payload[8..16].fill(0),
                1 => {
                    let first: Vec<u8> = payload[..16].to_vec();
                    payload.copy_within(16..32, 0);
                    payload[16..32].copy_from_slice(&first);
                }
                2 => payload[10 * 16] = 1,
                3 => payload[..8].copy_from_slice(&u64::MAX.to_le_bytes()),
                _ => unreachable!(),
            }
            let mut header = [0; HEADER_BYTES];
            header.copy_from_slice(&bytes[..HEADER_BYTES]);
            let checksum = page_checksum(
                &header,
                0,
                &bytes[HEADER_BYTES..HEADER_BYTES + PAYLOAD_BYTES],
            );
            bytes[HEADER_BYTES + PAYLOAD_BYTES..HEADER_BYTES + PAGE_BYTES]
                .copy_from_slice(&checksum);
            overwrite(path, 0, &bytes);
        })
        .err()
        .unwrap();
        assert!(format!("{error:#}").contains(expected), "{case}: {error:#}");
    }
}

#[test]
fn foreign_pages_and_post_open_corruption_are_errors_never_false_absence() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    let items = rows(600);
    let first = store(directory.path(), &budget, &items);
    let second = store(directory.path(), &budget, &items);
    let foreign = {
        let state = second.state.lock().unwrap();
        std::fs::read(state.file.path().unwrap()).unwrap()
    };
    {
        let mut state = first.state.lock().unwrap();
        state.cached_page = None;
        overwrite(
            state.file.path().unwrap(),
            HEADER_BYTES as u64,
            &foreign[HEADER_BYTES..HEADER_BYTES + PAGE_BYTES],
        );
    }
    let error = first.lookup(items[0].0).unwrap_err();
    assert!(error.to_string().contains("checksum"));
    assert!(first.lookup(items[599].0).is_err());
    let mut called = false;
    assert!(
        first
            .visit(&mut |_, _| {
                called = true;
                Ok(())
            })
            .is_err()
    );
    assert!(!called);
}

#[test]
fn concurrent_scans_read_one_page_per_scan_and_charge_each_live_buffer() {
    let directory = tempfile::tempdir().unwrap();
    let budget = test_budget(1 << 20);
    let items = rows(4 * PAGE_ENTRIES + 13);
    let store = store(directory.path(), &budget, &items);
    let pages = (items.len() as u64).div_ceil(PAGE_ENTRIES as u64);
    let before = store.state.lock().unwrap().page_reads;
    let barrier = std::sync::Barrier::new(2);
    std::thread::scope(|scope| {
        for slow in [false, true] {
            let store = &store;
            let items = &items;
            let barrier = &barrier;
            scope.spawn(move || {
                let mut observed = Vec::new();
                store
                    .visit(&mut |key, count| {
                        if observed.is_empty() {
                            barrier.wait();
                        }
                        if slow && key.bits() % 7 == 0 {
                            std::thread::yield_now();
                        }
                        observed.push((key, count));
                        Ok(())
                    })
                    .unwrap();
                assert_eq!(&observed, items);
            });
        }
    });
    assert_eq!(store.state.lock().unwrap().page_reads - before, 2 * pages);
    assert_eq!(store.peak_page_bytes(), 3 * PAGE_BYTES as u64);
    assert_eq!(budget.snapshot().reserved_memory_bytes, PAGE_BYTES as u64);
}
