use super::*;
use crate::kmer::budget::{ResourceExhausted, ResourceLimits};
use std::io::{Read, Seek, SeekFrom, Write};

fn budget(memory: u64, disk: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::ReplayBuffer as usize] = memory;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: memory,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk,
        temporary_open_files: 128,
    })
    .unwrap()
}

fn record(index: u64, paired: bool) -> ReplayRecord<'static> {
    ReplayRecord {
        index,
        fragment: if paired { index / 2 } else { index },
        mate: if paired {
            if index % 2 == 0 { Mate::R1 } else { Mate::R2 }
        } else {
            Mate::Unpaired
        },
        source: "original source.fastq",
        fields: [b"@read original identifier", b"ACNTAC", b"+read", b"!#I+@J"],
    }
}

fn complete(
    directory: &Path,
    budget: &Arc<ResourceBudget>,
    paired: bool,
    count: u64,
) -> (ReplayStore, [u8; 32]) {
    let mut writer = ReplayWriter::new(directory, budget, 1024, paired, count).unwrap();
    let mut identity = selected_hash();
    for index in 0..count {
        record(index, paired)
            .update_identity(&mut identity)
            .unwrap();
        writer.add(record(index, paired)).unwrap();
    }
    let digest = identity.finalize().into();
    (writer.finish(digest).unwrap(), digest)
}

#[test]
fn replay_preserves_original_fields_positions_mates_and_identity_with_bounded_visits() {
    let directory = tempfile::tempdir().unwrap();
    let resources = budget(4096, 1 << 20);
    for paired in [false, true] {
        let (store, _) = complete(directory.path(), &resources, paired, 101);
        assert_eq!(store.record_count(), 101);
        assert_eq!(store.paired(), paired);
        let mut next = 0;
        store
            .visit(&mut |observed| {
                let expected = record(next, paired);
                assert_eq!(observed.index, expected.index);
                assert_eq!(observed.fragment, expected.fragment);
                assert_eq!(observed.mate, expected.mate);
                assert_eq!(observed.source, expected.source);
                assert_eq!(observed.fields, expected.fields);
                if next == 0 {
                    let mut nested = 0;
                    store.visit(&mut |_| {
                        nested += 1;
                        Ok(())
                    })?;
                    assert_eq!(nested, 101);
                }
                next += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(next, 101);
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
        assert!(resources.snapshot().peak_reserved_memory_bytes <= 2 * 1024);
        store.close().unwrap();
        assert_eq!(resources.snapshot().reserved_temporary_disk_bytes, 0);
    }
    assert_eq!(directory.path().read_dir().unwrap().count(), 0);
}

#[test]
fn invalid_order_quality_cap_and_frame_size_cannot_publish_a_partial_replay() {
    let directory = tempfile::tempdir().unwrap();
    let resources = budget(4096, 1 << 20);
    for kind in 0..5 {
        let mut writer = ReplayWriter::new(
            directory.path(),
            &resources,
            if kind == 4 { 100 } else { 1024 },
            false,
            1,
        )
        .unwrap();
        let mut value = record(0, false);
        match kind {
            0 => value.index = 1,
            1 => value.fields[3] = b"I",
            2 => value.mate = Mate::R1,
            3 => {
                writer.add(record(0, false)).unwrap();
                value.index = 1;
                value.fragment = 1;
            }
            _ => {}
        }
        assert!(writer.add(value).is_err());
        assert!(writer.finish([0; 32]).is_err());
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
        assert_eq!(resources.snapshot().reserved_temporary_disk_bytes, 0);
    }
    let mut writer = ReplayWriter::new(directory.path(), &resources, 1024, false, 0).unwrap();
    writer.add(record(0, false)).unwrap();
    assert!(writer.finish([0; 32]).is_err());
}

fn mutate(path: &Path, offset: u64, bytes: &[u8]) {
    let mut file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
    file.seek(SeekFrom::Start(offset)).unwrap();
    file.write_all(bytes).unwrap();
    file.sync_all().unwrap();
}

#[test]
fn corruption_foreign_frames_truncation_and_oversized_lengths_are_errors_before_callbacks() {
    let directory = tempfile::tempdir().unwrap();
    let resources = budget(4096, 1 << 20);
    for kind in 0..6 {
        let (store, digest) = complete(directory.path(), &resources, false, 2);
        let file = store.file.into_inner().unwrap();
        let path = file.path().unwrap().to_path_buf();
        let before = resources.snapshot().peak_reserved_memory_bytes;
        match kind {
            0 => mutate(&path, HEADER as u64 + PREFIX as u64, b"x"),
            1 => mutate(&path, HEADER as u64 + 24, &u64::MAX.to_le_bytes()),
            2 => mutate(&path, 8, &99u32.to_le_bytes()),
            3 => {
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(&path)
                    .unwrap()
                    .set_len(file.len().unwrap() - 1)
                    .unwrap();
            }
            4 => {
                let (other, _) = complete(directory.path(), &resources, false, 2);
                let mut foreign =
                    std::fs::File::open(other.file.lock().unwrap().path().unwrap()).unwrap();
                foreign.seek(SeekFrom::Start(HEADER as u64)).unwrap();
                let mut frame = vec![
                    0;
                    PREFIX
                        + CHECKSUM
                        + record(0, false).source.len()
                        + record(0, false)
                            .fields
                            .iter()
                            .map(|field| field.len())
                            .sum::<usize>()
                ];
                foreign.read_exact(&mut frame).unwrap();
                mutate(&path, HEADER as u64, &frame);
                other.close().unwrap();
            }
            _ => {
                // The checksum depends on a random file ID. Writing a fixed
                // byte can leave it unchanged; always flip the observed byte.
                let offset = file.len().unwrap() - 1;
                let mut input = std::fs::File::open(&path).unwrap();
                input.seek(SeekFrom::Start(offset)).unwrap();
                let mut byte = [0];
                input.read_exact(&mut byte).unwrap();
                byte[0] ^= 1;
                mutate(&path, offset, &byte);
            }
        }
        assert!(
            ReplayStore::open(file, &resources, 1024, digest).is_err(),
            "corruption case {kind} must be rejected"
        );
        assert!(resources.snapshot().peak_reserved_memory_bytes <= before.max(1024));
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
        assert_eq!(resources.snapshot().reserved_temporary_disk_bytes, 0);
    }
}

#[test]
fn actual_append_budget_failure_is_typed_and_failed_replay_releases_its_owners() {
    let directory = tempfile::tempdir().unwrap();
    let resources = budget(4096, HEADER as u64 + 10);
    let mut writer = ReplayWriter::new(directory.path(), &resources, 1024, false, 0).unwrap();
    let error = writer.add(record(0, false)).unwrap_err();
    assert!(error.downcast_ref::<ResourceExhausted>().is_some());
    assert_eq!(
        resources.snapshot().reserved_temporary_disk_bytes,
        HEADER as u64
    );
    assert!(writer.finish([0; 32]).is_err());
    assert_eq!(resources.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn replay_and_count_files_share_one_descriptor_limit() {
    use crate::kmer::budget::TemporaryFilesExhausted;
    use crate::kmer::identity::{CanonicalKey, KmerSize};
    use crate::kmer::sorted_file::SortedCountWriter;
    use std::num::NonZeroU64;

    let directory = tempfile::tempdir().unwrap();
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::ReplayBuffer as usize] = 4096;
    classes[MemoryClass::LookupCache as usize] = 4096;
    let budget = ResourceBudget::new(ResourceLimits {
        memory_bytes: 8192,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 1 << 20,
        temporary_open_files: 1,
    })
    .unwrap();
    let (replay, _) = complete(directory.path(), &budget, false, 2);
    assert_eq!(budget.snapshot().reserved_temporary_open_files, 1);
    let k = KmerSize::new(3).unwrap();
    let error = SortedCountWriter::new(directory.path(), k, &budget, 1)
        .err()
        .unwrap();
    assert!(error.downcast_ref::<TemporaryFilesExhausted>().is_some());
    let mut visited = 0;
    replay
        .visit(&mut |_| {
            visited += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(visited, 2);
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
    replay.close().unwrap();
    let mut count = SortedCountWriter::new(directory.path(), k, &budget, 1).unwrap();
    count
        .push(
            CanonicalKey::new(1, k).unwrap(),
            NonZeroU64::new(1).unwrap(),
        )
        .unwrap();
    drop(count);
    assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
    assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().peak_reserved_temporary_open_files, 1);
    assert!(budget.snapshot().valid);
}
