use super::*;
use crate::kmer::budget::{ResourceExhausted, ResourceLimits};
use sha2::{Digest, Sha256};

fn resources(bytes: u64) -> Arc<ResourceBudget> {
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[MemoryClass::WorkerQueues as usize] = bytes;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: bytes.max(1),
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 0,
        temporary_open_files: 128,
    })
    .unwrap()
}

fn reference(hash: &mut Sha256, index: u64, mate: u8, source: &str, fields: [&[u8]; 4]) -> u64 {
    // Independently assemble a whole frame rather than using the chunk writer.
    let mut frame = Vec::new();
    frame.extend(index.to_le_bytes());
    frame.push(mate);
    for value in [
        source.as_bytes(),
        fields[0],
        fields[1],
        fields[2],
        fields[3],
    ] {
        frame.extend((value.len() as u64).to_le_bytes());
        frame.extend(value);
    }
    hash.update(&frame);
    frame.len() as u64
}

#[test]
fn empty_and_partial_streams_match_the_existing_serial_digest() {
    for count in [0, 1, 3, 127] {
        let budget = resources(IdentityWriter::reservation_bytes().unwrap());
        let mut writer = IdentityWriter::new(&budget).unwrap();
        let mut hash = Sha256::new();
        hash.update(DOMAIN);
        let mut bytes = 0;
        for index in 0..count {
            let fields = [b"@read original".as_slice(), b"ACNTG", b"+read", b"!@+#I"];
            writer
                .record(index, Mate::Unpaired, "source with spaces", fields)
                .unwrap();
            bytes += reference(&mut hash, index, 0, "source with spaces", fields);
        }
        let (digest, summary) = writer.finish().unwrap();
        assert_eq!(digest, <[u8; 32]>::from(hash.finalize()));
        assert_eq!(summary.selected_records, count);
        assert_eq!(summary.framed_bytes, bytes);
        assert_eq!(summary.buffer_bytes, BUFFER_BYTES);
        assert_eq!(summary.owned_buffers, 3);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
    }
}

#[test]
fn multiple_large_fields_and_mates_cross_every_chunk_without_reordering() {
    let budget = resources(IdentityWriter::reservation_bytes().unwrap());
    let mut writer = IdentityWriter::new(&budget).unwrap();
    let mut hash = Sha256::new();
    hash.update(DOMAIN);
    let sequence = [b"ACNT".as_slice(); BUFFER_BYTES].concat();
    let quality = vec![b'I'; sequence.len()];
    let header = "@α".repeat(BUFFER_BYTES / 2);
    let mut bytes = 0;
    for index in 0..7 {
        let mate = if index % 2 == 0 { Mate::R1 } else { Mate::R2 };
        let source = if index % 2 == 0 { "r1.gz" } else { "r2.gz" };
        let fields = [
            header.as_bytes(),
            sequence.as_slice(),
            b"+original",
            quality.as_slice(),
        ];
        writer.record(index, mate, source, fields).unwrap();
        bytes += reference(
            &mut hash,
            index,
            if index % 2 == 0 { 1 } else { 2 },
            source,
            fields,
        );
        assert_eq!(writer.current.capacity(), BUFFER_BYTES);
    }
    let (digest, summary) = writer.finish().unwrap();
    assert_eq!(digest, <[u8; 32]>::from(hash.finalize()));
    assert_eq!(summary.framed_bytes, bytes);
    assert_eq!(summary.selected_records, 7);
    assert_eq!(
        budget.snapshot().peak_reserved_memory_bytes,
        IdentityWriter::reservation_bytes().unwrap()
    );
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn rejected_order_preserves_the_previous_complete_prefix() {
    let budget = resources(IdentityWriter::reservation_bytes().unwrap());
    let mut writer = IdentityWriter::new(&budget).unwrap();
    let fields = [b"@id".as_slice(), b"NNAC", b"+", b"IIII"];
    assert!(writer.record(1, Mate::Unpaired, "stdin", fields).is_err());
    writer.record(0, Mate::Unpaired, "stdin", fields).unwrap();
    assert!(writer.record(0, Mate::Unpaired, "stdin", fields).is_err());
    let mut expected = Sha256::new();
    expected.update(DOMAIN);
    reference(&mut expected, 0, 0, "stdin", fields);
    let (digest, summary) = writer.finish().unwrap();
    assert_eq!(digest, <[u8; 32]>::from(expected.finalize()));
    assert_eq!(summary.selected_records, 1);
}

#[test]
fn reservation_failure_and_dropping_queued_bytes_release_all_owners() {
    let bytes = IdentityWriter::reservation_bytes().unwrap();
    let budget = resources(bytes - 1);
    let error = IdentityWriter::new(&budget).err().unwrap();
    assert!(error.downcast_ref::<ResourceExhausted>().is_some());
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 0);
    let budget = resources(bytes);
    let mut writer = IdentityWriter::new(&budget).unwrap();
    let data = vec![b'A'; 3 * BUFFER_BYTES + 17];
    for index in 0..9 {
        writer
            .record(index, Mate::Unpaired, "stdin", [b"@id", &data, b"+", &data])
            .unwrap();
    }
    drop(writer);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(Arc::strong_count(&budget), 1);
}

#[test]
fn a_stopped_or_panicking_worker_cannot_publish_an_input_identity() {
    let budget = resources(IdentityWriter::reservation_bytes().unwrap());
    let mut writer = IdentityWriter::new(&budget).unwrap();
    writer.close_channels();
    assert!(
        writer
            .record(0, Mate::Unpaired, "stdin", [b"@id", b"AC", b"+", b"II"])
            .is_err()
    );
    assert!(writer.finish().is_err());
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    let mut writer = IdentityWriter::new(&budget).unwrap();
    writer.close_channels();
    writer.handle.take().unwrap().join().unwrap().unwrap();
    let (sender, _receiver) = mpsc::sync_channel(1);
    writer.work = Some(sender);
    writer.handle = Some(std::thread::spawn(|| {
        panic!("injected identity worker panic")
    }));
    let error = writer.finish().unwrap_err();
    assert!(error.to_string().contains("panicked"));
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert_eq!(Arc::strong_count(&budget), 1);
}

#[test]
fn drained_worker_state_continues_the_identical_hash_stream_inline() {
    let budget = resources(IdentityWriter::reservation_bytes().unwrap());
    let mut writer = IdentityWriter::new(&budget).unwrap();
    let mut expected = Sha256::new();
    expected.update(DOMAIN);
    let bases = b"ACNT".repeat(BUFFER_BYTES);
    let quality = vec![b'I'; bases.len()];
    let fields = [
        b"@original".as_slice(),
        bases.as_slice(),
        b"+id",
        quality.as_slice(),
    ];
    for index in 0..3 {
        writer.record(index, Mate::R1, "source α", fields).unwrap();
        reference(&mut expected, index, 1, "source α", fields);
    }
    let (mut state, summary) = writer.finish_state().unwrap();
    assert_eq!(summary.selected_records, 3);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    for index in 3u64..7 {
        let mut frame = Vec::new();
        frame.extend(index.to_le_bytes());
        frame.push(2);
        for value in [
            "source β".as_bytes(),
            fields[0],
            fields[1],
            fields[2],
            fields[3],
        ] {
            frame.extend((value.len() as u64).to_le_bytes());
            frame.extend(value);
        }
        state.update(&frame).unwrap();
        reference(&mut expected, index, 2, "source β", fields);
    }
    assert_eq!(
        state.finalize().unwrap(),
        <[u8; 32]>::from(expected.finalize())
    );
}
