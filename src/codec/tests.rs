use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};
use std::io::Cursor;
use std::sync::atomic::{AtomicUsize, Ordering};

fn budget(memory: u64) -> Arc<ResourceBudget> {
    ResourceBudget::new(ResourceLimits {
        memory_bytes: memory.max(1),
        process_overhead_bytes: 0,
        class_bytes: [memory; MEMORY_CLASSES],
        temporary_disk_bytes: 0,
        temporary_open_files: 0,
    })
    .unwrap()
}

fn member(data: &[u8], fields: usize, header_crc: bool) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(data).unwrap();
    let compressed = encoder.finish().unwrap();
    if fields == 0 && !header_crc {
        return compressed;
    }
    let mut bytes = compressed[..10].to_vec();
    bytes[3] = 4 | 8 | 16 | if header_crc { 2 } else { 0 };
    let extra = fields.min(u16::MAX as usize);
    bytes.extend((extra as u16).to_le_bytes());
    bytes.extend(std::iter::repeat_n(0xfe, extra));
    for _ in 0..2 {
        bytes.extend(std::iter::repeat_n(0xfd, fields));
        bytes.push(0);
    }
    if header_crc {
        let mut crc = flate2::Crc::new();
        crc.update(&bytes);
        bytes.extend((crc.sum() as u16).to_le_bytes());
    }
    bytes.extend_from_slice(&compressed[10..]);
    bytes
}

struct ObservedRead {
    calls: Arc<AtomicUsize>,
    reader: Cursor<Vec<u8>>,
}
impl Read for ObservedRead {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        self.reader.read(bytes)
    }
}

struct ObservedWrite {
    calls: Arc<AtomicUsize>,
}
impl Write for ObservedWrite {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        Ok(())
    }
}

#[test]
fn admission_refuses_before_decoder_reads_or_encoder_writes() {
    for kind in 0..3 {
        let calls = Arc::new(AtomicUsize::new(0));
        let ledger = budget(0);
        let error = match kind {
            0 => Decoder::new(
                ObservedRead {
                    calls: calls.clone(),
                    reader: Cursor::new(member(b"reads", 64, true)),
                },
                Some(&ledger),
            )
            .err()
            .unwrap(),
            1 => Encoder::new(
                ObservedWrite {
                    calls: calls.clone(),
                },
                Some(&ledger),
            )
            .err()
            .unwrap(),
            _ => Buffered::new(
                ObservedRead {
                    calls: calls.clone(),
                    reader: Cursor::new(Vec::new()),
                },
                Some(&ledger),
            )
            .err()
            .unwrap(),
        };
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::StageMemory(MemoryClass::RetainedResults)
        );
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }
}

#[test]
fn gzip_members_and_maximum_header_fields_preserve_literal_bytes_and_late_errors() {
    let mut bytes = member(b"first\r\n", 65535, true);
    bytes.extend(member(b"", 65535, true));
    bytes.extend(member(b"second\n", 65535, false));
    let ledger = budget(2 << 20);
    let mut decoder = Decoder::new(Cursor::new(&bytes), Some(&ledger)).unwrap();
    let mut observed = Vec::new();
    decoder.read_to_end(&mut observed).unwrap();
    assert_eq!(observed, b"first\r\nsecond\n");
    assert!(ledger.snapshot().reserved_memory_bytes > 0);
    drop(decoder);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    let mut too_long = Decoder::new(Cursor::new(member(b"x", 65536, true)), Some(&ledger)).unwrap();
    assert!(io::copy(&mut too_long, &mut io::sink()).is_err());
    drop(too_long);
    bytes.extend(member(b"third\n", 0, false));
    *bytes.last_mut().unwrap() ^= 0xff;
    let mut corrupt = Decoder::new(Cursor::new(&bytes), Some(&ledger)).unwrap();
    let mut observed = Vec::new();
    assert!(corrupt.read_to_end(&mut observed).is_err());
    assert!(observed.starts_with(b"first\r\nsecond\n"));
    drop(corrupt);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn paired_decoders_share_the_unchanged_parent_limit_with_competing_count_memory() {
    let bytes = member(b"reads", 0, false);
    let one = Decoder::<Cursor<&Vec<u8>>>::allocation_bytes().unwrap() as u64;
    let ledger = budget(3 * one - 1);
    let slots = ledger.reserve(MemoryClass::CountSlots, one).unwrap();
    let first = Decoder::new(Cursor::new(&bytes), Some(&ledger)).unwrap();
    let error = Decoder::new(Cursor::new(&bytes), Some(&ledger))
        .err()
        .unwrap();
    assert_eq!(
        error.downcast_ref::<ResourceExhausted>().unwrap().resource,
        ResourceKind::RunMemory
    );
    drop(first);
    drop(slots);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    assert!(ledger.snapshot().valid);
}

#[test]
fn inner_reader_drops_before_its_codec_reservation() {
    struct Checked {
        ledger: Arc<ResourceBudget>,
        source: Cursor<Vec<u8>>,
        dropped: Arc<AtomicUsize>,
    }
    impl Read for Checked {
        fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
            self.source.read(bytes)
        }
    }
    impl Drop for Checked {
        fn drop(&mut self) {
            assert!(self.ledger.snapshot().reserved_memory_bytes > 0);
            self.dropped.fetch_add(1, Ordering::SeqCst);
        }
    }
    let ledger = budget(2 << 20);
    let dropped = Arc::new(AtomicUsize::new(0));
    let reader = Decoder::new(
        Checked {
            ledger: ledger.clone(),
            source: Cursor::new(member(b"reads", 0, false)),
            dropped: dropped.clone(),
        },
        Some(&ledger),
    )
    .unwrap();
    drop(reader);
    assert_eq!(dropped.load(Ordering::SeqCst), 1);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn encoder_guard_survives_trailer_flush_and_error_destruction() {
    struct Checked {
        ledger: Arc<ResourceBudget>,
        fail: bool,
        dropped: Arc<AtomicUsize>,
    }
    impl Write for Checked {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            assert!(self.ledger.snapshot().reserved_memory_bytes > 0);
            if self.fail {
                Err(io::Error::other("controlled output error"))
            } else {
                Ok(bytes.len())
            }
        }
        fn flush(&mut self) -> io::Result<()> {
            assert!(self.ledger.snapshot().reserved_memory_bytes > 0);
            Ok(())
        }
    }
    impl Drop for Checked {
        fn drop(&mut self) {
            if self.fail {
                assert!(self.ledger.snapshot().reserved_memory_bytes > 0);
            }
            self.dropped.fetch_add(1, Ordering::SeqCst);
        }
    }
    for fail in [false, true] {
        let ledger = budget(2 << 20);
        let dropped = Arc::new(AtomicUsize::new(0));
        let mut encoder = Encoder::new(
            Checked {
                ledger: ledger.clone(),
                fail,
                dropped: dropped.clone(),
            },
            Some(&ledger),
        )
        .unwrap();
        encoder.write_all(b"one small cache record\n").unwrap();
        let result = encoder.finish();
        assert_eq!(result.is_err(), fail);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        drop(result);
        assert_eq!(dropped.load(Ordering::SeqCst), 1);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }
}

#[test]
fn codec_heap_observations_fit_the_audited_envelope_and_release_every_allocation() {
    const ENV: &str = "SHARKMER_TEST_CODEC_ALLOCATIONS";
    if std::env::var_os(ENV).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "codec::tests::codec_heap_observations_fit_the_audited_envelope_and_release_every_allocation", "--nocapture", "--test-threads=1"])
            .env(ENV, "1").status().unwrap();
        assert!(status.success());
        return;
    }
    let mut random = Vec::new();
    let mut state = 17_u64;
    for _ in 0..(1 << 20) {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        random.push((state >> 32) as u8);
    }
    let mut corrupt = member(&random, 65535, true);
    *corrupt.last_mut().unwrap() ^= 1;
    let mut repeated = Vec::new();
    for _ in 0..8 {
        repeated.extend(member(b"reads\n", 65535, true));
    }
    let cases = [
        member(b"", 0, false),
        member(&random, 0, false),
        member(&random, 65535, true),
        member(b"x", 65536, true),
        corrupt,
        repeated,
    ];
    for (index, bytes) in cases.iter().enumerate() {
        let ledger = budget(2 << 20);
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let mut decoder = Decoder::new(Cursor::new(bytes), Some(&ledger)).unwrap();
        let decoded = io::copy(&mut decoder, &mut io::sink());
        let actual = crate::PEAK_ALLOC.peak_usage() - baseline;
        let reserved = ledger.snapshot().reserved_memory_bytes;
        assert!(
            actual as u64 <= reserved,
            "case={index} actual={actual} reserved={reserved}"
        );
        assert_eq!(decoded.is_err(), index == 3 || index == 4);
        drop(decoded);
        drop(decoder);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        println!("codec_allocator role=decode case={index} observed={actual} reserved={reserved}");
    }
    for bytes in [&b""[..], &b"reads\n"[..], &random[..]] {
        let ledger = budget(2 << 20);
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let mut encoder = Encoder::new(io::sink(), Some(&ledger)).unwrap();
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap();
        let actual = crate::PEAK_ALLOC.peak_usage() - baseline;
        let reserved = ledger.snapshot().peak_reserved_memory_bytes;
        assert!(
            actual as u64 <= reserved,
            "actual={actual} reserved={reserved}"
        );
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        println!(
            "codec_allocator role=encode input_bytes={} observed={actual} reserved={reserved}",
            bytes.len()
        );
    }
}
