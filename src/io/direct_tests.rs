use super::*;
use std::collections::BTreeMap;
use std::io::Cursor;

fn state() -> FastqReadState {
    FastqReadState {
        chunks: Vec::new(),
        chunk_index: 0,
        seqs: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
        direct: Some(DirectCounter::new(3, 1 << 20, 4096, 0, 0, 0).unwrap()),
    }
}

fn record(id: &str, sequence: &str) -> String {
    format!("@{id}\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()))
}

fn entries(counter: DirectCounter) -> (BTreeMap<u64, u64>, [u8; 32]) {
    let digest = counter.selected_digest.finalize().unwrap();
    let store = counter.builder.finalize(digest).unwrap();
    let mut counts = BTreeMap::new();
    store
        .as_store()
        .visit(&mut |key, count| {
            counts.insert(key.bits(), count.get());
            Ok(())
        })
        .unwrap();
    (counts, digest)
}

#[test]
fn bounded_remote_cache_miss_and_hit_preserve_exact_selection_identity_and_replay() {
    use clap::Parser;
    use sha2::{Digest, Sha256};
    use std::io::{BufRead, Read, Write};
    use std::net::TcpListener;

    let first = record("first", "ACGTACNAA");
    let second = record("second", "TTTTTT");
    for maximum in [0, 1] {
        let directory = tempfile::tempdir().unwrap();
        let cache = crate::cache::CacheConfig::new(Some(&directory.path().join("cache"))).unwrap();
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let url = format!("http://{}/reads.fastq.gz", listener.local_addr().unwrap());
        let mut body = Vec::new();
        for record in [&first, &second] {
            let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            gzip.write_all(record.replace('\n', "\r\n").as_bytes())
                .unwrap();
            body.extend(gzip.finish().unwrap());
        }
        // Exercise both original application gzip and an additional automatic
        // HTTP gzip layer without changing normalized count/replay evidence.
        let content_encoding = if maximum == 1 {
            let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            gzip.write_all(&body).unwrap();
            body = gzip.finish().unwrap();
            "Content-Encoding: gzip\r\n"
        } else {
            ""
        };
        let server = std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(std::time::Duration::from_secs(5)))
                .unwrap();
            let mut request = std::io::BufReader::new(stream.try_clone().unwrap());
            loop {
                let mut line = String::new();
                assert!(request.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
            }
            write!(
                stream,
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\n{content_encoding}Connection: close\r\n\r\n",
                body.len()
            )
            .unwrap();
            stream.write_all(&body).unwrap();
        });
        let args = crate::cli::Args::parse_from([
            "sharkmer",
            "--ena",
            "cached-metadata",
            "-s",
            "sample",
            "-k",
            "3",
            "--max-reads",
            &maximum.to_string(),
            "--read-threading",
            "--pcr-primers",
            "name=target,forward=ACG,reverse=ACG",
            "--fastq-record-bytes",
            "128",
            "--count-temp-dir",
            directory.path().to_str().unwrap(),
        ]);
        let selected = if maximum == 1 {
            vec![&first]
        } else {
            vec![&first, &second]
        };
        let mut framed = Sha256::new();
        framed.update(b"sharkmer-selected-fastq-v1\0");
        for (index, record) in selected.iter().enumerate() {
            framed.update((index as u64).to_le_bytes());
            framed.update([0]);
            for field in std::iter::once(url.as_bytes()).chain(record.lines().map(str::as_bytes)) {
                framed.update((field.len() as u64).to_le_bytes());
                framed.update(field);
            }
        }
        let digest: [u8; 32] = framed.finalize().into();
        let mut expected = BTreeMap::from([(6, 2), (44, 2)]);
        if maximum == 0 {
            expected.insert(0, 4);
        }

        // The second pass must use the verified cache after the sole HTTP server
        // has exited. Each pass still performs one exact ingestion and replay.
        for pass in 0..2 {
            let prepared = prepare_ingestion(&args).unwrap();
            let budget = Arc::clone(prepared.result_budget().unwrap());
            let (mut state, reads, bases, occurrences, plan) = ingest_reads(
                prepared,
                Some(ena::parse(&format!("fastq_ftp\n{url}\n"), "fixture", Some(&budget)).unwrap()),
                Some(&cache),
                false,
            )
            .unwrap();
            assert_eq!(reads, selected.len() as u64);
            assert_eq!(bases, if maximum == 1 { 8 } else { 14 });
            assert_eq!(occurrences, if maximum == 1 { 4 } else { 8 });
            assert!(Arc::ptr_eq(&budget, &state.direct.as_ref().unwrap().budget));
            assert!(state.chunks.is_empty() && state.seqs.is_empty());
            assert!(matches!(plan.source, ReadSourcePlan::CachedRemote(_)));
            let replay = plan.replay.as_ref().unwrap();
            assert_eq!(replay.store.record_count(), reads);
            let mut observed = Vec::new();
            replay
                .store
                .visit(&mut |record| {
                    assert_eq!(record.source, url);
                    observed.push(record.fields.map(<[u8]>::to_vec));
                    Ok(())
                })
                .unwrap();
            let wanted: Vec<_> = selected
                .iter()
                .map(|record| {
                    record
                        .lines()
                        .map(|line| line.as_bytes().to_vec())
                        .collect::<Vec<_>>()
                })
                .collect();
            assert_eq!(
                observed
                    .iter()
                    .map(|fields| fields.to_vec())
                    .collect::<Vec<_>>(),
                wanted
            );
            assert_eq!(
                entries(state.direct.take().unwrap()),
                (expected.clone(), digest)
            );
            drop(plan);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
            if pass == 0 {
                // Readback is normalized bytes, independent of gzip member and
                // original CRLF boundaries; no original records are retained.
                let path = cache.lookup(&url, maximum).unwrap().unwrap();
                let mut decoded = String::new();
                flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
                    .read_to_string(&mut decoded)
                    .unwrap();
                assert_eq!(
                    decoded,
                    selected.iter().map(|s| s.as_str()).collect::<String>()
                );
            }
        }
        server.join().unwrap();
    }
}

#[test]
fn direct_stream_matches_scalar_counts_and_never_retains_read_batches() {
    let sequences = ["ACGTACNAA", "NNACNGTNN", "NNN", "ANANAN", "TTTCGA"];
    let mut input = String::new();
    let mut expected = BTreeMap::new();
    let mut bases = 0;
    let mut raw_bases = 0;
    for index in 0..1005 {
        let sequence = sequences[index % sequences.len()];
        input += &record(&index.to_string(), sequence);
        bases += sequence.bytes().filter(|b| *b != b'N').count() as u64;
        raw_bases += sequence.len() as u64;
        for window in sequence.as_bytes().windows(3) {
            let encoded = window
                .iter()
                .map(|base| match base {
                    b'A' => Some(0u64),
                    b'C' => Some(1),
                    b'G' => Some(2),
                    b'T' => Some(3),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>();
            if let Some(encoded) = encoded {
                let forward = encoded.iter().fold(0, |bits, base| (bits << 2) | base);
                let reverse = encoded
                    .iter()
                    .rev()
                    .fold(0, |bits, base| (bits << 2) | (3 - base));
                *expected.entry(forward.min(reverse)).or_insert(0) += 1;
            }
        }
    }
    let mut state = state();
    assert!(
        !read_fastq(
            Cursor::new(input),
            &mut state,
            0,
            1,
            "test",
            &ProgressBar::hidden()
        )
        .unwrap()
    );
    assert_eq!(state.n_reads_read, 1005);
    assert_eq!(state.n_bases_read, raw_bases);
    assert!(state.chunks.is_empty());
    assert!(state.seqs.is_empty());
    let counter = state.direct.take().unwrap();
    assert_eq!(counter.valid_bases, bases);
    assert_eq!(counter.occurrences, expected.values().sum());
    assert_eq!(entries(counter).0, expected);
}

#[test]
fn reusable_key_buffer_preserves_long_read_windows_n_gaps_and_record_tails() {
    for k in [1, 3, 19, 31] {
        let mut counter = DirectCounter::new(k, 1 << 20, 16384, 0, 0, 0).unwrap();
        let pointer = counter.key_batch.as_ptr();
        let capacity = counter.key_batch.capacity();
        let mut random = 913_u64;
        let long = (0..3 * KMER_BATCH_CAPACITY + 31)
            .map(|position| {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                if position % 137 == 0 {
                    b'N'
                } else {
                    b"ACGT"[((random >> 32) & 3) as usize]
                }
            })
            .collect::<Vec<_>>();
        let short = vec![b'A'; k + 3];
        let mut expected = BTreeMap::new();
        let mut bases = 0;
        for (index, sequence) in [&long, &short].into_iter().enumerate() {
            bases += sequence.iter().filter(|&&byte| byte != b'N').count() as u64;
            for window in sequence.windows(k) {
                let encoded = window
                    .iter()
                    .map(|byte| match byte {
                        b'A' => Some(0_u64),
                        b'C' => Some(1),
                        b'G' => Some(2),
                        b'T' => Some(3),
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>();
                if let Some(encoded) = encoded {
                    let forward = encoded.iter().fold(0, |bits, base| (bits << 2) | base);
                    let reverse = encoded
                        .iter()
                        .rev()
                        .fold(0, |bits, base| (bits << 2) | (3 - base));
                    *expected.entry(forward.min(reverse)).or_insert(0_u64) += 1;
                }
            }
            let quality = vec![b'I'; sequence.len()];
            counter
                .ingest(
                    [b"@id", sequence, b"+", &quality],
                    "long",
                    index as u64,
                    Mate::Unpaired,
                )
                .unwrap();
            assert!(counter.key_batch.is_empty());
            assert_eq!(counter.next_key, 0);
            assert_eq!(counter.key_batch.as_ptr(), pointer);
            assert_eq!(counter.key_batch.capacity(), capacity);
        }
        assert_eq!(counter.valid_bases, bases);
        assert_eq!(counter.occurrences, expected.values().sum());
        assert_eq!(entries(counter).0, expected);
    }
}

#[test]
fn failed_batch_growth_retries_only_unaccepted_keys_and_releases_reservations() {
    let mut counter = DirectCounter::new(6, 1024, 256, 0, 0, 0).unwrap();
    let budget = Arc::clone(&counter.budget);
    let bytes = (counter.key_batch.capacity() * size_of::<CanonicalKey>()) as u64;
    assert_eq!(
        budget.snapshot().reserved_by_class[MemoryClass::WorkerQueues as usize],
        bytes
    );
    let keys = (0..4096)
        .map(|bits| CanonicalKey::new(bits, KmerSize::new(6).unwrap()).unwrap())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .take(15)
        .collect::<Vec<_>>();
    counter.key_batch.extend(&keys);
    let pressure = budget.reserve(MemoryClass::CountSlots, 512).unwrap();
    assert!(counter.flush_keys().is_err());
    assert_eq!(counter.builder.occurrences(), 14);
    assert_eq!(counter.next_key, 14);
    assert_eq!(counter.key_batch.len(), 15);
    assert!(
        counter
            .ingest(
                [b"@next", b"AAAAAAA", b"+", b"IIIIIII"],
                "test",
                0,
                Mate::Unpaired
            )
            .is_err()
    );
    assert_eq!(counter.builder.occurrences(), 14);
    drop(pressure);
    counter.flush_keys().unwrap();
    assert_eq!(counter.builder.occurrences(), 15);
    assert!(counter.key_batch.is_empty());
    let (counts, _) = entries(counter);
    assert_eq!(counts.len(), 15);
    assert!(counts.values().all(|count| *count == 1));
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    assert!(budget.snapshot().valid);
}

#[test]
fn selected_identity_retains_ids_qualities_mates_source_order_and_cap() {
    let input = record("a", "ACGT") + &record("b", "TGCA");
    let digest = |input: &str, cap, source: &str| {
        let mut state = state();
        read_fastq(
            Cursor::new(input),
            &mut state,
            cap,
            1,
            source,
            &ProgressBar::hidden(),
        )
        .unwrap();
        let n = state.n_reads_read;
        let (counts, sha) = entries(state.direct.take().unwrap());
        (n, counts, sha)
    };
    let baseline = digest(&input, 0, "source");
    assert_eq!(baseline, digest(&input.replace('\n', "\r\n"), 0, "source"));
    for changed in [
        input.replace("@a", "@changed"),
        input.replace("IIII", "JJJJ"),
        record("b", "TGCA") + &record("a", "ACGT"),
    ] {
        let observed = digest(&changed, 0, "source");
        assert_eq!(observed.1, baseline.1);
        assert_ne!(observed.2, baseline.2);
    }
    assert_ne!(digest(&input, 0, "other-source").2, baseline.2);
    let capped = digest(&input, 1, "source");
    assert_eq!(capped.0, 1);
    assert_ne!(capped.2, baseline.2);
    let mut state = state();
    read_fastq_paired(
        Cursor::new(record("a", "ACGT")),
        Cursor::new(record("b", "TGCA")),
        &mut state,
        0,
        1,
        "source",
        "source",
        &ProgressBar::hidden(),
    )
    .unwrap();
    let paired = entries(state.direct.take().unwrap());
    assert_eq!(paired.0, baseline.1);
    assert_ne!(paired.1, baseline.2);
}

#[test]
fn paired_selection_matches_legacy_at_cap_and_both_uneven_eof_directions() {
    for (r1, r2, cap) in [(3, 2, 0), (2, 3, 0), (3, 3, 4)] {
        let source1 = (0..r1)
            .map(|n| record(&format!("r1-{n}"), "AAAC"))
            .collect::<String>();
        let source2 = (0..r2)
            .map(|n| record(&format!("r2-{n}"), "TTGC"))
            .collect::<String>();
        let mut direct = state();
        let mut legacy = FastqReadState {
            chunks: vec![Chunk::new(&3)],
            chunk_index: 0,
            seqs: Vec::new(),
            n_reads_read: 0,
            n_bases_read: 0,
            direct: None,
        };
        let run = |state: &mut FastqReadState| {
            read_fastq_paired(
                Cursor::new(&source1),
                Cursor::new(&source2),
                state,
                cap,
                1,
                "r1",
                "r2",
                &ProgressBar::hidden(),
            )
            .unwrap()
        };
        assert_eq!(run(&mut direct), run(&mut legacy));
        assert_eq!(direct.n_reads_read, legacy.n_reads_read);
        assert_eq!(direct.n_bases_read, legacy.n_bases_read);
        drain_batch(&mut legacy, 1).unwrap();
        let expected = legacy.chunks[0]
            .get_kmer_counts()
            .iter()
            .map(|(&key, &count)| (key, u64::from(count)))
            .collect();
        assert_eq!(entries(direct.direct.take().unwrap()).0, expected);
    }
}

#[test]
fn direct_gzip_members_and_parse_errors_preserve_selection_and_bounds() {
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("reads.fastq.gz");
    let mut compressed = Vec::new();
    for sequence in ["ACGT", "NGCA"] {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder
            .write_all(record("id", sequence).as_bytes())
            .unwrap();
        compressed.extend(encoder.finish().unwrap());
    }
    std::fs::write(&path, &compressed).unwrap();
    let mut first = state();
    read_fastq(
        open_fastq_reader(&path).unwrap(),
        &mut first,
        0,
        1,
        "gzip",
        &ProgressBar::hidden(),
    )
    .unwrap();
    assert_eq!(first.n_reads_read, 2);
    assert_eq!(first.direct.as_ref().unwrap().valid_bases, 7);
    *compressed.last_mut().unwrap() ^= 0xff;
    std::fs::write(&path, &compressed).unwrap();
    assert!(
        read_fastq(
            open_fastq_reader(&path).unwrap(),
            &mut state(),
            0,
            1,
            "gzip",
            &ProgressBar::hidden()
        )
        .is_err()
    );
    for bad in ["@x\nACX\n+\nIII\n", "@x\nAC\n+\nI\n", "@x\nAC\n+\n"] {
        assert!(
            read_fastq(
                Cursor::new(bad),
                &mut state(),
                0,
                1,
                "bad",
                &ProgressBar::hidden()
            )
            .is_err()
        );
    }
}

#[test]
fn paired_replay_records_exactly_the_counted_selection_including_either_unmatched_tail() {
    for (r1, r2, cap, expected_count) in [(3, 2, 0, 5), (2, 3, 0, 5), (3, 3, 4, 4)] {
        let directory = tempfile::tempdir().unwrap();
        let mut state = state();
        let mut counter = DirectCounter::new(3, 1 << 20, 4096, 1 << 20, 8192, 0).unwrap();
        counter.replay = Some(
            crate::replay::ReplayWriter::new(directory.path(), &counter.budget, 4096, true, cap)
                .unwrap(),
        );
        state.direct = Some(counter);
        let source1 = (0..r1)
            .map(|n| record(&format!("r1-{n}"), "AANCG"))
            .collect::<String>();
        let source2 = (0..r2)
            .map(|n| record(&format!("r2-{n}"), "TTNGC"))
            .collect::<String>();
        read_fastq_paired(
            Cursor::new(&source1),
            Cursor::new(&source2),
            &mut state,
            cap,
            1,
            "original-r1",
            "original-r2",
            &ProgressBar::hidden(),
        )
        .unwrap();
        assert_eq!(state.n_reads_read, expected_count);
        let mut counter = state.direct.take().unwrap();
        let digest = counter.selected_digest.clone().finalize().unwrap();
        let replay = counter.replay.take().unwrap().finish(digest).unwrap();
        let mut observed = 0;
        replay
            .visit(&mut |row| {
                let first = row.mate == Mate::R1;
                assert_eq!(row.index, observed);
                assert_eq!(row.fragment, observed / 2);
                assert_eq!(
                    row.source,
                    if first { "original-r1" } else { "original-r2" }
                );
                assert_eq!(row.fields[1], if first { b"AANCG" } else { b"TTNGC" });
                assert_eq!(row.fields[3], b"IIIII");
                assert_eq!(
                    row.fields[0],
                    format!("@r{}-{}", if first { 1 } else { 2 }, observed / 2).as_bytes()
                );
                observed += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected_count);
        replay.close().unwrap();
        assert_eq!(counter.budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            counter.budget.snapshot().reserved_by_class[MemoryClass::ReplayBuffer as usize],
            0
        );
    }
}

#[test]
fn pipeline_digest_joins_once_matches_inline_and_seals_further_ingestion() {
    let directory = tempfile::tempdir().unwrap();
    let mut pipeline = DirectCounter::spillable_with_workers(
        3,
        1 << 20,
        4096,
        1024,
        1 << 24,
        directory.path(),
        0,
        0,
        true,
        0,
    )
    .unwrap();
    let mut inline = DirectCounter::new(3, 1 << 20, 1024, 1 << 24, 0, 0).unwrap();
    let fields = [b"@read id".as_slice(), b"ACNTACGT", b"+read", b"IIIIIIII"];
    for index in 0..3 {
        pipeline
            .ingest(fields, "stdin", index, Mate::Unpaired)
            .unwrap();
        inline
            .ingest(fields, "stdin", index, Mate::Unpaired)
            .unwrap();
    }
    let budget = Arc::clone(&pipeline.budget);
    let digest = pipeline.input_digest().unwrap();
    assert_eq!(digest, inline.input_digest().unwrap());
    assert_eq!(digest, pipeline.input_digest().unwrap());
    assert!(pipeline.identity_writer.is_none());
    assert_eq!(pipeline.identity_summary.unwrap().selected_records, 3);
    assert!(inline.identity_summary.is_none());
    let occurrences = pipeline.builder.occurrences();
    assert!(pipeline.ingest(fields, "stdin", 3, Mate::Unpaired).is_err());
    assert_eq!(pipeline.builder.occurrences(), occurrences);
    drop(pipeline);
    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn hybrid_transfer_preserves_record_tails_mates_hash_state_and_original_pass_replay() {
    for threshold in [256, usize::MAX] {
        let directory = tempfile::tempdir().unwrap();
        let mut hybrid = DirectCounter::spillable_with_policy(
            3,
            1 << 20,
            4096,
            8192,
            1 << 24,
            directory.path(),
            32768,
            0,
            true,
            Some(threshold),
            0,
        )
        .unwrap();
        hybrid.replay = Some(
            crate::replay::ReplayWriter::new(directory.path(), &hybrid.budget, 8192, true, 7)
                .unwrap(),
        );
        let mut inline = DirectCounter::new(3, 1 << 20, 8192, 0, 0, 0).unwrap();
        let sequence = format!("{}N{}", "ACGT".repeat(258), "TGCA".repeat(261));
        let quality = "I".repeat(sequence.len());
        let fields = [
            b"@original id".as_slice(),
            sequence.as_bytes(),
            b"+original",
            quality.as_bytes(),
        ];
        for index in 0..7 {
            let mate = if index % 2 == 0 { Mate::R1 } else { Mate::R2 };
            hybrid.ingest(fields, "paired α", index, mate).unwrap();
            inline.ingest(fields, "paired α", index, mate).unwrap();
        }
        let digest = hybrid.input_digest().unwrap();
        assert_eq!(digest, inline.input_digest().unwrap());
        assert_eq!(hybrid.valid_bases, inline.valid_bases);
        assert_eq!(hybrid.occurrences, inline.occurrences);
        let budget = Arc::clone(&hybrid.budget);
        let summary = hybrid.owner_transfer_summary;
        if threshold == 256 {
            let summary = summary.unwrap();
            assert_eq!(summary.selected_records_before_transfer, 1);
            assert_eq!(summary.threshold_bytes, 256);
            assert!(summary.current_table_bytes >= 256);
            assert_eq!(hybrid.identity_summary.unwrap().selected_records, 1);
        } else {
            assert!(summary.is_none());
            assert_eq!(hybrid.identity_summary.unwrap().selected_records, 7);
        }
        let replay = hybrid.replay.take().unwrap().finish(digest).unwrap();
        let mut visited = 0;
        replay
            .visit(&mut |record| {
                assert_eq!(record.fields, fields);
                assert_eq!(record.index, visited);
                visited += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(visited, 7);
        replay.close().unwrap();
        let expected = entries(inline).0;
        let store = hybrid.builder.finalize(digest).unwrap();
        let mut actual = BTreeMap::new();
        store
            .as_store()
            .visit(&mut |key, count| {
                actual.insert(key.bits(), count.get());
                Ok(())
            })
            .unwrap();
        assert_eq!(actual, expected);
        drop(store);
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
    }
}
