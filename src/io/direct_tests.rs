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
        direct: Some(DirectCounter::new(3, 1 << 20, 4096).unwrap()),
    }
}

fn record(id: &str, sequence: &str) -> String {
    format!("@{id}\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()))
}

fn entries(counter: DirectCounter) -> (BTreeMap<u64, u64>, [u8; 32]) {
    let digest = counter.selected_digest.finalize().into();
    let store = counter.builder.finalize(digest);
    let mut counts = BTreeMap::new();
    store
        .visit(&mut |key, count| {
            counts.insert(key.bits(), count.get());
            Ok(())
        })
        .unwrap();
    (counts, digest)
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
