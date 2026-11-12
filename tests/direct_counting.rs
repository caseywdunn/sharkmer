use serde_yaml_ng::Value;
use std::collections::BTreeMap;
use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::process::{Command, Stdio};

fn run(
    directory: &Path,
    threads: usize,
    backend: &str,
    stdin: Option<&[u8]>,
) -> std::process::Output {
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ERR571460_100k_R1.fastq.gz");
    let mut command = Command::new(env!("CARGO_BIN_EXE_sharkmer"));
    command
        .args([
            "-k",
            "31",
            "--chunks",
            "0",
            "-s",
            "counter",
            "--pcr-panel",
            "cnidaria",
            "--count-backend",
            backend,
            "--threads",
            &threads.to_string(),
            "-o",
        ])
        .arg(directory)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if backend == "legacy" || (threads == 2 && stdin.is_none()) {
        command.arg("--count-fingerprint");
    }
    if let Some(bytes) = stdin {
        let mut child = command.stdin(Stdio::piped()).spawn().unwrap();
        child.stdin.take().unwrap().write_all(bytes).unwrap();
        child.wait_with_output().unwrap()
    } else {
        command.arg(input).output().unwrap()
    }
}

fn outputs(directory: &Path) -> (BTreeMap<String, Vec<u8>>, Value) {
    let records = fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "fasta")
        })
        .map(|path| {
            (
                path.file_name().unwrap().to_string_lossy().to_string(),
                fs::read(path).unwrap(),
            )
        })
        .collect();
    let stats = serde_yaml_ng::from_slice(&fs::read(directory.join("counter.stats.yaml")).unwrap())
        .unwrap();
    (records, stats)
}

#[test]
fn exact_counter_matches_complete_fixture_at_all_supported_threads_and_one_shot_input() {
    let baseline = tempfile::tempdir().unwrap();
    let output = run(baseline.path(), 2, "legacy", None);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (records, stats) = outputs(baseline.path());
    assert!(!records.is_empty());
    let mut digest = None;
    for threads in [1, 2, 4, 8] {
        let directory = tempfile::tempdir().unwrap();
        let output = run(directory.path(), threads, "memory", None);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let (exact_records, exact_stats) = outputs(directory.path());
        assert_eq!(exact_records, records);
        if threads == 2 {
            assert_eq!(exact_stats["count_fingerprint"], stats["count_fingerprint"]);
        }
        for field in [
            "n_reads_read",
            "n_bases_read",
            "n_subreads_ingested",
            "n_bases_ingested",
            "n_kmers",
            "pcr_results",
        ] {
            assert_eq!(exact_stats[field], stats[field], "{field}");
        }
        assert_eq!(
            exact_stats["count_store"]["count_width_bits"].as_u64(),
            Some(64)
        );
        assert_eq!(
            exact_stats["count_store"]["occurrences"],
            exact_stats["n_kmers"]
        );
        let observed = exact_stats["count_store"]["selected_input_sha256"]
            .as_str()
            .unwrap()
            .to_string();
        if let Some(previous) = &digest {
            assert_eq!(&observed, previous);
        }
        digest = Some(observed);
        assert!(!directory.path().join("counter.histo").exists());
    }
    let fixture =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ERR571460_100k_R1.fastq.gz");
    let mut data = Vec::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(fixture).unwrap())
        .read_to_end(&mut data)
        .unwrap();
    let directory = tempfile::tempdir().unwrap();
    let output = run(directory.path(), 2, "memory", Some(&data));
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (stdin_records, stdin_stats) = outputs(directory.path());
    assert_eq!(stdin_records, records);
    assert_eq!(stdin_stats["n_kmers"], stats["n_kmers"]);
    assert_ne!(
        stdin_stats["count_store"]["selected_input_sha256"]
            .as_str()
            .unwrap(),
        digest.unwrap()
    );
}

#[test]
fn exact_counter_limits_fail_without_publishing_partial_sequences() {
    let input =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ERR571460_100k_R1.fastq.gz");
    for extra in [
        ["--count-slot-budget-mib", "1"],
        ["--fastq-record-bytes", "16"],
    ] {
        let directory = tempfile::tempdir().unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
            .args([
                "-k",
                "31",
                "--count-backend",
                "memory",
                "--pcr-panel",
                "cnidaria",
                "-s",
                "limits",
                "-o",
            ])
            .arg(directory.path())
            .args(extra)
            .arg(&input)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(!directory.path().join("limits.stats.yaml").exists());
        assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
            entry
                .unwrap()
                .path()
                .extension()
                .is_some_and(|extension| extension == "fasta")
        }));
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(
            error.contains("budget") || error.contains("byte buffer limit"),
            "{error}"
        );
    }
}
