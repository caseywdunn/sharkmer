use serde_yaml_ng::Value;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

fn dna(mut state: u64, length: usize) -> String {
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[((state >> 32) & 3) as usize] as char
        })
        .collect()
}

fn input() -> (Vec<u8>, Vec<String>) {
    let mut reads = String::new();
    let mut primers = Vec::new();
    for (name, sequence) in [
        ("repeat", dna(157, 60) + &"A".repeat(40) + &dna(158, 60)),
        ("clean", dna(159, 160)),
    ] {
        let reverse = sequence.as_bytes()[sequence.len() - 15..]
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => unreachable!(),
            })
            .collect::<String>();
        primers.extend(["--pcr-primers".to_string(), format!("name={name},forward={},reverse={reverse},trim=15,mismatches=0,min-length=100,max-length=220,dedup-edit-threshold=0", &sequence[..15])]);
        for index in 0..4 {
            reads += &format!(
                "@{name}_{index}\n{sequence}\n+\n{}\n",
                "I".repeat(sequence.len())
            );
        }
    }
    // Force multiple sort runs while keeping the relevant target graphs small.
    for index in 0..1000 {
        reads += &format!(
            "@noise_{index}\n{}\n+\n{}\n",
            dna(1000 + index, 150),
            "I".repeat(150)
        );
    }
    (reads.into_bytes(), primers)
}

fn run(
    output: &Path,
    temporary: &Path,
    backend: &str,
    threads: usize,
    bytes: &[u8],
    primers: &[String],
    extra: &[&str],
) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
        .args([
            "-k",
            "19",
            "--chunks",
            "0",
            "-s",
            "external",
            "--count-backend",
            backend,
            "--count-fingerprint",
            "--count-sort-buffer-mib",
            "1",
            "--threads",
            &threads.to_string(),
            "--count-temp-dir",
        ])
        .arg(temporary)
        .arg("-o")
        .arg(output)
        .args(primers)
        .args(extra)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let written = child.stdin.take().unwrap().write_all(bytes);
    if let Err(error) = &written {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    let output = child.wait_with_output().unwrap();
    if output.status.success() {
        assert!(written.is_ok());
    }
    output
}

fn result(directory: &Path) -> (BTreeMap<String, Vec<u8>>, Value) {
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
    let stats =
        serde_yaml_ng::from_slice(&fs::read(directory.join("external.stats.yaml")).unwrap())
            .unwrap();
    (records, stats)
}

#[test]
fn one_shot_replay_threads_original_reads_and_matches_file_backed_legacy_pcr() {
    let (bytes, primers) = input();
    let temporary = tempfile::tempdir().unwrap();
    let input_path = temporary.path().join("source.fastq");
    fs::write(&input_path, &bytes).unwrap();
    let baseline = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
        .args([
            "-k",
            "19",
            "--chunks",
            "0",
            "-s",
            "external",
            "--count-fingerprint",
            "--count-backend",
            "legacy",
            "--read-threading",
            "--threads",
            "2",
        ])
        .arg("-o")
        .arg(baseline.path())
        .args(&primers)
        .arg(&input_path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (expected, stats) = result(baseline.path());
    assert!(!expected.is_empty());
    fs::remove_file(&input_path).unwrap();
    for threads in [1, 2, 4, 8] {
        for backend in [
            "memory",
            "external",
            "parallel",
            "packed",
            "adaptive",
            "adaptive-parallel",
            "pipeline",
            "hybrid",
        ] {
            let directory = tempfile::tempdir().unwrap();
            let output = run(
                directory.path(),
                temporary.path(),
                backend,
                threads,
                &bytes,
                &primers,
                &["--read-threading", "-v"],
            );
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(output.status.success(), "{backend}/{threads}: {stderr}");
            assert!(
                stderr.contains("Threading reads through assembly graph"),
                "{stderr}"
            );
            assert!(!stderr.contains("Read threading unavailable"));
            let (observed, actual) = result(directory.path());
            assert_eq!(observed, expected, "{backend}/{threads}");
            assert_eq!(
                actual["pcr_results"], stats["pcr_results"],
                "{backend}/{threads}"
            );
            assert_eq!(actual["count_fingerprint"], stats["count_fingerprint"]);
            assert_eq!(
                actual["read_replay"]["selected_records"],
                actual["n_reads_read"]
            );
            assert_eq!(
                actual["read_replay"]["sequence_bases"],
                actual["n_bases_read"]
            );
            assert_eq!(
                actual["read_replay"]["selected_input_sha256"],
                actual["count_store"]["selected_input_sha256"]
            );
            assert_eq!(actual["read_replay"]["format_version"].as_u64(), Some(1));
            assert!(actual["read_replay"]["file_bytes"].as_u64().unwrap() > bytes.len() as u64);
            assert!(
                actual["read_replay"]["peak_replay_buffer_bytes"]
                    .as_u64()
                    .unwrap()
                    <= 64 << 20
            );
            assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn replay_disk_exhaustion_aborts_publication_and_preserves_foreign_files() {
    let (bytes, primers) = input();
    let temporary = tempfile::tempdir().unwrap();
    let directory = tempfile::tempdir().unwrap();
    let sentinel = directory.path().join("foreign.txt");
    fs::write(&sentinel, "keep").unwrap();
    let output = run(
        directory.path(),
        temporary.path(),
        "memory",
        2,
        &bytes.repeat(4),
        &primers,
        &["--read-threading", "--count-disk-budget-mib", "1"],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Temporary disk space"));
    assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
    assert!(fs::read_dir(directory.path()).unwrap().all(|entry| {
        let name = entry.unwrap().file_name().to_string_lossy().to_string();
        !name.ends_with(".fasta") && !name.ends_with(".stats.yaml")
    }));
    let manifest: Value = serde_yaml_ng::from_slice(
        &fs::read(directory.path().join("external.manifest.yaml")).unwrap(),
    )
    .unwrap();
    assert_eq!(manifest["status"].as_str(), Some("failed"));
    assert!(manifest["files"].as_sequence().unwrap().is_empty());
    assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
}

#[test]
fn exact_backends_match_complete_one_shot_pcr_at_all_thread_counts() {
    let (bytes, primers) = input();
    let baseline = tempfile::tempdir().unwrap();
    let temporary = tempfile::tempdir().unwrap();
    let output = run(
        baseline.path(),
        temporary.path(),
        "memory",
        2,
        &bytes,
        &primers,
        &[],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (expected, stats) = result(baseline.path());
    assert!(
        expected
            .values()
            .any(|bytes| String::from_utf8_lossy(bytes).contains("assembly_status=warning"))
    );
    assert!(expected.keys().any(|name| name.contains("clean")));
    assert!(stats["count_store"].get("external").is_none());
    for threads in [1, 2, 4, 8] {
        for backend in ["external", "parallel", "packed"] {
            let directory = tempfile::tempdir().unwrap();
            let output = run(
                directory.path(),
                temporary.path(),
                backend,
                threads,
                &bytes,
                &primers,
                &[],
            );
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            let (observed, external) = result(directory.path());
            assert_eq!(observed, expected);
            for field in [
                "n_reads_read",
                "n_bases_read",
                "n_subreads_ingested",
                "n_bases_ingested",
                "n_kmers",
                "count_fingerprint",
                "pcr_results",
            ] {
                assert_eq!(external[field], stats[field], "{field}");
            }
            assert_eq!(
                external["count_store"]["selected_input_sha256"],
                stats["count_store"]["selected_input_sha256"]
            );
            if backend == "external" {
                assert_eq!(
                    external["count_store"]["backend"].as_str(),
                    Some("checksummed_sorted_file")
                );
                assert_eq!(
                    external["count_store"]["owned_table_bytes"].as_u64(),
                    Some(0)
                );
                assert_eq!(external["count_table_capacity"].as_u64(), Some(0));
                let measurements = &external["count_store"]["external"];
                assert_eq!(
                    measurements["peak_sort_buffer_bytes"].as_u64(),
                    Some(1 << 20)
                );
                assert_eq!(measurements["lookup_cache_bytes"].as_u64(), Some(4096));
                let file_bytes = measurements["count_file_bytes"].as_u64().unwrap();
                assert!(file_bytes > 1 << 20);
                assert!(
                    measurements["peak_temporary_disk_bytes"].as_u64().unwrap() >= 2 * file_bytes
                );
                assert!(measurements["sorted_run_count"].as_u64().unwrap() > 1);
                assert_eq!(
                    measurements["merge_count"].as_u64().unwrap() + 1,
                    measurements["sorted_run_count"].as_u64().unwrap()
                );
                assert!(
                    external["count_store"]["peak_count_input_reserved_bytes"]
                        .as_u64()
                        .unwrap()
                        < 2 << 20
                );
            } else if backend == "packed" {
                assert_eq!(
                    external["count_store"]["backend"].as_str(),
                    Some("packed_columns_memory")
                );
                let measurements = &external["count_store"]["packed"];
                assert_eq!(measurements["key_bytes_per_bucket"].as_u64(), Some(5));
                assert_eq!(measurements["count_bytes_per_bucket"].as_u64(), Some(8));
                assert_eq!(
                    external["count_store"]["owned_table_bytes"]
                        .as_u64()
                        .unwrap(),
                    measurements["buckets"].as_u64().unwrap() * 13 + 7
                );
                assert!(external["count_store"].get("parallel").is_none());
            } else if threads > 1 {
                assert_eq!(
                    external["count_store"]["backend"].as_str(),
                    Some("sharded_hashbrown_memory")
                );
                assert_eq!(
                    external["count_store"]["parallel"]["counting_workers"].as_u64(),
                    Some((threads - 1) as u64)
                );
                assert_eq!(
                    external["count_store"]["parallel"]["reusable_key_buffers"].as_u64(),
                    Some((3 * (threads - 1)) as u64)
                );
                assert!(external["count_store"].get("external").is_none());
                assert!(
                    external["count_table_capacity"].as_u64().unwrap()
                        >= external["count_store"]["distinct_keys"].as_u64().unwrap()
                );
            } else {
                assert_eq!(
                    external["count_store"]["backend"].as_str(),
                    Some("hashbrown_memory")
                );
                assert!(external["count_store"].get("parallel").is_none());
            }
            assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn exact_backend_failures_clean_resources_and_never_publish_partial_sequences() {
    let (bytes, primers) = input();
    let temporary = tempfile::tempdir().unwrap();
    let foreign = temporary.path().join("keep.txt");
    fs::write(&foreign, b"preserve me").unwrap();
    let mut invalid = bytes.clone();
    invalid.extend_from_slice(b"@invalid\nACGTX\n+\nIIIII\n");
    for backend in ["external", "parallel", "packed"] {
        for (data, extra, expected_error) in [
            (
                bytes.as_slice(),
                if backend == "external" {
                    vec!["--count-disk-budget-mib", "1"]
                } else {
                    vec!["--count-table-budget-mib", "1"]
                },
                if backend == "external" {
                    "Temporary disk space reservation requires"
                } else {
                    "Count-table memory reservation requires"
                },
            ),
            (invalid.as_slice(), vec![], "Invalid character"),
        ] {
            let directory = tempfile::tempdir().unwrap();
            let output = run(
                directory.path(),
                temporary.path(),
                backend,
                2,
                data,
                &primers,
                &extra,
            );
            assert!(!output.status.success());
            assert!(
                String::from_utf8_lossy(&output.stderr).contains(expected_error),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(!directory.path().join("external.stats.yaml").exists());
            assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
                entry
                    .unwrap()
                    .path()
                    .extension()
                    .is_some_and(|extension| extension == "fasta")
            }));
            assert_eq!(fs::read(&foreign).unwrap(), b"preserve me");
            assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 1);
        }
    }
}

#[test]
fn adaptive_spill_preserves_one_shot_counts_products_and_threading_evidence() {
    let (bytes, primers) = input();
    let temporary = tempfile::tempdir().unwrap();
    let baseline = tempfile::tempdir().unwrap();
    let output = run(
        baseline.path(),
        temporary.path(),
        "memory",
        2,
        &bytes,
        &primers,
        &["--read-threading"],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let (expected, reference) = result(baseline.path());
    for backend in ["adaptive", "adaptive-parallel", "pipeline", "hybrid"] {
        for threads in [1, 2, 4, 8] {
            let directory = tempfile::tempdir().unwrap();
            let output = run(
                directory.path(),
                temporary.path(),
                backend,
                threads,
                &bytes,
                &primers,
                &["--count-table-budget-mib", "1", "--read-threading"],
            );
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            let (observed, stats) = result(directory.path());
            assert_eq!(observed, expected);
            assert_eq!(stats["pcr_results"], reference["pcr_results"]);
            assert_eq!(stats["count_fingerprint"], reference["count_fingerprint"]);
            assert_eq!(
                stats["count_store"]["selected_input_sha256"],
                reference["count_store"]["selected_input_sha256"]
            );
            assert_eq!(
                stats["read_replay"]["selected_records"],
                reference["read_replay"]["selected_records"]
            );
            assert_eq!(
                stats["read_replay"]["file_bytes"],
                reference["read_replay"]["file_bytes"]
            );
            if backend == "adaptive"
                || threads == 1
                || (matches!(backend, "pipeline" | "hybrid") && threads == 2)
            {
                assert_eq!(
                    stats["count_store"]["backend"].as_str(),
                    Some("checksummed_sorted_file")
                );
                assert_eq!(stats["count_store"]["owned_table_bytes"].as_u64(), Some(0));
                assert!(
                    stats["count_store"]["external"]["sorted_run_count"]
                        .as_u64()
                        .unwrap()
                        > 1
                );
                assert!(stats["count_store"].get("adaptive_parallel").is_none());
            } else {
                assert_eq!(
                    stats["count_store"]["backend"].as_str(),
                    Some("adaptive_sharded_exact")
                );
                let owners = threads
                    - if matches!(backend, "pipeline" | "hybrid") {
                        2
                    } else {
                        1
                    };
                let metrics = &stats["count_store"]["adaptive_parallel"];
                let shards = &metrics["shards"];
                assert!(shards["disk_shards"].as_u64().unwrap() > 0);
                assert_eq!(
                    shards["disk_shards"].as_u64().unwrap()
                        + shards["memory_shards"].as_u64().unwrap(),
                    owners as u64
                );
                assert!(metrics["peak_sort_buffer_bytes"].as_u64().unwrap() <= 1 << 20);
                assert!(metrics["peak_temporary_disk_bytes"].as_u64().unwrap() > 0);
                assert_eq!(
                    stats["count_store"]["parallel"]["counting_workers"].as_u64(),
                    Some(owners as u64)
                );
            }
            assert!(stats["count_store"]["peak_table_bytes"].as_u64().unwrap() > 0);
            assert!(stats["count_store"]["peak_table_bytes"].as_u64().unwrap() <= 1 << 20);
            if matches!(backend, "pipeline" | "hybrid") && threads > 1 {
                let identity = &stats["count_store"]["input_identity_worker"];
                let reads = stats["n_reads_read"].as_u64().unwrap();
                assert_eq!(identity["selected_records"].as_u64(), Some(reads));
                assert_eq!(
                    identity["framed_bytes"].as_u64(),
                    Some(bytes.len() as u64 + 50 * reads)
                );
                assert_eq!(identity["owned_buffers"].as_u64(), Some(3));
                assert_eq!(identity["buffer_bytes"].as_u64(), Some(64 << 10));
                assert_eq!(identity["queued_buffers"].as_u64(), Some(1));
            } else {
                assert!(stats["count_store"].get("input_identity_worker").is_none());
            }
            assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn adaptive_failure_after_spill_cannot_publish_a_partial_input() {
    let (bytes, primers) = input();
    let mut invalid = bytes.clone();
    invalid.extend_from_slice(b"@invalid\nACGTX\n+\nIIIII\n");
    let temporary = tempfile::tempdir().unwrap();
    fs::write(temporary.path().join("foreign"), b"keep").unwrap();
    for backend in ["adaptive", "adaptive-parallel", "pipeline", "hybrid"] {
        for threads in [1, 2, 4, 8] {
            for (data, extra, error) in [
                (
                    bytes.as_slice(),
                    vec!["--count-disk-budget-mib", "1"],
                    "Temporary disk space reservation requires",
                ),
                (invalid.as_slice(), vec![], "Invalid character"),
            ] {
                let directory = tempfile::tempdir().unwrap();
                fs::write(directory.path().join("foreign"), b"keep output").unwrap();
                let mut options = vec!["--count-table-budget-mib", "1", "--read-threading"];
                options.extend(extra);
                let output = run(
                    directory.path(),
                    temporary.path(),
                    backend,
                    threads,
                    data,
                    &primers,
                    &options,
                );
                assert!(!output.status.success());
                assert!(
                    String::from_utf8_lossy(&output.stderr).contains(error),
                    "{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                assert!(fs::read_dir(directory.path()).unwrap().all(|entry| {
                    let name = entry.unwrap().file_name().to_string_lossy().to_string();
                    !name.ends_with(".fasta") && !name.ends_with(".stats.yaml")
                }));
                let manifest: Value = serde_yaml_ng::from_slice(
                    &fs::read(directory.path().join("external.manifest.yaml")).unwrap(),
                )
                .unwrap();
                assert_eq!(manifest["status"].as_str(), Some("failed"));
                assert!(manifest["files"].as_sequence().unwrap().is_empty());
                assert_eq!(
                    fs::read(directory.path().join("foreign")).unwrap(),
                    b"keep output"
                );
                assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 1);
                assert_eq!(fs::read(temporary.path().join("foreign")).unwrap(), b"keep");
            }
        }
    }
}

#[test]
fn corrupt_gzip_after_a_valid_member_joins_workers_and_aborts_publication() {
    use flate2::{Compression, write::GzEncoder};

    let (bytes, primers) = input();
    let compress = |data: &[u8]| {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(data).unwrap();
        encoder.finish().unwrap()
    };
    // The valid first member exercises count spill and original-pass replay
    // before a decoder failure in a later concatenated member.
    let prefix = compress(&bytes);
    let tail = compress(b"@tail\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n");
    let mut crc_failure = tail.clone();
    let trailer = crc_failure.len() - 8;
    crc_failure[trailer] ^= 0x80;
    let truncated = &tail[..tail.len() - 4];
    let source = tempfile::tempdir().unwrap();
    let path = source.path().join("input.bin");
    for backend in ["adaptive", "adaptive-parallel", "pipeline", "hybrid"] {
        for threads in [1, 2, 4, 8] {
            for broken in [crc_failure.as_slice(), truncated] {
                let mut compressed = prefix.clone();
                compressed.extend_from_slice(broken);
                fs::write(&path, &compressed).unwrap();
                let temporary = tempfile::tempdir().unwrap();
                let directory = tempfile::tempdir().unwrap();
                fs::write(temporary.path().join("foreign"), b"keep temporary").unwrap();
                fs::write(directory.path().join("foreign"), b"keep output").unwrap();
                let output = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
                    .args([
                        "-k",
                        "19",
                        "--chunks",
                        "0",
                        "-s",
                        "external",
                        "--count-backend",
                        backend,
                        "--threads",
                        &threads.to_string(),
                        "--count-table-budget-mib",
                        "1",
                        "--count-sort-buffer-mib",
                        "1",
                        "--read-threading",
                        "--count-temp-dir",
                    ])
                    .arg(temporary.path())
                    .arg("-o")
                    .arg(directory.path())
                    .args(&primers)
                    .arg(&path)
                    .output()
                    .unwrap();
                assert!(!output.status.success(), "{backend}/{threads}");
                assert!(
                    String::from_utf8_lossy(&output.stderr).contains("Error"),
                    "{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                let manifest: Value = serde_yaml_ng::from_slice(
                    &fs::read(directory.path().join("external.manifest.yaml")).unwrap(),
                )
                .unwrap();
                assert_eq!(manifest["status"].as_str(), Some("failed"));
                assert!(manifest["files"].as_sequence().unwrap().is_empty());
                assert!(fs::read_dir(directory.path()).unwrap().all(|entry| {
                    let name = entry.unwrap().file_name().to_string_lossy().to_string();
                    !name.ends_with(".fasta") && !name.ends_with(".stats.yaml")
                }));
                assert_eq!(
                    fs::read(directory.path().join("foreign")).unwrap(),
                    b"keep output"
                );
                assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 1);
                assert_eq!(
                    fs::read(temporary.path().join("foreign")).unwrap(),
                    b"keep temporary"
                );
                assert_eq!(fs::read(&path).unwrap(), compressed);
            }
        }
    }
}

#[cfg(target_os = "linux")]
#[test]
fn parallel_ingestion_does_not_keep_an_idle_assembly_pool() {
    use std::time::{Duration, Instant};
    for backend in ["adaptive-parallel", "pipeline", "hybrid"] {
        let directory = tempfile::tempdir().unwrap();
        let mut child = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
            .args([
                "-k",
                "19",
                "--threads",
                "4",
                "--count-backend",
                backend,
                "--count-table-budget-mib",
                "1",
                "--count-sort-buffer-mib",
                "1",
                "-s",
                "pool_probe",
                "-o",
            ])
            .arg(directory.path())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut input = child.stdin.take().unwrap();
        let deadline = Instant::now() + Duration::from_secs(10);
        let mut observed = None;
        while Instant::now() < deadline {
            let names: Vec<_> = fs::read_dir(format!("/proc/{}/task", child.id()))
                .unwrap()
                .filter_map(|entry| fs::read_to_string(entry.ok()?.path().join("comm")).ok())
                .collect();
            if names
                .iter()
                .filter(|name| name.starts_with("count-owner-"))
                .count()
                == if matches!(backend, "pipeline" | "hybrid") {
                    2
                } else {
                    3
                }
                && names
                    .iter()
                    .filter(|name| name.trim_end() == "input-identity")
                    .count()
                    == usize::from(matches!(backend, "pipeline" | "hybrid"))
            {
                observed = Some(names);
                break;
            }
            if child.try_wait().unwrap().is_some() {
                break;
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        if observed.is_none() {
            let _ = child.kill();
            let output = child.wait_with_output().unwrap();
            panic!(
                "Counting owners never became ready: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
        input
            .write_all(b"@read\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
            .unwrap();
        drop(input);
        let output = child.wait_with_output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let names = observed.unwrap();
        assert_eq!(
            names.len(),
            4,
            "unexpected threads while ingestion is waiting for stdin: {names:?}"
        );
    }
}

#[cfg(target_os = "linux")]
#[test]
fn pipeline_reuses_coordinator_for_assembly_within_requested_os_threads() {
    use std::time::{Duration, Instant};

    let (bytes, primers) = input();
    for threads in [1, 2, 4, 8] {
        let directory = tempfile::tempdir().unwrap();
        let temporary = tempfile::tempdir().unwrap();
        let mut child = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
            .args([
                "-k",
                "19",
                "--threads",
                &threads.to_string(),
                "--count-backend",
                "pipeline",
                "--count-table-budget-mib",
                "1",
                "--count-sort-buffer-mib",
                "1",
                "--read-threading",
                "-s",
                "external",
                "--count-temp-dir",
            ])
            .arg(temporary.path())
            .arg("-o")
            .arg(directory.path())
            .args(&primers)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut input = child.stdin.take().unwrap();
        let data = bytes.clone();
        let producer = std::thread::spawn(move || input.write_all(&data));
        let deadline = Instant::now() + Duration::from_secs(30);
        let mut observed_workers = std::collections::BTreeSet::new();
        let mut unexpected = None;
        let mut timed_out = false;
        loop {
            if let Ok(entries) = fs::read_dir(format!("/proc/{}/task", child.id())) {
                let names: Vec<_> = entries
                    .filter_map(|entry| fs::read_to_string(entry.ok()?.path().join("comm")).ok())
                    .collect();
                if names.len() > threads {
                    unexpected = Some(names.clone());
                }
                for name in names {
                    if name.starts_with("assembly-") {
                        observed_workers.insert(name.trim_end().to_string());
                    }
                }
            }
            if child.try_wait().unwrap().is_some() {
                break;
            }
            if Instant::now() > deadline {
                timed_out = true;
                let _ = child.kill();
                break;
            }
            std::thread::sleep(Duration::from_millis(1));
        }
        producer.join().unwrap().unwrap();
        let output = child.wait_with_output().unwrap();
        assert!(!timed_out, "pipeline/{threads} timed out");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            unexpected.is_none(),
            "pipeline/{threads} exceeded OS thread allowance: {unexpected:?}"
        );
        assert_eq!(
            observed_workers.len(),
            threads - 1,
            "pipeline/{threads}: {observed_workers:?}"
        );
        let (records, _) = result(directory.path());
        assert!(!records.is_empty());
        assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
    }
}

#[cfg(target_os = "linux")]
#[test]
fn hybrid_transfers_worker_role_before_eof_then_spills_and_preserves_complete_pcr() {
    use std::time::{Duration, Instant};

    let (bytes, primers) = input();
    let directory = tempfile::tempdir().unwrap();
    let temporary = tempfile::tempdir().unwrap();
    let reference = run(
        directory.path(),
        temporary.path(),
        "adaptive",
        2,
        &bytes,
        &primers,
        &["--count-table-budget-mib", "4", "--read-threading"],
    );
    assert!(
        reference.status.success(),
        "{}",
        String::from_utf8_lossy(&reference.stderr)
    );
    let (expected_records, expected_stats) = result(directory.path());
    let directory = tempfile::tempdir().unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_sharkmer"))
        .args([
            "-k",
            "19",
            "-s",
            "external",
            "--threads",
            "2",
            "--count-backend",
            "hybrid",
            "--count-table-budget-mib",
            "4",
            "--count-sort-buffer-mib",
            "1",
            "--count-owner-threshold-mib",
            "1",
            "--count-fingerprint",
            "--read-threading",
            "--count-temp-dir",
        ])
        .arg(temporary.path())
        .arg("-o")
        .arg(directory.path())
        .args(&primers)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    input.write_all(&bytes).unwrap();
    let began = Instant::now();
    let mut observed = None;
    while began.elapsed() < Duration::from_secs(10) {
        let names: Vec<_> = fs::read_dir(format!("/proc/{}/task", child.id()))
            .unwrap()
            .filter_map(|entry| fs::read_to_string(entry.ok()?.path().join("comm")).ok())
            .collect();
        if names.iter().any(|name| name.trim_end() == "count-owner-0") {
            observed = Some(names);
            break;
        }
        if child.try_wait().unwrap().is_some() {
            break;
        }
        std::thread::sleep(Duration::from_millis(2));
    }
    drop(input);
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let names = observed.expect("count owner did not replace identity worker before EOF");
    assert_eq!(names.len(), 2, "{names:?}");
    assert!(!names.iter().any(|name| name.trim_end() == "input-identity"));
    let (actual_records, actual_stats) = result(directory.path());
    assert_eq!(actual_records, expected_records);
    assert_eq!(actual_stats["pcr_results"], expected_stats["pcr_results"]);
    assert_eq!(
        actual_stats["count_fingerprint"],
        expected_stats["count_fingerprint"]
    );
    assert_eq!(
        actual_stats["count_store"]["selected_input_sha256"],
        expected_stats["count_store"]["selected_input_sha256"]
    );
    assert_eq!(
        actual_stats["read_replay"]["selected_input_sha256"],
        actual_stats["count_store"]["selected_input_sha256"]
    );
    let metadata = &actual_stats["count_store"];
    let transfer = &metadata["count_owner_transfer"];
    let prefix = transfer["selected_records_before_transfer"]
        .as_u64()
        .unwrap();
    assert!(prefix > 0 && prefix < actual_stats["n_reads_read"].as_u64().unwrap());
    assert_eq!(
        metadata["input_identity_worker"]["selected_records"].as_u64(),
        Some(prefix)
    );
    assert_eq!(transfer["threshold_bytes"].as_u64(), Some(1 << 20));
    assert!(transfer["current_table_bytes"].as_u64().unwrap() >= 1 << 20);
    assert_eq!(
        metadata["adaptive_parallel"]["shards"]["disk_shards"].as_u64(),
        Some(1)
    );
    assert_eq!(metadata["parallel"]["counting_workers"].as_u64(), Some(1));
    assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);

    let directory = tempfile::tempdir().unwrap();
    fs::write(directory.path().join("foreign"), b"preserve").unwrap();
    let mut invalid = bytes;
    invalid.extend_from_slice(b"@invalid\nACGTX\n+\nIIIII\n");
    let failed = run(
        directory.path(),
        temporary.path(),
        "hybrid",
        2,
        &invalid,
        &primers,
        &[
            "--count-table-budget-mib",
            "4",
            "--count-owner-threshold-mib",
            "1",
            "--read-threading",
        ],
    );
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("Invalid character"));
    assert_eq!(
        fs::read(directory.path().join("foreign")).unwrap(),
        b"preserve"
    );
    assert!(!directory.path().join("external.stats.yaml").exists());
    assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
        entry
            .unwrap()
            .path()
            .extension()
            .is_some_and(|ext| ext == "fasta")
    }));
    let manifest: Value = serde_yaml_ng::from_slice(
        &fs::read(directory.path().join("external.manifest.yaml")).unwrap(),
    )
    .unwrap();
    assert_eq!(manifest["status"].as_str(), Some("failed"));
    assert!(manifest["files"].as_sequence().unwrap().is_empty());
    assert_eq!(fs::read_dir(temporary.path()).unwrap().count(), 0);
}

#[cfg(target_os = "linux")]
struct FileLimitRun<'a> {
    limit: u64,
    inherited: usize,
    backend: &'a str,
    replay: bool,
}

#[cfg(target_os = "linux")]
fn run_with_file_limit(
    output: &Path,
    temporary: &Path,
    threads: usize,
    bytes: &[u8],
    primers: &[String],
    settings: FileLimitRun<'_>,
) -> (std::process::Output, Option<u64>) {
    use std::os::unix::process::CommandExt;
    use std::time::{Duration, Instant};

    let inherited = (0..settings.inherited)
        .map(|_| fs::File::open("/dev/null").unwrap())
        .collect::<Vec<_>>();
    let mut command = Command::new(env!("CARGO_BIN_EXE_sharkmer"));
    command
        .args([
            "-k",
            "19",
            "--chunks",
            "0",
            "-s",
            "external",
            "--count-backend",
            settings.backend,
            "--count-fingerprint",
            "--count-sort-buffer-mib",
            "1",
            "--count-table-budget-mib",
            "1",
            "--threads",
            &threads.to_string(),
            "--count-temp-dir",
        ])
        .arg(temporary)
        .arg("-o")
        .arg(output)
        .args(primers)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if settings.replay {
        command.arg("--read-threading");
    }
    // SAFETY: after fork this closure performs only async-signal-safe fcntl and
    // setrlimit syscalls. It neither allocates nor takes userspace locks, and
    // changes descriptor flags/limits only in the test child before exec.
    unsafe {
        command.pre_exec(move || {
            for file in &inherited {
                rustix::io::fcntl_setfd(file, rustix::io::FdFlags::empty())
                    .map_err(std::io::Error::from)?;
            }
            rustix::process::setrlimit(
                rustix::process::Resource::Nofile,
                rustix::process::Rlimit {
                    current: Some(settings.limit),
                    maximum: Some(settings.limit),
                },
            )
            .map_err(std::io::Error::from)
        });
    }
    let mut child = command.spawn().unwrap();
    let proc = Path::new("/proc").join(child.id().to_string());
    let started = Instant::now();
    let observed = loop {
        if child.try_wait().unwrap().is_some() {
            break None;
        }
        let waiting = fs::read_to_string(proc.join("wchan")).unwrap_or_default();
        if waiting.trim().ends_with("pipe_read") {
            // The child is waiting for its first input. Owned replay files can
            // already exist, along with the later publication lease. Exclude
            // both known owners to recover the earlier startup census, without
            // the child's transient census directory.
            let mut non_owned = 0;
            let mut output_lock_seen = false;
            for entry in fs::read_dir(proc.join("fd")).unwrap() {
                let target = fs::read_link(entry.unwrap().path()).unwrap();
                if target == output.join("external.lock") {
                    output_lock_seen = true;
                    continue;
                }
                if !target.starts_with(temporary) {
                    non_owned += 1;
                }
            }
            assert!(
                output_lock_seen,
                "publication lease must open after startup preflight"
            );
            break Some(non_owned);
        }
        if started.elapsed() > Duration::from_secs(10) {
            child.kill().unwrap();
            let failed = child.wait_with_output().unwrap();
            panic!(
                "Child did not reach its input: {}",
                String::from_utf8_lossy(&failed.stderr)
            );
        }
        std::thread::sleep(Duration::from_millis(2));
    };
    let written = child.stdin.take().unwrap().write_all(bytes);
    let result = child.wait_with_output().unwrap();
    if result.status.success() {
        assert!(written.is_ok());
        assert!(observed.is_some());
    }
    (result, observed)
}

#[cfg(target_os = "linux")]
#[test]
fn process_file_limit_accounts_for_inherited_files_and_preserves_bounded_pcr() {
    let directory = tempfile::tempdir().unwrap();
    let (bytes, primers) = input();
    let reference_output = directory.path().join("reference");
    let reference_temporary = directory.path().join("reference-tmp");
    fs::create_dir(&reference_temporary).unwrap();
    let reference = run(
        &reference_output,
        &reference_temporary,
        "legacy",
        2,
        &bytes,
        &primers,
        &["--read-threading"],
    );
    assert!(
        reference.status.success(),
        "{}",
        String::from_utf8_lossy(&reference.stderr)
    );
    let (expected_records, expected) = result(&reference_output);
    for (index, (threads, inherited)) in [(1, 0), (2, 0), (4, 0), (8, 0), (2, 30)]
        .into_iter()
        .enumerate()
    {
        let output = directory.path().join(format!("accepted-{index}"));
        let temporary = directory.path().join(format!("accepted-tmp-{index}"));
        fs::create_dir(&temporary).unwrap();
        let (child, census) = run_with_file_limit(
            &output,
            &temporary,
            threads,
            &bytes,
            &primers,
            FileLimitRun {
                limit: 64,
                inherited,
                backend: "hybrid",
                replay: true,
            },
        );
        assert!(
            child.status.success(),
            "{}",
            String::from_utf8_lossy(&child.stderr)
        );
        let (records, stats) = result(&output);
        assert_eq!(records, expected_records);
        assert_eq!(stats["pcr_results"], expected["pcr_results"]);
        assert_eq!(stats["count_fingerprint"], expected["count_fingerprint"]);
        let limits = &stats["count_store"]["process_file_limits"];
        assert_eq!(limits["requested_temporary_file_limit"].as_u64(), Some(128));
        assert_eq!(limits["process_soft_limit"].as_u64(), Some(64));
        assert_eq!(limits["process_open_at_start"].as_u64(), census);
        assert_eq!(limits["other_file_allowance"].as_u64(), Some(16));
        assert_eq!(limits["output_lease_allowance"].as_u64(), Some(1));
        assert_eq!(
            limits["temporary_file_limit"].as_u64(),
            Some(64 - census.unwrap() - 16 - 1)
        );
        assert!(fs::read_dir(&temporary).unwrap().next().is_none());
    }
    // An exhausted temporary allowance does not prohibit an in-memory run
    // that opens no count/replay files. It still describes the actual preflight.
    let output = directory.path().join("memory-only");
    let temporary = directory.path().join("memory-only-tmp");
    fs::create_dir(&temporary).unwrap();
    let small_input = bytes
        .split_inclusive(|byte| *byte == b'\n')
        .take(32)
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    let (child, census) = run_with_file_limit(
        &output,
        &temporary,
        1,
        &small_input,
        &primers,
        FileLimitRun {
            limit: 16,
            inherited: 0,
            backend: "memory",
            replay: false,
        },
    );
    assert!(
        child.status.success(),
        "{}",
        String::from_utf8_lossy(&child.stderr)
    );
    let (_, stats) = result(&output);
    assert_eq!(
        stats["count_store"]["process_file_limits"]["temporary_file_limit"].as_u64(),
        Some(0)
    );
    assert_eq!(
        stats["count_store"]["process_file_limits"]["process_open_at_start"].as_u64(),
        census
    );
    for (index, (limit, inherited)) in [(16, 0), (64, 45)].into_iter().enumerate() {
        let output = directory.path().join(format!("refused-{index}"));
        let temporary = directory.path().join(format!("refused-tmp-{index}"));
        fs::create_dir_all(&output).unwrap();
        fs::create_dir(&temporary).unwrap();
        fs::write(temporary.join("foreign"), b"preserve").unwrap();
        let (child, _) = run_with_file_limit(
            &output,
            &temporary,
            2,
            &bytes,
            &primers,
            FileLimitRun {
                limit,
                inherited,
                backend: "hybrid",
                replay: true,
            },
        );
        assert!(!child.status.success());
        let error = String::from_utf8_lossy(&child.stderr);
        assert!(
            error.contains("Count/replay files require 1 open descriptors; configured budget is 0"),
            "{error}"
        );
        assert!(!error.contains("Too many open files"), "{error}");
        assert_eq!(fs::read(temporary.join("foreign")).unwrap(), b"preserve");
        assert_eq!(fs::read_dir(&temporary).unwrap().count(), 1);
        assert!(!output.join("external.stats.yaml").exists());
        assert!(fs::read_dir(&output).unwrap().all(|entry| {
            entry
                .unwrap()
                .path()
                .extension()
                .is_none_or(|extension| extension != "fasta")
        }));
        let manifest: Value =
            serde_yaml_ng::from_slice(&fs::read(output.join("external.manifest.yaml")).unwrap())
                .unwrap();
        assert_eq!(manifest["status"].as_str(), Some("failed"));
    }
}

fn pcr_without_graph_buffer_measurements(value: &Value) -> Value {
    let mut value = value.clone();
    for target in value.as_sequence_mut().unwrap() {
        target
            .as_mapping_mut()
            .unwrap()
            .remove(Value::String("graph_buffers".to_owned()));
    }
    value
}

#[test]
fn graph_buffer_limits_bound_target_admission_and_preserve_exact_products() {
    let directory = tempfile::tempdir().unwrap();
    let temporary = directory.path().join("temporary");
    fs::create_dir(&temporary).unwrap();
    let (bytes, primers) = input();
    let baseline = directory.path().join("baseline");
    let reference = run(
        &baseline,
        &temporary,
        "hybrid",
        2,
        &bytes,
        &primers,
        &["--read-threading"],
    );
    assert!(
        reference.status.success(),
        "{}",
        String::from_utf8_lossy(&reference.stderr)
    );
    let (expected, expected_stats) = result(&baseline);
    assert!(!expected.is_empty());
    assert!(expected_stats.get("assembly_graph_buffers").is_none());
    for threads in [1, 2, 4, 8] {
        for backend in [
            "memory",
            "external",
            "parallel",
            "packed",
            "adaptive",
            "adaptive-parallel",
            "pipeline",
            "hybrid",
        ] {
            let output = directory.path().join(format!("{backend}-{threads}"));
            let child = run(
                &output,
                &temporary,
                backend,
                threads,
                &bytes,
                &primers,
                &[
                    "--read-threading",
                    "--graph-buffer-budget-mib",
                    "2",
                    "--graph-buffer-target-mib",
                    "1",
                ],
            );
            assert!(
                child.status.success(),
                "{backend}/{threads}: {}",
                String::from_utf8_lossy(&child.stderr)
            );
            let (records, stats) = result(&output);
            assert_eq!(records, expected);
            assert_eq!(
                stats["count_fingerprint"],
                expected_stats["count_fingerprint"]
            );
            assert_eq!(
                pcr_without_graph_buffer_measurements(&stats["pcr_results"]),
                expected_stats["pcr_results"]
            );
            assert_eq!(
                stats["read_replay"]["selected_input_sha256"],
                stats["count_store"]["selected_input_sha256"]
            );
            let summary = &stats["assembly_graph_buffers"];
            assert_eq!(summary["shared_limit_bytes"].as_u64(), Some(2 << 20));
            assert_eq!(summary["per_target_limit_bytes"].as_u64(), Some(1 << 20));
            assert_eq!(summary["concurrency_limit"].as_u64(), Some(2));
            assert!(summary["peak_active_targets"].as_u64().unwrap() <= threads.min(2) as u64);
            assert_eq!(
                summary["peak_reserved_target_allowances_bytes"]
                    .as_u64()
                    .unwrap(),
                summary["peak_active_targets"].as_u64().unwrap() << 20
            );
            assert!(
                stats["pcr_results"]
                    .as_sequence()
                    .unwrap()
                    .iter()
                    .any(
                        |target| target["graph_buffers"]["peak_reserved_buffer_bytes"]
                            .as_u64()
                            .unwrap()
                            > 0
                    )
            );
            for target in stats["pcr_results"].as_sequence().unwrap() {
                assert_eq!(
                    target["graph_buffers"]["limit_bytes"].as_u64(),
                    Some(1 << 20)
                );
                assert!(
                    target["graph_buffers"]["peak_reserved_buffer_bytes"]
                        .as_u64()
                        .unwrap()
                        <= 1 << 20
                );
            }
            assert!(fs::read_dir(&temporary).unwrap().next().is_none());
        }
    }
    let output = directory.path().join("serial-admission");
    let child = run(
        &output,
        &temporary,
        "hybrid",
        8,
        &bytes,
        &primers,
        &[
            "--read-threading",
            "--graph-buffer-budget-mib",
            "1",
            "--graph-buffer-target-mib",
            "1",
        ],
    );
    assert!(
        child.status.success(),
        "{}",
        String::from_utf8_lossy(&child.stderr)
    );
    let (records, stats) = result(&output);
    assert_eq!(records, expected);
    assert_eq!(
        stats["assembly_graph_buffers"]["peak_active_targets"].as_u64(),
        Some(1)
    );
    assert_eq!(
        stats["assembly_graph_buffers"]["concurrency_limit"].as_u64(),
        Some(1)
    );
}

#[test]
fn graph_buffer_exhaustion_fails_without_partial_publication_or_temporary_leaks() {
    let directory = tempfile::tempdir().unwrap();
    let temporary = directory.path().join("temporary");
    let output = directory.path().join("outputs");
    fs::create_dir(&temporary).unwrap();
    fs::create_dir(&output).unwrap();
    fs::write(temporary.join("foreign"), b"preserve").unwrap();
    fs::write(output.join("foreign"), b"preserve").unwrap();
    let (bytes, primers) = long_graph_input();
    let reference = directory.path().join("roomy-reference");
    let successful = run(
        &reference,
        &temporary,
        "hybrid",
        4,
        &bytes,
        &primers,
        &[
            "--graph-buffer-budget-mib",
            "8",
            "--graph-buffer-target-mib",
            "8",
        ],
    );
    assert!(
        successful.status.success(),
        "{}",
        String::from_utf8_lossy(&successful.stderr)
    );
    assert!(!result(&reference).0.is_empty());
    let child = run(
        &output,
        &temporary,
        "hybrid",
        4,
        &bytes,
        &primers,
        &[
            "--graph-buffer-budget-mib",
            "1",
            "--graph-buffer-target-mib",
            "1",
        ],
    );
    assert!(!child.status.success());
    let error = String::from_utf8_lossy(&child.stderr);
    assert!(error.contains("Target-graph memory reservation"), "{error}");
    assert_eq!(fs::read(temporary.join("foreign")).unwrap(), b"preserve");
    assert_eq!(fs::read(output.join("foreign")).unwrap(), b"preserve");
    assert_eq!(fs::read_dir(&temporary).unwrap().count(), 1);
    assert!(!output.join("external.stats.yaml").exists());
    assert!(fs::read_dir(&output).unwrap().all(|entry| {
        entry
            .unwrap()
            .path()
            .extension()
            .is_none_or(|extension| extension != "fasta")
    }));
    let manifest: Value =
        serde_yaml_ng::from_slice(&fs::read(output.join("external.manifest.yaml")).unwrap())
            .unwrap();
    assert_eq!(manifest["status"].as_str(), Some("failed"));
}

fn long_graph_input() -> (Vec<u8>, Vec<String>) {
    let sequence = dna(91234, 25_000);
    let reverse = sequence.as_bytes()[sequence.len() - 15..]
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect::<String>();
    let primers = vec![
        "--pcr-primers".to_string(),
        format!(
            "name=long,forward={},reverse={reverse},trim=15,mismatches=0,min-length=100,max-length=30000",
            &sequence[..15]
        ),
    ];
    let mut bytes = Vec::new();
    for index in 0..3 {
        bytes.extend(
            format!(
                "@long_{index}\n{sequence}\n+\n{}\n",
                "I".repeat(sequence.len())
            )
            .bytes(),
        );
    }
    (bytes, primers)
}

#[test]
fn a_graph_budget_failure_preserves_the_previous_completed_generation() {
    let directory = tempfile::tempdir().unwrap();
    let temporary = directory.path().join("temporary");
    let output = directory.path().join("outputs");
    fs::create_dir(&temporary).unwrap();
    let (bytes, primers) = input();
    let successful = run(&output, &temporary, "hybrid", 2, &bytes, &primers, &[]);
    assert!(
        successful.status.success(),
        "{}",
        String::from_utf8_lossy(&successful.stderr)
    );
    let before: BTreeMap<_, _> = fs::read_dir(&output)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "fasta" || extension == "yaml")
        })
        .map(|path| (path.clone(), fs::read(path).unwrap()))
        .collect();
    let (long_bytes, long_primers) = long_graph_input();
    let failed = run(
        &output,
        &temporary,
        "hybrid",
        4,
        &long_bytes,
        &long_primers,
        &[
            "--graph-buffer-budget-mib",
            "1",
            "--graph-buffer-target-mib",
            "1",
        ],
    );
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("Target-graph memory reservation"));
    for (path, bytes) in before {
        assert_eq!(fs::read(path).unwrap(), bytes);
    }
    let manifest: Value =
        serde_yaml_ng::from_slice(&fs::read(output.join("external.manifest.yaml")).unwrap())
            .unwrap();
    assert_eq!(manifest["status"].as_str(), Some("complete"));
    assert!(!output.join("external_long.fasta").exists());
    assert!(fs::read_dir(&temporary).unwrap().next().is_none());
}

#[test]
fn an_empty_panel_closes_unused_replay_before_committing_stats() {
    let output = tempfile::tempdir().unwrap();
    let temporary = tempfile::tempdir().unwrap();
    let panel = temporary.path().join("empty-panel.yaml");
    fs::write(
        &panel,
        "name: empty\ndescription: Empty target set\nprimers: []\n",
    )
    .unwrap();
    let panel_string = panel.to_str().unwrap();
    let result = run(
        output.path(),
        temporary.path(),
        "external",
        2,
        b"@read\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
        &[],
        &[
            "--read-threading",
            "--pcr-panel-file",
            panel_string,
            "--graph-buffer-budget-mib",
            "2",
            "--graph-buffer-target-mib",
            "1",
        ],
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let stats: Value =
        serde_yaml_ng::from_slice(&fs::read(output.path().join("external.stats.yaml")).unwrap())
            .unwrap();
    assert_eq!(stats["run_status"].as_str(), Some("complete"));
    assert_eq!(stats["n_reads_read"].as_u64(), Some(1));
    assert!(stats["pcr_results"].is_null());
    let leftovers: Vec<_> = fs::read_dir(temporary.path())
        .unwrap()
        .map(|entry| entry.unwrap().file_name())
        .collect();
    assert_eq!(leftovers, [panel.file_name().unwrap()]);
    assert!(
        !fs::read_dir(output.path()).unwrap().any(|entry| entry
            .unwrap()
            .file_type()
            .unwrap()
            .is_dir())
    );
}
