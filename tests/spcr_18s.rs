use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn sharkmer_bin() -> String {
    env!("CARGO_BIN_EXE_sharkmer").to_string()
}

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("ERR571460_100k_R1.fastq.gz")
}

/// Integration test: verify that sharkmer recovers 18S from ERR571460 (Porites lutea)
/// using the cnidaria primer panel with 100k reads.
#[test]
fn test_18s_recovery_from_err571460() {
    let fixture = fixture_path();
    assert!(
        fixture.exists(),
        "Test fixture not found: {}",
        fixture.display()
    );

    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "ERR571460_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "31",
            "--pcr-panel",
            "cnidaria",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "100000",
            "--chunks",
            "1",
            "-v",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer exited with error: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // With reverse extension, forward and reverse graph extension
    // converge to produce complete amplicons. Verify the 18S FASTA is produced.
    let fasta_18s = outdir.path().join(format!("{}_cnidaria_18S.fasta", sample));
    assert!(
        fasta_18s.exists(),
        "Expected 18S FASTA to be produced with reverse extension"
    );
}

/// Verify that stats.yaml is produced and contains expected fields.
#[test]
fn test_stats_yaml_output() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "stats_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "--pcr-panel",
            "cnidaria",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "10000",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stats_path = outdir.path().join(format!("{}.stats.yaml", sample));
    assert!(stats_path.exists(), "stats.yaml not found");

    let stats_content = fs::read_to_string(&stats_path).expect("failed to read stats.yaml");

    // Check required fields are present
    assert!(stats_content.contains("sharkmer_version:"));
    assert!(stats_content.contains("sample:"));
    assert!(stats_content.contains("kmer_length:"));
    assert!(stats_content.contains("n_reads_read:"));
    assert!(stats_content.contains("n_bases_read:"));
    assert!(stats_content.contains("n_kmers:"));
    assert!(stats_content.contains("peak_memory_bytes:"));
    assert!(stats_content.contains("pcr_results:"));
}

/// Verify that incremental histograms are produced when --chunks > 0.
#[test]
fn test_histogram_output_with_chunks() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "histo_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "10000",
            "--chunks",
            "2",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let histo_path = outdir.path().join(format!("{}.histo", sample));
    let final_histo_path = outdir.path().join(format!("{}.final.histo", sample));

    assert!(histo_path.exists(), "Incremental histogram not found");
    assert!(final_histo_path.exists(), "Final histogram not found");

    // Histogram should have content
    let histo_content = fs::read_to_string(&histo_path).expect("failed to read histogram");
    assert!(!histo_content.is_empty(), "Histogram file is empty");
}

/// Verify that no histograms are produced when --chunks 0 (default).
#[test]
fn test_no_histogram_with_chunks_zero() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "nohisto_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "10000",
            "--chunks",
            "0",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let histo_path = outdir.path().join(format!("{}.histo", sample));
    assert!(
        !histo_path.exists(),
        "Histogram should not be produced with --chunks 0"
    );

    // Stats file should still be produced
    let stats_path = outdir.path().join(format!("{}.stats.yaml", sample));
    assert!(stats_path.exists(), "stats.yaml should always be produced");
}

/// Verify --dry-run exits successfully without producing output files.
#[test]
fn test_dry_run_no_output() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "dryrun_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "--pcr-panel",
            "cnidaria",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--dry-run",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer --dry-run failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // No output files should be created
    let stats_path = outdir.path().join(format!("{}.stats.yaml", sample));
    assert!(
        !stats_path.exists(),
        "stats.yaml should not be produced in dry-run mode"
    );
}

/// Verify --list-panels exits successfully and lists known panels.
#[test]
fn test_list_panels() {
    let output = Command::new(sharkmer_bin())
        .args(["--list-panels"])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer --list-panels failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("cnidaria"), "Should list cnidaria panel");
    assert!(stdout.contains("human"), "Should list human panel");
    assert!(stdout.contains("bacteria"), "Should list bacteria panel");
    assert!(stdout.contains("teleostei"), "Should list teleostei panel");
}

/// Verify --cite exits successfully and outputs citation info.
#[test]
fn test_cite_flag() {
    let output = Command::new(sharkmer_bin())
        .args(["--cite"])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer --cite failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(!stdout.is_empty(), "Citation output should not be empty");
}

/// Verify sharkmer errors when --sample is missing and no --ena to derive from.
#[test]
fn test_error_missing_sample() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-o",
            outdir.path().to_str().unwrap(),
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        !output.status.success(),
        "sharkmer should fail without --sample"
    );
}

/// Verify sharkmer errors on nonexistent input file.
#[test]
fn test_error_nonexistent_input() {
    let outdir = tempfile::tempdir().expect("failed to create temp dir");

    let output = Command::new(sharkmer_bin())
        .args([
            "-s",
            "test_sample",
            "-o",
            outdir.path().to_str().unwrap(),
            "/nonexistent/file.fastq",
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        !output.status.success(),
        "sharkmer should fail on nonexistent input"
    );
}

/// Verify gene names in output files are prefixed with panel name.
#[test]
fn test_gene_name_panel_prefix() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "prefix_test";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "31",
            "--pcr-panel",
            "cnidaria",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "100000",
            "-v",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // With reverse extension, amplicons are produced. Verify the run
    // succeeds and output files use panel-prefixed gene names.
    let stats_path = outdir.path().join(format!("{}.stats.yaml", sample));
    let stats_content = fs::read_to_string(&stats_path).expect("failed to read stats.yaml");
    assert!(
        stats_content.contains("cnidaria_18S"),
        "stats.yaml should contain panel-prefixed gene name cnidaria_18S"
    );

    // FASTA output should exist and use panel-prefixed gene name
    let fasta_18s = outdir.path().join(format!("{}_cnidaria_18S.fasta", sample));
    assert!(
        fasta_18s.exists(),
        "Expected cnidaria_18S FASTA with panel-prefixed gene name"
    );
}

/// Verify kmer histogram correctness against jellyfish-verified reference values.
/// Reference values were validated by running both sharkmer and jellyfish on the
/// fixture file and confirming identical output (see scripts/compare_jellyfish.sh).
#[test]
fn test_histogram_correctness() {
    let fixture = fixture_path();
    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "histo_correct";

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--chunks",
            "1",
            "--histo-max",
            "10000",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Parse the final histogram into a count -> frequency map
    let final_histo_path = outdir.path().join(format!("{}.final.histo", sample));
    assert!(final_histo_path.exists(), "Final histogram not found");
    let content = fs::read_to_string(&final_histo_path).expect("failed to read final histogram");

    let mut histo: std::collections::HashMap<u64, u64> = std::collections::HashMap::new();
    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("count") {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            let count: u64 = parts[0].parse().expect("invalid count");
            let freq: u64 = parts[1].parse().expect("invalid frequency");
            histo.insert(count, freq);
        }
    }

    // Spot-check values verified against jellyfish output
    assert_eq!(histo[&1], 10146199, "singleton count mismatch");
    assert_eq!(histo[&2], 389066, "count-2 frequency mismatch");
    assert_eq!(histo[&10], 4010, "count-10 frequency mismatch");
    assert_eq!(histo[&100], 7, "count-100 frequency mismatch");
    assert_eq!(histo[&1000], 0, "count-1000 frequency mismatch");
    assert_eq!(histo[&10001], 13, "overflow bucket mismatch");

    // Verify stats match
    let stats_path = outdir.path().join(format!("{}.stats.yaml", sample));
    let stats = fs::read_to_string(&stats_path).expect("failed to read stats.yaml");
    assert!(
        stats.contains("n_kmers: 12997261"),
        "total kmer count mismatch"
    );
    assert!(
        stats.contains("n_singleton_kmers: 10146199"),
        "singleton kmer count mismatch"
    );
}

/// Verify incremental counting produces the same final histogram regardless
/// of chunk count. Runs with 20 chunks and checks that the final histogram
/// matches the single-chunk reference.
#[test]
fn test_incremental_histogram_consistency() {
    let fixture = fixture_path();
    let outdir_1 = tempfile::tempdir().expect("failed to create temp dir");
    let outdir_20 = tempfile::tempdir().expect("failed to create temp dir");

    // Run with 1 chunk
    let output_1 = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            "chunk1",
            "-o",
            outdir_1.path().to_str().unwrap(),
            "--chunks",
            "1",
            "--histo-max",
            "10000",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output_1.status.success(),
        "sharkmer (1 chunk) failed: {}",
        String::from_utf8_lossy(&output_1.stderr)
    );

    // Run with 20 chunks
    let output_20 = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            "chunk20",
            "-o",
            outdir_20.path().to_str().unwrap(),
            "--chunks",
            "20",
            "--histo-max",
            "10000",
            fixture.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output_20.status.success(),
        "sharkmer (20 chunks) failed: {}",
        String::from_utf8_lossy(&output_20.stderr)
    );

    // Parse both final histograms and compare
    let parse_histo = |path: std::path::PathBuf| -> Vec<(u64, u64)> {
        let content = fs::read_to_string(&path)
            .unwrap_or_else(|_| panic!("failed to read {}", path.display()));
        let mut entries = Vec::new();
        for line in content.lines() {
            if line.starts_with('#') || line.starts_with("count") {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                let count: u64 = parts[0].parse().expect("invalid count");
                let freq: u64 = parts[1].parse().expect("invalid frequency");
                entries.push((count, freq));
            }
        }
        entries
    };

    let histo_1 = parse_histo(outdir_1.path().join("chunk1.final.histo"));
    let histo_20 = parse_histo(outdir_20.path().join("chunk20.final.histo"));

    assert_eq!(
        histo_1.len(),
        histo_20.len(),
        "Histogram row counts differ: 1-chunk={} vs 20-chunk={}",
        histo_1.len(),
        histo_20.len()
    );

    for (a, b) in histo_1.iter().zip(histo_20.iter()) {
        assert_eq!(
            a, b,
            "Histogram mismatch at count {}: 1-chunk freq={} vs 20-chunk freq={}",
            a.0, a.1, b.1
        );
    }
}

/// Verify --export-panel outputs valid YAML for a known panel.
#[test]
fn test_export_panel() {
    let output = Command::new(sharkmer_bin())
        .args(["--export-panel", "cnidaria"])
        .output()
        .expect("failed to run sharkmer");

    assert!(
        output.status.success(),
        "sharkmer --export-panel failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("name: cnidaria"),
        "Should contain panel name"
    );
    assert!(
        stdout.contains("primers:"),
        "Should contain primers section"
    );
    assert!(
        stdout.contains("forward_seq:"),
        "Should contain primer sequences"
    );
}

// --- Malformed FASTQ tests ---

/// Helper to write a temporary FASTQ file and run sharkmer on it,
/// returning the Command output.
fn run_with_fastq_content(content: &str) -> std::process::Output {
    let tmpdir = tempfile::tempdir().expect("failed to create temp dir");
    let fastq_path = tmpdir.path().join("test.fastq");
    let mut f = fs::File::create(&fastq_path).expect("failed to create temp fastq");
    f.write_all(content.as_bytes())
        .expect("failed to write temp fastq");

    let outdir = tempfile::tempdir().expect("failed to create temp dir");

    let output = Command::new(sharkmer_bin())
        .args([
            "-k",
            "21",
            "-s",
            "malformed_test",
            "-o",
            outdir.path().to_str().unwrap(),
            "--chunks",
            "0",
            fastq_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");
    output
}

/// Verify sharkmer rejects FASTA input with a clear error message.
/// Uses 4 lines so the reader doesn't hit truncation before validation.
#[test]
fn test_error_fasta_input() {
    let content = ">seq1\nACGTACGTACGTACGTACGTACGT\n>seq2\nACGTACGTACGTACGTACGTACGT\n";
    let output = run_with_fastq_content(content);

    assert!(
        !output.status.success(),
        "sharkmer should fail on FASTA input"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("FASTA") || stderr.contains("fasta"),
        "Error should mention FASTA format, got: {}",
        stderr
    );
}

/// Verify sharkmer rejects FASTQ with mismatched sequence/quality lengths.
#[test]
fn test_error_mismatched_quality_length() {
    let content = "@read1\nACGTACGTACGTACGTACGTACGT\n+\nIIII\n";
    let output = run_with_fastq_content(content);

    assert!(
        !output.status.success(),
        "sharkmer should fail on mismatched seq/qual lengths"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("mismatched") || stderr.contains("length"),
        "Error should mention length mismatch, got: {}",
        stderr
    );
}

/// Verify sharkmer rejects FASTQ with invalid header (no @ prefix).
#[test]
fn test_error_invalid_fastq_header() {
    let content = "BAD_HEADER\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n";
    let output = run_with_fastq_content(content);

    assert!(
        !output.status.success(),
        "sharkmer should fail on invalid FASTQ header"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("header") || stderr.contains("'@'"),
        "Error should mention invalid header, got: {}",
        stderr
    );
}

/// Verify sharkmer rejects an empty input file.
#[test]
fn test_error_empty_input() {
    let output = run_with_fastq_content("");

    assert!(
        !output.status.success(),
        "sharkmer should fail on empty input"
    );
}

// --- Panel validation tests ---

/// Helper to write a temporary YAML panel file and run --validate-panels on it.
fn run_validate_panel(yaml: &str) -> std::process::Output {
    let tmpdir = tempfile::tempdir().expect("failed to create temp dir");
    let panel_path = tmpdir.path().join("test_panel.yaml");
    let mut f = fs::File::create(&panel_path).expect("failed to create temp panel");
    f.write_all(yaml.as_bytes())
        .expect("failed to write temp panel");

    let output = Command::new(sharkmer_bin())
        .args([
            "--validate-panels",
            "--pcr-panel-file",
            panel_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run sharkmer");
    output
}

/// Verify --validate-panels accepts a well-formed user panel.
#[test]
fn test_validate_panel_valid() {
    let yaml = "\
name: test_panel
description: A test panel
primers:
  - gene_name: test_gene
    forward_seq: ACGTACGTACGT
    reverse_seq: TGCATGCATGCA
    min_length: 100
    max_length: 500
    min_count: 2
    mismatches: 2
    trim: 15
";
    let output = run_validate_panel(yaml);

    assert!(
        output.status.success(),
        "Valid panel should pass validation: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("All primers valid"));
}

/// Verify --validate-panels rejects a primer with invalid nucleotide characters.
#[test]
fn test_validate_panel_invalid_nucleotide() {
    let yaml = "\
name: bad_panel
description: Panel with invalid bases
primers:
  - gene_name: bad_gene
    forward_seq: ACGTXZQP
    reverse_seq: TGCATGCA
    min_length: 100
    max_length: 500
    min_count: 2
    mismatches: 2
    trim: 15
";
    let output = run_validate_panel(yaml);

    assert!(
        !output.status.success(),
        "Panel with invalid nucleotides should fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Invalid nucleotide"),
        "Error should mention invalid nucleotide, got: {}",
        stderr
    );
}

/// Verify --validate-panels rejects min_length > max_length.
#[test]
fn test_validate_panel_min_exceeds_max() {
    let yaml = "\
name: bad_panel
description: Panel with inverted lengths
primers:
  - gene_name: bad_gene
    forward_seq: ACGTACGTACGT
    reverse_seq: TGCATGCATGCA
    min_length: 1000
    max_length: 100
    min_count: 2
    mismatches: 2
    trim: 15
";
    let output = run_validate_panel(yaml);

    assert!(
        !output.status.success(),
        "Panel with min > max length should fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("min-length") || stderr.contains("greater than max"),
        "Error should mention length constraint, got: {}",
        stderr
    );
}

/// Verify --pcr-panel-file rejects malformed YAML.
#[test]
fn test_validate_panel_malformed_yaml() {
    let yaml = "this is not: valid: yaml: [[[";
    let output = run_validate_panel(yaml);

    assert!(!output.status.success(), "Malformed YAML should fail");
}

/// Verify --pcr-panel-file rejects YAML missing required fields.
#[test]
fn test_validate_panel_missing_fields() {
    let yaml = "\
name: incomplete_panel
description: Missing primer fields
primers:
  - gene_name: no_seqs
";
    let output = run_validate_panel(yaml);

    assert!(
        !output.status.success(),
        "Panel missing required fields should fail"
    );
}
