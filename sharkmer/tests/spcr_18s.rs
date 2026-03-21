use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// Integration test: verify that sharkmer recovers 18S from ERR571460 (Porites lutea)
/// using the cnidaria primer panel with 100k reads.
#[test]
fn test_18s_recovery_from_err571460() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("ERR571460_100k_R1.fastq.gz");
    assert!(
        fixture.exists(),
        "Test fixture not found: {}",
        fixture.display()
    );

    let outdir = tempfile::tempdir().expect("failed to create temp dir");
    let sample = "ERR571460_test";

    let sharkmer_bin = env!("CARGO_BIN_EXE_sharkmer");

    // Pass the .gz file directly — sharkmer handles gzip natively
    let output = Command::new(sharkmer_bin)
        .args([
            "-k",
            "31",
            "--pcr",
            "cnidaria",
            "-s",
            sample,
            "-o",
            outdir.path().to_str().unwrap(),
            "--max-reads",
            "100000",
            "-n",
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

    // Check that 18S FASTA was produced
    let fasta_18s = outdir.path().join(format!("{}_cnidaria_18S.fasta", sample));
    assert!(
        fasta_18s.exists(),
        "18S FASTA not found at {}",
        fasta_18s.display()
    );

    let content = fs::read_to_string(&fasta_18s).expect("failed to read 18S FASTA");

    // Should have at least one sequence
    assert!(
        content.contains(">"),
        "18S FASTA contains no sequence records"
    );

    // 18S amplicon should be roughly 1700-1900 bp
    let seq: String = content.lines().filter(|l| !l.starts_with('>')).collect();
    assert!(
        seq.len() > 1700 && seq.len() < 1900,
        "18S sequence length {} outside expected range 1700-1900",
        seq.len()
    );
}
