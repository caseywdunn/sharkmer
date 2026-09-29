//! Inclusive classified FASTA, exact provenance, and generation ownership.

use bio::io::fasta;
use serde_yaml_ng::Value;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SAMPLE: &str = "uncertain_test";

fn deterministic_dna(mut state: u64, length: usize) -> String {
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[((state >> 32) & 3) as usize] as char
        })
        .collect()
}

struct Fixture {
    _directory: tempfile::TempDir,
    input: PathBuf,
    primers: Vec<String>,
    clean_target: String,
}

impl Fixture {
    fn new() -> Self {
        // Keep the repeat sequence identical to the existing diagnostic CLI
        // fixture. An independent clean locus makes standard parity nonempty.
        let left = deterministic_dna(157, 60);
        let right = deterministic_dna(158, 60);
        let repeat = format!("{}{}{}", left, "A".repeat(40), right);
        let clean_target = deterministic_dna(159, 160);
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("reads.fastq");
        let mut raw = String::new();
        let mut primers = Vec::new();
        for (gene, sequence) in [("repeat", &repeat), ("clean", &clean_target)] {
            for index in 0..4 {
                raw.push_str(&format!(
                    "@{gene}_{index}\n{sequence}\n+\n{}\n",
                    "I".repeat(sequence.len())
                ));
            }
            let reverse = String::from_utf8(bio::alphabets::dna::revcomp(
                &sequence.as_bytes()[sequence.len() - 15..],
            ))
            .unwrap();
            primers.push(format!(
                "name={gene},forward={},reverse={reverse},trim=15,mismatches=0,min-length=100,max-length=220,dedup-edit-threshold=0",
                &sequence[..15]
            ));
        }
        fs::write(&input, raw).unwrap();
        Self {
            _directory: directory,
            input,
            primers,
            clean_target,
        }
    }

    fn run(&self, directory: &Path, emit: bool, diagnose: bool) -> Output {
        let mut command = Command::new(env!("CARGO_BIN_EXE_sharkmer"));
        command.args(["-k", "19", "--threads", "1", "-s", SAMPLE, "-o"]);
        command.arg(directory);
        for specification in &self.primers {
            command.args(["--pcr-primers", specification]);
        }
        if !emit {
            command.arg("--no-emit-uncertain");
        }
        if diagnose {
            command.arg("--diagnose-withheld-paths");
        }
        command.arg(&self.input).output().unwrap()
    }

    fn run_success(&self, directory: &Path, emit: bool, diagnose: bool) -> Value {
        let output = self.run(directory, emit, diagnose);
        assert!(
            output.status.success(),
            "sharkmer failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        read_yaml(&directory.join(format!("{SAMPLE}.stats.yaml")))
    }
}

fn read_yaml(path: &Path) -> Value {
    serde_yaml_ng::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

fn gene<'a>(stats: &'a Value, name: &str) -> &'a Value {
    stats["pcr_results"]
        .as_sequence()
        .unwrap()
        .iter()
        .find(|value| value["gene_name"].as_str() == Some(name))
        .unwrap()
}

fn standard_gene(value: &Value) -> Value {
    let mut normalized = value.clone();
    for field in [
        "status",
        "n_products",
        "n_uncertain_products",
        "product_lengths",
        "output_file",
    ] {
        normalized
            .as_mapping_mut()
            .unwrap()
            .remove(Value::String(field.to_string()));
    }
    normalized
        .as_mapping_mut()
        .unwrap()
        .remove(Value::String("uncertain_products".to_string()));
    if let Some(thresholds) = normalized["threshold_diagnostics"].as_sequence_mut() {
        for threshold in thresholds {
            threshold
                .as_mapping_mut()
                .unwrap()
                .remove(Value::String("withheld_path_diagnostics".to_string()));
        }
    }
    normalized
}

fn sha256(contents: &[u8]) -> String {
    format!("{:x}", Sha256::digest(contents))
}

fn assert_manifest(directory: &Path, stats: &Value, expected_files: &[String]) {
    let manifest = read_yaml(&directory.join(format!("{SAMPLE}.manifest.yaml")));
    assert_eq!(manifest["schema_version"].as_u64(), Some(2));
    assert_eq!(stats["product_output_schema_version"].as_u64(), Some(2));
    assert_eq!(
        stats["product_output_layout"].as_str(),
        Some("combined_classified")
    );
    assert_eq!(manifest["status"].as_str(), Some("complete"));
    assert_eq!(stats["run_status"].as_str(), Some("complete"));
    assert_eq!(manifest["run_id"], stats["run_id"]);
    let receipts = manifest["files"].as_sequence().unwrap();
    let names: BTreeSet<String> = receipts
        .iter()
        .map(|receipt| receipt["path"].as_str().unwrap().to_string())
        .collect();
    assert_eq!(names.len(), receipts.len());
    assert_eq!(
        names,
        expected_files.iter().cloned().collect::<BTreeSet<_>>()
    );
    for receipt in receipts {
        let bytes = fs::read(directory.join(receipt["path"].as_str().unwrap())).unwrap();
        assert_eq!(receipt["sha256"].as_str().unwrap(), sha256(&bytes));
    }
}

fn expected_files(emit: bool) -> Vec<String> {
    let mut files = vec![
        format!("{SAMPLE}.stats.yaml"),
        format!("{SAMPLE}_clean.fasta"),
    ];
    if emit {
        files.push(format!("{SAMPLE}_repeat.fasta"));
    }
    files
}

fn assert_uncertain_agreement(directory: &Path, stats: &Value) {
    let repeat = gene(stats, "repeat");
    assert_eq!(repeat["status"].as_str(), Some("uncertain"));
    assert_eq!(repeat["n_standard_products"].as_u64(), Some(0));
    assert!(repeat["failure_reason"].as_str().is_some());
    let uncertainty = &repeat["uncertain_products"];
    assert_eq!(
        uncertainty["read_support_evaluation"].as_str(),
        Some("not_evaluated")
    );
    assert_eq!(
        uncertainty["repeat_length_evaluation"].as_str(),
        Some("unresolved")
    );
    let products = uncertainty["products"].as_sequence().unwrap();
    assert!(!products.is_empty());
    assert_eq!(
        uncertainty["n_products"].as_u64().unwrap() as usize,
        products.len()
    );
    assert_eq!(uncertainty["unique_products"], uncertainty["n_products"]);
    assert_eq!(repeat["n_products"], uncertainty["n_products"]);
    assert_eq!(repeat["n_uncertain_products"], uncertainty["n_products"]);
    assert_eq!(repeat["product_lengths"], uncertainty["product_lengths"]);
    let basename = uncertainty["output_file"].as_str().unwrap();
    assert_eq!(basename, format!("{SAMPLE}_repeat.fasta"));
    assert_eq!(repeat["output_file"], uncertainty["output_file"]);
    assert!(
        !directory
            .join(format!("{SAMPLE}_repeat.uncertain.fasta"))
            .exists()
    );
    let records: Vec<_> = fasta::Reader::from_file(directory.join(basename))
        .unwrap()
        .records()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(records.len(), products.len());
    let mut oriented_sequences = BTreeSet::new();
    let mut source_count = 0_u64;
    for (index, (record, product)) in records.iter().zip(products).enumerate() {
        assert_eq!(record.id(), product["id"].as_str().unwrap());
        assert!(record.id().ends_with(".uncertain"));
        assert_eq!(record.desc(), product["description"].as_str());
        let description = record.desc().unwrap();
        for label in [
            "assembly_status=uncertain",
            "uncertainty=repeat_structure",
            "repeat_length=unresolved",
            "read_support=not_evaluated",
            "enumeration=bounded",
            "limit_scope=gene",
        ] {
            assert!(description.split_whitespace().any(|field| field == label));
        }
        assert_eq!(
            record.seq(),
            product["oriented_sequence"].as_str().unwrap().as_bytes()
        );
        assert_eq!(
            sha256(record.seq()),
            product["sequence_sha256"].as_str().unwrap()
        );
        assert_eq!(
            record.seq().len(),
            product["sequence_length"].as_u64().unwrap() as usize
        );
        assert_eq!(
            record.seq().len(),
            uncertainty["product_lengths"][index].as_u64().unwrap() as usize
        );
        assert!(oriented_sequences.insert(record.seq().to_vec()));
        let sources = product["sources"].as_sequence().unwrap();
        assert!(!sources.is_empty());
        source_count += sources.len() as u64;
        for source in sources {
            assert!(source["threshold"].as_u64().unwrap() > 0);
            let markers = source["marker_occurrences"].as_sequence().unwrap();
            assert!(!markers.is_empty());
            for marker in markers {
                let start = marker["start"].as_u64().unwrap() as usize;
                let end = marker["end"].as_u64().unwrap() as usize;
                let identity = marker["identity_sequence"].as_str().unwrap().as_bytes();
                assert_eq!(record.seq().get(start..end), Some(identity));
                assert_eq!(
                    sha256(identity),
                    marker["identity_sha256"].as_str().unwrap()
                );
            }
        }
    }
    let retained = uncertainty["retained_path_records"].as_u64().unwrap();
    let excluded = uncertainty["excluded_standard_equivalent_path_records"]
        .as_u64()
        .unwrap();
    assert_eq!(source_count, retained - excluded);
    assert_eq!(
        retained,
        excluded
            + uncertainty["exact_duplicate_path_records"]
                .as_u64()
                .unwrap()
            + products.len() as u64
    );

    let clean_uncertainty = &gene(stats, "clean")["uncertain_products"];
    assert_eq!(gene(stats, "clean")["n_products"].as_u64(), Some(1));
    assert_eq!(
        gene(stats, "clean")["n_standard_products"].as_u64(),
        Some(1)
    );
    assert_eq!(
        gene(stats, "clean")["n_uncertain_products"].as_u64(),
        Some(0)
    );
    assert_eq!(clean_uncertainty["n_products"].as_u64(), Some(0));
    assert!(clean_uncertainty["output_file"].is_null());
    assert!(
        !directory
            .join(format!("{SAMPLE}_clean.uncertain.fasta"))
            .exists()
    );
}

#[test]
fn uncertain_cli_preserves_standard_outputs_and_agrees_with_stats_and_manifest() {
    let fixture = Fixture::new();
    let off = tempfile::tempdir().unwrap();
    let on = tempfile::tempdir().unwrap();
    let both = tempfile::tempdir().unwrap();
    let diagnose = tempfile::tempdir().unwrap();
    let default_stats = fixture.run_success(off.path(), false, false);
    let emitted_stats = fixture.run_success(on.path(), true, false);
    let both_stats = fixture.run_success(both.path(), true, true);
    let diagnostic_stats = fixture.run_success(diagnose.path(), false, true);

    let standard_name = format!("{SAMPLE}_clean.fasta");
    let standard_bytes = fs::read(off.path().join(&standard_name)).unwrap();
    let clean_records: Vec<_> = fasta::Reader::new(standard_bytes.as_slice())
        .records()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(clean_records.len(), 1);
    assert_eq!(clean_records[0].seq(), fixture.clean_target.as_bytes());
    assert!(
        clean_records[0]
            .desc()
            .unwrap()
            .split_whitespace()
            .any(|field| field == "assembly_status=standard")
    );
    for (directory, stats, emit, hidden) in [
        (off.path(), &default_stats, false, false),
        (on.path(), &emitted_stats, true, false),
        (both.path(), &both_stats, true, true),
        (diagnose.path(), &diagnostic_stats, false, true),
    ] {
        assert_eq!(
            fs::read(directory.join(&standard_name)).unwrap(),
            standard_bytes
        );
        for field in [
            "n_reads_read",
            "n_bases_read",
            "n_subreads_ingested",
            "n_bases_ingested",
            "n_kmers",
            "n_multi_kmers",
            "n_singleton_kmers",
        ] {
            assert_eq!(
                stats[field], default_stats[field],
                "count field {field} changed"
            );
        }
        for name in ["repeat", "clean"] {
            assert_eq!(
                standard_gene(gene(stats, name)),
                standard_gene(gene(&default_stats, name))
            );
            assert_eq!(!gene(stats, name)["uncertain_products"].is_null(), emit);
        }
        let thresholds = gene(stats, "repeat")["threshold_diagnostics"]
            .as_sequence()
            .unwrap();
        assert_eq!(
            thresholds
                .iter()
                .any(|threshold| !threshold["withheld_path_diagnostics"].is_null()),
            hidden
        );
        assert_manifest(directory, stats, &expected_files(emit));
        if emit {
            assert_uncertain_agreement(directory, stats);
        } else {
            assert!(!directory.join(format!("{SAMPLE}_repeat.fasta")).exists());
        }
    }
    assert_eq!(
        gene(&emitted_stats, "repeat")["uncertain_products"],
        gene(&both_stats, "repeat")["uncertain_products"]
    );
    assert_eq!(
        fs::read(on.path().join(format!("{SAMPLE}_repeat.fasta"))).unwrap(),
        fs::read(both.path().join(format!("{SAMPLE}_repeat.fasta"))).unwrap()
    );
}

#[test]
fn turning_uncertain_emission_off_removes_only_receipted_uncertain_output() {
    let fixture = Fixture::new();
    let directory = tempfile::tempdir().unwrap();
    let first = fixture.run_success(directory.path(), true, false);
    assert_uncertain_agreement(directory.path(), &first);
    let foreign = directory
        .path()
        .join(format!("{SAMPLE}_unowned.uncertain.fasta"));
    fs::write(&foreign, ">user_record\nACGT\n").unwrap();
    let standard = fs::read(directory.path().join(format!("{SAMPLE}_clean.fasta"))).unwrap();
    let second = fixture.run_success(directory.path(), false, false);
    assert_ne!(first["run_id"], second["run_id"]);
    assert!(
        !directory
            .path()
            .join(format!("{SAMPLE}_repeat.fasta"))
            .exists()
    );
    assert_eq!(fs::read_to_string(foreign).unwrap(), ">user_record\nACGT\n");
    assert_eq!(
        fs::read(directory.path().join(format!("{SAMPLE}_clean.fasta"))).unwrap(),
        standard
    );
    assert_manifest(directory.path(), &second, &expected_files(false));
    assert!(gene(&second, "repeat")["uncertain_products"].is_null());
}

#[test]
fn modified_owned_uncertain_output_refuses_rerun_and_preserves_prior_generation() {
    let fixture = Fixture::new();
    let directory = tempfile::tempdir().unwrap();
    let first = fixture.run_success(directory.path(), true, false);
    assert_manifest(directory.path(), &first, &expected_files(true));
    let uncertain_path = directory.path().join(format!("{SAMPLE}_repeat.fasta"));
    let mut modified = fs::read(&uncertain_path).unwrap();
    modified.extend_from_slice(b"\n>user_appended\nACGT\n");
    fs::write(&uncertain_path, &modified).unwrap();
    let preserved = [
        format!("{SAMPLE}.manifest.yaml"),
        format!("{SAMPLE}.stats.yaml"),
        format!("{SAMPLE}_clean.fasta"),
    ]
    .into_iter()
    .map(|name| (name.clone(), fs::read(directory.path().join(name)).unwrap()))
    .collect::<Vec<_>>();

    let failed = fixture.run(directory.path(), false, false);
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("modified prior output"));
    assert_eq!(fs::read(&uncertain_path).unwrap(), modified);
    for (name, contents) in preserved {
        assert_eq!(fs::read(directory.path().join(name)).unwrap(), contents);
    }
}

#[test]
fn mixed_target_emits_standard_then_uncertain_and_opt_out_preserves_standard() {
    let fixture = Fixture::new();
    let mixed = format!(
        "{}{}{}",
        deterministic_dna(157, 60),
        deterministic_dna(160, 40),
        deterministic_dna(158, 60)
    );
    let mut raw = fs::read_to_string(&fixture.input).unwrap();
    for index in 0..4 {
        raw.push_str(&format!(
            "@mixed_{index}\n{mixed}\n+\n{}\n",
            "I".repeat(mixed.len())
        ));
    }
    fs::write(&fixture.input, raw).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let first = fixture.run_success(directory.path(), true, false);
    let result = gene(&first, "repeat");
    assert_eq!(result["status"].as_str(), Some("success"));
    assert_eq!(result["n_standard_products"].as_u64(), Some(1));
    let uncertain = result["n_uncertain_products"].as_u64().unwrap();
    assert!(uncertain > 0);
    assert_eq!(result["n_products"].as_u64(), Some(1 + uncertain));
    let path = directory.path().join(format!("{SAMPLE}_repeat.fasta"));
    let records = fasta::Reader::from_file(&path)
        .unwrap()
        .records()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records.len(), 1 + uncertain as usize);
    assert_eq!(records[0].seq(), mixed.as_bytes());
    assert!(
        records[0]
            .desc()
            .unwrap()
            .split_whitespace()
            .any(|field| field == "assembly_status=standard")
    );
    assert!(records[1..].iter().all(|record| {
        record
            .desc()
            .unwrap()
            .split_whitespace()
            .any(|field| field == "assembly_status=uncertain")
    }));
    let lengths = result["product_lengths"].as_sequence().unwrap();
    assert_eq!(lengths.len(), records.len());
    for (length, record) in lengths.iter().zip(&records) {
        assert_eq!(length.as_u64(), Some(record.seq().len() as u64));
    }
    assert_eq!(
        result["output_file"],
        result["uncertain_products"]["output_file"]
    );
    let second = fixture.run_success(directory.path(), false, false);
    let standard_only = fasta::Reader::from_file(path)
        .unwrap()
        .records()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(standard_only.len(), 1);
    assert_eq!(standard_only[0].id(), records[0].id());
    assert_eq!(standard_only[0].desc(), records[0].desc());
    assert_eq!(standard_only[0].seq(), records[0].seq());
    assert_eq!(gene(&second, "repeat")["n_products"].as_u64(), Some(1));
    assert_eq!(
        gene(&second, "repeat")["n_uncertain_products"].as_u64(),
        Some(0)
    );
}

#[test]
fn inclusive_console_reports_uncertain_only_output_and_zero_targets_remain_explicit() {
    let mut fixture = Fixture::new();
    fixture.primers.push(format!(
        "name=missing,forward={},reverse={},trim=15,mismatches=0,min-length=100,max-length=220",
        deterministic_dna(170, 20),
        deterministic_dna(171, 20)
    ));
    let directory = tempfile::tempdir().unwrap();
    let output = fixture.run(directory.path(), true, false);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("0 standard,") && stderr.contains("uncertain"));
    assert!(stderr.contains("targets with output"));
    assert!(!stderr.contains("repeat (no products"));
    let stats = read_yaml(&directory.path().join(format!("{SAMPLE}.stats.yaml")));
    let missing = gene(&stats, "missing");
    assert_eq!(missing["status"].as_str(), Some("fail"));
    for field in ["n_products", "n_standard_products", "n_uncertain_products"] {
        assert_eq!(missing[field].as_u64(), Some(0));
    }
    assert!(missing["output_file"].is_null());
    assert!(missing["failure_reason"].as_str().is_some());
    assert!(
        !directory
            .path()
            .join(format!("{SAMPLE}_missing.fasta"))
            .exists()
    );
}

#[test]
fn legacy_gene_name_suffix_cannot_collide_with_uncertain_record_ids() {
    let mut fixture = Fixture::new();
    fixture.primers[1] = fixture.primers[1].replace("name=clean,", "name=repeat.uncertain,");
    let directory = tempfile::tempdir().unwrap();
    let stats = fixture.run_success(directory.path(), true, false);
    let mut ids = BTreeSet::new();
    for name in ["repeat", "repeat.uncertain"] {
        let path = directory
            .path()
            .join(gene(&stats, name)["output_file"].as_str().unwrap());
        let records = fasta::Reader::from_file(path)
            .unwrap()
            .records()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!records.is_empty());
        for record in records {
            assert!(ids.insert(record.id().to_string()));
        }
    }
}
