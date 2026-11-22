use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};
use std::path::PathBuf;

fn snapshot(directory: &Path) -> std::collections::BTreeMap<PathBuf, Vec<u8>> {
    fn visit(
        base: &Path,
        directory: &Path,
        rows: &mut std::collections::BTreeMap<PathBuf, Vec<u8>>,
    ) {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(base, &path, rows);
            } else {
                rows.insert(
                    path.strip_prefix(base).unwrap().to_path_buf(),
                    fs::read(path).unwrap(),
                );
            }
        }
    }
    let mut rows = std::collections::BTreeMap::new();
    visit(directory, directory, &mut rows);
    rows
}

#[test]
fn markerless_recovery_refuses_its_original_disk_ceiling_before_creating_grants_then_retries_exactly()
 {
    let directory = tempfile::tempdir().unwrap();
    let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
    for (name, bytes) in [
        ("sample_a.fasta", b"old sequence".as_slice()),
        ("sample.stats.yaml", b"old stats".as_slice()),
    ] {
        first
            .stage_file(name, |file| {
                file.write_all(bytes)?;
                Ok(())
            })
            .unwrap();
    }
    first.commit().unwrap();
    let public = snapshot(directory.path());
    let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
    for (name, bytes) in [
        ("sample_a.fasta", b"new sequence".as_slice()),
        ("sample.stats.yaml", b"new stats".as_slice()),
    ] {
        next.stage_file(name, |file| {
            file.write_all(bytes)?;
            Ok(())
        })
        .unwrap();
    }
    assert!(
        next.commit_inner(|step| {
            ensure!(
                step != PublicationStep::BackupReady,
                "interrupted after journal"
            );
            Ok(())
        })
        .is_err()
    );
    let staging = next.staging_directory.to_path_buf();
    next.committed = true;
    drop(next);
    for entry in fs::read_dir(&staging).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap();
        if name.starts_with(".sharkmer-owner-v1-") || name == ".sharkmer-preparation-v1" {
            fs::remove_file(path).unwrap();
        }
    }
    let before = snapshot(directory.path());
    let required: u64 = before
        .iter()
        .filter(|(name, _)| {
            name.starts_with(staging.file_name().unwrap())
                || name.as_os_str() == ".sharkmer-recovery-sample.yaml"
        })
        .map(|(_, bytes)| bytes.len() as u64)
        .sum();
    assert!(required > 1);
    for allowance in [required - 1, required] {
        let ledger = ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: allowance,
            temporary_open_files: 128,
        })
        .unwrap();
        let root = BudgetPath::from_path(
            directory.path(),
            Some(&ledger),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let result = recover(&root, "sample");
        if allowance < required {
            let error = result.err().unwrap();
            assert_eq!(
                error.downcast_ref::<ResourceExhausted>().unwrap().resource,
                ResourceKind::TemporaryDisk
            );
            assert_eq!(snapshot(directory.path()), before);
        } else {
            assert!(matches!(result.unwrap(), Some(RecoveryOutcome::Restored)));
            assert_eq!(snapshot(directory.path()), public);
            assert!(!staging.exists());
            assert_eq!(
                ledger.snapshot().peak_reserved_temporary_disk_bytes,
                required
            );
        }
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        drop(root);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}
