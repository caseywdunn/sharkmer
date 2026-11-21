use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};
use std::collections::BTreeMap;

const URL: &str = "http://example.com/reads.fastq.gz";

fn budget(disk: u64) -> Arc<ResourceBudget> {
    budget_with_memory(disk, 8 << 20)
}

fn budget_with_memory(disk: u64, memory: u64) -> Arc<ResourceBudget> {
    let mut classes = [8 << 20; MEMORY_CLASSES];
    classes[crate::kmer::budget::MemoryClass::RetainedResults as usize] = memory;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: memory,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk,
        temporary_open_files: 64,
    })
    .unwrap()
}

fn pair(root: &Path, data: &[u8]) -> [Vec<u8>; 2] {
    let temporary = root.join("fixture");
    fs::write(&temporary, data).unwrap();
    let meta = CacheMeta {
        url: URL.to_owned(),
        sha256: compute_sha256(&temporary).unwrap(),
        complete: true,
        n_reads: 2,
    };
    fs::remove_file(temporary).unwrap();
    [
        data.to_vec(),
        serde_yaml_ng::to_string(&meta).unwrap().into_bytes(),
    ]
}

fn prepared(prior: bool) -> (tempfile::TempDir, Transaction, [Vec<u8>; 2], [Vec<u8>; 2]) {
    let root = tempfile::tempdir().unwrap();
    let key = cache_key(URL);
    let locations = Locations::new(root.path(), None).unwrap();
    let public = locations.public(key.as_bytes()).unwrap();
    let old = pair(root.path(), b"previous reads");
    let new = pair(root.path(), b"new reads");
    if prior {
        for i in 0..2 {
            fs::write(&public[i], &old[i]).unwrap();
        }
    }
    fs::write(root.path().join("foreign.txt"), b"foreign cache data").unwrap();
    let mut transaction = Transaction::begin(root.path(), &key, prior, None).unwrap();
    for (i, bytes) in new.iter().enumerate() {
        let mut temporary = CacheTemporary::create(
            transaction.directory(),
            if i == 0 {
                ".new-data-"
            } else {
                ".sharkmer-cache-meta-"
            },
            None,
        )
        .unwrap();
        temporary.writer().write_all(bytes).unwrap();
        temporary.finish(Ok(())).unwrap();
        if i == 0 {
            transaction.data = Some(temporary);
        } else {
            transaction.metadata = Some(temporary);
        }
    }
    transaction.prepare().unwrap();
    (root, transaction, old, new)
}

fn publish_one(transaction: &mut Transaction, index: usize) {
    let file = if index == 0 {
        transaction.data.as_mut().unwrap()
    } else {
        transaction.metadata.as_mut().unwrap()
    };
    file.begin_link_publication().unwrap();
    fs::hard_link(file.path(), &transaction.locations.link[index]).unwrap();
    fs::rename(
        &transaction.locations.link[index],
        &transaction.public[index],
    )
    .unwrap();
}

fn snapshot(root: &Path) -> BTreeMap<PathBuf, Vec<u8>> {
    fn visit(root: &Path, directory: &Path, files: &mut BTreeMap<PathBuf, Vec<u8>>) {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if fs::symlink_metadata(&path).unwrap().file_type().is_dir() {
                visit(root, &path, files);
            } else {
                files.insert(
                    path.strip_prefix(root).unwrap().to_path_buf(),
                    fs::read(path).unwrap(),
                );
            }
        }
    }
    let mut files = BTreeMap::new();
    visit(root, root, &mut files);
    files
}

fn check(root: &Path, wanted: Option<&[Vec<u8>; 2]>) {
    let locations = Locations::new(root, None).unwrap();
    let public = locations.public(cache_key(URL).as_bytes()).unwrap();
    for i in 0..2 {
        if let Some(wanted) = wanted {
            assert_eq!(fs::read(&public[i]).unwrap(), wanted[i]);
        } else {
            assert!(!public[i].exists());
        }
    }
    assert!(!locations.root.exists());
    assert!(!locations.owner.exists());
    assert_eq!(
        fs::read(root.join("foreign.txt")).unwrap(),
        b"foreign cache data"
    );
}

#[test]
fn interrupted_fresh_and_replacement_publication_restore_only_the_receipted_generation() {
    for prior in [false, true] {
        for published in 0..=2 {
            let (root, mut transaction, old, _) = prepared(prior);
            for i in 0..published {
                publish_one(&mut transaction, i);
            }
            // Model process death: destructors must not clean the abandoned scope.
            std::mem::forget(transaction);
            recover(root.path(), None).unwrap();
            check(root.path(), prior.then_some(&old));
            recover(root.path(), None).unwrap();
        }
    }
}

#[test]
fn accepted_pair_survives_interrupted_payload_and_control_cleanup() {
    for prior in [false, true] {
        for removed in [0, 1, 2] {
            let (root, mut transaction, _, new) = prepared(prior);
            publish_one(&mut transaction, 0);
            publish_one(&mut transaction, 1);
            sync_directory(root.path()).unwrap();
            File::create(&transaction.locations.accepted)
                .unwrap()
                .sync_all()
                .unwrap();
            sync_directory(&transaction.locations.root).unwrap();
            if removed >= 1 && prior {
                fs::remove_file(&transaction.locations.old[0]).unwrap();
            }
            if removed == 2 {
                for entry in fs::read_dir(transaction.directory()).unwrap() {
                    let path = entry.unwrap().path();
                    if path != transaction.locations.context.as_ref()
                        && path != transaction.locations.accepted.as_ref()
                    {
                        fs::remove_file(path).unwrap();
                    }
                }
            }
            std::mem::forget(transaction);
            recover(root.path(), None).unwrap();
            check(root.path(), Some(&new));
        }
    }
}

#[test]
fn missing_or_modified_receipt_refuses_without_discarding_backup_evidence() {
    for missing in [false, true] {
        let (root, mut transaction, _, _) = prepared(true);
        publish_one(&mut transaction, 0);
        if missing {
            fs::remove_file(&transaction.locations.journal).unwrap();
        } else {
            fs::write(&transaction.locations.journal, b"changed receipt").unwrap();
        }
        std::mem::forget(transaction);
        let before = snapshot(root.path());
        assert!(recover(root.path(), None).is_err());
        assert_eq!(snapshot(root.path()), before);
    }
}

#[test]
fn failed_restoration_preserves_both_backups_and_retries_after_one_public_file_was_restored() {
    let (root, mut transaction, old, _) = prepared(true);
    publish_one(&mut transaction, 0);
    publish_one(&mut transaction, 1);
    fs::hard_link(
        &transaction.locations.old[0],
        &transaction.locations.restore[0],
    )
    .unwrap();
    fs::rename(&transaction.locations.restore[0], &transaction.public[0]).unwrap();
    // A foreign non-file blocks restoration of the remaining receipt.
    fs::create_dir(&transaction.locations.restore[1]).unwrap();
    std::mem::forget(transaction);
    let locations = Locations::new(root.path(), None).unwrap();
    let before = snapshot(root.path());
    assert!(recover(root.path(), None).is_err());
    assert_eq!(snapshot(root.path()), before);
    for (i, expected) in old.iter().enumerate() {
        assert_eq!(fs::read(&locations.old[i]).unwrap(), *expected);
    }
    fs::remove_dir(&locations.restore[1]).unwrap();
    recover(root.path(), None).unwrap();
    check(root.path(), Some(&old));
}

#[test]
fn cold_recovery_quota_refuses_before_mutation_then_recovers_under_the_original_allowance() {
    let (root, mut transaction, old, _) = prepared(true);
    publish_one(&mut transaction, 0);
    std::mem::forget(transaction);
    let before = snapshot(root.path());
    let refused = budget(1);
    let error = recover(root.path(), Some(&refused)).unwrap_err();
    assert_eq!(
        error.downcast_ref::<ResourceExhausted>().unwrap().resource,
        ResourceKind::TemporaryDisk
    );
    assert_eq!(snapshot(root.path()), before);
    assert_eq!(refused.snapshot().reserved_memory_bytes, 0);
    assert_eq!(refused.snapshot().reserved_temporary_disk_bytes, 0);
    let admitted = budget(1 << 20);
    recover(root.path(), Some(&admitted)).unwrap();
    check(root.path(), Some(&old));
    assert_eq!(admitted.snapshot().reserved_memory_bytes, 0);
    assert_eq!(admitted.snapshot().reserved_temporary_disk_bytes, 0);
}

#[test]
fn unregistered_root_and_foreign_legacy_temporary_names_are_never_cleaned_by_prefix() {
    let root = tempfile::tempdir().unwrap();
    let locations = Locations::new(root.path(), None).unwrap();
    fs::create_dir(&locations.root).unwrap();
    fs::write(locations.root.try_join("foreign").unwrap(), b"unowned").unwrap();
    let before = snapshot(root.path());
    assert!(recover(root.path(), None).is_err());
    assert_eq!(snapshot(root.path()), before);
    fs::remove_file(locations.root.try_join("foreign").unwrap()).unwrap();
    fs::remove_dir(&locations.root).unwrap();
    let legacy = root.path().join(".sharkmer-cache-data-backup-foreign");
    fs::write(&legacy, b"unregistered old recovery bytes").unwrap();
    recover(root.path(), None).unwrap();
    assert_eq!(
        fs::read(legacy).unwrap(),
        b"unregistered old recovery bytes"
    );
}

#[test]
fn preparatory_scope_and_final_empty_scope_can_resume_cleanup() {
    for context in [false, true] {
        let root = tempfile::tempdir().unwrap();
        let mut transaction =
            Transaction::begin(root.path(), &cache_key(URL), false, None).unwrap();
        if !context {
            for entry in fs::read_dir(transaction.directory()).unwrap() {
                fs::remove_file(entry.unwrap().path()).unwrap();
            }
        }
        transaction.registered = false;
        // Keep the precise scope; its abandoned temporary aliases may disappear.
        drop(transaction);
        recover(root.path(), None).unwrap();
        assert!(!root.path().join(ROOT).exists());
        assert!(!root.path().join(OWNER).exists());
    }
}

#[test]
fn modified_backup_or_public_file_refuses_before_any_restoration() {
    for private in [false, true] {
        let (root, mut transaction, _, _) = prepared(true);
        publish_one(&mut transaction, 0);
        let changed = if private {
            &transaction.locations.old[0]
        } else {
            &transaction.public[0]
        };
        fs::write(changed, b"foreign modification").unwrap();
        std::mem::forget(transaction);
        let before = snapshot(root.path());
        assert!(recover(root.path(), None).is_err());
        assert_eq!(snapshot(root.path()), before);
    }
}

// Native tests launch only this fixture in a separate process. The environment
// configures a test fixture, never a fault hook in the production application.
#[test]
#[ignore = "external native HTTP/fault driver supplies the fixture directory and URL"]
fn native_cache_fixture() {
    let root = PathBuf::from(std::env::var_os("SHARKMER_TEST_CACHE_DIRECTORY").unwrap());
    let url = std::env::var("SHARKMER_TEST_CACHE_URL").unwrap();
    let action = std::env::var("SHARKMER_TEST_CACHE_ACTION").unwrap();
    let maximum = std::env::var("SHARKMER_TEST_CACHE_MAXIMUM")
        .unwrap()
        .parse()
        .unwrap();
    let config = CacheConfig::new(Some(&root)).unwrap();
    // Keep the original application workspace and disk/descriptor limits.
    // HTTP client internals are outside the application ledger.
    let ledger = budget(64 << 30);
    if action == "fill" {
        config
            .download_to_cache_bounded(&url, maximum, 128, &ledger)
            .unwrap();
    } else {
        assert_eq!(action, "recover");
        config.lookup_bounded(&url, maximum, &ledger).unwrap();
    }
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
    assert!(ledger.snapshot().valid);
    assert!(!root.join(ROOT).exists());
    assert!(!root.join(OWNER).exists());
}
