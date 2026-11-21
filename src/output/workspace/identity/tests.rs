use super::*;
#[test]
fn identity_backing_stays_charged_and_reuses_its_allocation_for_verification() {
    let parent = tempfile::tempdir().unwrap();
    let ledger = super::super::tests::budget(0);
    let path =
        BudgetPath::from_path(parent.path(), Some(&ledger), MemoryClass::RetainedResults).unwrap();
    let baseline = ledger.snapshot().reserved_memory_bytes;
    let identity = Identity::capture(&path).unwrap();
    assert_eq!(
        ledger.snapshot().reserved_memory_bytes,
        baseline + parent.path().as_os_str().len() as u64 + 1
    );
    let peak = ledger.snapshot().peak_reserved_memory_bytes;
    identity.verify(&path).unwrap();
    assert_eq!(ledger.snapshot().peak_reserved_memory_bytes, peak);
    drop(identity);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, baseline);
    drop(path);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn identity_heap_observations_include_scan_and_long_path_workspace_and_release() {
    const ENV: &str = "SHARKMER_TEST_WORKSPACE_IDENTITY_ALLOCATIONS";
    if std::env::var_os(ENV).is_none() {
        let status=std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact","output::workspace::identity::tests::identity_heap_observations_include_scan_and_long_path_workspace_and_release","--nocapture","--test-threads=1"]).env(ENV,"1").status().unwrap();
        assert!(status.success());
        return;
    }
    for (index, shared, long) in [(0, false, false), (1, true, false), (2, false, true)] {
        let parent = if shared {
            tempfile::tempdir_in(".").unwrap()
        } else {
            tempfile::tempdir().unwrap()
        };
        let mut directory = parent.path().to_path_buf();
        if long {
            for _ in 0..30 {
                directory.push("identity-path-workspace");
            }
            fs::create_dir_all(&directory).unwrap();
        }
        let ledger = super::super::tests::budget(0);
        let path =
            BudgetPath::from_path(&directory, Some(&ledger), MemoryClass::RetainedResults).unwrap();
        let before = ledger.snapshot().reserved_memory_bytes;
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let identity = Identity::capture(&path).unwrap();
        identity.verify(&path).unwrap();
        let observed = crate::PEAK_ALLOC.peak_usage() - baseline;
        let reserved = ledger.snapshot().reserved_memory_bytes - before;
        assert!(
            observed as u64 <= reserved,
            "case={index} observed={observed} reserved={reserved}"
        );
        drop(identity);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, before);
        println!(
            "workspace_identity_allocator case={index} observed={observed} reserved={reserved}"
        );
    }
}

#[test]
fn creating_and_removing_children_keeps_parent_identity_stable() {
    for shared in [false, true] {
        let directory = if shared {
            tempfile::tempdir_in(".").unwrap()
        } else {
            tempfile::tempdir().unwrap()
        };
        let parent =
            BudgetPath::from_path(directory.path(), None, MemoryClass::RetainedResults).unwrap();
        let identity = Identity::capture(&parent).unwrap();
        fs::create_dir(parent.join("child")).unwrap();
        fs::write(parent.join("child/data"), b"payload").unwrap();
        identity.verify(&parent).unwrap();
        fs::remove_file(parent.join("child/data")).unwrap();
        fs::remove_dir(parent.join("child")).unwrap();
        identity.verify(&parent).unwrap();
    }
}

#[test]
fn pathname_workspace_refuses_before_inspecting_an_unavailable_parent() {
    use crate::kmer::budget::{ResourceExhausted, ResourceKind};
    let directory = tempfile::tempdir().unwrap();
    let ledger = super::super::tests::budget(0);
    let parent = BudgetPath::from_path(
        &directory.path().join("unavailable"),
        Some(&ledger),
        MemoryClass::RetainedResults,
    )
    .unwrap();
    let used = ledger.snapshot().reserved_memory_bytes;
    let competing = ledger
        .reserve(MemoryClass::RetainedResults, (8 << 20) - used)
        .unwrap();
    let error = Identity::capture(&parent).err().unwrap();
    assert_eq!(
        error.downcast_ref::<ResourceExhausted>().unwrap().resource,
        ResourceKind::StageMemory(MemoryClass::RetainedResults)
    );
    assert!(!parent.exists());
    drop(competing);
    drop(parent);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}
