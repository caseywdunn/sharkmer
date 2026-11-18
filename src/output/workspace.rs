//! Registered private storage outside the output directory.
//!
//! The durable registry precedes directory creation, and a synced marker
//! precedes payloads. Cleanup owns only that exact child, never its parent.

use super::*;

mod identity;

pub(super) const REGISTRY: &str = ".sharkmer-workspace-v1";
const TEMPORARY: &str = ".sharkmer-workspace-";
const MAX_PATH_BYTES: usize = 128 << 10;
const CANONICALIZATION_BYTES: u64 = 256 << 10;
#[cfg(unix)]
const MAGIC: &[u8; 8] = b"SHKWRK3U";
#[cfg(windows)]
const MAGIC: &[u8; 8] = b"SHKWRK3W";
#[cfg(not(any(unix, windows)))]
const MAGIC: &[u8; 8] = b"SHKWRK3F";

fn canonical(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetPath> {
    // Retain a bounded allowance across the platform library's temporary path
    // buffers and the separately charged surviving path owner.
    let _workspace = budget
        .map(|budget| budget.reserve(MemoryClass::RetainedResults, CANONICALIZATION_BYTES))
        .transpose()?;
    let canonical = fs::canonicalize(path)?;
    BudgetPath::from_path(&canonical, budget, MemoryClass::RetainedResults)
}

fn root_path(parent: &BudgetPath, run_id: &str) -> Result<BudgetPath> {
    let name = BudgetString::format(
        format_args!(".sharkmer-workspace-{run_id}"),
        parent.budget(),
        MemoryClass::RetainedResults,
    )?;
    parent.try_join(name.as_str())
}

fn owner_path(root: &BudgetPath, run_id: &str) -> Result<BudgetPath> {
    let name = BudgetString::format(
        format_args!(".sharkmer-workspace-owner-v1-{run_id}"),
        root.budget(),
        MemoryClass::RetainedResults,
    )?;
    root.try_join(name.as_str())
}

#[cfg(unix)]
fn encode(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetVec<u8>> {
    use std::os::unix::ffi::OsStrExt;
    BudgetVec::collect_in_class(
        path.as_os_str().as_bytes().iter().copied(),
        budget,
        MemoryClass::RetainedResults,
    )
}
#[cfg(windows)]
fn encode(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetVec<u8>> {
    use std::os::windows::ffi::OsStrExt;
    BudgetVec::collect_in_class(
        path.as_os_str().encode_wide().flat_map(u16::to_le_bytes),
        budget,
        MemoryClass::RetainedResults,
    )
}
#[cfg(not(any(unix, windows)))]
fn encode(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetVec<u8>> {
    BudgetVec::collect_in_class(
        path.to_str()
            .context("Unsupported workspace path encoding")?
            .bytes(),
        budget,
        MemoryClass::RetainedResults,
    )
}

#[cfg(unix)]
fn decode(bytes: &[u8], budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetPath> {
    use std::os::unix::ffi::OsStrExt;
    BudgetPath::from_path(
        Path::new(std::ffi::OsStr::from_bytes(bytes)),
        budget,
        MemoryClass::RetainedResults,
    )
}
#[cfg(windows)]
fn decode(bytes: &[u8], budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetPath> {
    use std::os::windows::ffi::OsStringExt;
    ensure!(
        bytes.len().is_multiple_of(2),
        "Invalid workspace path encoding"
    );
    let units = BudgetVec::collect_in_class(
        bytes
            .as_chunks::<2>()
            .0
            .iter()
            .map(|pair| u16::from_le_bytes(*pair)),
        budget,
        MemoryClass::RetainedResults,
    )?;
    let _temporary = budget
        .map(|b| {
            b.reserve(
                MemoryClass::RetainedResults,
                u64::try_from(bytes.len())?
                    .checked_mul(3)
                    .context("Workspace path allocation overflow")?,
            )
        })
        .transpose()?;
    let path = std::ffi::OsString::from_wide(&units);
    BudgetPath::from_path(Path::new(&path), budget, MemoryClass::RetainedResults)
}
#[cfg(not(any(unix, windows)))]
fn decode(bytes: &[u8], budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetPath> {
    BudgetPath::from_path(
        Path::new(std::str::from_utf8(bytes)?),
        budget,
        MemoryClass::RetainedResults,
    )
}

fn marker(path: &Path) -> Result<()> {
    let metadata = fs::symlink_metadata(path)?;
    ensure!(
        metadata.file_type().is_file() && metadata.len() == 0,
        "Invalid workspace owner marker {}",
        path.display()
    );
    Ok(())
}

/// Live immutable registry allowance and preallocated error-cleanup paths.
pub(super) struct Lease {
    parent: BudgetPath,
    root: BudgetPath,
    owner: BudgetPath,
    registry: BudgetPath,
    registry_digest: [u8; 32],
    identity: identity::Identity,
    charge: disk::RecoveryCharge,
}

impl Lease {
    /// Publish the registration before the external child can exist.
    pub(super) fn register(staging: &BudgetPath, run_id: &str, parent: &Path) -> Result<Self> {
        let budget = staging.budget();
        let parent = canonical(parent, budget)?;
        let identity = identity::Identity::capture(&parent)?;
        let root = root_path(&parent, run_id)?;
        let owner = owner_path(&root, run_id)?;
        let registry = staging.try_join(REGISTRY)?;
        let encoded = encode(&parent, budget)?;
        ensure!(
            encoded.len() <= MAX_PATH_BYTES,
            "Workspace parent path is too long"
        );
        let mut header = [0; 44];
        header[..8].copy_from_slice(MAGIC);
        header[8..12].copy_from_slice(&u32::try_from(encoded.len())?.to_le_bytes());
        header[12..].copy_from_slice(identity.digest());
        let mut digest = Sha256::new();
        digest.update(header);
        digest.update(&encoded[..]);
        let checksum = digest.finalize();
        let mut temporary = disk::MetadataTemporary::create(staging, TEMPORARY, budget)?;
        let written = (|| {
            temporary.writer().write_all(&header)?;
            temporary.writer().write_all(&encoded)?;
            temporary.writer().write_all(&checksum)?;
            Ok(())
        })();
        temporary.finish(written)?;
        let registry_digest = sha256_bytes(temporary.path())?;
        let charge = temporary.publish(&registry, true)?;
        sync_directory(staging)?;
        Ok(Self {
            parent,
            root,
            owner,
            registry,
            registry_digest,
            identity,
            charge,
        })
    }

    /// Create only the registered child; no payload precedes its synced marker.
    pub(super) fn create(&self) -> Result<()> {
        self.identity.verify(&self.parent)?;
        fs::create_dir(&self.root)?;
        let owner = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&self.owner)?;
        owner.sync_all()?;
        drop(owner);
        sync_directory(&self.root)?;
        sync_directory(&self.parent)
    }

    /// Copy the bounded root for ingestion without modifying command-line args.
    pub(super) fn root(&self) -> Result<BudgetPath> {
        BudgetPath::from_path(&self.root, self.root.budget(), MemoryClass::RetainedResults)
    }

    /// Registry bytes whose allowance belongs to this transaction's ledger.
    pub(super) fn bytes(&self) -> u64 {
        self.charge.bytes()
    }

    /// Release only after the registry and external private data are gone.
    pub(super) fn release(&mut self) {
        self.charge.release();
    }

    /// Refuse edits before publication or live cleanup.
    pub(super) fn verify(&self) -> Result<()> {
        ensure!(
            sha256_bytes(&self.registry)? == self.registry_digest,
            "Modified live workspace registry"
        );
        self.identity.verify(&self.parent)
    }

    /// Normal/error cleanup needs no new allocation after store owners close.
    /// A nonempty root falls back to a fully inventoried recovery plan.
    pub(super) fn close_empty(&self, staging: &Path) -> Result<bool> {
        self.verify()?;
        match fs::symlink_metadata(&self.root) {
            Ok(metadata) => {
                ensure!(
                    metadata.file_type().is_dir(),
                    "Non-directory live workspace root"
                );
                for entry in fs::read_dir(&self.root)? {
                    if Some(entry?.file_name().as_os_str()) != self.owner.file_name() {
                        return Ok(false);
                    }
                }
                match fs::symlink_metadata(&self.owner) {
                    Ok(_) => {
                        marker(&self.owner)?;
                        sync_directory(&self.root)?;
                        fs::remove_file(&self.owner)?;
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                    Err(error) => {
                        return Err(error).context("Failed to inspect live workspace marker");
                    }
                }
                fs::remove_dir(&self.root)?;
                sync_directory(&self.parent)?;
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).context("Failed to inspect live workspace"),
        }
        fs::remove_file(&self.registry)?;
        sync_directory(staging)?;
        Ok(true)
    }
}

/// Exact removal paths for an interrupted external root, admitted before I/O.
pub(super) struct Cleanup {
    parent: BudgetPath,
    identity: identity::Identity,
    directories: BudgetVec<BudgetPath>,
    files: BudgetVec<BudgetPath>,
    owner: Option<BudgetPath>,
}

impl Cleanup {
    /// Decode only a bounded, checksummed registry in the granted staging scope.
    pub(super) fn read(staging: &BudgetPath, run_id: &str) -> Result<Option<Self>> {
        let budget = staging.budget();
        let registry = staging.try_join(REGISTRY)?;
        let metadata = match fs::symlink_metadata(&registry) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(error).context("Failed to inspect workspace registry"),
        };
        ensure!(
            metadata.file_type().is_file()
                && (76..=MAX_PATH_BYTES as u64 + 76).contains(&metadata.len()),
            "Invalid workspace registry size/type"
        );
        let length = usize::try_from(metadata.len())?;
        let mut bytes = BudgetVec::collect_in_class(
            std::iter::repeat_n(0u8, length),
            budget,
            MemoryClass::RetainedResults,
        )?;
        File::open(&registry)?.read_exact(&mut bytes)?;
        ensure!(
            &bytes[..8] == MAGIC,
            "Unsupported workspace registry version/path encoding; registration and private data preserved"
        );
        let path_bytes = usize::try_from(u32::from_le_bytes(bytes[8..12].try_into()?))?;
        ensure!(
            path_bytes.checked_add(76) == Some(length)
                && Sha256::digest(&bytes[..length - 32])[..] == bytes[length - 32..],
            "Damaged workspace registry"
        );
        let parent = decode(&bytes[44..length - 32], budget)?;
        ensure!(
            parent.is_absolute(),
            "Workspace parent must be canonical and absolute"
        );
        // Validate before even accepting a missing child: the pathname can
        // resolve to another host's local disk or a replacement mount/parent.
        let identity = identity::Identity::recorded(&parent, bytes[12..44].try_into()?)?;
        let root = root_path(&parent, run_id)?;
        let owner = owner_path(&root, run_id)?;
        let mut directories = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
        let mut files = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
        match fs::symlink_metadata(&root) {
            Ok(metadata) => ensure!(
                metadata.file_type().is_dir(),
                "Refusing non-directory workspace root"
            ),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return Ok(Some(Self {
                    parent,
                    identity,
                    directories,
                    files,
                    owner: None,
                }));
            }
            Err(error) => return Err(error).context("Failed to inspect interrupted workspace"),
        }
        ensure!(
            canonical(&parent, budget)?.as_ref() == AsRef::<Path>::as_ref(&parent),
            "Workspace parent changed; private data preserved"
        );
        let owner = match fs::symlink_metadata(&owner) {
            Ok(_) => {
                marker(&owner)?;
                Some(owner)
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                ensure!(
                    fs::read_dir(&root)?.next().transpose()?.is_none(),
                    "Unmarked nonempty workspace; private data preserved"
                );
                None
            }
            Err(error) => {
                return Err(error).context("Failed to inspect interrupted workspace owner");
            }
        };
        directories.push(root)?;
        let mut next = 0;
        while next < directories.len() {
            for entry in fs::read_dir(&directories[next])? {
                let path = directories[next].try_join(entry?.file_name())?;
                if owner
                    .as_ref()
                    .is_some_and(|owner| owner.as_ref() == AsRef::<Path>::as_ref(&path))
                {
                    continue;
                }
                let metadata = fs::symlink_metadata(&path)?;
                if metadata.file_type().is_dir() {
                    directories.push(path)?;
                } else {
                    ensure!(
                        metadata.file_type().is_file(),
                        "Refusing symbolic link or special file in workspace"
                    );
                    files.push(path)?;
                }
            }
            next += 1;
        }
        Ok(Some(Self {
            parent,
            identity,
            directories,
            files,
            owner,
        }))
    }

    /// Include every regular extent in the caller's common inode inventory.
    pub(super) fn files(&self) -> impl Iterator<Item = &BudgetPath> {
        self.files.iter().chain(self.owner.iter())
    }

    /// Remove payloads before revoking their external ownership marker.
    pub(super) fn discard(&self) -> Result<()> {
        self.identity.verify(&self.parent)?;
        for path in &self.files {
            fs::remove_file(path)?;
        }
        for path in self.directories.iter().skip(1).rev() {
            sync_directory(path)?;
            fs::remove_dir(path)?;
        }
        if let Some(root) = self.directories.first() {
            sync_directory(root)?;
            if let Some(owner) = &self.owner {
                fs::remove_file(owner)?;
            }
            fs::remove_dir(root)?;
            sync_directory(&self.parent)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};

    fn crash(transaction: &mut OutputTransaction) {
        // Model process death after closing handles, retaining durable paths.
        transaction.committed = true;
    }

    pub(super) fn budget(disk: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: disk,
            temporary_open_files: 128,
        })
        .unwrap()
    }

    #[test]
    fn namespace_replacement_or_unavailability_preserves_registration_until_original_storage_returns()
     {
        for replacement in [None, Some(false), Some(true)] {
            let out = tempfile::tempdir().unwrap();
            let base = tempfile::tempdir().unwrap();
            let parent = base.path().join("parent");
            let displaced = base.path().join("displaced");
            fs::create_dir(&parent).unwrap();
            fs::write(parent.join("foreign"), b"original foreign bytes").unwrap();
            let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
            let root = transaction.workspace_root(Some(&parent)).unwrap();
            fs::write(root.join("counts"), b"original private counts").unwrap();
            let staging = transaction.staging_directory.to_path_buf();
            let registration = fs::read(staging.join(REGISTRY)).unwrap();
            let basename = root.file_name().unwrap().to_owned();
            let owner = owner_path(&root, &transaction.run_id).unwrap();
            let owner_basename = owner.file_name().unwrap().to_owned();
            crash(&mut transaction);
            drop(transaction);
            fs::rename(&parent, &displaced).unwrap();
            if let Some(copy_child) = replacement {
                fs::create_dir(&parent).unwrap();
                fs::write(parent.join("foreign"), b"replacement foreign bytes").unwrap();
                if copy_child {
                    let foreign_root = parent.join(&basename);
                    fs::create_dir(&foreign_root).unwrap();
                    fs::write(foreign_root.join(&owner_basename), b"").unwrap();
                    fs::write(foreign_root.join("counts"), b"replacement private bytes").unwrap();
                }
            }
            let retry = OutputTransaction::begin(out.path(), "sample");
            assert!(
                retry.is_err(),
                "namespace change {replacement:?} was accepted"
            );
            assert_eq!(fs::read(staging.join(REGISTRY)).unwrap(), registration);
            assert_eq!(
                fs::read(displaced.join(&basename).join("counts")).unwrap(),
                b"original private counts"
            );
            if let Some(copy_child) = replacement {
                assert_eq!(
                    fs::read(parent.join("foreign")).unwrap(),
                    b"replacement foreign bytes"
                );
                if copy_child {
                    assert_eq!(
                        fs::read(parent.join(&basename).join("counts")).unwrap(),
                        b"replacement private bytes"
                    );
                }
                fs::remove_dir_all(&parent).unwrap();
            }
            fs::rename(&displaced, &parent).unwrap();
            let retry = OutputTransaction::begin(out.path(), "sample").unwrap();
            assert!(!root.exists());
            assert_eq!(
                fs::read(parent.join("foreign")).unwrap(),
                b"original foreign bytes"
            );
            drop(retry);
        }
    }

    #[test]
    fn interrupted_workspace_is_reclaimed_without_touching_its_parent_or_foreign_files() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        fs::write(base.path().join("foreign"), b"keep").unwrap();
        let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        let run = root.join(".sharkmer-count-run-example/shard-0");
        fs::create_dir_all(&run).unwrap();
        fs::write(run.join("counts"), b"interrupted count bytes").unwrap();
        fs::create_dir(root.join(".sharkmer-replay-run-example")).unwrap();
        fs::write(
            root.join(".sharkmer-replay-run-example/replay"),
            b"original records",
        )
        .unwrap();
        crash(&mut transaction);
        drop(transaction);
        let retry = OutputTransaction::begin(out.path(), "sample").unwrap();
        assert!(!root.exists());
        assert_eq!(fs::read(base.path().join("foreign")).unwrap(), b"keep");
        assert_eq!(fs::read_dir(base.path()).unwrap().count(), 1);
        drop(retry);
    }

    #[cfg(unix)]
    #[test]
    fn workspace_and_stage_share_one_inode_inventory_and_refuse_before_mutation() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        fs::write(root.join("counts"), [7; 2048]).unwrap();
        fs::hard_link(root.join("counts"), root.join("alias")).unwrap();
        let registry_bytes = fs::metadata(transaction.staging_directory.join(REGISTRY))
            .unwrap()
            .len();
        let required = 73 + registry_bytes + 2048;
        let staging = transaction.staging_directory.to_path_buf();
        crash(&mut transaction);
        drop(transaction);
        let ledger = budget(required);
        let competing = ledger.reserve_disk(1).unwrap();
        let error = OutputTransaction::begin_with_budget(out.path(), "sample", Some(&ledger))
            .err()
            .unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert_eq!(fs::read(root.join("counts")).unwrap(), [7; 2048]);
        assert_eq!(fs::read(root.join("alias")).unwrap(), [7; 2048]);
        assert_eq!(
            fs::metadata(staging.join(REGISTRY)).unwrap().len(),
            registry_bytes
        );
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 1);
        drop(competing);
        let retry =
            OutputTransaction::begin_with_budget(out.path(), "sample", Some(&ledger)).unwrap();
        assert!(!root.exists());
        assert_eq!(
            ledger.snapshot().peak_reserved_temporary_disk_bytes,
            required
        );
        drop(retry);
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
        assert!(ledger.snapshot().valid);
    }

    #[test]
    fn bootstrap_without_external_marker_reclaims_only_an_empty_registered_child() {
        for payload in [false, true] {
            let out = tempfile::tempdir().unwrap();
            let base = tempfile::tempdir().unwrap();
            let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
            let root = transaction.workspace_root(Some(base.path())).unwrap();
            let owner = owner_path(&root, &transaction.run_id).unwrap();
            fs::remove_file(owner).unwrap();
            if payload {
                fs::write(root.join("foreign"), b"keep").unwrap();
            }
            crash(&mut transaction);
            drop(transaction);
            let retry = OutputTransaction::begin(out.path(), "sample");
            if payload {
                assert!(retry.is_err());
                assert_eq!(fs::read(root.join("foreign")).unwrap(), b"keep");
            } else {
                drop(retry.unwrap());
                assert!(!root.exists());
            }
        }
    }

    #[cfg(unix)]
    #[test]
    fn workspace_cleanup_never_follows_a_symbolic_link() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let foreign = tempfile::tempdir().unwrap();
        fs::write(foreign.path().join("important"), b"keep").unwrap();
        let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        std::os::unix::fs::symlink(foreign.path(), root.join("linked-directory")).unwrap();
        crash(&mut transaction);
        drop(transaction);
        assert!(OutputTransaction::begin(out.path(), "sample").is_err());
        assert_eq!(fs::read(foreign.path().join("important")).unwrap(), b"keep");
        assert!(root.join("linked-directory").is_symlink());
    }

    #[cfg(unix)]
    #[test]
    fn registry_preserves_non_utf8_parent_paths() {
        use std::os::unix::ffi::OsStringExt;
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let parent = base
            .path()
            .join(std::ffi::OsString::from_vec(b"parent-\xff".to_vec()));
        fs::create_dir(&parent).unwrap();
        let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
        let root = transaction.workspace_root(Some(&parent)).unwrap();
        fs::write(root.join("counts"), b"payload").unwrap();
        crash(&mut transaction);
        drop(transaction);
        let retry = OutputTransaction::begin(out.path(), "sample").unwrap();
        assert!(!root.exists() && parent.is_dir());
        drop(retry);
    }

    #[test]
    fn changed_registry_refuses_before_deleting_any_external_data() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let mut transaction = OutputTransaction::begin(out.path(), "sample").unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        fs::write(root.join("counts"), b"keep").unwrap();
        let registry = transaction.staging_directory.join(REGISTRY);
        let mut bytes = fs::read(&registry).unwrap();
        bytes[12] ^= 1;
        fs::write(&registry, bytes).unwrap();
        crash(&mut transaction);
        drop(transaction);
        assert!(OutputTransaction::begin(out.path(), "sample").is_err());
        assert_eq!(fs::read(root.join("counts")).unwrap(), b"keep");
    }

    #[test]
    fn success_and_live_error_cleanup_remove_closed_workspaces_and_release_the_ledger() {
        for success in [false, true] {
            let out = tempfile::tempdir().unwrap();
            let base = tempfile::tempdir().unwrap();
            let ledger = budget(1 << 20);
            let mut transaction =
                OutputTransaction::begin_with_budget(out.path(), "sample", Some(&ledger)).unwrap();
            let root = transaction.workspace_root(Some(base.path())).unwrap();
            assert!(ledger.snapshot().reserved_temporary_disk_bytes > 73);
            transaction
                .stage_file("sample.stats.yaml", |file| {
                    file.write_all(b"stats")?;
                    Ok(())
                })
                .unwrap();
            if success {
                transaction.commit().unwrap();
            } else {
                drop(transaction);
            }
            assert!(!root.exists());
            drop(root);
            assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
            assert_eq!(Arc::strong_count(&ledger), 1);
            assert!(ledger.snapshot().valid);
        }
    }
    #[test]
    fn restored_generation_finishes_workspace_ownership_before_transaction_drop() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(out.path(), "sample").unwrap();
        first
            .stage_file("sample.stats.yaml", |file| {
                file.write_all(b"old stats")?;
                Ok(())
            })
            .unwrap();
        first.commit().unwrap();
        let original = fs::read(out.path().join("sample.manifest.yaml")).unwrap();
        let ledger = budget(1 << 20);
        let mut transaction =
            OutputTransaction::begin_with_budget(out.path(), "sample", Some(&ledger)).unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        transaction
            .stage_file("sample.stats.yaml", |file| {
                file.write_all(b"new stats")?;
                Ok(())
            })
            .unwrap();
        assert!(
            transaction
                .commit_inner(|step| {
                    ensure!(
                        step != PublicationStep::BackupReady,
                        "interrupt replacement"
                    );
                    Ok(())
                })
                .is_err()
        );
        transaction.recover_publication().unwrap();
        assert!(transaction.private_cleaned);
        assert!(!root.exists() && !transaction.staging_directory.exists());
        drop(transaction);
        drop(root);
        assert_eq!(
            fs::read(out.path().join("sample.manifest.yaml")).unwrap(),
            original
        );
        assert_eq!(
            fs::read(out.path().join("sample.stats.yaml")).unwrap(),
            b"old stats"
        );
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }
    #[test]
    fn accepted_cleanup_reuses_the_live_workspace_paths_under_memory_pressure() {
        let out = tempfile::tempdir().unwrap();
        let base = tempfile::tempdir().unwrap();
        let ledger = budget(1 << 20);
        let mut transaction =
            OutputTransaction::begin_with_budget(out.path(), "sample", Some(&ledger)).unwrap();
        let root = transaction.workspace_root(Some(base.path())).unwrap();
        transaction
            .stage_file("sample.stats.yaml", |file| {
                file.write_all(b"stats")?;
                Ok(())
            })
            .unwrap();
        let mut competing = None;
        transaction
            .commit_inner(|step| {
                if step == PublicationStep::CompleteManifest {
                    let current = ledger.snapshot().reserved_memory_bytes;
                    competing = Some(
                        ledger.reserve(MemoryClass::RetainedResults, (8 << 20) - current - 4096)?,
                    );
                }
                Ok(())
            })
            .unwrap();
        assert!(!root.exists() && !transaction.staging_directory.exists());
        assert!(transaction.workspace.is_none());
        assert_eq!(ledger.snapshot().reserved_temporary_disk_bytes, 0);
        drop(competing);
        drop(root);
        drop(transaction);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
    }
}
