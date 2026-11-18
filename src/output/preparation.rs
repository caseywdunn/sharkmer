//! Durable ownership of private data before a recovery journal exists.
//!
//! A directory name alone never authorizes recursive deletion. A synced empty
//! marker grants the private scope; an immutable context binds it to the public
//! generation that existed before any payload was written. Control files stay
//! until the payload is gone, so interrupted cleanup remains recoverable.

use super::*;

const OWNER: &str = ".sharkmer-owner-v1-";
const CONTEXT: &str = ".sharkmer-preparation-v1";
const TEMPORARY: &str = ".sharkmer-preparation-";
const ACCEPTED: &str = ".sharkmer-accepted-v1-";
const MAGIC: &[u8; 8] = b"SHKPREP1";
const CONTEXT_BYTES: usize = 73;

/// Live ownership paths are retained so error cleanup needs no new allocation.
pub(super) struct Lease {
    charge: disk::RecoveryCharge,
    owner: BudgetPath,
    context: BudgetPath,
    context_digest: [u8; 32],
}

impl Lease {
    /// Extent retained by this live control record.
    pub(super) fn bytes(&self) -> u64 {
        self.charge.bytes()
    }
    /// Release only after its private path has been removed.
    pub(super) fn release(&mut self) {
        self.charge.release();
    }

    /// Clean known live outputs without requiring another memory allowance.
    pub(super) fn discard(
        &self,
        staging: &BudgetPath,
        sample: &str,
        run_id: &str,
        outputs: &[StagedOutput],
    ) -> Result<()> {
        let mut simple = true;
        for entry in fs::read_dir(staging)? {
            let name = entry?.file_name();
            if Some(name.as_os_str()) != self.owner.file_name()
                && Some(name.as_os_str()) != self.context.file_name()
                && !outputs
                    .iter()
                    .any(|o| Some(name.as_os_str()) == o.staging_file.path().file_name())
            {
                simple = false;
            }
        }
        if !simple {
            return discard_owned(staging, sample, run_id);
        }
        marker(&self.owner)?;
        ensure!(
            regular(&self.context)?.len() == CONTEXT_BYTES as u64
                && sha256_bytes(&self.context)? == self.context_digest,
            "Modified live preparation context"
        );
        for output in outputs {
            match fs::remove_file(output.staging_file.path()) {
                Ok(()) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => return Err(error).context("Failed to remove live private output"),
            }
        }
        sync_directory(staging)?;
        fs::remove_file(&self.context)?;
        sync_directory(staging)?;
        fs::remove_file(&self.owner)?;
        fs::remove_dir(staging)?;
        Ok(())
    }
}

/// Create a durable scope before writing any payload or backup.
pub(super) fn create(
    directory: &BudgetPath,
    staging: &BudgetPath,
    run_id: &str,
    previous_manifest: &Path,
) -> Result<Lease> {
    let mut context = [0; CONTEXT_BYTES];
    context[..8].copy_from_slice(MAGIC);
    match fs::symlink_metadata(previous_manifest) {
        Ok(metadata) => {
            ensure!(
                metadata.file_type().is_file(),
                "Non-regular preparation manifest"
            );
            context[8] = 1;
            context[9..41].copy_from_slice(&sha256_bytes(previous_manifest)?);
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error).context("Failed to inspect preparation manifest"),
    }
    let checksum = Sha256::digest(&context[..41]);
    context[41..].copy_from_slice(&checksum);
    let owner = owner_path(staging, run_id)?;
    let destination = staging.try_join(CONTEXT)?;
    let reservation = staging
        .budget()
        .map(|budget| budget.reserve_disk(CONTEXT_BYTES as u64))
        .transpose()?;
    fs::create_dir(staging).context("Failed to create private output scope")?;
    let mut published = None;
    let result = (|| {
        create_marker(&owner)?;
        sync_directory(staging)?;
        sync_directory(directory)?;
        let mut temporary = disk::MetadataTemporary::create_reserved(
            staging,
            TEMPORARY,
            staging.budget(),
            reservation,
        )?;
        let written = temporary.writer().write_all(&context).map_err(Into::into);
        temporary.finish(written)?;
        published = Some(temporary.publish(&destination, true)?);
        sync_directory(staging)?;
        Ok(())
    })();
    if let Err(error) = result {
        // No payload has been granted yet. All known handles have closed.
        // Preserve the owner if a failed atomic operation left an unknown path.
        let removed_context = match fs::remove_file(&destination) {
            Ok(()) => true,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => true,
            Err(_) => false,
        };
        if removed_context {
            if let Some(charge) = &mut published {
                charge.release();
            }
            let only_owner = fs::read_dir(staging).is_ok_and(|entries| {
                entries.into_iter().all(|entry| {
                    entry
                        .is_ok_and(|entry| Some(entry.file_name().as_os_str()) == owner.file_name())
                })
            });
            if only_owner {
                let _ = fs::remove_file(&owner);
                let _ = fs::remove_dir(staging);
                let _ = sync_directory(directory);
            }
        }
        return Err(error);
    }
    Ok(Lease {
        charge: published.context("Missing preparation allowance")?,
        owner,
        context: destination,
        context_digest: Sha256::digest(context).into(),
    })
}

/// Preplan ownership before mutation, including validated historical journals.
pub(super) fn owner_path(staging: &BudgetPath, run_id: &str) -> Result<BudgetPath> {
    let name = BudgetString::format(
        format_args!("{OWNER}{run_id}"),
        staging.budget(),
        MemoryClass::RetainedResults,
    )?;
    staging.try_join(name.as_str())
}

/// The caller has validated the journal's entire private scope and admitted its
/// extents. Persist an empty owner before removing that journal; a later seal
/// binds payload cleanup to the accepted public manifest. Never replace a path
/// that appeared after planning.
pub(super) fn adopt_owner(staging: &Path, owner: &Path) -> Result<()> {
    create_marker(owner)?;
    sync_directory(staging)
}

fn create_marker(path: &Path) -> Result<()> {
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.sync_all()?;
    drop(file);
    Ok(())
}

fn regular(path: &Path) -> Result<fs::Metadata> {
    let metadata = fs::symlink_metadata(path)?;
    ensure!(
        metadata.file_type().is_file(),
        "Non-regular private ownership path {}",
        path.display()
    );
    Ok(metadata)
}

fn marker(path: &Path) -> Result<()> {
    ensure!(
        regular(path)?.len() == 0,
        "Modified private ownership marker {}",
        path.display()
    );
    Ok(())
}

fn generated_run_id(run_id: &str) -> bool {
    let mut parts = run_id.split('-');
    let hexadecimal = |part: &str| {
        !part.is_empty()
            && part
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    };
    let Some(time) = parts.next() else {
        return false;
    };
    let Some(pid) = parts.next() else {
        return false;
    };
    let Some(nonce) = parts.next() else {
        return false;
    };
    time.len() <= 32
        && pid.len() <= 8
        && nonce.len() == 16
        && [time, pid, nonce].into_iter().all(hexadecimal)
        && parts.next().is_none()
}

/// Immutable control paths in their removal order, with the owner last.
pub(super) struct Controls {
    context: Option<BudgetPath>,
    accepted: BudgetVec<(BudgetPath, [u8; 32])>,
    owner: BudgetPath,
    previous: Option<Option<[u8; 32]>>,
    registry: Option<BudgetPath>,
}

impl Controls {
    pub(super) fn read(staging: &BudgetPath, run_id: &str) -> Result<Option<Self>> {
        let owner = owner_path(staging, run_id)?;
        match fs::symlink_metadata(&owner) {
            Ok(_) => marker(&owner)?,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(error).context("Failed to inspect private owner"),
        }
        let path = staging.try_join(CONTEXT)?;
        let (context, previous) = match fs::symlink_metadata(&path) {
            Ok(_) => {
                ensure!(
                    regular(&path)?.len() == CONTEXT_BYTES as u64,
                    "Invalid preparation context size"
                );
                let mut bytes = [0; CONTEXT_BYTES];
                File::open(&path)?.read_exact(&mut bytes)?;
                ensure!(
                    &bytes[..8] == MAGIC
                        && bytes[8] <= 1
                        && Sha256::digest(&bytes[..41])[..] == bytes[41..],
                    "Invalid preparation context"
                );
                let digest: [u8; 32] = bytes[9..41].try_into()?;
                ensure!(
                    bytes[8] == 1 || digest == [0; 32],
                    "Invalid absent-manifest context"
                );
                (Some(path), Some((bytes[8] == 1).then_some(digest)))
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => (None, None),
            Err(error) => return Err(error).context("Failed to inspect preparation context"),
        };
        let mut accepted = BudgetVec::in_class(staging.budget(), MemoryClass::RetainedResults);
        for entry in fs::read_dir(staging)? {
            let entry = entry?;
            let name = entry.file_name();
            if let Some(digest) = name.to_str().and_then(|name| name.strip_prefix(ACCEPTED)) {
                ensure!(accepted.len() < 2, "Too many private acceptance markers");
                let mut decoded = [0; 32];
                ensure!(digest.len() == 64, "Invalid private acceptance marker");
                for (pair, byte) in digest
                    .as_bytes()
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .zip(&mut decoded)
                {
                    let digit = |c: u8| match c {
                        b'0'..=b'9' => Some(c - b'0'),
                        b'a'..=b'f' => Some(c - b'a' + 10),
                        _ => None,
                    };
                    *byte = (digit(pair[0]).context("Invalid acceptance digest")? << 4)
                        | digit(pair[1]).context("Invalid acceptance digest")?;
                }
                let path = staging.try_join(&name)?;
                marker(&path)?;
                accepted.push((path, decoded))?;
            }
        }
        let registry = staging.try_join(workspace::REGISTRY)?;
        let registry = match fs::symlink_metadata(&registry) {
            Ok(_) => {
                regular(&registry)?;
                Some(registry)
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
            Err(error) => {
                return Err(error).context("Failed to inspect private workspace registration");
            }
        };
        Ok(Some(Self {
            context,
            accepted,
            owner,
            previous,
            registry,
        }))
    }

    pub(super) fn contains(&self, path: &Path) -> bool {
        AsRef::<Path>::as_ref(&self.owner) == path
            || self.context.as_deref() == Some(path)
            || self.registry.as_deref() == Some(path)
            || self
                .accepted
                .iter()
                .any(|(p, _)| AsRef::<Path>::as_ref(p) == path)
    }

    pub(super) fn append(self, files: &mut BudgetVec<BudgetPath>) -> Result<()> {
        if let Some(path) = self.registry {
            files.push(path)?;
        }
        if let Some(path) = self.context {
            files.push(path)?;
        }
        for (path, _) in self.accepted {
            files.push(path)?;
        }
        files.push(self.owner)
    }

    /// Allocate the acceptance path before mutation; creating it is zero bytes.
    pub(super) fn acceptance_path(staging: &BudgetPath, digest: &[u8; 32]) -> Result<BudgetPath> {
        let mut name = BudgetString::with_capacity(
            ACCEPTED.len() + 64,
            staging.budget(),
            MemoryClass::RetainedResults,
        )?;
        name.push_str(ACCEPTED)?;
        name.push_hex(digest)?;
        staging.try_join(name.as_str())
    }
}

/// Persist acceptance before removing a journal or the first private payload.
pub(super) fn seal(staging: &Path, acceptance: &Path) -> Result<()> {
    match fs::symlink_metadata(acceptance) {
        Ok(_) => marker(acceptance)?,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => create_marker(acceptance)?,
        Err(error) => return Err(error).context("Failed to inspect private acceptance marker"),
    }
    sync_directory(staging)
}

pub(super) fn temporary(name: &str) -> bool {
    [
        TEMPORARY,
        ".sharkmer-manifest-",
        ".sharkmer-recovery-",
        ".sharkmer-workspace-",
    ]
    .into_iter()
    .any(|prefix| {
        name.strip_prefix(prefix).is_some_and(|suffix| {
            !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_alphanumeric())
        })
    })
}

/// Restore copies live inside the granted backup scope, including after an
/// earlier interrupted recovery. Root files with this spelling are not owned.
pub(super) fn restore_temporary_name(name: &str) -> bool {
    name.strip_prefix(".sharkmer-restore-")
        .is_some_and(|suffix| {
            !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_alphanumeric())
        })
}

/// A fixed-depth removal plan; all paths are allocated before mutation.
struct ScopeCleanup {
    staging: BudgetPath,
    backup: Option<BudgetPath>,
    payload: BudgetVec<BudgetPath>,
    controls: BudgetVec<BudgetPath>,
    workspace: Option<workspace::Cleanup>,
}

impl ScopeCleanup {
    fn read(
        staging: BudgetPath,
        sample: &str,
        controls: Controls,
        payload_granted: bool,
        run_id: &str,
    ) -> Result<Self> {
        let mut payload = BudgetVec::in_class(staging.budget(), MemoryClass::RetainedResults);
        let mut backup = None;
        ensure!(
            payload_granted || controls.registry.is_none(),
            "Workspace registry has no payload grant"
        );
        let workspace = workspace::Cleanup::read(&staging, run_id)?;
        for entry in fs::read_dir(&staging)? {
            let name = entry?.file_name();
            let path = staging.try_join(&name)?;
            if controls.contains(&path) {
                continue;
            }
            let name = name.to_str().context("Non-UTF8 private path")?;
            if name == ".previous" && payload_granted {
                ensure!(
                    fs::symlink_metadata(&path)?.file_type().is_dir(),
                    "Non-directory preparation backup"
                );
                for entry in fs::read_dir(&path)? {
                    let name = entry?.file_name();
                    let name = name.to_str().context("Non-UTF8 private backup path")?;
                    ensure!(
                        name == "manifest.yaml"
                            || restore_temporary_name(name)
                            || ensure_owned_output_basename(sample, name).is_ok(),
                        "Unexpected private backup path"
                    );
                    let path = path.try_join(name)?;
                    regular(&path)?;
                    payload.push(path)?;
                }
                backup = Some(path);
            } else {
                ensure!(
                    temporary(name)
                        || (payload_granted && ensure_owned_output_basename(sample, name).is_ok()),
                    "Unregistered private preparation path"
                );
                regular(&path)?;
                payload.push(path)?;
            }
        }
        let mut control_files = BudgetVec::in_class(staging.budget(), MemoryClass::RetainedResults);
        controls.append(&mut control_files)?;
        Ok(Self {
            staging,
            backup,
            payload,
            controls: control_files,
            workspace,
        })
    }

    fn files(&self) -> impl Iterator<Item = &BudgetPath> {
        self.payload
            .iter()
            .chain(self.controls.iter())
            .chain(self.workspace.iter().flat_map(workspace::Cleanup::files))
    }

    fn discard(self) -> Result<()> {
        if let Some(workspace) = &self.workspace {
            workspace.discard()?;
        }
        for path in self.payload {
            fs::remove_file(&path)?;
        }
        if let Some(backup) = self.backup {
            fs::remove_dir(&backup)?;
        }
        sync_directory(&self.staging)?;
        if let Some((owner, controls)) = self.controls.split_last() {
            for path in controls {
                fs::remove_file(path)?;
            }
            sync_directory(&self.staging)?;
            fs::remove_file(owner)?;
        }
        fs::remove_dir(&self.staging)?;
        Ok(())
    }
}

/// Reclaim granted scopes bound to the unchanged or receipted public generation.
/// Unregistered nonempty legacy directories and root temporaries are preserved.
pub(super) fn reclaim(
    directory: &BudgetPath,
    sample: &str,
    manifest_path: &BudgetPath,
    current: Option<&OutputManifest>,
) -> Result<()> {
    let prefix = BudgetString::format(
        format_args!(".sharkmer-stage-{sample}-"),
        directory.budget(),
        MemoryClass::RetainedResults,
    )?;
    let current_digest = current.map(|_| sha256_bytes(manifest_path)).transpose()?;
    let mut plans = BudgetVec::in_class(directory.budget(), MemoryClass::RetainedResults);
    let mut inventory = recovery::storage::Inventory::private(directory.budget());
    for entry in fs::read_dir(directory)? {
        let entry = entry?;
        let name = entry.file_name();
        let Some(run_id) = name.to_str().and_then(|n| n.strip_prefix(prefix.as_str())) else {
            continue;
        };
        let generated = generated_run_id(run_id);
        // Older validated journals can carry a different hexadecimal run-ID
        // spelling. Only an explicit owner grants those scopes; unmarked
        // nonstandard names remain untouched, including empty directories.
        if !generated
            && (run_id.is_empty()
                || !run_id
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit() || byte == b'-'))
        {
            continue;
        }
        let staging = directory.try_join(&name)?;
        if !fs::symlink_metadata(&staging)?.file_type().is_dir() {
            continue;
        }
        let recorded = current.is_some_and(|m| {
            m.status != "complete" && m.staging_directory.as_deref() == name.to_str()
        });
        let plan = if let Some(controls) = Controls::read(&staging, run_id)? {
            let accepted = current_digest.is_some_and(|digest| {
                controls
                    .accepted
                    .iter()
                    .any(|(_, expected)| *expected == digest)
            });
            if accepted {
                ensure!(
                    current.is_some_and(|m| m.status == "complete"),
                    "Incomplete accepted generation; private scope preserved"
                );
                for receipt in &current.context("Missing accepted generation")?.files {
                    recovery::verify_file(directory, receipt)?;
                }
            } else if !recorded {
                if let Some(previous) = controls.previous {
                    ensure!(
                        current_digest == previous,
                        "Public generation changed; preparation scope preserved"
                    );
                } else {
                    ensure!(
                        controls.accepted.is_empty(),
                        "Modified accepted generation; private scope preserved"
                    );
                }
            }
            let granted = recorded || controls.previous.is_some() || accepted;
            ScopeCleanup::read(staging, sample, controls, granted, run_id)?
        } else {
            if !generated || recorded || fs::read_dir(&staging)?.next().transpose()?.is_some() {
                continue;
            }
            ScopeCleanup {
                staging,
                backup: None,
                payload: BudgetVec::in_class(directory.budget(), MemoryClass::RetainedResults),
                controls: BudgetVec::in_class(directory.budget(), MemoryClass::RetainedResults),
                workspace: None,
            }
        };
        for path in plan.files() {
            inventory.include(path)?;
        }
        plans.push(plan)?;
    }
    let mut charge = inventory.admit(0, 0, 0, directory.budget())?;
    for plan in plans {
        plan.discard()?;
    }
    sync_directory(directory)?;
    charge.release();
    Ok(())
}

/// A live transaction already owns these extents; do not adopt them twice.
pub(super) fn discard_owned(staging: &BudgetPath, sample: &str, run_id: &str) -> Result<()> {
    let controls = Controls::read(staging, run_id)?.context("Missing private output owner")?;
    let granted = controls.previous.is_some() || !controls.accepted.is_empty();
    ScopeCleanup::read(
        BudgetPath::from_path(staging, staging.budget(), MemoryClass::RetainedResults)?,
        sample,
        controls,
        granted,
        run_id,
    )?
    .discard()
}

/// Remove an accepted first generation after all staged names were published.
pub(super) fn finish(staging: &BudgetPath, run_id: &str, digest: &[u8; 32]) -> Result<()> {
    let controls = Controls::read(staging, run_id)?.context("Missing private output owner")?;
    ensure!(
        controls
            .accepted
            .iter()
            .any(|(_, expected)| expected == digest),
        "Missing accepted output anchor"
    );
    let workspace = workspace::Cleanup::read(staging, run_id)?;
    for entry in fs::read_dir(staging)? {
        let path = staging.try_join(entry?.file_name())?;
        ensure!(
            controls.contains(&path),
            "Unpublished data in accepted private scope"
        );
    }
    let mut files = BudgetVec::in_class(staging.budget(), MemoryClass::RetainedResults);
    controls.append(&mut files)?;
    if let Some(workspace) = workspace {
        workspace.discard()?;
    }
    sync_directory(staging)?;
    let (owner, controls) = files
        .split_last()
        .context("Missing accepted private owner")?;
    for path in controls {
        fs::remove_file(path)?;
    }
    sync_directory(staging)?;
    fs::remove_file(owner)?;
    fs::remove_dir(staging)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};

    const RUN: &str = "1234-56-0123456789abcdef";

    fn path(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> BudgetPath {
        BudgetPath::from_path(path, budget, MemoryClass::RetainedResults).unwrap()
    }

    fn setup(directory: &Path) -> (BudgetPath, BudgetPath, BudgetPath) {
        let directory = path(directory, None);
        let staging = directory
            .try_join(format!(".sharkmer-stage-sample-{RUN}"))
            .unwrap();
        let manifest = directory.try_join("sample.manifest.yaml").unwrap();
        (directory, staging, manifest)
    }

    fn complete(directory: &Path) {
        let mut transaction = OutputTransaction::begin(directory, "sample").unwrap();
        for name in ["sample_a.fasta", "sample.stats.yaml"] {
            transaction
                .stage_file(name, |file| {
                    file.write_all(b"accepted bytes")?;
                    Ok(())
                })
                .unwrap();
        }
        transaction.commit().unwrap();
    }

    #[test]
    fn fresh_retry_reclaims_premanifest_data_and_partial_backups() {
        for prior_complete in [false, true] {
            let root = tempfile::tempdir().unwrap();
            if prior_complete {
                complete(root.path());
            }
            let (directory, staging, manifest) = setup(root.path());
            create(&directory, &staging, RUN, &manifest).unwrap();
            fs::write(staging.join("sample_a.fasta"), b"partial sequence").unwrap();
            fs::create_dir(staging.join(".previous")).unwrap();
            fs::write(staging.join(".previous/manifest.yaml"), b"partial backup").unwrap();
            let previous = fs::read(&manifest).ok();
            let next = OutputTransaction::begin(root.path(), "sample").unwrap();
            assert!(!staging.exists());
            if prior_complete {
                assert_eq!(fs::read(&manifest).unwrap(), previous.unwrap());
            }
            drop(next);
        }
    }

    #[test]
    fn bootstrap_retry_reclaims_only_empty_or_explicitly_granted_scopes() {
        for boundary in 0..3 {
            let root = tempfile::tempdir().unwrap();
            let (directory, staging, manifest) = setup(root.path());
            fs::create_dir(&staging).unwrap();
            if boundary >= 1 {
                create_marker(&owner_path(&staging, RUN).unwrap()).unwrap();
            }
            if boundary >= 2 {
                fs::write(
                    staging.join(".sharkmer-preparation-Abc123"),
                    b"partial context",
                )
                .unwrap();
            }
            reclaim(&directory, "sample", &manifest, None).unwrap();
            assert!(!staging.exists());
        }
        let root = tempfile::tempdir().unwrap();
        let (directory, staging, manifest) = setup(root.path());
        fs::create_dir(&staging).unwrap();
        fs::write(staging.join("user-data"), b"keep").unwrap();
        fs::write(
            directory.join(".sharkmer-restore-Abc123"),
            b"legacy evidence",
        )
        .unwrap();
        reclaim(&directory, "sample", &manifest, None).unwrap();
        assert_eq!(fs::read(staging.join("user-data")).unwrap(), b"keep");
        assert_eq!(
            fs::read(directory.join(".sharkmer-restore-Abc123")).unwrap(),
            b"legacy evidence"
        );
    }

    #[test]
    fn damaged_context_or_ungranted_payload_is_preserved() {
        for corrupted in [false, true] {
            let root = tempfile::tempdir().unwrap();
            let (directory, staging, manifest) = setup(root.path());
            create(&directory, &staging, RUN, &manifest).unwrap();
            fs::write(staging.join("sample_a.fasta"), b"evidence").unwrap();
            if corrupted {
                fs::write(staging.join(CONTEXT), [0; CONTEXT_BYTES]).unwrap();
            } else {
                fs::remove_file(staging.join(CONTEXT)).unwrap();
            }
            assert!(reclaim(&directory, "sample", &manifest, None).is_err());
            assert_eq!(
                fs::read(staging.join("sample_a.fasta")).unwrap(),
                b"evidence"
            );
            assert!(owner_path(&staging, RUN).unwrap().exists());
        }
    }

    #[test]
    fn refusal_admits_all_scopes_before_deleting_any_evidence() {
        let root = tempfile::tempdir().unwrap();
        let (directory, staging, manifest) = setup(root.path());
        create(&directory, &staging, RUN, &manifest).unwrap();
        fs::write(staging.join("sample_a.fasta"), [7; 100]).unwrap();
        let other = directory
            .try_join(".sharkmer-stage-sample-1234-56-fedcba9876543210")
            .unwrap();
        create(&directory, &other, "1234-56-fedcba9876543210", &manifest).unwrap();
        fs::write(other.join("sample_a.fasta"), [8; 100]).unwrap();
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1 << 20,
            process_overhead_bytes: 0,
            class_bytes: [1 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: 345,
            temporary_open_files: 128,
        })
        .unwrap();
        let limited = path(root.path(), Some(&budget));
        let error = reclaim(&limited, "sample", &manifest, None).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert_eq!(fs::read(staging.join("sample_a.fasta")).unwrap(), [7; 100]);
        assert_eq!(fs::read(other.join("sample_a.fasta")).unwrap(), [8; 100]);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        drop(limited);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        let unlimited = path(root.path(), None);
        reclaim(&unlimited, "sample", &manifest, None).unwrap();
        assert!(!staging.exists() && !other.exists());
    }

    #[test]
    fn accepted_anchor_survives_each_private_cleanup_boundary() {
        for boundary in 0..7 {
            let root = tempfile::tempdir().unwrap();
            let (directory, staging, manifest) = setup(root.path());
            create(&directory, &staging, RUN, &manifest).unwrap();
            complete(root.path());
            // complete's startup legitimately reclaimed the old empty scope.
            create(&directory, &staging, RUN, &manifest).unwrap();
            let anchor =
                Controls::acceptance_path(&staging, &sha256_bytes(&manifest).unwrap()).unwrap();
            seal(&staging, &anchor).unwrap();
            let backup = staging.join(".previous");
            fs::create_dir(&backup).unwrap();
            fs::write(backup.join("manifest.yaml"), b"obsolete backup").unwrap();
            fs::write(
                backup.join(".sharkmer-restore-Ab12Cd"),
                b"earlier interrupted restore",
            )
            .unwrap();
            fs::write(staging.join("sample_a.fasta"), b"obsolete private link").unwrap();
            if boundary >= 1 {
                fs::remove_file(staging.join("sample_a.fasta")).unwrap();
            }
            if boundary >= 2 {
                fs::remove_file(backup.join("manifest.yaml")).unwrap();
            }
            if boundary >= 3 {
                fs::remove_file(backup.join(".sharkmer-restore-Ab12Cd")).unwrap();
            }
            if boundary >= 4 {
                fs::remove_dir(&backup).unwrap();
            }
            if boundary >= 5 {
                fs::remove_file(staging.join(CONTEXT)).unwrap();
            }
            if boundary >= 6 {
                fs::remove_file(&anchor).unwrap();
            }
            let current = read_previous_manifest(&manifest, "sample")
                .unwrap()
                .unwrap();
            reclaim(&directory, "sample", &manifest, Some(&current)).unwrap();
            assert!(!staging.exists());
            assert_eq!(
                fs::read(directory.join("sample_a.fasta")).unwrap(),
                b"accepted bytes"
            );
        }
    }

    #[test]
    fn accepted_producer_scopes_allow_historical_hex_run_ids_and_preserve_unmarked_names() {
        for run_id in ["a1", "A1-b2", RUN] {
            let root = tempfile::tempdir().unwrap();
            complete(root.path());
            let directory = path(root.path(), None);
            let manifest = directory.try_join("sample.manifest.yaml").unwrap();
            let staging = directory
                .try_join(format!(".sharkmer-stage-sample-{run_id}"))
                .unwrap();
            fs::create_dir(&staging).unwrap();
            let owner = owner_path(&staging, run_id).unwrap();
            adopt_owner(&staging, &owner).unwrap();
            let anchor =
                Controls::acceptance_path(&staging, &sha256_bytes(&manifest).unwrap()).unwrap();
            seal(&staging, &anchor).unwrap();
            fs::write(staging.join("sample_a.fasta"), b"obsolete private payload").unwrap();
            let unmarked = directory.join(".sharkmer-stage-sample-a2");
            fs::create_dir(&unmarked).unwrap();
            fs::write(unmarked.join("foreign"), b"keep").unwrap();
            let current = read_previous_manifest(&manifest, "sample")
                .unwrap()
                .unwrap();
            reclaim(&directory, "sample", &manifest, Some(&current)).unwrap();
            assert!(!staging.exists());
            assert_eq!(
                fs::read(directory.join("sample_a.fasta")).unwrap(),
                b"accepted bytes"
            );
            assert_eq!(fs::read(unmarked.join("foreign")).unwrap(), b"keep");
        }
    }

    #[test]
    fn edited_accepted_manifest_preserves_private_cleanup_evidence() {
        let root = tempfile::tempdir().unwrap();
        complete(root.path());
        let (directory, staging, manifest) = setup(root.path());
        create(&directory, &staging, RUN, &manifest).unwrap();
        let anchor =
            Controls::acceptance_path(&staging, &sha256_bytes(&manifest).unwrap()).unwrap();
        seal(&staging, &anchor).unwrap();
        fs::write(staging.join("sample_a.fasta"), b"evidence").unwrap();
        OpenOptions::new()
            .append(true)
            .open(&manifest)
            .unwrap()
            .write_all(b"\n")
            .unwrap();
        let current = read_previous_manifest(&manifest, "sample")
            .unwrap()
            .unwrap();
        assert!(reclaim(&directory, "sample", &manifest, Some(&current)).is_err());
        assert_eq!(
            fs::read(staging.join("sample_a.fasta")).unwrap(),
            b"evidence"
        );
    }

    #[cfg(unix)]
    #[test]
    fn granted_scope_never_follows_symbolic_links() {
        let root = tempfile::tempdir().unwrap();
        let outside = tempfile::tempdir().unwrap();
        let (directory, staging, manifest) = setup(root.path());
        create(&directory, &staging, RUN, &manifest).unwrap();
        let foreign = outside.path().join("evidence");
        fs::write(&foreign, b"keep").unwrap();
        std::os::unix::fs::symlink(&foreign, staging.join("sample_a.fasta")).unwrap();
        assert!(reclaim(&directory, "sample", &manifest, None).is_err());
        assert_eq!(fs::read(foreign).unwrap(), b"keep");
        assert!(owner_path(&staging, RUN).unwrap().exists());
    }
}
