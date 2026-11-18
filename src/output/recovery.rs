//! Durable rollback for replacement of a completed output generation.
//!
//! The sample lock is held by the caller throughout preparation and recovery.
//! Backups are independent copies, so editing a public file cannot edit its
//! recovery copy through a shared inode. The journal is synced before the
//! first public change, and removed durably before its backups are discarded.

use super::*;

mod copy;
use copy::CopyWorkspace;
mod accepted;
pub(super) use accepted::accept;
pub(super) mod storage;

const RECOVERY_SCHEMA: u32 = 1;
const BACKUP_DIRECTORY: &str = ".previous";
const BACKUP_MANIFEST: &str = "manifest.yaml";

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RecoveryJournal {
    pub(super) schema_version: u32,
    pub(super) producer: String,
    pub(super) sample: String,
    pub(super) run_id: String,
    pub(super) staging_directory: String,
    pub(super) previous_manifest_sha256: String,
    pub(super) backup_files: Vec<OutputReceipt>,
    pub(super) replacement_files: Vec<OutputReceipt>,
}

/// Immutable new-journal view; historical parsing remains a separate owner.
#[derive(Serialize)]
pub(crate) struct RecoveryJournalView<'a> {
    schema_version: u32,
    producer: &'a str,
    sample: &'a str,
    run_id: &'a str,
    staging_directory: &'a str,
    previous_manifest_sha256: &'a str,
    backup_files: &'a [&'a OutputReceipt],
    replacement_files: ReceiptList<'a>,
}

pub(super) enum RecoveryOutcome {
    Restored,
    Accepted,
}

fn journal_path(directory: &BudgetPath, sample: &str) -> Result<BudgetPath> {
    let name = BudgetString::format(
        format_args!(".sharkmer-recovery-{sample}.yaml"),
        directory.budget(),
        MemoryClass::RetainedResults,
    )?;
    directory.try_join(name.as_str())
}

pub(super) fn prepare(transaction: &mut OutputTransaction) -> Result<()> {
    let mut copy = CopyWorkspace::new(transaction.budget.as_ref())?;
    let previous = transaction
        .previous_complete
        .as_ref()
        .context("Missing previous generation")?;
    let expected_manifest = transaction
        .previous_manifest_sha256
        .as_ref()
        .context("Missing previous manifest receipt")?;
    let backup = transaction.staging_directory.try_join(BACKUP_DIRECTORY)?;
    let backup_manifest = backup.try_join(BACKUP_MANIFEST)?;
    // Admit the complete backup before creating its directory or any copy.
    // Public input remains outside the temporary allowance. Copies are capped
    // at this admitted extent even if a source changes during preparation.
    require_regular(&transaction.manifest_path)?;
    let mut remaining = fs::metadata(&transaction.manifest_path)?.len();
    for receipt in &previous.files {
        let source = transaction
            .output_directory
            .try_join(receipt.path.as_str())?;
        match fs::symlink_metadata(&source) {
            Ok(metadata) => {
                ensure!(
                    metadata.file_type().is_file(),
                    "Refusing non-regular prior output"
                );
                remaining = remaining
                    .checked_add(metadata.len())
                    .context("Output recovery backup length overflow")?;
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).context("Failed to inspect prior backup extent"),
        }
    }
    let mut charge = disk::RecoveryCharge::reserve(remaining, transaction.budget.as_ref())
        .context("Unable to admit output recovery backup bytes")?;
    if let Err(error) = fs::create_dir(&backup) {
        charge.release();
        return Err(error).context("Failed to create output recovery directory");
    }
    transaction.backup_disk = charge;
    remaining -= copy_synced(
        &transaction.manifest_path,
        &backup_manifest,
        remaining,
        &mut copy,
    )?;
    ensure!(
        sha256_file_owned(&backup_manifest, transaction.budget.as_ref())? == *expected_manifest,
        "Prior output manifest changed while preparing recovery"
    );
    let mut backup_files =
        BudgetVec::in_class(transaction.budget.as_ref(), MemoryClass::RetainedResults);
    let mut inventory = backup_inventory();
    for receipt in &previous.files {
        let source = transaction
            .output_directory
            .try_join(receipt.path.as_str())?;
        match fs::symlink_metadata(&source) {
            Ok(_) => {
                verify_file(&transaction.output_directory, receipt)?;
                remaining -= copy_synced(
                    &source,
                    &backup.try_join(receipt.path.as_str())?,
                    remaining,
                    &mut copy,
                )?;
                verify_file(&backup, receipt)?;
                inventory_entry(&mut inventory, receipt)?;
                backup_files.push(receipt)?;
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).context("Failed to inspect prior output for backup"),
        }
    }
    sync_directory(&backup)?;
    sync_directory(&transaction.staging_directory)?;
    transaction.backup_inventory_sha256 = Some(inventory.finalize().into());
    drop(copy);
    let journal = RecoveryJournalView {
        schema_version: RECOVERY_SCHEMA,
        producer: MANIFEST_PRODUCER,
        sample: &transaction.sample,
        run_id: &transaction.run_id,
        staging_directory: &transaction.staging_basename,
        previous_manifest_sha256: expected_manifest,
        backup_files: &backup_files,
        replacement_files: ReceiptList::Staged(&transaction.staged_outputs),
    };
    let yaml = crate::stats::yaml::prepare_recovery(&journal, transaction.budget.as_ref())?;
    let journal_basename = BudgetString::format(
        format_args!(".sharkmer-recovery-{}.yaml", transaction.sample),
        transaction.budget.as_ref(),
        MemoryClass::RetainedResults,
    )?;
    let journal_path = transaction
        .output_directory
        .try_join(journal_basename.as_str())?;
    let mut temporary = disk::MetadataTemporary::create(
        &transaction.staging_directory,
        ".sharkmer-recovery-",
        transaction.budget.as_ref(),
    )?;
    let written = yaml.write(temporary.writer());
    temporary
        .finish(written)
        .context("Failed to write output recovery receipt")?;
    transaction.recovery_journal_sha256 = Some(sha256_bytes(temporary.path())?);
    transaction.journal_disk = temporary.publish(&journal_path, true).context(
        "Failed to publish output recovery receipt without overwriting an existing file",
    )?;
    sync_directory(&transaction.staging_directory)?;
    sync_directory(&transaction.output_directory)
}

fn backup_inventory() -> Sha256 {
    let mut digest = Sha256::new();
    digest.update(b"sharkmer-recovery-backups-v1\0");
    digest
}

fn inventory_entry(digest: &mut Sha256, receipt: &OutputReceipt) -> Result<()> {
    digest.update(u64::try_from(receipt.path.len())?.to_le_bytes());
    digest.update(receipt.path.as_bytes());
    digest.update(receipt.sha256.as_bytes());
    Ok(())
}

/// Recover before reading the public manifest or cleaning any incomplete run.
/// Refuse the entire recovery before mutation if any relevant file has changed.
pub(super) fn recover(directory: &BudgetPath, sample: &str) -> Result<Option<RecoveryOutcome>> {
    recover_with_credit(directory, sample, 0, 0)
}

/// The credit is owned by this transaction's live staged/backup guards in the
/// same ledger. Fresh-process recovery must always pass zero.
pub(super) fn recover_with_credit(
    directory: &BudgetPath,
    sample: &str,
    credit: u64,
    journal_credit: u64,
) -> Result<Option<RecoveryOutcome>> {
    let budget = directory.budget();
    let path = journal_path(directory, sample)?;
    match fs::symlink_metadata(&path) {
        Ok(metadata) => ensure!(
            metadata.file_type().is_file(),
            "Refusing non-regular output recovery receipt {}",
            path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error).context("Failed to inspect output recovery receipt"),
    }
    let journal = parsing::recovery(&path)
        .context("Failed to parse output recovery receipt; backups preserved")?;
    ensure!(
        journal.schema_version == RECOVERY_SCHEMA
            && journal.producer == MANIFEST_PRODUCER
            && journal.sample == sample,
        "Refusing foreign or unsupported output recovery receipt"
    );
    ensure!(
        !journal.run_id.is_empty()
            && journal
                .run_id
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() || byte == b'-'),
        "Invalid recovery run ID"
    );
    ensure_safe_basename(&journal.staging_directory)?;
    ensure!(
        journal.staging_directory.as_str()
            == BudgetString::format(
                format_args!(".sharkmer-stage-{sample}-{}", journal.run_id),
                budget,
                MemoryClass::RetainedResults,
            )?
            .as_str(),
        "Invalid recovery staging directory"
    );
    let staging = directory.try_join(journal.staging_directory.as_str())?;
    require_directory(&staging)?;
    let backup = staging.try_join(BACKUP_DIRECTORY)?;
    require_directory(&backup)?;
    let backup_manifest = backup.try_join(BACKUP_MANIFEST)?;
    require_regular(&backup_manifest)?;
    ensure!(
        sha256_file_owned(&backup_manifest, budget)?.as_str()
            == journal.previous_manifest_sha256.as_str(),
        "Modified output recovery manifest; backups preserved"
    );
    let previous =
        read_previous_manifest(&backup_manifest, sample)?.context("Missing recovery manifest")?;
    ensure!(
        previous.status == "complete",
        "Recovery manifest must describe a completed generation"
    );
    // This validates the full historical receipt dictionary even for files that
    // were already missing when the rerun began.
    verify_previous_outputs(&backup, &previous)?;
    let old = receipt_map(sample, &journal.backup_files, budget)?;
    let new = receipt_map(sample, &journal.replacement_files, budget)?;
    ensure!(
        new.contains_key(
            BudgetString::format(
                format_args!("{sample}.stats.yaml"),
                budget,
                MemoryClass::RetainedResults,
            )?
            .as_str()
        ),
        "Recovery replacement is missing its stats receipt"
    );
    for receipt in &journal.backup_files {
        ensure!(
            previous.files.contains(receipt),
            "Recovery backup receipt is absent from the prior manifest"
        );
        verify_file(&backup, receipt)?;
    }
    // Backups must cover every physically present prior file, not an arbitrary
    // subset of the manifest's receipts.
    for receipt in &previous.files {
        if fs::symlink_metadata(backup.try_join(receipt.path.as_str())?).is_ok() {
            ensure!(
                old.contains_key(receipt.path.as_str()),
                "Unlisted prior backup {}",
                receipt.path
            );
        }
    }
    let manifest_basename = BudgetString::format(
        format_args!("{sample}.manifest.yaml"),
        budget,
        MemoryClass::RetainedResults,
    )?;
    let manifest_path = directory.try_join(manifest_basename.as_str())?;
    require_regular(&manifest_path)?;
    let public_sha = sha256_file_owned(&manifest_path, budget)?;
    let current = read_previous_manifest(&manifest_path, sample)?
        .context("Missing public manifest during recovery")?;
    if public_sha.as_str() != journal.previous_manifest_sha256.as_str() {
        ensure!(
            current.run_id == journal.run_id
                && matches!(current.status.as_str(), "publishing" | "complete")
                && current.files == journal.replacement_files,
            "Refusing modified output manifest during recovery; backups preserved"
        );
    }
    if public_sha.as_str() != journal.previous_manifest_sha256.as_str()
        && current.status == "complete"
    {
        // Completion is the acceptance point. A crash during subsequent cleanup
        // must retain the new generation, never restore the previous one.
        for receipt in &journal.replacement_files {
            verify_file(directory, receipt)?;
        }
        let cleanup = storage::Cleanup::read(
            &staging,
            &backup,
            &old,
            &new,
            &journal.run_id,
            &sha256_bytes(&manifest_path)?,
        )?;
        let inventory = storage::Inventory::read(&path, &cleanup)?;
        let mut adopted_disk = inventory.admit(0, credit, journal_credit, budget)?;
        sync_directory(directory)?;
        cleanup.discard(directory, &path, &staging, &backup)?;
        adopted_disk.release();
        return Ok(Some(RecoveryOutcome::Accepted));
    }

    // Validate the complete affected namespace before removing even one file.
    // Missing entries are recoverable; edits, symlinks and foreign files are not.
    let mut plan = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
    for name in old
        .keys()
        .chain(new.keys().filter(|name| !old.contains_key(name)))
    {
        let destination = directory.try_join(name)?;
        let actual = match fs::symlink_metadata(&destination) {
            Ok(_) => {
                require_regular(&destination)?;
                let actual = sha256_file_owned(&destination, budget)?;
                ensure!(
                    old.get(name)
                        .is_some_and(|expected| expected == actual.as_str())
                        || new
                            .get(name)
                            .is_some_and(|expected| expected == actual.as_str()),
                    "Refusing modified or unowned output during recovery: {}; backups preserved",
                    destination.display()
                );
                Some(actual)
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
            Err(error) => return Err(error).context("Failed to inspect output during recovery"),
        };
        let already_previous = old
            .get(name)
            .is_some_and(|expected| actual.as_deref() == Some(expected));
        let restore_source = if !already_previous && old.contains_key(name) {
            let source = backup.try_join(name)?;
            let extent = CopyExtent {
                bytes: fs::metadata(&source)?.len(),
                sha256: old.get(name).context("Missing prior output receipt")?,
            };
            Some(RestoreSource { source, extent })
        } else {
            None
        };
        // Own every destination, restore path and plan slot before any public
        // mutation. Only the historical journal strings are borrowed here.
        plan.push(RecoveryFile {
            name,
            destination,
            restore_source,
            present: actual.is_some(),
            already_previous,
        })?;
    }
    // Preserve the original BTreeMap's lexical mutation order without nodes or
    // sort scratch. Each basename occurs exactly once in the merged namespace.
    plan.sort_unstable_by(|left, right| left.name.cmp(right.name));
    let restore_manifest_extent =
        if public_sha.as_str() != journal.previous_manifest_sha256.as_str() {
            Some(CopyExtent {
                bytes: fs::metadata(&backup_manifest)?.len(),
                sha256: journal.previous_manifest_sha256.as_str(),
            })
        } else {
            None
        };
    let restore_bytes = plan.iter().try_fold(
        restore_manifest_extent
            .as_ref()
            .map_or(0, |extent| extent.bytes),
        |bytes, file| {
            bytes
                .checked_add(
                    file.restore_source
                        .as_ref()
                        .map_or(0, |source| source.extent.bytes),
                )
                .context("Output restoration length overflow")
        },
    )?;
    // Admit the shared copy buffer before deleting or restoring any public file.
    let mut copy = CopyWorkspace::new(budget)?;
    let cleanup = storage::Cleanup::read(
        &staging,
        &backup,
        &old,
        &new,
        &journal.run_id,
        &sha256_bytes(&backup_manifest)?,
    )?;
    let mut inventory = storage::Inventory::read(&path, &cleanup)?;
    for file in &plan {
        // An old hash can also be a restored copy whose private name is gone.
        // Until the old complete manifest is durable, conservatively retain
        // every affected public payload in a publishing generation.
        if file.present
            && (!file.already_previous
                || public_sha.as_str() != journal.previous_manifest_sha256.as_str())
        {
            inventory.include(&file.destination)?;
        }
    }
    // This is the last admission before mutation: a refusal drops all memory
    // owners and leaves both the filesystem and disk ledger unchanged.
    let mut recovery_disk = inventory.admit(restore_bytes, credit, journal_credit, budget)?;
    for file in &plan {
        if file.already_previous {
            continue;
        }
        if file.present {
            fs::remove_file(&file.destination).with_context(|| {
                format!(
                    "Failed to remove interrupted replacement {}; backups preserved",
                    file.name
                )
            })?;
        }
        if let Some(source) = &file.restore_source {
            restore_no_clobber(source, &file.destination, &mut copy)?;
        }
    }
    sync_directory(directory)?;
    if let Some(extent) = &restore_manifest_extent {
        restore_manifest(
            &backup_manifest,
            &manifest_path,
            directory,
            &mut copy,
            extent,
        )?;
    }
    // Restored files become persistent only when their complete prior manifest
    // is durable. Check the accepted bytes before releasing the copy allowance.
    for file in &plan {
        if let Some(expected) = old.get(file.name) {
            require_regular(&file.destination)?;
            ensure!(
                digest_matches(&file.destination, expected)?,
                "Modified restored output {}; backups preserved",
                file.destination.display()
            );
        }
    }
    ensure!(
        digest_matches(&manifest_path, journal.previous_manifest_sha256.as_str())?,
        "Modified restored output manifest; backups preserved"
    );
    sync_directory(directory)?;
    drop(copy);
    cleanup.discard(directory, &path, &staging, &backup)?;
    recovery_disk.release();
    Ok(Some(RecoveryOutcome::Restored))
}

/// A sorted borrowed receipt dictionary with owned, precharged backing.
struct ReceiptIndex<'a> {
    entries: BudgetVec<(&'a str, &'a str)>,
}

impl<'a> ReceiptIndex<'a> {
    fn get(&self, name: &str) -> Option<&'a str> {
        self.entries
            .binary_search_by(|(key, _)| key.cmp(&name))
            .ok()
            .map(|index| self.entries[index].1)
    }

    fn contains_key(&self, name: &str) -> bool {
        self.get(name).is_some()
    }

    fn keys(&self) -> impl Iterator<Item = &'a str> + '_ {
        self.entries.iter().map(|(name, _)| *name)
    }
}

/// Immutable recovery decision and admitted paths survive through restoration.
struct RecoveryFile<'a> {
    name: &'a str,
    destination: BudgetPath,
    restore_source: Option<RestoreSource<'a>>,
    present: bool,
    already_previous: bool,
}

struct RestoreSource<'a> {
    source: BudgetPath,
    extent: CopyExtent<'a>,
}

struct CopyExtent<'a> {
    bytes: u64,
    sha256: &'a str,
}

/// Field order closes the writer before unlinking its path on every error.
/// Successful copies transfer only a closed path to the publication helpers.
struct RestoreTemporary {
    file: File,
    path: tempfile::TempPath,
}

fn digest_matches(path: &Path, expected: &str) -> Result<bool> {
    let digest = sha256_bytes(path)?;
    let mut hex = [0; 64];
    for (index, byte) in digest.into_iter().enumerate() {
        hex[2 * index] = b"0123456789abcdef"[usize::from(byte >> 4)];
        hex[2 * index + 1] = b"0123456789abcdef"[usize::from(byte & 15)];
    }
    Ok(hex == expected.as_bytes())
}

fn receipt_map<'a>(
    sample: &str,
    receipts: &'a [OutputReceipt],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ReceiptIndex<'a>> {
    receipt_index(
        sample,
        receipts
            .iter()
            .map(|receipt| (receipt.path.as_str(), receipt.sha256.as_str())),
        receipts.len(),
        budget,
    )
}

fn receipt_index<'a>(
    sample: &str,
    receipts: impl Iterator<Item = (&'a str, &'a str)>,
    length: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ReceiptIndex<'a>> {
    let mut entries =
        BudgetVec::with_capacity_in_class(length, budget, MemoryClass::RetainedResults)?;
    for (path, sha256) in receipts {
        ensure_owned_output_basename(sample, path)?;
        ensure!(
            sha256.len() == 64
                && sha256
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()),
            "Invalid recovery output hash"
        );
        entries.push((path, sha256))?;
    }
    entries.sort_unstable_by(|left, right| left.0.cmp(right.0));
    ensure!(
        entries.windows(2).all(|pair| pair[0].0 != pair[1].0),
        "Duplicate recovery output receipt"
    );
    Ok(ReceiptIndex { entries })
}

fn require_regular(path: &Path) -> Result<()> {
    ensure!(
        fs::symlink_metadata(path)?.file_type().is_file(),
        "Refusing non-regular recovery path {}",
        path.display()
    );
    Ok(())
}

fn require_directory(path: &Path) -> Result<()> {
    ensure!(
        fs::symlink_metadata(path)?.file_type().is_dir(),
        "Refusing non-directory recovery path {}",
        path.display()
    );
    Ok(())
}

pub(super) fn verify_file(directory: &BudgetPath, receipt: &OutputReceipt) -> Result<()> {
    let path = directory.try_join(receipt.path.as_str())?;
    require_regular(&path)?;
    ensure!(
        sha256_file_owned(&path, directory.budget())?.as_str() == receipt.sha256.as_str(),
        "Modified or incomplete recovery file {}; backups preserved",
        path.display()
    );
    Ok(())
}

fn copy_synced(
    source: &Path,
    destination: &Path,
    allowance: u64,
    copy: &mut CopyWorkspace,
) -> Result<u64> {
    require_regular(source)?;
    let mut reader = File::open(source)?;
    let length = reader.metadata()?.len();
    ensure!(
        length <= allowance,
        "Prior output grew beyond its admitted backup extent"
    );
    let mut writer = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(destination)?;
    copy.copy(&mut Read::by_ref(&mut reader).take(length), &mut writer)?;
    ensure!(
        writer.metadata()?.len() == length && reader.metadata()?.len() == length,
        "Prior output extent changed during backup"
    );
    writer.flush()?;
    writer.sync_all()?;
    Ok(length)
}

fn restore_temporary(
    source: &Path,
    copy: &mut CopyWorkspace,
    extent: &CopyExtent<'_>,
) -> Result<tempfile::TempPath> {
    require_regular(source)?;
    let mut reader = File::open(source)?;
    ensure!(
        reader.metadata()?.len() == extent.bytes,
        "Recovery source extent changed before restoration; backups preserved"
    );
    // The journal owns the source's backup directory through its staging
    // receipt. A killed copy remains there for restart cleanup to reclaim.
    let directory = source
        .parent()
        .context("Recovery source has no parent directory")?;
    let (file, path) = tempfile::Builder::new()
        .prefix(".sharkmer-restore-")
        .tempfile_in(directory)?
        .into_parts();
    let mut temporary = RestoreTemporary { file, path };
    copy.copy(
        &mut Read::by_ref(&mut reader).take(extent.bytes),
        &mut temporary.file,
    )?;
    ensure!(
        temporary.file.metadata()?.len() == extent.bytes
            && reader.metadata()?.len() == extent.bytes,
        "Recovery source extent changed during restoration; backups preserved"
    );
    temporary.file.flush()?;
    temporary.file.sync_all()?;
    let RestoreTemporary { file, path } = temporary;
    drop(file);
    ensure!(
        digest_matches(&path, extent.sha256)?,
        "Recovery source content changed during restoration; backups preserved"
    );
    Ok(path)
}

fn restore_no_clobber(
    source: &RestoreSource<'_>,
    destination: &Path,
    copy: &mut CopyWorkspace,
) -> Result<()> {
    restore_temporary(&source.source, copy, &source.extent)?
        .persist_noclobber(destination)
        .map_err(|error| error.error)
        .context(
            "Failed to restore prior output without overwriting a concurrently created file",
        )?;
    Ok(())
}

fn restore_manifest(
    source: &Path,
    destination: &Path,
    directory: &Path,
    copy: &mut CopyWorkspace,
    extent: &CopyExtent<'_>,
) -> Result<()> {
    restore_temporary(source, copy, extent)?
        .persist(destination)
        .map_err(|error| error.error)
        .context("Failed to restore prior output manifest")?;
    sync_directory(directory)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restoration_rejects_changed_extents_hashes_and_foreign_destinations() {
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("backup");
        let destination = directory.path().join("public");
        fs::write(&source, b"abc").unwrap();
        fs::write(&destination, b"unowned data").unwrap();
        let digest = sha256_file(&source).unwrap();
        let wrong = "f".repeat(64);
        let mut copy = CopyWorkspace::new(None).unwrap();
        for (bytes, sha256) in [
            (2, digest.as_str()),
            (3, wrong.as_str()),
            (3, digest.as_str()),
        ] {
            let restore = RestoreSource {
                source: BudgetPath::from_path(&source, None, MemoryClass::RetainedResults).unwrap(),
                extent: CopyExtent { bytes, sha256 },
            };
            assert!(restore_no_clobber(&restore, &destination, &mut copy).is_err());
            assert_eq!(fs::read(&destination).unwrap(), b"unowned data");
            assert_eq!(fs::read(&source).unwrap(), b"abc");
            assert!(!fs::read_dir(directory.path()).unwrap().any(|entry| {
                entry
                    .unwrap()
                    .file_name()
                    .to_string_lossy()
                    .starts_with(".sharkmer-restore-")
            }));
        }
        let extent = CopyExtent {
            bytes: 3,
            sha256: &digest,
        };
        let temporary = restore_temporary(&source, &mut copy, &extent).unwrap();
        assert_eq!(fs::read(&temporary).unwrap(), b"abc");
        let path = temporary.to_path_buf();
        drop(temporary);
        assert!(!path.exists());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn restored_copy_is_closed_before_path_publication() {
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("backup");
        fs::write(&source, b"abc").unwrap();
        let digest = sha256_file(&source).unwrap();
        let extent = CopyExtent {
            bytes: 3,
            sha256: &digest,
        };
        let mut copy = CopyWorkspace::new(None).unwrap();
        let temporary = restore_temporary(&source, &mut copy, &extent).unwrap();
        let expected: &Path = temporary.as_ref();
        for entry in fs::read_dir("/proc/self/fd").unwrap() {
            if let Ok(path) = fs::read_link(entry.unwrap().path()) {
                assert_ne!(path.as_path(), expected);
            }
        }
        let output = directory.path().join("accepted");
        temporary.persist_noclobber(&output).unwrap();
        assert_eq!(fs::read(output).unwrap(), b"abc");
    }

    #[test]
    fn restart_reclaims_unpublished_restore_copies_inside_receipted_staging() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        for (name, bytes) in [
            ("sample_gene.fasta", b"old fasta".as_slice()),
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
        let previous_manifest = fs::read(directory.path().join("sample.manifest.yaml")).unwrap();
        fs::write(directory.path().join("foreign.txt"), b"unowned data").unwrap();
        let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
        for (name, bytes) in [
            ("sample_gene.fasta", b"new fasta".as_slice()),
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
                    step != PublicationStep::OutputPublished,
                    "interrupt publication"
                );
                Ok(())
            })
            .is_err()
        );
        let staging = next.staging_directory.to_path_buf();
        next.committed = true;
        drop(next);
        let source = staging.join(".previous/sample_gene.fasta");
        let digest = sha256_file(&source).unwrap();
        let extent = CopyExtent {
            bytes: 9,
            sha256: &digest,
        };
        let mut copy = CopyWorkspace::new(None).unwrap();
        let temporary = restore_temporary(&source, &mut copy, &extent).unwrap();
        assert_eq!(temporary.parent(), source.parent());
        let orphan = temporary.keep().unwrap();
        assert!(orphan.starts_with(&staging));
        assert!(!directory.path().read_dir().unwrap().any(|entry| {
            entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".sharkmer-restore-")
        }));
        let root =
            BudgetPath::from_path(directory.path(), None, MemoryClass::RetainedResults).unwrap();
        assert!(matches!(
            recover(&root, "sample").unwrap(),
            Some(RecoveryOutcome::Restored)
        ));
        assert!(!orphan.exists() && !staging.exists());
        assert_eq!(
            fs::read(directory.path().join("sample.manifest.yaml")).unwrap(),
            previous_manifest
        );
        assert_eq!(
            fs::read(directory.path().join("sample_gene.fasta")).unwrap(),
            b"old fasta"
        );
        assert_eq!(
            fs::read(directory.path().join("sample.stats.yaml")).unwrap(),
            b"old stats"
        );
        assert_eq!(
            fs::read(directory.path().join("foreign.txt")).unwrap(),
            b"unowned data"
        );
    }

    #[test]
    fn backup_copy_refuses_an_extent_above_admission_before_creating_the_file() {
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("source");
        let destination = directory.path().join("backup");
        fs::write(&source, b"abc").unwrap();
        let mut copy = CopyWorkspace::new(None).unwrap();
        assert!(copy_synced(&source, &destination, 2, &mut copy).is_err());
        assert!(!destination.exists());
        assert_eq!(copy_synced(&source, &destination, 3, &mut copy).unwrap(), 3);
        assert_eq!(fs::read(&destination).unwrap(), b"abc");
    }

    fn index_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn receipt_index_preserves_order_aliases_and_rejects_invalid_entries() {
        let mut receipts: Vec<_> = (0..2000)
            .rev()
            .map(|i| OutputReceipt {
                path: format!("sample_λ_{i}.fasta"),
                sha256: format!("{i:064x}"),
            })
            .collect();
        let expected: std::collections::BTreeMap<_, _> = receipts
            .iter()
            .map(|r| (r.path.as_str(), r.sha256.as_str()))
            .collect();
        let budget = index_budget(1 << 20);
        let index = receipt_map("sample", &receipts, Some(&budget)).unwrap();
        assert!(index.keys().eq(expected.keys().copied()));
        for receipt in &receipts {
            let borrowed = index.get(&receipt.path).unwrap();
            assert_eq!(borrowed, receipt.sha256);
            assert_eq!(borrowed.as_ptr(), receipt.sha256.as_ptr());
        }
        assert!(index.get("sample_missing.fasta").is_none());
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            2000 * size_of::<(&str, &str)>() as u64
        );
        drop(index);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        let small = index_budget(1);
        let error = receipt_map("sample", &receipts, Some(&small))
            .err()
            .unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(small.snapshot().reserved_memory_bytes, 0);
        receipts.push(receipts[0].clone());
        assert!(
            receipt_map("sample", &receipts, Some(&budget))
                .err()
                .unwrap()
                .to_string()
                .contains("Duplicate")
        );
        receipts.pop();
        receipts[0].sha256 = "F".repeat(64);
        assert!(
            receipt_map("sample", &receipts, Some(&budget))
                .err()
                .unwrap()
                .to_string()
                .contains("hash")
        );
        receipts[0].sha256 = "f".repeat(64);
        receipts[0].path = "foreign.fasta".to_owned();
        assert!(receipt_map("sample", &receipts, Some(&budget)).is_err());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn receipt_index_reservations_cover_isolated_allocator_peaks() {
        const ENV: &str = "SHARKMER_RECOVERY_INDEX_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "output::recovery::tests::receipt_index_reservations_cover_isolated_allocator_peaks", "--nocapture", "--test-threads=1"])
                .env(ENV, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        for rows in [0, 1, 13, 2000] {
            let receipts: Vec<_> = (0..rows)
                .rev()
                .map(|i| OutputReceipt {
                    path: format!("sample_{i}.fasta"),
                    sha256: format!("{i:064x}"),
                })
                .collect();
            let budget = index_budget(1 << 20);
            let baseline = crate::PEAK_ALLOC.current_usage();
            crate::PEAK_ALLOC.reset_peak_usage();
            let index = receipt_map("sample", &receipts, Some(&budget)).unwrap();
            let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
            let reserved = budget.snapshot().peak_reserved_memory_bytes;
            assert!(
                peak as u64 <= reserved,
                "recovery index heap {peak} exceeded {reserved}"
            );
            drop(index);
            assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
            println!(
                "recovery_index_allocator rows={rows} peak_bytes={peak} reserved_bytes={reserved}"
            );
        }
    }

    #[test]
    fn borrowed_recovery_journal_preserves_prior_bytes_and_owned_reference_lifetimes() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        for (name, data) in [
            ("sample_gene.fasta", "old fasta"),
            ("sample.stats.yaml", "old stats"),
        ] {
            first
                .stage_file(name, |file| {
                    file.write_all(data.as_bytes())?;
                    Ok(())
                })
                .unwrap();
        }
        first.commit().unwrap();
        let budget = ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; crate::kmer::budget::MEMORY_CLASSES],
            // This fixture measures journal memory while writing real payloads.
            temporary_disk_bytes: 1 << 20,
            temporary_open_files: 0,
        })
        .unwrap();
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        for (name, data) in [
            ("sample_gene.fasta", "new fasta"),
            ("sample.stats.yaml", "new stats"),
        ] {
            next.stage_file(name, |file| {
                file.write_all(data.as_bytes())?;
                Ok(())
            })
            .unwrap();
        }
        next.staged_outputs
            .sort_unstable_by(|a, b| a.receipt.path.cmp(&b.receipt.path));
        let expected = RecoveryJournal {
            schema_version: RECOVERY_SCHEMA,
            producer: MANIFEST_PRODUCER.to_owned(),
            sample: next.sample.to_string(),
            run_id: next.run_id.to_string(),
            staging_directory: next.staging_basename.to_string(),
            previous_manifest_sha256: next.previous_manifest_sha256.as_ref().unwrap().to_string(),
            backup_files: next.previous_complete.as_ref().unwrap().files.clone(),
            replacement_files: next
                .staged_outputs
                .iter()
                .map(|s| OutputReceipt {
                    path: s.receipt.path.to_string(),
                    sha256: s.receipt.sha256.to_string(),
                })
                .collect(),
        };
        let before = budget.snapshot().reserved_memory_bytes;
        prepare(&mut next).unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        assert_eq!(
            fs::read(directory.path().join(".sharkmer-recovery-sample.yaml")).unwrap(),
            serde_yaml_ng::to_string(&expected).unwrap().as_bytes()
        );
        next.recovery_pending = true;
        next.recover_publication().unwrap();
        drop(next);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
            "old fasta"
        );
    }
}
