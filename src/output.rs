use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use crate::kmer::budget::path::BudgetPath;
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use anyhow::{Context, Result, bail, ensure};
use fs2::FileExt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Component, Path};
use std::sync::Arc;

mod disk;
pub(crate) mod parsing;
mod preparation;
mod recovery;
mod workspace;

pub(crate) use disk::{MetadataTemporary, RecoveryCharge, StagedFile};
pub(crate) use recovery::RecoveryJournalView;
use recovery::RecoveryOutcome;
pub(crate) use recovery::storage::Inventory as StorageInventory;

const MANIFEST_SCHEMA_VERSION: u32 = 4;
const MANIFEST_PRODUCER: &str = "sharkmer";

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub(crate) struct OutputReceipt {
    pub(crate) path: String,
    pub(crate) sha256: String,
}

#[derive(Debug, Deserialize, Serialize)]
pub(crate) struct OutputManifest {
    schema_version: u32,
    producer: String,
    sample: String,
    run_id: String,
    status: String,
    files: Vec<OutputReceipt>,
    #[serde(skip_serializing_if = "Option::is_none")]
    failure_reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    staging_directory: Option<String>,
}

#[derive(Serialize)]
struct OwnedReceipt {
    path: BudgetString,
    sha256: BudgetString,
}

pub(crate) struct StagedOutput {
    receipt: OwnedReceipt,
    staging_file: StagedFile,
}

/// Immutable views serialize the existing receipt owners without copying text.
pub(crate) enum ReceiptList<'a> {
    Historical(&'a [OutputReceipt]),
    Staged(&'a [StagedOutput]),
    Empty,
}
impl Serialize for ReceiptList<'_> {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;
        match self {
            Self::Historical(receipts) => receipts.serialize(serializer),
            Self::Empty => serializer.serialize_seq(Some(0))?.end(),
            Self::Staged(outputs) => {
                let mut sequence = serializer.serialize_seq(Some(outputs.len()))?;
                for staged in *outputs {
                    sequence.serialize_element(&staged.receipt)?;
                }
                sequence.end()
            }
        }
    }
}

#[derive(Serialize)]
pub(crate) struct ManifestView<'a> {
    schema_version: u32,
    producer: &'a str,
    sample: &'a str,
    run_id: &'a str,
    status: &'a str,
    files: ReceiptList<'a>,
    #[serde(skip_serializing_if = "Option::is_none")]
    failure_reason: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    staging_directory: Option<&'a str>,
}

mod manifest_model {
    pub trait Sealed {}
}
/// Only these immutable manifest models may enter the preflight serializer.
pub(crate) trait ManifestModel: Serialize + manifest_model::Sealed {}
impl manifest_model::Sealed for OutputManifest {}
impl ManifestModel for OutputManifest {}
impl manifest_model::Sealed for ManifestView<'_> {}
impl ManifestModel for ManifestView<'_> {}

pub(crate) struct OutputTransaction {
    output_directory: BudgetPath,
    sample: BudgetString,
    run_id: BudgetString,
    manifest_basename: BudgetString,
    manifest_path: BudgetPath,
    staging_basename: BudgetString,
    staging_directory: BudgetPath,
    staged_outputs: BudgetVec<StagedOutput>,
    backup_disk: disk::RecoveryCharge,
    journal_disk: disk::RecoveryCharge,
    preparation_disk: preparation::Lease,
    workspace: Option<workspace::Lease>,
    backup_inventory_sha256: Option<[u8; 32]>,
    recovery_journal_sha256: Option<[u8; 32]>,
    previous_incomplete: Option<parsing::Parsed<OutputManifest>>,
    previous_complete: Option<parsing::Parsed<OutputManifest>>,
    previous_manifest_sha256: Option<BudgetString>,
    replacement_receipts_active: bool,
    publication_started: bool,
    recovery_pending: bool,
    private_cleaned: bool,
    _lock_file: File,
    committed: bool,
    budget: Option<std::sync::Arc<crate::kmer::budget::ResourceBudget>>,
}

impl OutputTransaction {
    #[cfg(test)]
    pub(crate) fn begin(output_directory: &Path, sample: &str) -> Result<Self> {
        Self::begin_with_budget(output_directory, sample, None)
    }

    /// Borrow the run ledger from startup, before publication owner creation.
    pub(crate) fn begin_with_budget(
        output_directory: &Path,
        sample: &str,
        budget: Option<&std::sync::Arc<crate::kmer::budget::ResourceBudget>>,
    ) -> Result<Self> {
        let output_directory =
            BudgetPath::from_path(output_directory, budget, MemoryClass::RetainedResults)?;
        let sample = BudgetString::from_str(sample, budget, MemoryClass::RetainedResults)?;
        fs::create_dir_all(&output_directory).with_context(|| {
            format!(
                "Failed to create output directory {}",
                output_directory.display()
            )
        })?;
        let lock_basename = BudgetString::format(
            format_args!("{sample}.lock"),
            budget,
            MemoryClass::RetainedResults,
        )?;
        let lock_path = output_directory.try_join(lock_basename.as_str())?;
        reject_non_regular_existing_path(&lock_path, "output lock")?;
        let lock_file = OpenOptions::new()
            .create(true)
            .truncate(false)
            .read(true)
            .write(true)
            .open(&lock_path)
            .with_context(|| format!("Failed to open output lock {}", lock_path.display()))?;
        lock_file
            .lock_exclusive()
            .with_context(|| format!("Failed to lock output sample {sample}"))?;

        let manifest_basename = BudgetString::format(
            format_args!("{sample}.manifest.yaml"),
            budget,
            MemoryClass::RetainedResults,
        )?;
        let manifest_path = output_directory.try_join(manifest_basename.as_str())?;
        recovery::recover(&output_directory, &sample)?;
        let mut previous_manifest = read_previous_manifest(&manifest_path, &sample)?;
        if previous_manifest.is_none() {
            reject_unowned_stats_output(&output_directory, &sample, budget)?;
        }
        let cleanup = previous_manifest
            .as_ref()
            .map(|manifest| verify_previous_outputs(&output_directory, manifest))
            .transpose()?;
        if let Some(staging) = cleanup
            .as_ref()
            .and_then(|cleanup| cleanup.staging_directory.as_ref())
        {
            let backup = staging.try_join(".previous")?;
            match fs::symlink_metadata(&backup) {
                Ok(_) => bail!(
                    "Refusing incomplete replacement without its recovery receipt; backups preserved in {}",
                    staging.display()
                ),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => {
                    return Err(error).context("Failed to inspect incomplete replacement backups");
                }
            }
        }
        preparation::reclaim(
            &output_directory,
            &sample,
            &manifest_path,
            previous_manifest.as_deref(),
        )?;
        let preserve_complete = previous_manifest
            .as_ref()
            .is_some_and(|manifest| manifest.status == "complete");

        let previous_manifest_sha256 = preserve_complete
            .then(|| sha256_file_owned(&manifest_path, budget))
            .transpose()?;

        if !preserve_complete {
            if let Some(cleanup) = cleanup {
                cleanup.remove()?;
            }
            previous_manifest = None;
        }

        let run_id = new_run_id(budget)?;
        let staging_basename = BudgetString::format(
            format_args!(".sharkmer-stage-{sample}-{run_id}"),
            budget,
            MemoryClass::RetainedResults,
        )?;
        let staging_directory = output_directory.try_join(staging_basename.as_str())?;
        let preparation_disk = preparation::create(
            &output_directory,
            &staging_directory,
            &run_id,
            &manifest_path,
        )?;

        let transaction = Self {
            output_directory,
            sample,
            run_id,
            manifest_basename,
            manifest_path,
            staging_basename,
            staging_directory,
            staged_outputs: BudgetVec::in_class(budget, MemoryClass::RetainedResults),
            backup_disk: disk::RecoveryCharge::default(),
            journal_disk: disk::RecoveryCharge::default(),
            preparation_disk,
            workspace: None,
            backup_inventory_sha256: None,
            recovery_journal_sha256: None,
            previous_incomplete: if preserve_complete {
                None
            } else {
                previous_manifest.take()
            },
            previous_complete: if preserve_complete {
                previous_manifest
            } else {
                None
            },
            previous_manifest_sha256,
            replacement_receipts_active: false,
            publication_started: false,
            recovery_pending: false,
            private_cleaned: false,
            _lock_file: lock_file,
            committed: false,
            budget: budget.cloned(),
        };
        if !preserve_complete {
            transaction.write_manifest("in_progress", ReceiptList::Empty, None)?;
        }

        Ok(transaction)
    }

    /// Register the one private count/replay root before ingestion opens files.
    pub(crate) fn workspace_root(&mut self, parent: Option<&Path>) -> Result<BudgetPath> {
        ensure!(
            self.workspace.is_none(),
            "Run workspace was already registered"
        );
        let fallback;
        let parent = match parent {
            Some(parent) => parent,
            None => {
                fallback = std::env::temp_dir();
                &fallback
            }
        };
        self.workspace = Some(workspace::Lease::register(
            &self.staging_directory,
            &self.run_id,
            parent,
        )?);
        let lease = self
            .workspace
            .as_ref()
            .context("Missing registered workspace")?;
        let root = lease.root()?;
        lease.create()?;
        Ok(root)
    }

    fn release_preparation(&mut self) {
        self.preparation_disk.release();
        if let Some(workspace) = &mut self.workspace {
            workspace.release();
        }
    }

    pub(crate) fn run_id(&self) -> &str {
        &self.run_id
    }

    pub(crate) fn manifest_basename(&self) -> &str {
        &self.manifest_basename
    }

    pub(crate) fn stage_file(
        &mut self,
        basename: &str,
        write_file: impl FnOnce(&mut StagedFile) -> Result<()>,
    ) -> Result<()> {
        ensure_owned_output_basename(&self.sample, basename)?;
        ensure!(
            !self
                .staged_outputs
                .iter()
                .any(|staged| staged.receipt.path == basename),
            "Output file {} was staged more than once",
            basename
        );

        self.staged_outputs.reserve_one()?;
        let mut sha256 =
            BudgetString::with_capacity(64, self.budget.as_ref(), MemoryClass::RetainedResults)?;
        let staging_path = self.staging_directory.try_join(basename)?;
        let receipt_path =
            BudgetString::from_str(basename, self.budget.as_ref(), MemoryClass::RetainedResults)?;
        let mut file = StagedFile::create(staging_path, self.budget.as_ref())?;
        let written = write_file(&mut file)
            .with_context(|| format!("Failed to write staged output {basename}"));
        file.finish(written)?;
        sha256.push_hex(&sha256_bytes(file.path())?)?;
        self.staged_outputs.push(StagedOutput {
            receipt: OwnedReceipt {
                path: receipt_path,
                sha256,
            },
            staging_file: file,
        })?;
        Ok(())
    }

    pub(crate) fn commit(mut self) -> Result<()> {
        let result = self.commit_inner(|_| Ok(()));
        if result.is_err() && self.recovery_pending {
            self.recover_publication().with_context(|| {
                format!(
                    "Output publication failed: {:#}. Recovery remains pending; keep the recovery receipt and staging directory in {}",
                    result.as_ref().err().map(ToString::to_string).unwrap_or_default(),
                    self.output_directory.display()
                )
            })?;
        }
        result
    }

    fn commit_inner(
        &mut self,
        mut observe: impl FnMut(PublicationStep) -> Result<()>,
    ) -> Result<()> {
        if let Some(workspace) = &mut self.workspace {
            ensure!(
                workspace.close_empty(&self.staging_directory)?,
                "Count/replay owners must close before output commit"
            );
            workspace.release();
        }
        self.workspace = None;
        // Basenames are unique, so unstable sorting preserves the exact order
        // while avoiding a hidden allocation for sorting scratch.
        self.staged_outputs
            .sort_unstable_by(|left, right| left.receipt.path.cmp(&right.receipt.path));
        self.verify_replacement()?;
        for staged in &self.staged_outputs {
            let owned = self.previous_complete.as_ref().is_some_and(|manifest| {
                manifest
                    .files
                    .iter()
                    .any(|previous| previous.path.as_str() == staged.receipt.path.as_str())
            });
            if !owned {
                reject_existing_destination(
                    &self
                        .output_directory
                        .try_join(staged.receipt.path.as_str())?,
                    &staged.receipt.path,
                )?;
            }
        }

        if self.previous_complete.is_some() {
            // Set this before journal publication: an error syncing the journal
            // can still leave a durable recovery receipt behind.
            self.recovery_pending = true;
            recovery::prepare(self)?;
            observe(PublicationStep::BackupReady)?;
        }
        // Backup copying can take time. Recheck both receipts and the public
        // manifest before changing the accepted generation.
        let cleanup = self.verify_replacement()?;
        self.publication_started = true;
        self.write_manifest(
            "publishing",
            ReceiptList::Staged(&self.staged_outputs),
            None,
        )?;
        observe(PublicationStep::PublishingManifest)?;
        if let Some(cleanup) = cleanup {
            for path in cleanup.files {
                fs::remove_file(&path).with_context(|| {
                    format!("Failed to remove prior owned output {}", path.display())
                })?;
                observe(PublicationStep::PreviousFileRemoved)?;
            }
        }
        self.replacement_receipts_active = true;
        for staged in self.staged_outputs.iter_mut() {
            staged.staging_file.publishing();
            publish_no_clobber(
                staged.staging_file.path(),
                &self
                    .output_directory
                    .try_join(staged.receipt.path.as_str())?,
            )?;
            observe(PublicationStep::OutputPublished)?;
        }
        sync_directory(&self.output_directory)?;
        let accepted_manifest =
            self.write_manifest("complete", ReceiptList::Staged(&self.staged_outputs), None)?;
        self.committed = true;
        for staged in self.staged_outputs.iter_mut() {
            staged.staging_file.release();
        }
        observe(PublicationStep::CompleteManifest)?;
        if self.recovery_pending {
            // A durable complete manifest accepts the new generation. Cleanup
            // must never turn it into a failed generation or roll it back.
            if let Err(error) = recovery::accept(self, &accepted_manifest) {
                eprintln!(
                    "Warning: completed outputs retained; recovery cleanup pending: {error:#}"
                );
            } else {
                self.release_preparation();
                self.backup_disk.release();
                self.journal_disk.release();
                self.recovery_pending = false;
            }
        } else {
            preparation::finish(&self.staging_directory, &self.run_id, &accepted_manifest)?;
            self.release_preparation();
            sync_directory(&self.output_directory)?;
        }
        Ok(())
    }

    fn verify_replacement(&self) -> Result<Option<VerifiedCleanup>> {
        if let Some(expected) = &self.previous_manifest_sha256 {
            reject_non_regular_existing_path(&self.manifest_path, "output manifest")?;
            ensure!(
                sha256_file_owned(&self.manifest_path, self.budget.as_ref())? == *expected,
                "Refusing modified prior output manifest {}",
                self.manifest_path.display()
            );
        }
        self.previous_complete
            .as_ref()
            .map(|manifest| verify_previous_outputs(&self.output_directory, manifest))
            .transpose()
    }

    fn recover_publication(&mut self) -> Result<()> {
        let credit = self.staged_outputs.iter().try_fold(
            self.backup_disk
                .bytes()
                .checked_add(self.preparation_disk.bytes())
                .and_then(|bytes| {
                    bytes.checked_add(self.workspace.as_ref().map_or(0, workspace::Lease::bytes))
                })
                .context("Live preparation allowance overflow")?,
            |bytes, staged| {
                bytes
                    .checked_add(staged.staging_file.reserved_bytes())
                    .context("Live output disk allowance overflow")
            },
        )?;
        match recovery::recover_with_credit(
            &self.output_directory,
            &self.sample,
            credit,
            self.journal_disk.bytes(),
        )? {
            Some(RecoveryOutcome::Accepted) => {
                self.private_cleaned = true;
                self.committed = true;
                self.release_preparation();
                self.backup_disk.release();
                self.journal_disk.release();
                for staged in self.staged_outputs.iter_mut() {
                    staged.staging_file.release();
                }
            }
            Some(RecoveryOutcome::Restored) => {
                self.private_cleaned = true;
                self.publication_started = false;
                self.release_preparation();
                self.backup_disk.release();
                self.journal_disk.release();
                for staged in self.staged_outputs.iter_mut() {
                    staged.staging_file.release();
                }
            }
            None => {
                // No journal is safe only while the accepted prior generation
                // is unchanged. A removed receipt cannot authorize discarding
                // backups after publication has already changed public files.
                self.verify_replacement()?;
                self.publication_started = false;
                self.journal_disk.release();
            }
        }
        self.recovery_pending = false;
        Ok(())
    }

    fn write_manifest(
        &self,
        status: &str,
        files: ReceiptList<'_>,
        failure_reason: Option<&str>,
    ) -> Result<[u8; 32]> {
        let manifest = ManifestView {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER,
            sample: &self.sample,
            run_id: &self.run_id,
            status,
            files,
            failure_reason,
            staging_directory: (status != "complete").then_some(self.staging_basename.as_str()),
        };
        if self.previous_complete.is_some() && !self.publication_started && status == "failed" {
            let basename = BudgetString::format(
                format_args!(".sharkmer-failed-{}-{}.yaml", self.sample, self.run_id),
                self.budget.as_ref(),
                MemoryClass::RetainedResults,
            )?;
            let path = self.output_directory.try_join(basename.as_str())?;
            write_manifest_atomic_in(
                &self.staging_directory,
                &self.output_directory,
                &path,
                &manifest,
                true,
                false,
                self.budget.as_ref(),
            )
        } else {
            write_manifest_atomic_in(
                &self.staging_directory,
                &self.output_directory,
                &self.manifest_path,
                &manifest,
                false,
                status == "complete",
                self.budget.as_ref(),
            )
        }
    }
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        if self.committed || self.private_cleaned {
            // The transaction has handed its paths to accepted publication (or
            // a crash-state fixture). Private-file destructors must not reclaim
            // paths after this ownership transfer.
            for staged in self.staged_outputs.iter_mut() {
                staged.staging_file.publishing();
            }
            return;
        }
        if self.recovery_pending {
            if let Err(error) = self.recover_publication() {
                eprintln!("Warning: output recovery pending; backups preserved: {error:#}");
                for staged in self.staged_outputs.iter_mut() {
                    staged.staging_file.publishing();
                }
                return;
            }
            if self.committed || self.private_cleaned {
                return;
            }
        }
        let _ = self.write_manifest(
            "failed",
            if self.previous_complete.is_some() && !self.publication_started {
                ReceiptList::Empty
            } else if self.replacement_receipts_active {
                ReceiptList::Staged(&self.staged_outputs)
            } else {
                self.previous_incomplete
                    .as_ref()
                    .map_or(ReceiptList::Empty, |m| ReceiptList::Historical(&m.files))
            },
            Some("run ended before output publication completed"),
        );
        if let Some(workspace) = &mut self.workspace {
            match workspace.close_empty(&self.staging_directory) {
                Ok(true) => workspace.release(),
                Ok(false) => {}
                Err(error) => {
                    for staged in self.staged_outputs.iter_mut() {
                        staged.staging_file.publishing();
                    }
                    log::error!("Unable to verify private workspace; data preserved: {error:#}");
                    return;
                }
            }
        }
        match self.preparation_disk.discard(
            &self.staging_directory,
            &self.sample,
            &self.run_id,
            &self.staged_outputs,
        ) {
            Ok(()) => {
                self.backup_disk.release();
                self.release_preparation();
            }
            Err(error) => {
                for staged in self.staged_outputs.iter_mut() {
                    staged.staging_file.publishing();
                }
                log::error!(
                    "Unable to remove private output scope {}; its disk charge is retained: {error:#}",
                    self.staging_directory.display()
                );
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PublicationStep {
    BackupReady,
    PublishingManifest,
    PreviousFileRemoved,
    OutputPublished,
    CompleteManifest,
}

struct VerifiedCleanup {
    files: BudgetVec<BudgetPath>,
    staging_directory: Option<BudgetPath>,
}

impl VerifiedCleanup {
    fn remove(self) -> Result<()> {
        for path in self.files {
            fs::remove_file(&path).with_context(|| {
                format!("Failed to remove prior owned output {}", path.display())
            })?;
        }
        if let Some(path) = self.staging_directory {
            match fs::remove_dir_all(&path) {
                Ok(()) => Ok(()),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
                Err(error) => Err(error),
            }
            .with_context(|| {
                format!(
                    "Failed to remove prior staging directory {}",
                    path.display()
                )
            })?;
        }
        Ok(())
    }
}

fn read_previous_manifest(
    path: &BudgetPath,
    sample: &str,
) -> Result<Option<parsing::Parsed<OutputManifest>>> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("Failed to inspect output manifest {}", path.display()));
        }
    };
    ensure!(
        metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
        "Refusing non-regular output manifest {}",
        path.display()
    );
    let manifest = parsing::manifest(path)
        .with_context(|| format!("Failed to parse output manifest {}", path.display()))?;
    ensure!(
        matches!(manifest.schema_version, 1 | 2 | 3 | MANIFEST_SCHEMA_VERSION)
            && manifest.producer == MANIFEST_PRODUCER
            && manifest.sample == sample,
        "Refusing foreign or unsupported output manifest {}",
        path.display()
    );
    ensure!(
        matches!(
            manifest.status.as_str(),
            "in_progress" | "publishing" | "complete" | "failed"
        ),
        "Output manifest {} has unknown status {}",
        path.display(),
        manifest.status
    );
    ensure!(
        !manifest.run_id.is_empty()
            && manifest
                .run_id
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() || byte == b'-'),
        "Output manifest {} has an invalid run ID",
        path.display()
    );
    if manifest.status == "complete" {
        ensure!(
            manifest.staging_directory.is_none() && manifest.failure_reason.is_none(),
            "Complete output manifest {} contains incomplete-run metadata",
            path.display()
        );
    } else {
        ensure!(
            manifest.staging_directory.is_some(),
            "Noncomplete output manifest {} is missing its staging directory",
            path.display()
        );
    }
    Ok(Some(manifest))
}

fn verify_previous_outputs(
    output_directory: &BudgetPath,
    manifest: &OutputManifest,
) -> Result<VerifiedCleanup> {
    let budget = output_directory.budget();
    let mut files = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
    let mut seen = BudgetMap::in_class(budget, MemoryClass::RetainedResults);
    let stats_basename = BudgetString::format(
        format_args!("{}.stats.yaml", manifest.sample),
        budget,
        MemoryClass::RetainedResults,
    )?;
    if matches!(manifest.status.as_str(), "publishing" | "complete") {
        ensure!(
            manifest
                .files
                .iter()
                .filter(|receipt| receipt.path.as_str() == stats_basename.as_str())
                .count()
                == 1,
            "{} output manifest must contain exactly one stats receipt",
            manifest.status
        );
    }
    for receipt in &manifest.files {
        ensure_owned_output_basename(&manifest.sample, &receipt.path)?;
        ensure!(
            receipt.sha256.len() == 64
                && receipt
                    .sha256
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()),
            "Output manifest receipt for {} has an invalid SHA-256",
            receipt.path
        );
        ensure!(
            seen.get(receipt.path.as_str()).is_none(),
            "Output manifest contains duplicate receipt {}",
            receipt.path
        );
        seen.insert(receipt.path.as_str(), ())?;
        let path = output_directory.try_join(receipt.path.as_str())?;
        let metadata = match fs::symlink_metadata(&path) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("Failed to inspect prior output {}", path.display()));
            }
        };
        ensure!(
            metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
            "Refusing to replace non-regular prior output {}",
            path.display()
        );
        let actual_sha256 = sha256_file_owned(&path, budget)?;
        ensure!(
            actual_sha256.as_str() == receipt.sha256.as_str(),
            "Refusing to replace modified prior output {}",
            path.display()
        );
        files.push(path)?;
    }

    let staging_directory = match &manifest.staging_directory {
        Some(basename) => {
            ensure_safe_basename(basename)?;
            ensure!(
                basename.as_str()
                    == BudgetString::format(
                        format_args!(".sharkmer-stage-{}-{}", manifest.sample, manifest.run_id),
                        budget,
                        MemoryClass::RetainedResults,
                    )?
                    .as_str(),
                "Output manifest has an invalid staging directory {}",
                basename
            );
            let path = output_directory.try_join(basename.as_str())?;
            match fs::symlink_metadata(&path) {
                Ok(metadata) => {
                    ensure!(
                        metadata.file_type().is_dir() && !metadata.file_type().is_symlink(),
                        "Refusing non-directory prior staging path {}",
                        path.display()
                    );
                    Some(path)
                }
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!("Failed to inspect prior staging path {}", path.display())
                    });
                }
            }
        }
        None => None,
    };
    Ok(VerifiedCleanup {
        files,
        staging_directory,
    })
}

fn reject_non_regular_existing_path(path: &Path, label: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            ensure!(
                metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
                "Refusing non-regular {} {}",
                label,
                path.display()
            );
            Ok(())
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect {}", path.display())),
    }
}

fn reject_existing_destination(path: &Path, basename: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(_) => bail!(
            "Refusing to overwrite unowned or concurrently created output {}",
            basename
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect output {basename}")),
    }
}

fn reject_unowned_stats_output(
    output_directory: &BudgetPath,
    sample: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<()> {
    let basename = BudgetString::format(
        format_args!("{sample}.stats.yaml"),
        budget,
        MemoryClass::RetainedResults,
    )?;
    let path = output_directory.try_join(basename.as_str())?;
    match fs::symlink_metadata(&path) {
        Ok(_) => bail!(
            "Refusing existing unowned stats output {}. Use a fresh output directory or move the legacy file.",
            path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect output {basename}")),
    }
}

fn publish_no_clobber(staging_path: &Path, destination: &Path) -> Result<()> {
    fs::hard_link(staging_path, destination).with_context(|| {
        format!(
            "Failed to publish output {} without overwriting an existing file",
            destination.display()
        )
    })?;
    fs::remove_file(staging_path).with_context(|| {
        format!(
            "Failed to remove staged output after publishing {}",
            destination.display()
        )
    })?;
    Ok(())
}

#[cfg(test)]
fn write_manifest_atomic(
    output_directory: &Path,
    manifest_path: &Path,
    manifest: &OutputManifest,
) -> Result<()> {
    write_manifest_atomic_impl(output_directory, manifest_path, manifest, false, None).map(|_| ())
}

#[cfg(test)]
fn write_manifest_atomic_impl<M: ManifestModel>(
    output_directory: &Path,
    manifest_path: &Path,
    manifest: &M,
    no_clobber: bool,
    budget: Option<&std::sync::Arc<crate::kmer::budget::ResourceBudget>>,
) -> Result<[u8; 32]> {
    write_manifest_atomic_in(
        output_directory,
        output_directory,
        manifest_path,
        manifest,
        no_clobber,
        false,
        budget,
    )
}

fn write_manifest_atomic_in<M: ManifestModel>(
    temporary_directory: &Path,
    output_directory: &Path,
    manifest_path: &Path,
    manifest: &M,
    no_clobber: bool,
    seal_complete: bool,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<[u8; 32]> {
    let yaml = crate::stats::yaml::prepare_manifest(manifest, budget)?;
    let mut temporary =
        disk::MetadataTemporary::create(temporary_directory, ".sharkmer-manifest-", budget)
            .context("Failed to create temporary output manifest")?;
    let written = yaml.write(temporary.writer());
    temporary
        .finish(written)
        .context("Failed to write output manifest")?;
    let digest = sha256_bytes(temporary.path())?;
    if seal_complete {
        let staging =
            BudgetPath::from_path(temporary_directory, budget, MemoryClass::RetainedResults)?;
        let acceptance = preparation::Controls::acceptance_path(&staging, &digest)?;
        preparation::seal(temporary_directory, &acceptance)?;
    }
    let mut published = temporary
        .publish(manifest_path, no_clobber)
        .with_context(|| {
            format!(
                "Failed to atomically publish output manifest {}",
                manifest_path.display()
            )
        })?;
    sync_directory(temporary_directory)?;
    sync_directory(output_directory)?;
    published.release();
    Ok(digest)
}

fn ensure_owned_output_basename(sample: &str, basename: &str) -> Result<()> {
    ensure_safe_basename(basename)?;
    let suffix = basename.strip_prefix(sample);
    let stats = suffix == Some(".stats.yaml");
    let gene_name = suffix
        .and_then(|rest| rest.strip_prefix('_'))
        .and_then(|rest| rest.strip_suffix(".fasta"));
    ensure!(
        stats || gene_name.is_some_and(|gene| !gene.is_empty()),
        "Output receipt path is outside the owned sample namespace: {}",
        basename
    );
    Ok(())
}

fn ensure_safe_basename(basename: &str) -> Result<()> {
    let path = Path::new(basename);
    let mut components = path.components();
    ensure!(
        matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none(),
        "Output manifest path must be a basename: {}",
        basename
    );
    Ok(())
}

fn sha256_file_owned(path: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetString> {
    BudgetString::hex(&sha256_bytes(path)?, budget, MemoryClass::RetainedResults)
}

#[cfg(test)]
fn sha256_file(path: &Path) -> Result<String> {
    let bytes = sha256_bytes(path)?;
    let mut encoded = [0u8; 64];
    for (i, byte) in bytes.into_iter().enumerate() {
        encoded[2 * i] = b"0123456789abcdef"[usize::from(byte >> 4)];
        encoded[2 * i + 1] = b"0123456789abcdef"[usize::from(byte & 15)];
    }
    Ok(std::str::from_utf8(&encoded)
        .context("Invalid SHA-256 encoding")?
        .to_owned())
}

fn sha256_bytes(path: &Path) -> Result<[u8; 32]> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open {} for SHA-256", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let bytes_read = file
            .read(&mut buffer)
            .with_context(|| format!("Failed to read {} for SHA-256", path.display()))?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }
    Ok(hasher.finalize().into())
}

/// Make completed metadata publication/removal durable within its directory.
#[cfg(unix)]
pub(crate) fn sync_directory(path: &Path) -> Result<()> {
    File::open(path)
        .with_context(|| format!("Failed to open directory {} for sync", path.display()))?
        .sync_all()
        .with_context(|| format!("Failed to sync directory {}", path.display()))
}

/// Directory sync is unavailable through the portable non-Unix file API.
#[cfg(not(unix))]
pub(crate) fn sync_directory(_path: &Path) -> Result<()> {
    Ok(())
}

fn new_run_id(budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetString> {
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    BudgetString::format(
        format_args!(
            "{timestamp:x}-{:x}-{:016x}",
            std::process::id(),
            rand::random::<u64>()
        ),
        budget,
        MemoryClass::RetainedResults,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stage_text(transaction: &mut OutputTransaction, basename: &str, value: &str) {
        transaction
            .stage_file(basename, |file| {
                file.write_all(value.as_bytes())?;
                Ok(())
            })
            .unwrap();
    }

    fn read_manifest(directory: &Path, sample: &str) -> OutputManifest {
        serde_yaml_ng::from_reader(
            File::open(directory.join(format!("{sample}.manifest.yaml"))).unwrap(),
        )
        .unwrap()
    }

    fn publication_budget(bytes: u64) -> std::sync::Arc<crate::kmer::budget::ResourceBudget> {
        use crate::kmer::budget::{MEMORY_CLASSES, ResourceBudget, ResourceLimits};
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; MEMORY_CLASSES],
            // These fixtures isolate memory refusal; staged disk admission is
            // exercised separately with exact byte ceilings below.
            temporary_disk_bytes: 1 << 20,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn staged_disk_refusal_preserves_completed_outputs_and_allows_a_clean_retry() {
        use crate::kmer::budget::{
            MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits,
        };
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: 73 + 3,
            temporary_open_files: 128,
        })
        .unwrap();
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        let error = next
            .stage_file("sample_a.fasta", |file| {
                file.write_all(b"too large")?;
                Ok(())
            })
            .unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert!(next.staged_outputs.is_empty());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 73);
        assert_old_generation(directory.path(), &original);
        stage_text(&mut next, "sample_a.fasta", "ok");
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 73 + 2);
        assert_old_generation(directory.path(), &original);
        drop(next);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_old_generation(directory.path(), &original);

        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 73 + 14);
        next.commit().unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            fs::read(directory.path().join("sample_a.fasta")).unwrap(),
            b"new a"
        );
        assert_eq!(
            fs::read_to_string(directory.path().join("foreign.txt")).unwrap(),
            "keep"
        );
    }

    #[test]
    fn previous_output_verification_owns_cleanup_paths_and_refuses_without_mutation() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let manifest = read_manifest(directory.path(), "sample");
        let mut refused = 0;
        let mut succeeded = 0;
        for available in [0, 1, 8, 64, 128, 256, 512, 1024, 2048, 4096] {
            let budget = publication_budget(1 << 20);
            let root = BudgetPath::from_path(
                directory.path(),
                Some(&budget),
                MemoryClass::RetainedResults,
            )
            .unwrap();
            let baseline = budget.snapshot().reserved_memory_bytes;
            let competing = budget
                .reserve(
                    MemoryClass::RetainedResults,
                    (1 << 20) - baseline - available,
                )
                .unwrap();
            let before = budget.snapshot().reserved_memory_bytes;
            match verify_previous_outputs(&root, &manifest) {
                Ok(cleanup) => {
                    succeeded += 1;
                    assert_eq!(
                        cleanup
                            .files
                            .iter()
                            .map(|p| p.file_name().unwrap().to_str().unwrap())
                            .collect::<Vec<_>>(),
                        manifest
                            .files
                            .iter()
                            .map(|r| r.path.as_str())
                            .collect::<Vec<_>>()
                    );
                    assert!(budget.snapshot().reserved_memory_bytes > before);
                    drop(cleanup);
                }
                Err(error) => {
                    refused += 1;
                    assert!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some(),
                        "{error:#}"
                    );
                }
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, before);
            assert_old_generation(directory.path(), &original);
            drop(competing);
            let cleanup = verify_previous_outputs(&root, &manifest).unwrap();
            assert_eq!(cleanup.files.len(), 3);
            drop(cleanup);
            assert_eq!(budget.snapshot().reserved_memory_bytes, baseline);
            assert_eq!(
                budget.snapshot().peak_reserved_by_class[MemoryClass::TargetGraphs as usize],
                0
            );
            drop(root);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
        assert!(refused > 0 && succeeded > 0);
    }

    #[test]
    fn backup_copy_admission_preserves_previous_generation_and_can_retry() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        let baseline = budget.snapshot().reserved_memory_bytes;
        let competing = budget
            .reserve(MemoryClass::RetainedResults, (8 << 20) - baseline - 65535)
            .unwrap();
        let error = recovery::prepare(&mut next).unwrap_err();
        assert!(
            error.to_string().contains("output copy workspace"),
            "{error:#}"
        );
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_old_generation(directory.path(), &original);
        assert!(!next.staging_directory.join(".previous").exists());
        assert!(
            !directory
                .path()
                .join(".sharkmer-recovery-sample.yaml")
                .exists()
        );
        drop(competing);
        assert_eq!(budget.snapshot().reserved_memory_bytes, baseline);
        next.commit().unwrap();
        assert_eq!(
            fs::read(directory.path().join("sample_a.fasta")).unwrap(),
            b"new a"
        );
        assert_eq!(
            fs::read(directory.path().join("foreign.txt")).unwrap(),
            b"keep"
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn backup_disk_admission_precedes_copying_and_retries_without_public_changes() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let backup_bytes = fs::metadata(directory.path().join("sample.manifest.yaml"))
            .unwrap()
            .len()
            + original
                .iter()
                .filter(|(name, _)| name.as_str() != "sample.manifest.yaml")
                .filter(|(name, _)| name.as_str() != "foreign.txt")
                .map(|(_, bytes)| bytes.len() as u64)
                .sum::<u64>();
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        let baseline = budget.snapshot().reserved_temporary_disk_bytes;
        assert_eq!(baseline, 73 + 14);
        let competing = budget
            .reserve_disk((1 << 20) - baseline - backup_bytes + 1)
            .unwrap();
        let before = budget.snapshot().reserved_temporary_disk_bytes;
        let error = recovery::prepare(&mut next).unwrap_err();
        assert_eq!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .unwrap()
                .resource,
            crate::kmer::budget::ResourceKind::TemporaryDisk,
        );
        assert!(error.to_string().contains("recovery backup bytes"));
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, before);
        assert!(!next.staging_directory.join(".previous").exists());
        assert!(
            !directory
                .path()
                .join(".sharkmer-recovery-sample.yaml")
                .exists()
        );
        assert_old_generation(directory.path(), &original);
        drop(competing);
        next.recovery_pending = true;
        recovery::prepare(&mut next).unwrap();
        let journal_bytes = fs::metadata(directory.path().join(".sharkmer-recovery-sample.yaml"))
            .unwrap()
            .len();
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            baseline + backup_bytes + journal_bytes
        );
        assert_eq!(next.journal_disk.bytes(), journal_bytes);
        next.recover_publication().unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_old_generation(directory.path(), &original);
        drop(next);
        let mut retry =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut retry, "sample_a.fasta", "new a");
        stage_text(&mut retry, "sample.stats.yaml", "new stats");
        retry.commit().unwrap();
        assert_eq!(
            fs::read(directory.path().join("sample_a.fasta")).unwrap(),
            b"new a"
        );
        assert_eq!(
            fs::read(directory.path().join("foreign.txt")).unwrap(),
            b"keep"
        );
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn backup_directory_collision_releases_the_unwritten_disk_allowance() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        fs::create_dir(next.staging_directory.join(".previous")).unwrap();
        let before = budget.snapshot().reserved_temporary_disk_bytes;
        assert!(recovery::prepare(&mut next).is_err());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, before);
        assert_old_generation(directory.path(), &original);
        drop(next);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    fn recovery_tree_snapshot(
        root: &Path,
    ) -> std::collections::BTreeMap<std::path::PathBuf, Vec<u8>> {
        fn walk(
            root: &Path,
            directory: &Path,
            files: &mut std::collections::BTreeMap<std::path::PathBuf, Vec<u8>>,
        ) {
            for entry in fs::read_dir(directory).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    walk(root, &path, files);
                } else {
                    files.insert(
                        path.strip_prefix(root).unwrap().to_path_buf(),
                        fs::read(&path).unwrap(),
                    );
                }
            }
        }
        let mut files = std::collections::BTreeMap::new();
        walk(root, root, &mut files);
        files
    }
    #[test]
    fn recovery_workspace_refusal_preserves_all_files_and_can_retry() {
        let mut refused = 0;
        let mut succeeded = 0;
        for available in [
            0, 1, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 32768, 262144, 524288, 786432, 950000,
        ] {
            let directory = tempfile::tempdir().unwrap();
            let original = previous_generation(directory.path());
            let mut next = replacement_generation(directory.path());
            assert!(
                next.commit_inner(|step| {
                    ensure!(
                        step != PublicationStep::OutputPublished,
                        "simulated interruption"
                    );
                    Ok(())
                })
                .is_err()
            );
            next.committed = true;
            drop(next);
            let interrupted = recovery_tree_snapshot(directory.path());
            let budget = publication_budget(1 << 20);
            let root = BudgetPath::from_path(
                directory.path(),
                Some(&budget),
                MemoryClass::RetainedResults,
            )
            .unwrap();
            let baseline = budget.snapshot().reserved_memory_bytes;
            let competing = budget
                .reserve(
                    MemoryClass::RetainedResults,
                    (1 << 20) - baseline - available,
                )
                .unwrap();
            let before = budget.snapshot().reserved_memory_bytes;
            match recovery::recover(&root, "sample") {
                Ok(Some(RecoveryOutcome::Restored)) => {
                    succeeded += 1;
                    assert_old_generation(directory.path(), &original);
                }
                Err(error) => {
                    refused += 1;
                    assert!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some(),
                        "{error:#}"
                    );
                    assert_eq!(recovery_tree_snapshot(directory.path()), interrupted);
                    assert_eq!(budget.snapshot().reserved_memory_bytes, before);
                    drop(competing);
                    assert!(matches!(
                        recovery::recover(&root, "sample").unwrap(),
                        Some(RecoveryOutcome::Restored)
                    ));
                    assert_old_generation(directory.path(), &original);
                    drop(root);
                    assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                    assert!(budget.snapshot().valid);
                    continue;
                }
                _ => panic!("Interrupted generation was not restored"),
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, before);
            drop(competing);
            drop(root);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
        assert!(refused > 0 && succeeded > 0);
    }

    #[test]
    fn restore_disk_admission_preserves_interrupted_files_and_allows_exact_budget_retry() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let restore_bytes: u64 = original.iter().map(|(_, bytes)| bytes.len() as u64).sum();
        let mut interrupted = replacement_generation(directory.path());
        assert!(
            interrupted
                .commit_inner(|step| {
                    ensure!(
                        step != PublicationStep::OutputPublished,
                        "interrupt publication"
                    );
                    Ok(())
                })
                .is_err()
        );
        let staging = interrupted.staging_directory.to_path_buf();
        interrupted.committed = true;
        drop(interrupted);
        // A crash after linking can leave two names for the same payload, and
        // an earlier failed restore can leave a separate partial private copy.
        let public_stats = directory.path().join("sample.stats.yaml");
        fs::hard_link(&public_stats, staging.join("sample.stats.yaml")).unwrap();
        fs::write(
            staging.join(".previous/.sharkmer-restore-Ab12Cd"),
            b"partial old stats",
        )
        .unwrap();
        let before = recovery_tree_snapshot(directory.path());
        let private_bytes: u64 = before
            .iter()
            .filter(|(path, _)| {
                path.components().next().is_some_and(|component| {
                    component
                        .as_os_str()
                        .to_str()
                        .is_some_and(|name| name.starts_with(".sharkmer-stage-"))
                }) || path.as_path() == Path::new(".sharkmer-recovery-sample.yaml")
            })
            .map(|(_, bytes)| bytes.len() as u64)
            .sum();
        let expected_bytes =
            private_bytes + restore_bytes + fs::metadata(&public_stats).unwrap().len();
        #[cfg(unix)]
        let expected_bytes = expected_bytes - fs::metadata(&public_stats).unwrap().len();
        let budget = ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: expected_bytes,
            temporary_open_files: 128,
        })
        .unwrap();
        let root = BudgetPath::from_path(
            directory.path(),
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let competing = budget.reserve_disk(1).unwrap();
        let error = recovery::recover(&root, "sample").err().unwrap();
        assert_eq!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .unwrap()
                .resource,
            crate::kmer::budget::ResourceKind::TemporaryDisk
        );
        assert!(error.to_string().contains("output restoration bytes"));
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 1);
        drop(competing);
        assert!(matches!(
            recovery::recover(&root, "sample").unwrap(),
            Some(RecoveryOutcome::Restored)
        ));
        assert_old_generation(directory.path(), &original);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            expected_bytes
        );
        drop(root);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn restart_accounts_for_restored_public_files_until_the_old_manifest_is_accepted() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let mut next = replacement_generation(directory.path());
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
        let restored = directory.path().join("sample.stats.yaml");
        fs::remove_file(&restored).unwrap();
        fs::copy(staging.join(".previous/sample.stats.yaml"), &restored).unwrap();
        assert!(!staging.join("sample.stats.yaml").exists());
        let before = recovery_tree_snapshot(directory.path());
        // This public copy has its old hash, but no private alias and no old
        // public manifest yet. It still occupies recovery space on restart.
        let existing: u64 = before
            .iter()
            .filter(|(path, _)| {
                path.components().count() > 1
                    || path.as_path() == Path::new(".sharkmer-recovery-sample.yaml")
                    || path.as_path() == Path::new("sample.stats.yaml")
            })
            .map(|(_, bytes)| bytes.len() as u64)
            .sum();
        let restore: u64 = original
            .iter()
            .filter(|(name, _)| name != "sample.stats.yaml")
            .map(|(_, bytes)| bytes.len() as u64)
            .sum();
        let required = existing + restore;
        let budget = ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: required,
            temporary_open_files: 128,
        })
        .unwrap();
        let root = BudgetPath::from_path(
            directory.path(),
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let competing = budget.reserve_disk(1).unwrap();
        let error = recovery::recover(&root, "sample").err().unwrap();
        let exhausted = error
            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
            .unwrap();
        assert_eq!(
            exhausted.resource,
            crate::kmer::budget::ResourceKind::TemporaryDisk
        );
        assert_eq!(exhausted.required_bytes, u128::from(required + 1));
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 1);
        drop(competing);
        assert!(matches!(
            recovery::recover(&root, "sample").unwrap(),
            Some(RecoveryOutcome::Restored)
        ));
        assert_old_generation(directory.path(), &original);
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            required
        );
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        drop(root);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn manifest_disk_refusal_preserves_accepted_bytes_and_retries_on_the_same_budget() {
        use crate::kmer::budget::{ResourceExhausted, ResourceKind};
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let before = recovery_tree_snapshot(directory.path());
        let path = directory.path().join("sample.manifest.yaml");
        let mut manifest = read_manifest(directory.path(), "sample");
        manifest.run_id.push_str("-a");
        let expected = serde_yaml_ng::to_string(&manifest).unwrap();
        let budget = publication_budget(8 << 20);
        let competing = budget
            .reserve_disk((1 << 20) - expected.len() as u64 + 1)
            .unwrap();
        let error =
            write_manifest_atomic_impl(directory.path(), &path, &manifest, false, Some(&budget))
                .unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            competing.bytes()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        drop(competing);
        let digest =
            write_manifest_atomic_impl(directory.path(), &path, &manifest, false, Some(&budget))
                .unwrap();
        assert_eq!(fs::read(&path).unwrap(), expected.as_bytes());
        assert_eq!(digest, sha256_bytes(&path).unwrap());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn zero_metadata_allowance_refuses_first_manifest_without_leaving_private_files() {
        use crate::kmer::budget::{
            MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits,
        };
        let directory = tempfile::tempdir().unwrap();
        fs::write(directory.path().join("foreign.txt"), b"keep").unwrap();
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 128,
        })
        .unwrap();
        let error = OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
            .err()
            .unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert!(!directory.path().join("sample.manifest.yaml").exists());
        assert!(!directory.path().read_dir().unwrap().any(|entry| {
            entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".sharkmer-")
        }));
        assert_eq!(
            fs::read(directory.path().join("foreign.txt")).unwrap(),
            b"keep"
        );
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn journal_disk_refusal_preserves_completed_outputs_and_releases_attempt_owners() {
        use crate::kmer::budget::{ResourceExhausted, ResourceKind};
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let backup_bytes: u64 = original.iter().map(|(_, bytes)| bytes.len() as u64).sum();
        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        let baseline = budget.snapshot().reserved_temporary_disk_bytes;
        let competing = budget
            .reserve_disk((1 << 20) - baseline - backup_bytes)
            .unwrap();
        next.recovery_pending = true;
        let error = recovery::prepare(&mut next).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert!(error.to_string().contains("recovery receipt"), "{error:#}");
        assert_eq!(next.journal_disk.bytes(), 0);
        assert_eq!(next.backup_disk.bytes(), backup_bytes);
        assert!(
            !directory
                .path()
                .join(".sharkmer-recovery-sample.yaml")
                .exists()
        );
        assert!(!directory.path().read_dir().unwrap().any(|entry| {
            entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".sharkmer-recovery-")
        }));
        assert_old_generation(directory.path(), &original);
        drop(next);
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            competing.bytes()
        );
        assert_old_generation(directory.path(), &original);
        drop(competing);
        let mut retry =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut retry, "sample_a.fasta", "new a");
        stage_text(&mut retry, "sample.stats.yaml", "new stats");
        retry.commit().unwrap();
        assert_eq!(
            fs::read(directory.path().join("sample_a.fasta")).unwrap(),
            b"new a"
        );
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn in_process_recovery_credits_both_live_payloads_and_the_journal() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        next.recovery_pending = true;
        recovery::prepare(&mut next).unwrap();
        let journal_bytes = fs::metadata(directory.path().join(".sharkmer-recovery-sample.yaml"))
            .unwrap()
            .len();
        let baseline = budget.snapshot().reserved_temporary_disk_bytes;
        assert_eq!(next.journal_disk.bytes(), journal_bytes);
        let competing = budget.reserve_disk((1 << 20) - baseline).unwrap();
        next.recover_publication().unwrap();
        assert_old_generation(directory.path(), &original);
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            competing.bytes()
        );
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            1 << 20
        );
        drop(competing);
        drop(next);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn accepted_restart_adopts_cleanup_storage_without_charging_persistent_outputs() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let mut next = replacement_generation(directory.path());
        assert!(
            next.commit_inner(|step| {
                ensure!(
                    step != PublicationStep::CompleteManifest,
                    "stop after acceptance"
                );
                Ok(())
            })
            .is_err()
        );
        drop(next);
        let before = recovery_tree_snapshot(directory.path());
        let accepted: std::collections::BTreeMap<_, _> = before
            .iter()
            .filter(|(path, _)| {
                path.components().count() == 1
                    && !path
                        .as_os_str()
                        .as_encoded_bytes()
                        .starts_with(b".sharkmer-")
            })
            .map(|(path, bytes)| (path.clone(), bytes.clone()))
            .collect();
        let cleanup_bytes = before
            .iter()
            .filter(|(path, _)| !accepted.contains_key(*path))
            .map(|(_, bytes)| bytes.len() as u64)
            .sum();
        let budget = ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: 8 << 20,
            process_overhead_bytes: 0,
            class_bytes: [8 << 20; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: cleanup_bytes,
            temporary_open_files: 128,
        })
        .unwrap();
        let root = BudgetPath::from_path(
            directory.path(),
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let competing = budget.reserve_disk(1).unwrap();
        let error = recovery::recover(&root, "sample").err().unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 1);
        drop(competing);
        assert!(matches!(
            recovery::recover(&root, "sample").unwrap(),
            Some(RecoveryOutcome::Accepted)
        ));
        assert_eq!(recovery_tree_snapshot(directory.path()), accepted);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_disk_bytes,
            cleanup_bytes
        );
        drop(root);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn recovery_refuses_unlisted_private_paths_before_changing_outputs() {
        for in_backup in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let original = previous_generation(directory.path());
            let mut next = replacement_generation(directory.path());
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
            let unexpected = if in_backup {
                staging.join(".previous/foreign.txt")
            } else {
                staging.join("foreign.txt")
            };
            fs::write(&unexpected, b"foreign private file").unwrap();
            let before = recovery_tree_snapshot(directory.path());
            let budget = publication_budget(8 << 20);
            let root = BudgetPath::from_path(
                directory.path(),
                Some(&budget),
                MemoryClass::RetainedResults,
            )
            .unwrap();
            let error = recovery::recover(&root, "sample").err().unwrap();
            assert!(error.to_string().contains("Unlisted"), "{error:#}");
            assert_eq!(recovery_tree_snapshot(directory.path()), before);
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
            fs::remove_file(unexpected).unwrap();
            assert!(matches!(
                recovery::recover(&root, "sample").unwrap(),
                Some(RecoveryOutcome::Restored)
            ));
            assert_old_generation(directory.path(), &original);
            drop(root);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
            assert!(budget.snapshot().valid);
        }
    }

    #[test]
    fn historical_parse_refusal_preserves_completed_generation_and_can_retry() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        let budget = publication_budget(1 << 20);
        let competing = budget
            .reserve(MemoryClass::RetainedResults, (1 << 20) - 700)
            .unwrap();
        let before = budget.snapshot().reserved_memory_bytes;
        let result =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget));
        let error = result.err().unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some(),
            "{error:#}"
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        assert_old_generation(directory.path(), &original);
        assert!(
            !fs::read_dir(directory.path())
                .unwrap()
                .any(|e| e.unwrap().path().is_dir())
        );
        drop(competing);
        let next = OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
            .unwrap();
        assert!(next.previous_complete.is_some());
        assert_old_generation(directory.path(), &original);
        drop(next);
        assert_old_generation(directory.path(), &original);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn manifest_workspace_refuses_before_temporary_creation_and_preserves_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("sample.manifest.yaml");
        let manifest = OutputManifest {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER.to_owned(),
            sample: "sample".to_owned(),
            run_id: "run".to_owned(),
            status: "in_progress".to_owned(),
            files: vec![OutputReceipt {
                path: "sample_gene.fasta".to_owned(),
                sha256: "f".repeat(64),
            }],
            failure_reason: Some("line one\nline two λ🧬".to_owned()),
            staging_directory: Some("stage".to_owned()),
        };
        let small = publication_budget(1);
        let error =
            write_manifest_atomic_impl(directory.path(), &path, &manifest, false, Some(&small))
                .unwrap_err();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        assert_eq!(small.snapshot().reserved_memory_bytes, 0);
        let enough = publication_budget(1 << 20);
        write_manifest_atomic_impl(directory.path(), &path, &manifest, false, Some(&enough))
            .unwrap();
        assert_eq!(
            fs::read(&path).unwrap(),
            serde_yaml_ng::to_string(&manifest).unwrap().as_bytes()
        );
        assert_eq!(enough.snapshot().reserved_memory_bytes, 0);
        assert!(enough.snapshot().peak_reserved_memory_bytes > 0);
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn publication_workspace_refusal_cleans_new_staging_and_preserves_completed_generation() {
        let directory = tempfile::tempdir().unwrap();
        let small = publication_budget(1);
        let err = OutputTransaction::begin_with_budget(directory.path(), "new", Some(&small))
            .err()
            .unwrap();
        assert!(
            err.downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 0);
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "old fasta");
        stage_text(&mut first, "sample.stats.yaml", "old stats");
        first.commit().unwrap();
        let names = [
            "sample_gene.fasta",
            "sample.stats.yaml",
            "sample.manifest.yaml",
        ];
        let expected: Vec<_> = names
            .iter()
            .map(|name| fs::read(directory.path().join(name)).unwrap())
            .collect();
        let budget = publication_budget(1 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_gene.fasta", "replacement fasta");
        stage_text(&mut next, "sample.stats.yaml", "replacement stats");
        let remaining = (1 << 20) - budget.snapshot().reserved_memory_bytes;
        let competing = budget
            .reserve(MemoryClass::RetainedResults, remaining - 1)
            .unwrap();
        let err = next.commit().unwrap_err();
        assert!(
            err.downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        for (name, bytes) in names.into_iter().zip(expected) {
            assert_eq!(
                fs::read(directory.path().join(name)).unwrap(),
                bytes,
                "{name}"
            );
        }
        drop(competing);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
        assert!(!fs::read_dir(directory.path()).unwrap().any(|p| {
            p.unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".sharkmer-stage-")
        }));
    }

    #[test]
    fn completed_generation_survives_an_aborted_rerun_before_publication() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "old sequence");
        stage_text(&mut first, "sample.stats.yaml", "old stats");
        first.commit().unwrap();
        let original_manifest = fs::read(directory.path().join("sample.manifest.yaml")).unwrap();
        {
            let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
            assert_eq!(
                fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
                "old sequence"
            );
            assert_eq!(
                fs::read(directory.path().join("sample.manifest.yaml")).unwrap(),
                original_manifest
            );
            stage_text(&mut next, "sample_gene.fasta", "unfinished replacement");
        }
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
            "old sequence"
        );
        assert_eq!(
            fs::read_to_string(directory.path().join("sample.stats.yaml")).unwrap(),
            "old stats"
        );
        assert_eq!(
            fs::read(directory.path().join("sample.manifest.yaml")).unwrap(),
            original_manifest
        );
        let failures: Vec<_> = fs::read_dir(directory.path())
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| {
                path.file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(".sharkmer-failed-sample-")
            })
            .collect();
        assert_eq!(failures.len(), 1);
        let failed: OutputManifest =
            serde_yaml_ng::from_reader(File::open(&failures[0]).unwrap()).unwrap();
        assert_eq!(failed.status, "failed");
        assert!(failed.files.is_empty());
        assert!(
            !fs::read_dir(directory.path())
                .unwrap()
                .any(|entry| entry.unwrap().path().is_dir())
        );
    }

    #[test]
    fn rerun_rechecks_owned_edits_and_unowned_destinations_before_cleanup() {
        for edit_owned in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
            stage_text(&mut first, "sample_gene.fasta", "old sequence");
            stage_text(&mut first, "sample.stats.yaml", "old stats");
            first.commit().unwrap();
            let original_manifest =
                fs::read(directory.path().join("sample.manifest.yaml")).unwrap();
            let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
            let basename = if edit_owned {
                "sample_gene.fasta"
            } else {
                "sample_foreign.fasta"
            };
            fs::write(directory.path().join(basename), "user edit").unwrap();
            stage_text(&mut next, "sample_gene.fasta", "replacement");
            if !edit_owned {
                stage_text(&mut next, basename, "unowned replacement");
            }
            stage_text(&mut next, "sample.stats.yaml", "replacement stats");
            let error = next.commit().unwrap_err().to_string();
            assert!(
                error.contains(if edit_owned {
                    "modified prior output"
                } else {
                    "unowned or concurrently created"
                }),
                "{error}"
            );
            assert_eq!(
                fs::read_to_string(directory.path().join(basename)).unwrap(),
                "user edit"
            );
            assert_eq!(
                fs::read_to_string(directory.path().join("sample.stats.yaml")).unwrap(),
                "old stats"
            );
            assert_eq!(
                fs::read(directory.path().join("sample.manifest.yaml")).unwrap(),
                original_manifest
            );
        }
    }

    fn previous_generation(directory: &Path) -> Vec<(String, Vec<u8>)> {
        let mut first = OutputTransaction::begin(directory, "sample").unwrap();
        stage_text(&mut first, "sample_a.fasta", "old a");
        stage_text(&mut first, "sample_b.fasta", "old b");
        stage_text(&mut first, "sample.stats.yaml", "old stats");
        first.commit().unwrap();
        fs::write(directory.join("foreign.txt"), "keep").unwrap();
        [
            "sample_a.fasta",
            "sample_b.fasta",
            "sample.stats.yaml",
            "sample.manifest.yaml",
        ]
        .into_iter()
        .map(|name| (name.to_string(), fs::read(directory.join(name)).unwrap()))
        .collect()
    }

    fn replacement_generation(directory: &Path) -> OutputTransaction {
        let mut next = OutputTransaction::begin(directory, "sample").unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample_c.fasta", "new c");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        next
    }

    fn assert_old_generation(directory: &Path, original: &[(String, Vec<u8>)]) {
        for (name, bytes) in original {
            assert_eq!(fs::read(directory.join(name)).unwrap(), *bytes, "{name}");
        }
        assert!(!directory.join("sample_c.fasta").exists());
        assert_eq!(
            fs::read_to_string(directory.join("foreign.txt")).unwrap(),
            "keep"
        );
        assert!(!directory.join(".sharkmer-recovery-sample.yaml").exists());
    }

    #[test]
    fn publication_errors_restore_every_previous_file_and_exact_manifest() {
        // Include each unlink/link position, not only the first publication.
        for stop in 1..=8 {
            let directory = tempfile::tempdir().unwrap();
            let original = previous_generation(directory.path());
            let mut next = replacement_generation(directory.path());
            let mut step = 0;
            let result = next.commit_inner(|_| {
                step += 1;
                ensure!(step != stop, "injected publication I/O failure");
                Ok(())
            });
            assert!(result.is_err());
            drop(next);
            assert_old_generation(directory.path(), &original);
            assert!(
                !fs::read_dir(directory.path())
                    .unwrap()
                    .any(|entry| entry.unwrap().path().is_dir())
            );
        }
    }

    #[test]
    fn restart_recovers_each_interrupted_publication_and_accepts_only_complete() {
        for stop in 1..=9 {
            let directory = tempfile::tempdir().unwrap();
            let original = previous_generation(directory.path());
            let mut next = replacement_generation(directory.path());
            let mut step = 0;
            assert!(
                next.commit_inner(|_| {
                    step += 1;
                    ensure!(step != stop, "simulated process termination");
                    Ok(())
                })
                .is_err()
            );
            // Simulate termination without running the destructor's rollback.
            next.committed = true;
            drop(next);
            assert!(
                directory
                    .path()
                    .join(".sharkmer-recovery-sample.yaml")
                    .exists()
            );
            let restarted = OutputTransaction::begin(directory.path(), "sample").unwrap();
            if stop < 9 {
                assert_old_generation(directory.path(), &original);
            } else {
                assert_eq!(
                    fs::read_to_string(directory.path().join("sample_a.fasta")).unwrap(),
                    "new a"
                );
                assert_eq!(
                    fs::read_to_string(directory.path().join("sample_c.fasta")).unwrap(),
                    "new c"
                );
                assert!(!directory.path().join("sample_b.fasta").exists());
                assert_eq!(read_manifest(directory.path(), "sample").status, "complete");
            }
            drop(restarted);
        }
    }

    #[test]
    fn failed_rollback_keeps_private_payloads_and_their_disk_charges() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        assert!(
            next.commit_inner(|step| {
                ensure!(
                    step != PublicationStep::BackupReady,
                    "stop before public changes"
                );
                Ok(())
            })
            .is_err()
        );
        let staging = next.staging_directory.to_path_buf();
        let backup_bytes: u64 = fs::read_dir(staging.join(".previous"))
            .unwrap()
            .map(|entry| entry.unwrap().metadata().unwrap().len())
            .sum();
        let backup = staging.join(".previous/sample_a.fasta");
        let journal_bytes = next.journal_disk.bytes();
        assert_eq!(
            journal_bytes,
            fs::metadata(directory.path().join(".sharkmer-recovery-sample.yaml"))
                .unwrap()
                .len()
        );
        fs::write(&backup, "user edit").unwrap();
        drop(next);
        assert_eq!(fs::read(staging.join("sample_a.fasta")).unwrap(), b"new a");
        assert_eq!(
            fs::read(staging.join("sample.stats.yaml")).unwrap(),
            b"new stats"
        );
        assert_eq!(fs::read(&backup).unwrap(), b"user edit");
        assert_eq!(
            budget.snapshot().reserved_temporary_disk_bytes,
            73 + 14 + backup_bytes + journal_bytes
        );
        assert_eq!(
            fs::read(directory.path().join("sample_a.fasta")).unwrap(),
            b"old a"
        );
        assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
    }

    #[test]
    fn accepted_cleanup_uses_existing_receipts_when_yaml_decoding_would_refuse() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let budget = publication_budget(8 << 20);
        let journal = BudgetPath::from_path(
            &directory.path().join(".sharkmer-recovery-sample.yaml"),
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let mut next =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        stage_text(&mut next, "sample_a.fasta", "new a");
        stage_text(&mut next, "sample.stats.yaml", "new stats");
        let mut competing = None;
        let mut parser_refused = false;
        next.commit_inner(|step| {
            if step == PublicationStep::CompleteManifest {
                let current = budget.snapshot().reserved_memory_bytes;
                competing =
                    Some(budget.reserve(MemoryClass::RetainedResults, (8 << 20) - current - 4096)?);
                let error = parsing::recovery(&journal).err().unwrap();
                assert!(
                    error
                        .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                        .is_some()
                );
                parser_refused = true;
            }
            Ok(())
        })
        .unwrap();
        assert!(parser_refused && next.committed && !next.recovery_pending);
        assert!(!journal.exists());
        assert!(!next.staging_directory.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        drop(competing);
        drop(next);
        drop(journal);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn accepted_cleanup_preserves_modified_or_missing_recovery_evidence() {
        for change in [
            "journal",
            "backup_manifest",
            "backup_file",
            "missing_backup",
            "manifest",
            "output",
        ] {
            let directory = tempfile::tempdir().unwrap();
            previous_generation(directory.path());
            let budget = publication_budget(8 << 20);
            let mut next =
                OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                    .unwrap();
            stage_text(&mut next, "sample_a.fasta", "new a");
            stage_text(&mut next, "sample.stats.yaml", "new stats");
            let staging = next.staging_directory.to_path_buf();
            let journal = directory.path().join(".sharkmer-recovery-sample.yaml");
            let modified = match change {
                "journal" => journal.clone(),
                "backup_manifest" => staging.join(".previous/manifest.yaml"),
                "backup_file" | "missing_backup" => staging.join(".previous/sample_a.fasta"),
                "manifest" => directory.path().join("sample.manifest.yaml"),
                "output" => directory.path().join("sample_a.fasta"),
                _ => unreachable!(),
            };
            next.commit_inner(|step| {
                if step == PublicationStep::CompleteManifest {
                    if change == "missing_backup" {
                        fs::remove_file(&modified)?;
                    } else {
                        fs::write(&modified, b"user edit")?;
                    }
                }
                Ok(())
            })
            .unwrap();
            assert!(next.committed && next.recovery_pending, "{change}");
            let charged = budget.snapshot().reserved_temporary_disk_bytes;
            assert!(charged > 0);
            assert!(journal.exists() && staging.exists());
            drop(next);
            assert!(journal.exists() && staging.exists(), "{change}");
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, charged);
            if change == "missing_backup" {
                assert!(!modified.exists());
            } else {
                assert_eq!(fs::read(&modified).unwrap(), b"user edit");
            }
            assert!(
                OutputTransaction::begin(directory.path(), "sample").is_err(),
                "{change}"
            );
            assert!(budget.snapshot().valid);
        }
    }

    #[test]
    fn recovery_refuses_user_edits_before_modifying_any_other_output() {
        for modify_backup in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            previous_generation(directory.path());
            let mut next = replacement_generation(directory.path());
            assert!(
                next.commit_inner(|step| {
                    ensure!(
                        step != PublicationStep::OutputPublished,
                        "simulated interruption"
                    );
                    Ok(())
                })
                .is_err()
            );
            let backup = next.staging_directory.join(".previous/sample_a.fasta");
            let edited = if modify_backup {
                backup.clone()
            } else {
                directory.path().join("sample_a.fasta")
            };
            fs::write(&edited, "user edit").unwrap();
            // A public edit must not mutate the independent backup inode.
            if !modify_backup {
                assert_eq!(fs::read_to_string(&backup).unwrap(), "old a");
            }
            let manifest_before = fs::read(&next.manifest_path).unwrap();
            next.committed = true;
            drop(next);
            assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
            assert_eq!(fs::read_to_string(&edited).unwrap(), "user edit");
            assert_eq!(
                fs::read(directory.path().join("sample.manifest.yaml")).unwrap(),
                manifest_before
            );
            assert!(
                directory
                    .path()
                    .join(".sharkmer-recovery-sample.yaml")
                    .exists()
            );
            assert!(backup.exists());
        }
    }

    #[test]
    fn missing_old_outputs_and_partially_restored_recovery_are_idempotent() {
        let directory = tempfile::tempdir().unwrap();
        let original = previous_generation(directory.path());
        fs::remove_file(directory.path().join("sample_b.fasta")).unwrap();
        let mut next = replacement_generation(directory.path());
        let mut published = 0;
        assert!(
            next.commit_inner(|step| {
                if step == PublicationStep::OutputPublished {
                    published += 1;
                }
                ensure!(published != 3, "simulated interruption");
                Ok(())
            })
            .is_err()
        );
        // Represent a second interruption during rollback: one restored file,
        // one removed new-only file, publishing manifest still in place.
        fs::write(directory.path().join("sample_a.fasta"), "old a").unwrap();
        fs::remove_file(directory.path().join("sample_c.fasta")).unwrap();
        next.committed = true;
        drop(next);
        let restarted = OutputTransaction::begin(directory.path(), "sample").unwrap();
        for (name, bytes) in original {
            if name != "sample_b.fasta" {
                assert_eq!(fs::read(directory.path().join(name)).unwrap(), bytes);
            }
        }
        assert!(!directory.path().join("sample_b.fasta").exists());
        assert!(!directory.path().join("sample_c.fasta").exists());
        drop(restarted);
    }

    #[test]
    fn rerun_refuses_manifest_edits_even_when_prior_file_receipts_match() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let next = replacement_generation(directory.path());
        let manifest = directory.path().join("sample.manifest.yaml");
        let mut modified = fs::read(&manifest).unwrap();
        modified.extend_from_slice(b"# user annotation\n");
        fs::write(&manifest, &modified).unwrap();
        assert!(
            next.commit()
                .unwrap_err()
                .to_string()
                .contains("modified prior output manifest")
        );
        assert_eq!(fs::read(manifest).unwrap(), modified);
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_a.fasta")).unwrap(),
            "old a"
        );
    }

    #[test]
    fn changed_output_set_removes_only_receipted_prior_files() {
        let directory = tempfile::tempdir().unwrap();
        let sentinel = directory.path().join("sentinel.txt");
        fs::write(&sentinel, "keep").unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_a.fasta", "a");
        stage_text(&mut first, "sample_b.fasta", "b");
        stage_text(&mut first, "sample.stats.yaml", "first");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "sample").unwrap();
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_a.fasta")).unwrap(),
            "a"
        );
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_b.fasta")).unwrap(),
            "b"
        );
        stage_text(&mut second, "sample_a.fasta", "new-a");
        stage_text(&mut second, "sample.stats.yaml", "second");
        second.commit().unwrap();

        assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_a.fasta")).unwrap(),
            "new-a"
        );
        assert!(!directory.path().join("sample_b.fasta").exists());
        let manifest = read_manifest(directory.path(), "sample");
        assert_eq!(manifest.status, "complete");
        assert_eq!(manifest.files.len(), 2);
    }

    #[test]
    fn success_followed_by_no_products_removes_prior_fasta() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "product");
        stage_text(&mut first, "sample.stats.yaml", "success");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut second, "sample.stats.yaml", "no product");
        second.commit().unwrap();

        assert!(!directory.path().join("sample_gene.fasta").exists());
        assert!(directory.path().join("sample.stats.yaml").exists());
    }

    #[test]
    fn schema_one_split_outputs_migrate_and_modified_owned_sidecars_are_protected() {
        for modified in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
            stage_text(&mut first, "sample_gene.fasta", "standard");
            stage_text(&mut first, "sample_gene.uncertain.fasta", "uncertain");
            stage_text(&mut first, "sample.stats.yaml", "legacy stats");
            first.commit().unwrap();
            let manifest_path = directory.path().join("sample.manifest.yaml");
            let mut legacy = read_manifest(directory.path(), "sample");
            legacy.schema_version = 1;
            write_manifest_atomic(directory.path(), &manifest_path, &legacy).unwrap();
            let manifest_before = fs::read(&manifest_path).unwrap();
            let sidecar = directory.path().join("sample_gene.uncertain.fasta");
            let unowned = directory.path().join("sample_unowned.uncertain.fasta");
            fs::write(&unowned, "user file").unwrap();
            if modified {
                fs::write(&sidecar, "user edit").unwrap();
                let error = OutputTransaction::begin(directory.path(), "sample")
                    .err()
                    .unwrap();
                assert!(error.to_string().contains("modified prior output"));
                assert_eq!(fs::read(&manifest_path).unwrap(), manifest_before);
                assert_eq!(fs::read_to_string(&sidecar).unwrap(), "user edit");
                assert_eq!(
                    fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
                    "standard"
                );
            } else {
                let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
                stage_text(
                    &mut next,
                    "sample_gene.fasta",
                    "combined classified records",
                );
                stage_text(&mut next, "sample.stats.yaml", "new stats");
                next.commit().unwrap();
                assert!(!sidecar.exists());
                let current = read_manifest(directory.path(), "sample");
                assert_eq!(current.schema_version, 4);
                assert_eq!(current.files.len(), 2);
            }
            assert_eq!(fs::read_to_string(unowned).unwrap(), "user file");
        }
    }

    #[test]
    fn schemas_two_and_three_combined_outputs_migrate_without_touching_unowned_files() {
        for version in [2, 3] {
            let directory = tempfile::tempdir().unwrap();
            let old_header = if version == 2 {
                ">sample_gene_0 assembly_status=standard\nACGT\n"
            } else {
                ">sample_gene_0 sample=sample gene=gene product=0 length=4\nACGT\n"
            };
            let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
            stage_text(&mut first, "sample_gene.fasta", old_header);
            stage_text(&mut first, "sample.stats.yaml", "old stats");
            first.commit().unwrap();

            let manifest_path = directory.path().join("sample.manifest.yaml");
            let mut previous = read_manifest(directory.path(), "sample");
            previous.schema_version = version;
            write_manifest_atomic(directory.path(), &manifest_path, &previous).unwrap();
            let user_file = directory.path().join("user.fasta");
            fs::write(&user_file, old_header).unwrap();

            let mut next = OutputTransaction::begin(directory.path(), "sample").unwrap();
            let new_header = ">sample_gene_0\nACGT\n";
            stage_text(&mut next, "sample_gene.fasta", new_header);
            stage_text(&mut next, "sample.stats.yaml", "new stats");
            next.commit().unwrap();

            let current = read_manifest(directory.path(), "sample");
            assert_eq!(current.schema_version, 4);
            assert_ne!(current.run_id, previous.run_id);
            assert_eq!(current.files.len(), 2);
            assert_eq!(
                fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
                new_header
            );
            assert_eq!(fs::read_to_string(user_file).unwrap(), old_header);
        }
    }

    #[test]
    fn modified_owned_file_is_preserved_and_refused() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "owned");
        stage_text(&mut first, "sample.stats.yaml", "stats");
        first.commit().unwrap();
        fs::write(directory.path().join("sample_gene.fasta"), "modified").unwrap();

        let error = OutputTransaction::begin(directory.path(), "sample")
            .err()
            .expect("modified file must be refused");
        assert!(error.to_string().contains("modified prior output"));
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
            "modified"
        );
    }

    #[test]
    fn legacy_collision_is_not_overwritten() {
        let directory = tempfile::tempdir().unwrap();
        let stats_path = directory.path().join("sample.stats.yaml");
        fs::write(&stats_path, "legacy").unwrap();

        let error = OutputTransaction::begin(directory.path(), "sample")
            .err()
            .expect("legacy stats must be refused during begin");
        assert!(error.to_string().contains("existing unowned stats"));
        assert_eq!(fs::read_to_string(stats_path).unwrap(), "legacy");
        assert!(!directory.path().join("sample.manifest.yaml").exists());
    }

    #[test]
    fn staged_metadata_refusal_precedes_file_creation_and_can_retry_after_release() {
        for available in [0, 1, 8, 64, 128, 256, 512, 1024, 8192] {
            let directory = tempfile::tempdir().unwrap();
            let budget = publication_budget(1 << 20);
            let mut transaction =
                OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                    .unwrap();
            let path = transaction.staging_directory.join("sample_gene.fasta");
            let held = (1 << 20) - budget.snapshot().reserved_memory_bytes - available;
            let competing = budget.reserve(MemoryClass::RetainedResults, held).unwrap();
            let called = std::cell::Cell::new(false);
            let result = transaction.stage_file("sample_gene.fasta", |file| {
                called.set(true);
                file.write_all(b">seq\nACGT\n")?;
                Ok(())
            });
            if let Err(error) = result {
                assert!(
                    error
                        .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                        .is_some()
                );
                assert!(!called.get());
                assert!(!path.exists());
                assert!(transaction.staged_outputs.is_empty());
                drop(competing);
                stage_text(&mut transaction, "sample_gene.fasta", ">seq\nACGT\n");
            } else {
                assert!(called.get());
                drop(competing);
            }
            assert_eq!(transaction.staged_outputs.len(), 1);
            assert_eq!(
                transaction.staged_outputs[0].receipt.sha256.as_str(),
                sha256_file(&path).unwrap()
            );
            drop(transaction);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
    }

    #[test]
    fn borrowed_publication_views_preserve_bytes_without_metadata_copy_allocations() {
        const ENV: &str = "SHARKMER_PUBLICATION_METADATA_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "output::tests::borrowed_publication_views_preserve_bytes_without_metadata_copy_allocations", "--nocapture", "--test-threads=1"])
                .env(ENV, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        let directory = tempfile::tempdir().unwrap();
        let budget = publication_budget(8 << 20);
        let mut transaction =
            OutputTransaction::begin_with_budget(directory.path(), "sample", Some(&budget))
                .unwrap();
        for index in 0..2000 {
            let path = BudgetString::format(
                format_args!("sample_gene_{index}.fasta"),
                Some(&budget),
                MemoryClass::RetainedResults,
            )
            .unwrap();
            let staging_path = transaction
                .staging_directory
                .try_join(path.as_str())
                .unwrap();
            transaction
                .staged_outputs
                .push(StagedOutput {
                    receipt: OwnedReceipt {
                        path,
                        sha256: BudgetString::from_str(
                            &"f".repeat(64),
                            Some(&budget),
                            MemoryClass::RetainedResults,
                        )
                        .unwrap(),
                    },
                    staging_file: StagedFile::fixture(staging_path),
                })
                .unwrap();
        }
        let view = ManifestView {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER,
            sample: &transaction.sample,
            run_id: &transaction.run_id,
            status: "publishing",
            files: ReceiptList::Staged(&transaction.staged_outputs),
            failure_reason: Some("line one\nλ🧬"),
            staging_directory: Some(&transaction.staging_basename),
        };
        let historical = OutputManifest {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER.to_owned(),
            sample: transaction.sample.to_string(),
            run_id: transaction.run_id.to_string(),
            status: "publishing".to_owned(),
            files: transaction
                .staged_outputs
                .iter()
                .map(|s| OutputReceipt {
                    path: s.receipt.path.to_string(),
                    sha256: s.receipt.sha256.to_string(),
                })
                .collect(),
            failure_reason: Some("line one\nλ🧬".to_owned()),
            staging_directory: Some(transaction.staging_basename.to_string()),
        };
        assert_eq!(
            serde_yaml_ng::to_string(&view).unwrap(),
            serde_yaml_ng::to_string(&historical).unwrap()
        );
        let owned = budget.snapshot().reserved_memory_bytes;
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let prepared = crate::stats::yaml::prepare_manifest(&view, Some(&budget)).unwrap();
        let workspace = budget.snapshot().reserved_memory_bytes - owned;
        prepared.write(&mut std::io::sink()).unwrap();
        let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
        assert!(
            peak as u64 <= workspace,
            "borrowed metadata allocated beyond emitter workspace"
        );
        drop(prepared);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(budget.snapshot().reserved_memory_bytes, owned);
        println!(
            "publication_metadata_allocator rows=2000 peak_bytes={peak} reserved_workspace_bytes={workspace}"
        );
        transaction.committed = true;
        drop(transaction);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn publishing_transaction_with_partial_files_is_invalidated_on_restart() {
        let directory = tempfile::tempdir().unwrap();
        let mut interrupted = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut interrupted, "sample_gene.fasta", "partial");
        stage_text(&mut interrupted, "sample.stats.yaml", "partial stats");
        interrupted
            .staged_outputs
            .sort_by(|left, right| left.receipt.path.cmp(&right.receipt.path));
        interrupted.replacement_receipts_active = true;
        interrupted
            .write_manifest(
                "publishing",
                ReceiptList::Staged(&interrupted.staged_outputs),
                None,
            )
            .unwrap();
        let published_basename = interrupted.staged_outputs[0].receipt.path.to_string();
        publish_no_clobber(
            interrupted.staged_outputs[0].staging_file.path(),
            &directory.path().join(&published_basename),
        )
        .unwrap();
        interrupted.committed = true;
        drop(interrupted);
        assert_eq!(
            read_manifest(directory.path(), "sample").status,
            "publishing"
        );
        assert!(directory.path().join(&published_basename).exists());

        let mut restarted = OutputTransaction::begin(directory.path(), "sample").unwrap();
        assert!(!directory.path().join("sample_gene.fasta").exists());
        assert!(!directory.path().join("sample.stats.yaml").exists());
        stage_text(&mut restarted, "sample.stats.yaml", "restarted");
        restarted.commit().unwrap();
        assert_eq!(read_manifest(directory.path(), "sample").status, "complete");
    }

    #[test]
    fn overlapping_sample_and_gene_names_cannot_clobber_outputs() {
        let directory = tempfile::tempdir().unwrap();
        let collision_basename = "a_b_c.fasta";
        let mut first = OutputTransaction::begin(directory.path(), "a").unwrap();
        stage_text(&mut first, collision_basename, "first product");
        stage_text(&mut first, "a.stats.yaml", "first stats");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "a_b").unwrap();
        stage_text(&mut second, collision_basename, "second product");
        stage_text(&mut second, "a_b.stats.yaml", "second stats");
        assert!(second.commit().is_err());

        assert_eq!(
            fs::read_to_string(directory.path().join(collision_basename)).unwrap(),
            "first product"
        );
        assert_eq!(read_manifest(directory.path(), "a").status, "complete");
    }

    #[test]
    fn sample_lock_is_held_for_transaction_lifetime() {
        let directory = tempfile::tempdir().unwrap();
        let transaction = OutputTransaction::begin(directory.path(), "sample").unwrap();
        let lock_path = directory.path().join("sample.lock");
        let second_handle = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&lock_path)
            .unwrap();
        assert!(second_handle.try_lock_exclusive().is_err());
        drop(transaction);
        second_handle.try_lock_exclusive().unwrap();
        second_handle.unlock().unwrap();
        assert!(lock_path.exists());
    }

    #[test]
    fn tampered_receipt_cannot_claim_unrelated_or_lock_files() {
        assert!(ensure_owned_output_basename("sample", "sample.manifest.yaml").is_err());
        for basename in ["sentinel.txt", "sample.lock"] {
            let directory = tempfile::tempdir().unwrap();
            let claimed_path = directory.path().join(basename);
            fs::write(&claimed_path, "keep").unwrap();
            let manifest = OutputManifest {
                schema_version: MANIFEST_SCHEMA_VERSION,
                producer: MANIFEST_PRODUCER.to_string(),
                sample: "sample".to_string(),
                run_id: "abc-123".to_string(),
                status: "complete".to_string(),
                files: vec![OutputReceipt {
                    path: basename.to_string(),
                    sha256: sha256_file(&claimed_path).unwrap(),
                }],
                failure_reason: None,
                staging_directory: None,
            };
            let manifest_path = directory.path().join("sample.manifest.yaml");
            serde_yaml_ng::to_writer(File::create(&manifest_path).unwrap(), &manifest).unwrap();

            assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
            assert_eq!(fs::read_to_string(claimed_path).unwrap(), "keep");
        }
    }

    #[test]
    fn malformed_complete_manifest_is_refused() {
        let directory = tempfile::tempdir().unwrap();
        let manifest = OutputManifest {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER.to_string(),
            sample: "sample".to_string(),
            run_id: "abc-123".to_string(),
            status: "complete".to_string(),
            files: Vec::new(),
            failure_reason: None,
            staging_directory: None,
        };
        let manifest_path = directory.path().join("sample.manifest.yaml");
        serde_yaml_ng::to_writer(File::create(&manifest_path).unwrap(), &manifest).unwrap();

        assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
        assert_eq!(read_manifest(directory.path(), "sample").status, "complete");
    }

    #[cfg(unix)]
    #[test]
    fn symlink_collision_is_preserved_and_refused() {
        use std::os::unix::fs::symlink;

        let directory = tempfile::tempdir().unwrap();
        let sentinel = directory.path().join("sentinel");
        fs::write(&sentinel, "keep").unwrap();
        symlink(&sentinel, directory.path().join("sample.stats.yaml")).unwrap();

        assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
        assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        assert!(!directory.path().join("sample.manifest.yaml").exists());
    }

    #[cfg(unix)]
    #[test]
    fn symlink_lock_and_manifest_are_refused() {
        use std::os::unix::fs::symlink;

        for basename in ["sample.lock", "sample.manifest.yaml"] {
            let directory = tempfile::tempdir().unwrap();
            let sentinel = directory.path().join("sentinel");
            fs::write(&sentinel, "keep").unwrap();
            symlink(&sentinel, directory.path().join(basename)).unwrap();

            assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
            assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        }
    }
    #[test]
    fn missing_recovery_receipt_preserves_incomplete_replacement_and_backups_on_restart() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let mut next = replacement_generation(directory.path());
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
        fs::remove_file(directory.path().join(".sharkmer-recovery-sample.yaml")).unwrap();
        let before = recovery_tree_snapshot(directory.path());
        let error = OutputTransaction::begin(directory.path(), "sample")
            .err()
            .unwrap();
        assert!(
            error.to_string().contains("without its recovery receipt"),
            "{error:#}"
        );
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        assert_eq!(
            fs::read(staging.join(".previous/sample_a.fasta")).unwrap(),
            b"old a"
        );
    }

    #[test]
    fn live_recovery_cannot_treat_a_removed_journal_as_unstarted_publication() {
        let directory = tempfile::tempdir().unwrap();
        previous_generation(directory.path());
        let mut next = replacement_generation(directory.path());
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
        fs::remove_file(directory.path().join(".sharkmer-recovery-sample.yaml")).unwrap();
        let before = recovery_tree_snapshot(directory.path());
        assert!(next.recover_publication().is_err());
        assert!(next.recovery_pending);
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
        next.committed = true;
        drop(next);
        assert_eq!(recovery_tree_snapshot(directory.path()), before);
    }
}
