//! Adopt interrupted payloads before making more recovery data or cleaning up.

use super::*;
#[cfg(unix)]
use std::os::unix::fs::MetadataExt;

/// Own every private removal path before mutation. Enumerating a directory
/// while deleting its entries can skip names on network filesystems.
pub(super) struct Cleanup {
    files: BudgetVec<BudgetPath>,
    payload_files: usize,
    acceptance: Option<BudgetPath>,
    adoption_owner: Option<BudgetPath>,
    workspace: Option<workspace::Cleanup>,
}

impl Cleanup {
    pub(super) fn read(
        staging: &BudgetPath,
        backup: &BudgetPath,
        old: &ReceiptIndex<'_>,
        new: &ReceiptIndex<'_>,
        run_id: &str,
        accepted_digest: &[u8; 32],
    ) -> Result<Self> {
        let mut files = BudgetVec::in_class(staging.budget(), MemoryClass::RetainedResults);
        let controls = preparation::Controls::read(staging, run_id)?;
        let workspace = if controls.is_some() {
            workspace::Cleanup::read(staging, run_id)?
        } else {
            None
        };
        // The validated journal owns exactly this private, fixed-depth tree.
        // Refuse unexpected entries rather than following links or removing
        // unrecorded subdirectories as a side effect of recovery cleanup.
        for entry in fs::read_dir(staging)? {
            let entry = entry?;
            let name = entry.file_name();
            if name == BACKUP_DIRECTORY {
                require_directory(backup)?;
                continue;
            }
            let path = staging.try_join(&name)?;
            if controls
                .as_ref()
                .is_some_and(|controls| controls.contains(&path))
            {
                continue;
            }
            ensure!(
                name.to_str().is_some_and(|name| new.contains_key(name)
                    || (controls.is_some() && preparation::temporary(name))),
                "Unlisted staged recovery path; backups preserved"
            );
            let path = staging.try_join(&name)?;
            require_regular(&path)?;
            files.push(path)?;
        }
        for entry in fs::read_dir(backup)? {
            let entry = entry?;
            let name = entry.file_name();
            ensure!(
                name.to_str().is_some_and(|name| {
                    name == BACKUP_MANIFEST
                        || old.contains_key(name)
                        || preparation::restore_temporary_name(name)
                }),
                "Unlisted backup recovery path; backups preserved"
            );
            let path = backup.try_join(&name)?;
            require_regular(&path)?;
            files.push(path)?;
        }
        let payload_files = files.len();
        let acceptance = Some(preparation::Controls::acceptance_path(
            staging,
            accepted_digest,
        )?);
        // This path does not exist yet and has no logical extent. Keep it out
        // of the existing-inode inventory, but allocate it before admission or
        // public mutation. Creation is deferred until cleanup accepts the
        // verified complete generation while the journal still owns the scope.
        let adoption_owner = controls
            .is_none()
            .then(|| preparation::owner_path(staging, run_id))
            .transpose()?;
        if let Some(controls) = controls {
            controls.append(&mut files)?;
        }
        Ok(Self {
            files,
            payload_files,
            acceptance,
            adoption_owner,
            workspace,
        })
    }

    pub(super) fn discard(
        &self,
        directory: &Path,
        journal: &Path,
        staging: &Path,
        backup: &Path,
    ) -> Result<()> {
        if let Some(owner) = &self.adoption_owner {
            preparation::adopt_owner(staging, owner)?;
        }
        if let Some(acceptance) = &self.acceptance {
            preparation::seal(staging, acceptance)?;
        }
        if let Some(workspace) = &self.workspace {
            workspace.discard()?;
        }
        fs::remove_file(journal).context("Failed to remove finished output recovery receipt")?;
        sync_directory(directory)?;
        for path in &self.files[..self.payload_files] {
            fs::remove_file(path).with_context(|| {
                format!("Failed to remove planned recovery file {}", path.display())
            })?;
        }
        fs::remove_dir(backup).context("Failed to remove finished output backup directory")?;
        sync_directory(staging)?;
        if let Some(acceptance) = &self.acceptance {
            // Keep the owner last, including when sealing added a second anchor.
            let recorded = &self.files[self.payload_files..];
            let (owner, controls) = match &self.adoption_owner {
                Some(owner) => (owner, recorded),
                None => recorded
                    .split_last()
                    .context("Missing private owner cleanup")?,
            };
            for path in controls {
                fs::remove_file(path)?;
            }
            if !controls
                .iter()
                .any(|path| path.as_ref() == AsRef::<Path>::as_ref(acceptance))
            {
                fs::remove_file(acceptance)?;
            }
            sync_directory(staging)?;
            fs::remove_file(owner)?;
        }
        fs::remove_dir(staging).context("Failed to remove finished output staging directory")?;
        sync_directory(directory)
    }
}

/// One logical extent per inode on Unix. Other platforms conservatively count
/// every pathname, since portable metadata does not identify hard-link aliases.
pub(crate) struct Inventory {
    #[cfg(unix)]
    inodes: BudgetMap<(u64, u64), u64>,
    payload_bytes: u64,
    journal_bytes: u64,
}

impl Inventory {
    /// Start one common inode inventory for unjournaled private scopes.
    pub(crate) fn private(budget: Option<&Arc<ResourceBudget>>) -> Self {
        Self {
            #[cfg(unix)]
            inodes: BudgetMap::in_class(budget, MemoryClass::RetainedResults),
            payload_bytes: 0,
            journal_bytes: 0,
        }
    }

    pub(super) fn read(journal: &BudgetPath, cleanup: &Cleanup) -> Result<Self> {
        let mut inventory = Self {
            #[cfg(unix)]
            inodes: BudgetMap::in_class(journal.budget(), MemoryClass::RetainedResults),
            payload_bytes: 0,
            journal_bytes: fs::symlink_metadata(journal)?.len(),
        };
        for path in &cleanup.files {
            inventory.include(path)?;
        }
        if let Some(workspace) = &cleanup.workspace {
            for path in workspace.files() {
                inventory.include(path)?;
            }
        }
        Ok(inventory)
    }

    pub(crate) fn include(&mut self, path: &BudgetPath) -> Result<()> {
        let metadata = fs::symlink_metadata(path)?;
        ensure!(
            metadata.file_type().is_file(),
            "Refusing non-regular recovery storage {}; backups preserved",
            path.display()
        );
        #[cfg(unix)]
        {
            let key = (metadata.dev(), metadata.ino());
            if let Some(bytes) = self.inodes.get(&key) {
                ensure!(
                    *bytes == metadata.len(),
                    "Recovery extent changed while inventorying; backups preserved"
                );
                return Ok(());
            }
            self.inodes.insert(key, metadata.len())?;
        }
        self.payload_bytes = self
            .payload_bytes
            .checked_add(metadata.len())
            .context("Existing recovery storage length overflow")?;
        Ok(())
    }

    pub(crate) fn admit(
        &self,
        restore_bytes: u64,
        credit: u64,
        journal_credit: u64,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<disk::RecoveryCharge> {
        // Credit only matching live owners on this ledger. Journal credit is
        // separate so unused payload allowance cannot hide unowned metadata.
        let required = self
            .payload_bytes
            .saturating_sub(credit)
            .checked_add(self.journal_bytes.saturating_sub(journal_credit))
            .and_then(|bytes| bytes.checked_add(restore_bytes))
            .context("Combined recovery storage length overflow")?;
        disk::RecoveryCharge::reserve(required, budget)
            .context("Unable to admit existing recovery storage and output restoration bytes")
    }
}

#[cfg(test)]
mod tests;
