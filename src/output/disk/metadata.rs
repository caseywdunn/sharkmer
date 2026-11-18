//! Closed, charged temporary metadata before atomic publication.

use super::*;

/// Keep the file before its path so allocation failures close before unlinking.
struct TemporaryParts {
    file: File,
    path: tempfile::TempPath,
}

/// Reuse staged-data write accounting without exposing an uncharged handle.
/// Field order also closes and removes an unpublished file before its path drops.
pub(crate) struct MetadataTemporary {
    file: StagedFile,
    path: Option<tempfile::TempPath>,
}

impl MetadataTemporary {
    /// Establish ownership before creating the empty metadata file.
    pub(crate) fn create(
        directory: &Path,
        prefix: &str,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Self> {
        let reservation = budget.map(|budget| budget.reserve_disk(0)).transpose()?;
        Self::create_reserved(directory, prefix, budget, reservation)
    }

    /// Use an allowance admitted before creating the surrounding private scope.
    pub(crate) fn create_reserved(
        directory: &Path,
        prefix: &str,
        budget: Option<&Arc<ResourceBudget>>,
        reservation: Option<DiskReservation>,
    ) -> Result<Self> {
        let temporary = tempfile::Builder::new()
            .prefix(prefix)
            .tempfile_in(directory)?;
        let (file, path) = temporary.into_parts();
        let parts = TemporaryParts { file, path };
        let owned_path = BudgetPath::from_path(&parts.path, budget, MemoryClass::RetainedResults)?;
        let TemporaryParts { file, path } = parts;
        Ok(Self {
            file: StagedFile {
                path: owned_path,
                file: Some(file),
                reservation,
                length: 0,
                failure: None,
                poisoned: false,
                publishing: false,
            },
            path: Some(path),
        })
    }

    /// Borrow only the write interface whose growth is admitted before I/O.
    pub(crate) fn writer(&mut self) -> &mut StagedFile {
        &mut self.file
    }

    /// Preserve typed admission failures, then flush, verify, sync and close.
    pub(crate) fn finish<T>(&mut self, written: Result<T>) -> Result<T> {
        self.file.finish(written)
    }

    /// Hash only the closed private file before publishing its name.
    pub(crate) fn path(&self) -> &Path {
        self.file.path()
    }

    /// Live extent credit for recovery on this temporary's existing ledger.
    pub(crate) fn reserved_bytes(&self) -> u64 {
        self.file.reserved_bytes()
    }

    /// A public hard link may survive until durable acceptance or rollback.
    pub(crate) fn begin_link_publication(&mut self) -> Result<()> {
        ensure!(
            self.file.file.is_none() && !self.file.poisoned && self.file.failure.is_none(),
            "Temporary must be finished and closed before publication"
        );
        self.file.publishing();
        Ok(())
    }

    /// Release only after private links are gone and public acceptance or
    /// rollback has been synced. The caller proves both transitions.
    pub(crate) fn confirm_removed_or_accepted(&mut self) {
        self.file.publishing = false;
        self.file.release();
    }

    /// Transfer the allowance to the published name until confirmed acceptance
    /// or removal. A failed atomic operation can still have created that name.
    pub(crate) fn publish(
        mut self,
        destination: &Path,
        no_clobber: bool,
    ) -> Result<RecoveryCharge> {
        ensure!(
            self.file.file.is_none() && !self.file.poisoned && self.file.failure.is_none(),
            "Metadata must be finished and closed before publication"
        );
        let path = self.path.take().context("Metadata path is unavailable")?;
        self.file.publishing();
        let published = if no_clobber {
            path.persist_noclobber(destination)
        } else {
            path.persist(destination)
        };
        published.map_err(|error| error.error).with_context(|| {
            format!(
                "Failed to atomically publish metadata {}",
                destination.display()
            )
        })?;
        Ok(RecoveryCharge(self.file.reservation.take()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 1 << 20,
            process_overhead_bytes: 0,
            class_bytes: [1 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: bytes,
            temporary_open_files: 128,
        })
        .unwrap()
    }

    #[test]
    fn metadata_write_refusal_keeps_the_original_error_and_releases_private_storage() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(0);
        let mut temporary =
            MetadataTemporary::create(directory.path(), ".metadata-", Some(&budget)).unwrap();
        let path = temporary.path().to_path_buf();
        assert!(temporary.writer().write_all(b"abc").is_err());
        let error = temporary.finish(Ok(())).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        assert!(fs::read(&path).unwrap().is_empty());
        drop(temporary);
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn published_closed_metadata_keeps_its_charge_until_durable_acceptance() {
        for accepted in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let budget = budget(3);
            let mut temporary =
                MetadataTemporary::create(directory.path(), ".metadata-", Some(&budget)).unwrap();
            temporary.writer().write_all(b"abc").unwrap();
            temporary.finish(Ok(())).unwrap();
            let path = temporary.path().to_path_buf();
            #[cfg(target_os = "linux")]
            for entry in fs::read_dir("/proc/self/fd").unwrap() {
                if let Ok(open) = fs::read_link(entry.unwrap().path()) {
                    assert_ne!(open, path);
                }
            }
            let destination = directory.path().join("published.yaml");
            let mut charge = temporary.publish(&destination, true).unwrap();
            assert!(!path.exists());
            assert_eq!(fs::read(&destination).unwrap(), b"abc");
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 3);
            if accepted {
                sync_directory(directory.path()).unwrap();
                charge.release();
            }
            drop(charge);
            assert_eq!(
                budget.snapshot().reserved_temporary_disk_bytes,
                if accepted { 0 } else { 3 }
            );
            assert!(budget.snapshot().valid);
        }
    }
}
