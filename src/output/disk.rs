//! Staged output bytes remain charged until removal or accepted publication.

use super::*;
use crate::kmer::budget::{DiskReservation, ResourceBudget};
use std::sync::Arc;

mod metadata;
pub(crate) use metadata::MetadataTemporary;

/// Recovery copies keep their allowance until removal or durable acceptance.
/// The caller owns that transition; dropping this guard alone cannot prove it.
#[derive(Default)]
pub(crate) struct RecoveryCharge(Option<DiskReservation>);

impl RecoveryCharge {
    /// Live allowance belonging to this transaction's recovery payloads.
    pub(crate) fn bytes(&self) -> u64 {
        self.0.as_ref().map_or(0, DiskReservation::bytes)
    }

    /// Admit the entire copy set before its first filesystem mutation.
    pub(crate) fn reserve(bytes: u64, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        Ok(Self(budget.map(|b| b.reserve_disk(bytes)).transpose()?))
    }

    /// Release after confirmed removal or acceptance as persistent output.
    pub(crate) fn release(&mut self) {
        drop(self.0.take());
    }
}

impl Drop for RecoveryCharge {
    fn drop(&mut self) {
        if let Some(reservation) = self.0.take() {
            // Recovery data can still need these bytes. The caller releases
            // only after confirmed cleanup or durable output acceptance.
            std::mem::forget(reservation);
        }
    }
}

/// Write-only staged data. Callers cannot seek or access an uncharged file handle.
pub(crate) struct StagedFile {
    path: BudgetPath,
    file: Option<File>,
    reservation: Option<DiskReservation>,
    length: u64,
    failure: Option<anyhow::Error>,
    poisoned: bool,
    publishing: bool,
}

impl StagedFile {
    /// The allowance remains live even after its private name is published.
    pub(super) fn reserved_bytes(&self) -> u64 {
        self.reservation.as_ref().map_or(0, DiskReservation::bytes)
    }

    /// Establish byte ownership before creating a private output file.
    pub(super) fn create(path: BudgetPath, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let reservation = budget.map(|b| b.reserve_disk(0)).transpose()?;
        let file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&path)
            .with_context(|| format!("Failed to create staged output {}", path.display()))?;
        Ok(Self {
            path,
            file: Some(file),
            reservation,
            length: 0,
            failure: None,
            poisoned: false,
            publishing: false,
        })
    }

    /// Borrow the private path without exposing a mutable file handle.
    pub(super) fn path(&self) -> &Path {
        &self.path
    }

    /// Reject swallowed failures and changed extents, then sync and close.
    pub(super) fn finish<T>(&mut self, written: Result<T>) -> Result<T> {
        if let Some(error) = self.failure.take() {
            return Err(error).context("Unable to admit staged output bytes");
        }
        let value = written?;
        ensure!(
            !self.poisoned,
            "Staged output has an incomplete I/O operation"
        );
        let file = self.file.as_mut().context("Staged output is closed")?;
        file.flush().context("Failed to flush staged output")?;
        file.sync_all().context("Failed to sync staged output")?;
        ensure!(
            file.metadata()
                .context("Failed to inspect staged output")?
                .len()
                == self.length,
            "Staged output extent changed outside its writer"
        );
        drop(self.file.take());
        Ok(value)
    }

    /// A hard-linked public name can outlive the private path on failure.
    pub(super) fn publishing(&mut self) {
        self.publishing = true;
    }

    /// Called only after durable acceptance or confirmed rollback cleanup.
    pub(super) fn release(&mut self) {
        drop(self.reservation.take());
    }

    #[cfg(test)]
    /// Construct metadata fixtures without creating or charging file data.
    pub(super) fn fixture(path: BudgetPath) -> Self {
        Self {
            path,
            file: None,
            reservation: None,
            length: 0,
            failure: None,
            poisoned: false,
            publishing: false,
        }
    }
}

impl Write for StagedFile {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self.poisoned || self.failure.is_some() {
            return Err(std::io::Error::other(
                "Staged output has an incomplete I/O operation",
            ));
        }
        if self.file.is_none() {
            return Err(std::io::Error::other("Staged output is closed"));
        }
        let length = self
            .length
            .checked_add(u64::try_from(bytes.len()).map_err(std::io::Error::other)?)
            .ok_or_else(|| std::io::Error::other("Staged output length overflow"))?;
        if let Some(reservation) = &mut self.reservation
            && let Err(error) = reservation.resize(length)
        {
            self.failure = Some(error);
            self.poisoned = true;
            return Err(std::io::Error::other("Unable to admit staged output bytes"));
        }
        let result = self
            .file
            .as_mut()
            .ok_or_else(|| std::io::Error::other("Staged output is closed"))?
            .write(bytes);
        match result {
            Ok(written) => {
                self.length += written as u64;
                if let Some(reservation) = &mut self.reservation
                    && let Err(error) = reservation.resize(self.length)
                {
                    self.failure = Some(error);
                    self.poisoned = true;
                    return Err(std::io::Error::other(
                        "Unable to account for staged output bytes",
                    ));
                }
                Ok(written)
            }
            // Write::write_all retries EINTR without advancing the file offset.
            // Keep the attempted extent charged while allowing that retry.
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => Err(error),
            Err(error) => {
                self.poisoned = true;
                Err(error)
            }
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        if self.poisoned || self.failure.is_some() {
            return Err(std::io::Error::other(
                "Staged output has an incomplete I/O operation",
            ));
        }
        let result = self
            .file
            .as_mut()
            .ok_or_else(|| std::io::Error::other("Staged output is closed"))?
            .flush();
        if result.is_err() {
            self.poisoned = true;
        }
        result
    }
}

impl Drop for StagedFile {
    fn drop(&mut self) {
        drop(self.file.take());
        if self.publishing {
            // Recovery failed or publication was interrupted. The public link
            // can still own these bytes; releasing would undercount that data.
            if let Some(reservation) = self.reservation.take() {
                std::mem::forget(reservation);
            }
            return;
        }
        match fs::remove_file(&self.path) {
            Ok(()) => self.release(),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => self.release(),
            Err(error) => {
                if let Some(reservation) = self.reservation.take() {
                    std::mem::forget(reservation);
                }
                log::error!(
                    "Unable to remove staged output {}; its disk charge is retained: {error}",
                    self.path.display()
                );
            }
        }
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

    fn staged(path: &Path, budget: &Arc<ResourceBudget>) -> StagedFile {
        StagedFile::create(
            BudgetPath::from_path(path, Some(budget), MemoryClass::RetainedResults).unwrap(),
            Some(budget),
        )
        .unwrap()
    }

    #[test]
    fn refused_growth_does_not_write_and_cannot_be_hidden_by_the_callback() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("stage");
        let budget = budget(6);
        let mut file = staged(&path, &budget);
        file.write_all(b"abc").unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 3);
        assert!(file.write_all(b"defg").is_err());
        assert_eq!(fs::read(&path).unwrap(), b"abc");
        let error = file.finish(Ok(())).unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::TemporaryDisk
        );
        drop(file);
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 3);
    }

    #[test]
    fn linked_output_stays_charged_until_acceptance_and_survives_owner_drop() {
        for accepted in [false, true] {
            let directory = tempfile::tempdir().unwrap();
            let path = directory.path().join("stage");
            let public = directory.path().join("public");
            let budget = budget(3);
            let mut file = staged(&path, &budget);
            file.write_all(b"abc").unwrap();
            file.finish(Ok(())).unwrap();
            assert!(file.write_all(b"x").is_err());
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 3);
            file.publishing();
            fs::hard_link(&path, &public).unwrap();
            fs::remove_file(&path).unwrap();
            if accepted {
                file.release();
            }
            drop(file);
            assert_eq!(fs::read(&public).unwrap(), b"abc");
            assert_eq!(
                budget.snapshot().reserved_temporary_disk_bytes,
                if accepted { 0 } else { 3 }
            );
        }
    }

    #[test]
    fn cleanup_failure_retains_the_charge_and_preserves_the_remaining_path() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("stage");
        let budget = budget(3);
        let mut file = staged(&path, &budget);
        file.write_all(b"abc").unwrap();
        file.finish(Ok(())).unwrap();
        fs::remove_file(&path).unwrap();
        fs::create_dir(&path).unwrap();
        fs::write(path.join("foreign"), b"keep").unwrap();
        drop(file);
        assert_eq!(fs::read(path.join("foreign")).unwrap(), b"keep");
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 3);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn unwinding_closes_and_removes_private_data_before_releasing_its_charge() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("stage");
        let budget = budget(3);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut file = staged(&path, &budget);
            file.write_all(b"abc").unwrap();
            panic!("stop during staged output");
        }));
        assert!(result.is_err());
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }

    #[test]
    fn external_extent_changes_refuse_finalization_and_clean_private_data() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("stage");
        let budget = budget(3);
        let mut file = staged(&path, &budget);
        file.write_all(b"abc").unwrap();
        OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap()
            .write_all(b"x")
            .unwrap();
        assert!(
            file.finish(Ok(()))
                .unwrap_err()
                .to_string()
                .contains("extent changed")
        );
        drop(file);
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }
}
