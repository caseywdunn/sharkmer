//! Owned temporary-file extents for external counts and replay candidates.
//!
//! Charges cover logical file bytes. Filesystem metadata, allocation rounding
//! and process overhead still need measured headroom in a production budget.

use super::budget::{DiskReservation, ResourceBudget, TemporaryFileReservation};
use anyhow::{Context, Result, ensure};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;
use tempfile::NamedTempFile;

/// A temporary file whose disk charge outlives its data, including on errors.
/// No mutable file handle escapes this owner to bypass extent accounting.
pub(crate) struct OwnedDiskFile {
    file: Option<NamedTempFile>,
    reservation: Option<DiskReservation>,
    descriptor: Option<TemporaryFileReservation>,
    length: u64,
    poisoned: bool,
}

impl OwnedDiskFile {
    /// Create an identifiable incomplete file in the configured directory.
    pub(crate) fn new(directory: &Path, budget: &Arc<ResourceBudget>) -> Result<Self> {
        let reservation = budget.reserve_disk(0)?;
        let descriptor = budget.reserve_temporary_file()?;
        let file = tempfile::Builder::new()
            .prefix(".sharkmer-incomplete-")
            .tempfile_in(directory)
            .context("Unable to create owned temporary file")?;
        Ok(Self {
            file: Some(file),
            reservation: Some(reservation),
            descriptor: Some(descriptor),
            length: 0,
            poisoned: false,
        })
    }

    /// Locate this owned file for diagnostics; callers must not modify it.
    pub(crate) fn path(&self) -> Result<&Path> {
        Ok(self
            .file
            .as_ref()
            .context("Temporary file is closed")?
            .path())
    }

    /// Report the accepted extent, rejecting a possibly partial I/O write.
    pub(crate) fn len(&self) -> Result<u64> {
        ensure!(
            !self.poisoned,
            "Temporary file has an incomplete I/O operation"
        );
        ensure!(self.file.is_some(), "Temporary file is closed");
        Ok(self.length)
    }

    /// Reject unexpected physical truncation or extension before format access.
    pub(crate) fn verify_extent(&self) -> Result<()> {
        let expected = self.len()?;
        let actual = self
            .file
            .as_ref()
            .context("Temporary file is closed")?
            .as_file()
            .metadata()
            .context("Unable to inspect owned temporary-file extent")?
            .len();
        ensure!(actual == expected, "Owned temporary-file extent changed");
        Ok(())
    }

    /// Reserve the new extent before writing. A configured-limit failure can
    /// be retried after other owners release space; a partial I/O failure cannot.
    pub(crate) fn append(&mut self, bytes: &[u8]) -> Result<()> {
        let length = self
            .len()?
            .checked_add(u64::try_from(bytes.len())?)
            .context("Temporary file length overflow")?;
        self.reservation
            .as_mut()
            .context("Temporary disk reservation is missing")?
            .resize(length)?;
        let file = self.file.as_mut().context("Temporary file is closed")?;
        let result = file
            .as_file_mut()
            .seek(SeekFrom::Start(self.length))
            .and_then(|_| file.as_file_mut().write_all(bytes));
        if let Err(error) = result {
            // Keep the entire attempted extent charged until confirmed removal.
            self.poisoned = true;
            return Err(error).context("Failed to extend owned temporary file");
        }
        self.length = length;
        Ok(())
    }

    /// Read a bounded caller-owned buffer; failures never become key absence.
    pub(crate) fn read_exact_at(&mut self, offset: u64, output: &mut [u8]) -> Result<()> {
        let end = offset
            .checked_add(u64::try_from(output.len())?)
            .context("Temporary-file read range overflow")?;
        ensure!(
            end <= self.len()?,
            "Read exceeds owned temporary-file extent"
        );
        let file = self.file.as_mut().context("Temporary file is closed")?;
        let result = file
            .as_file_mut()
            .seek(SeekFrom::Start(offset))
            .and_then(|_| file.as_file_mut().read_exact(output));
        if let Err(error) = result {
            self.poisoned = true;
            return Err(error).context("Failed to read owned temporary file");
        }
        Ok(())
    }

    /// Flush file contents and metadata before a format's completion step.
    pub(crate) fn sync(&self) -> Result<()> {
        self.len()?;
        self.file
            .as_ref()
            .context("Temporary file is closed")?
            .as_file()
            .sync_all()
            .context("Failed to sync owned temporary file")
    }

    fn close_inner(&mut self) -> Result<()> {
        let result = self.file.take().map(|file| file.close()).transpose();
        // NamedTempFile::close closes the descriptor even if unlinking fails.
        // Disk space and descriptor ownership therefore end independently.
        drop(self.descriptor.take());
        if let Err(error) = result {
            // Releasing this guard would pretend a surviving file is gone.
            // Retain its charge conservatively for the rest of the run.
            if let Some(reservation) = self.reservation.take() {
                std::mem::forget(reservation);
            }
            return Err(error).context(
                "Unable to remove owned temporary file; its disk reservation remains charged",
            );
        }
        drop(self.reservation.take());
        Ok(())
    }

    /// Remove the file before releasing its charge, reporting cleanup failures.
    pub(crate) fn close(mut self) -> Result<()> {
        self.close_inner()
    }
}

impl Drop for OwnedDiskFile {
    fn drop(&mut self) {
        if let Err(error) = self.close_inner() {
            log::error!("{error:#}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{ResourceExhausted, ResourceKind, ResourceLimits};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 1,
            process_overhead_bytes: 0,
            class_bytes: [0; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: bytes,
            temporary_open_files: 128,
        })
        .unwrap()
    }

    #[test]
    fn shared_extent_limit_rejects_before_write_and_retry_preserves_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(12);
        let mut first = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        let mut second = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        first.append(b"ACGT").unwrap();
        second.append(b"IJKL").unwrap();
        let first_path = first.path().unwrap().to_path_buf();
        let second_path = second.path().unwrap().to_path_buf();
        let error = first.append(b"EFGHI").unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::TemporaryDisk,
                required_bytes: 13,
                limit_bytes: 12,
            })
        );
        assert_eq!(first.len().unwrap(), 4);
        assert_eq!(std::fs::read(&first_path).unwrap(), b"ACGT");
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 8);
        second.close().unwrap();
        assert!(!second_path.exists());
        first.append(b"EFGHI").unwrap();
        first.sync().unwrap();
        assert_eq!(std::fs::read(&first_path).unwrap(), b"ACGTEFGHI");
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 9);
        first.close().unwrap();
        assert!(!first_path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 9);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn bounded_reads_do_not_overwrite_the_append_position_or_hide_range_errors() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(100);
        let mut file = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        file.append(b"AABB").unwrap();
        let mut bytes = [0; 2];
        file.read_exact_at(1, &mut bytes).unwrap();
        assert_eq!(&bytes, b"AB");
        let error = file.read_exact_at(4, &mut bytes).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
        file.append(b"CC").unwrap();
        assert_eq!(std::fs::read(file.path().unwrap()).unwrap(), b"AABBCC");
        file.read_exact_at(6, &mut []).unwrap();
        let path = file.path().unwrap().to_path_buf();
        drop(file);
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }

    #[test]
    fn truncated_file_read_poisoning_prevents_reuse_and_keeps_extent_charged() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(100);
        let mut file = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        file.append(b"complete").unwrap();
        let path = file.path().unwrap().to_path_buf();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(3)
            .unwrap();
        let error = file.read_exact_at(0, &mut [0; 8]).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
        assert_eq!(
            error.downcast_ref::<std::io::Error>().unwrap().kind(),
            std::io::ErrorKind::UnexpectedEof
        );
        assert!(file.len().is_err());
        assert!(file.append(b"more").is_err());
        assert!(file.sync().is_err());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 8);
        drop(file);
        assert!(!path.exists());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }

    #[cfg(unix)]
    #[test]
    fn failed_removal_keeps_surviving_data_charged_instead_of_reusing_its_space() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(12);
        let mut file = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        file.append(b"owned").unwrap();
        let original = file.path().unwrap().to_path_buf();
        let moved = directory.path().join("surviving-data");
        // Force a real remove-file error without relying on user permissions.
        std::fs::rename(&original, &moved).unwrap();
        std::fs::create_dir(&original).unwrap();
        let error = file.close().unwrap_err();
        assert!(error.to_string().contains("remains charged"));
        assert_eq!(std::fs::read(&moved).unwrap(), b"owned");
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 5);
        let mut other = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        assert!(other.append(b"12345678").is_err());
        assert_eq!(other.len().unwrap(), 0);
        other.close().unwrap();
        std::fs::remove_dir(original).unwrap();
        std::fs::remove_file(moved).unwrap();
    }
    #[test]
    fn descriptor_refusal_precedes_file_creation_and_preserves_existing_file() {
        let directory = tempfile::tempdir().unwrap();
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1,
            process_overhead_bytes: 0,
            class_bytes: [0; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 100,
            temporary_open_files: 1,
        })
        .unwrap();
        let mut file = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        file.append(b"accepted").unwrap();
        let path = file.path().unwrap().to_path_buf();
        let error = OwnedDiskFile::new(directory.path(), &budget).err().unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::TemporaryFilesExhausted>()
                .is_some()
        );
        assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
        assert_eq!(std::fs::read(&path).unwrap(), b"accepted");
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 1);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 8);
        file.close().unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
        let replacement = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        drop(replacement);
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }

    #[test]
    fn failed_file_creation_returns_descriptor_reservation_for_retry() {
        let directory = tempfile::tempdir().unwrap();
        let budget = budget(100);
        let error = OwnedDiskFile::new(&directory.path().join("missing"), &budget)
            .err()
            .unwrap();
        assert_eq!(
            error.downcast_ref::<std::io::Error>().unwrap().kind(),
            std::io::ErrorKind::NotFound
        );
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        let file = OwnedDiskFile::new(directory.path(), &budget).unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 1);
        drop(file);
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
        assert!(budget.snapshot().valid);
    }
}
