//! One owned copy buffer shared by all backup or restore files in a transaction.
use super::*;

const COPY_BYTES: usize = 64 << 10;

/// Reserve before any recovery mutation and retain through the final synced copy.
pub(super) struct CopyWorkspace {
    buffer: Option<BudgetVec<u8>>,
}

impl CopyWorkspace {
    pub(super) fn new(budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let buffer = budget
            .map(|budget| {
                BudgetVec::collect_in_class(
                    std::iter::repeat_n(0, COPY_BYTES),
                    Some(budget),
                    MemoryClass::RetainedResults,
                )
                .context("Unable to admit output copy workspace")
            })
            .transpose()?;
        Ok(Self { buffer })
    }

    /// Preserve the existing legacy copy path; bounded copies use no library buffer.
    pub(super) fn copy(
        &mut self,
        source: &mut impl Read,
        destination: &mut impl Write,
    ) -> Result<()> {
        let Some(buffer) = &mut self.buffer else {
            std::io::copy(source, destination)?;
            return Ok(());
        };
        loop {
            let count = match source.read(buffer) {
                Ok(0) => return Ok(()),
                Ok(count) => count,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error).context("Failed to read output recovery copy"),
            };
            destination
                .write_all(&buffer[..count])
                .context("Failed to write output recovery copy")?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{self, Cursor};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = bytes;
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn copy_retries_interrupted_reads_and_short_writes_without_losing_bytes() {
        struct Interrupted<'a>(bool, Cursor<&'a [u8]>);
        impl Read for Interrupted<'_> {
            fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
                if self.0 {
                    self.0 = false;
                    return Err(io::ErrorKind::Interrupted.into());
                }
                self.1.read(bytes)
            }
        }
        struct ShortWriter {
            bytes: Vec<u8>,
            refuse: bool,
        }
        impl Write for ShortWriter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                if self.refuse {
                    return Ok(0);
                }
                let count = bytes.len().min(7);
                self.bytes.extend_from_slice(&bytes[..count]);
                Ok(count)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let data: Vec<u8> = (0..3 * COPY_BYTES + 19).map(|i| (i % 251) as u8).collect();
        for enabled in [false, true] {
            let budget = budget(COPY_BYTES as u64);
            let mut workspace = CopyWorkspace::new(enabled.then_some(&budget)).unwrap();
            let mut output = ShortWriter {
                bytes: Vec::new(),
                refuse: false,
            };
            workspace
                .copy(&mut Interrupted(true, Cursor::new(&data)), &mut output)
                .unwrap();
            assert_eq!(output.bytes, data);
            output.refuse = true;
            let error = workspace
                .copy(&mut Cursor::new(&data), &mut output)
                .unwrap_err();
            assert_eq!(
                error.downcast_ref::<io::Error>().unwrap().kind(),
                io::ErrorKind::WriteZero
            );
            assert_eq!(
                budget.snapshot().reserved_memory_bytes,
                if enabled { COPY_BYTES as u64 } else { 0 }
            );
            drop(workspace);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
    }

    #[test]
    fn copy_workspace_covers_isolated_heap_and_releases_on_panic() {
        const ENV: &str = "SHARKMER_OUTPUT_COPY_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "output::recovery::copy::tests::copy_workspace_covers_isolated_heap_and_releases_on_panic", "--nocapture", "--test-threads=1"])
                .env(ENV, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        for bytes in [0, 1, COPY_BYTES - 1, COPY_BYTES, COPY_BYTES + 1, 4 << 20] {
            let data = vec![b'A'; bytes];
            let budget = budget(COPY_BYTES as u64);
            let _ = budget.snapshot();
            let baseline = crate::PEAK_ALLOC.current_usage();
            crate::PEAK_ALLOC.reset_peak_usage();
            let mut workspace = CopyWorkspace::new(Some(&budget)).unwrap();
            workspace
                .copy(&mut Cursor::new(&data), &mut io::sink())
                .unwrap();
            let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
            let reserved = budget.snapshot().peak_reserved_memory_bytes;
            assert_eq!(peak as u64, reserved);
            assert_eq!(reserved, COPY_BYTES as u64);
            drop(workspace);
            assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            println!(
                "publication_copy_allocator bytes={bytes} peak_bytes={peak} reserved_bytes={reserved}"
            );
        }
        struct Panicking;
        impl Read for Panicking {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                panic!("injected copy-reader panic");
            }
        }
        let budget = budget(COPY_BYTES as u64);
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let mut workspace = CopyWorkspace::new(Some(&budget)).unwrap();
                assert_eq!(budget.snapshot().reserved_memory_bytes, COPY_BYTES as u64);
                workspace.copy(&mut Panicking, &mut io::sink()).unwrap();
            }))
            .is_err()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&budget), 1);
        assert!(budget.snapshot().valid);
    }
}
