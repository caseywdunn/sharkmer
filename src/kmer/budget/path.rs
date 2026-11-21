//! Opaque OS paths whose backing allocation is reserved before construction.
use super::collections::Charge;
use super::{MemoryClass, ResourceBudget};
use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Immutable path data drops before the reservation for its backing buffer.
pub(crate) struct BudgetPath {
    data: PathBuf,
    charge: Charge,
}

impl BudgetPath {
    /// Borrow the original ledger for metadata derived from this path.
    pub(crate) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.charge.budget()
    }

    /// Copy encoded OS text without a lossy Unicode conversion.
    pub(crate) fn from_path(
        path: &Path,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let bound = path.as_os_str().len();
        let mut charge = Charge::in_class(budget, class);
        let reservation = charge.reserve(u64::try_from(bound)?)?;
        let mut data = PathBuf::new();
        data.try_reserve_exact(bound)
            .context("Unable to allocate output path")?;
        data.as_mut_os_string().push(path.as_os_str());
        charge.replace(
            reservation,
            u64::try_from(data.capacity())?,
            u64::try_from(bound)?,
        )?;
        Ok(Self { data, charge })
    }

    /// Preallocate a fresh joined path while keeping the original owner alive.
    pub(crate) fn try_join(&self, tail: impl AsRef<Path>) -> Result<Self> {
        let tail = tail.as_ref();
        let length = self
            .data
            .as_os_str()
            .len()
            .checked_add(tail.as_os_str().len())
            .and_then(|n| n.checked_add(1))
            .context("Output path length overflow")?;
        let (bound, workspace) = join_bounds(&self.data, tail, length)?;
        let mut charge = Charge::in_class(self.charge.budget(), self.charge.class());
        let reservation = charge.reserve(u64::try_from(bound)?)?;
        let transient = charge.reserve(u64::try_from(workspace)?)?;
        let mut data = PathBuf::new();
        data.try_reserve_exact(length)
            .context("Unable to allocate joined output path")?;
        data.as_mut_os_string().push(self.data.as_os_str());
        data.push(tail);
        charge.replace(
            reservation,
            u64::try_from(data.capacity())?,
            u64::try_from(bound)?,
        )?;
        drop(transient);
        Ok(Self { data, charge })
    }

    /// Copy independently under the original ledger and class.
    #[cfg(test)]
    pub(crate) fn try_clone(&self) -> Result<Self> {
        Self::from_path(&self.data, self.charge.budget(), self.charge.class())
    }
}

fn join_bounds(_base: &Path, _tail: &Path, length: usize) -> Result<(usize, usize)> {
    #[cfg(windows)]
    if matches!(_base.components().next(), Some(std::path::Component::Prefix(prefix)) if prefix.kind().is_verbatim())
        && !_tail.as_os_str().is_empty()
        && !_tail.is_absolute()
        && !matches!(
            _tail.components().next(),
            Some(std::path::Component::Prefix(_))
        )
    {
        // PathBuf::push's verbatim-prefix branch collects borrowed Components
        // then constructs a second OsString. Precharge Vec growth overlap and
        // both old/replacement text buffers before invoking the standard API.
        // Other joins use their encoded text bound. No normalization rule
        // or encoded-path behavior is reimplemented here.
        let text = length
            .checked_mul(2)
            .context("Output path capacity overflow")?
            .max(8);
        let components = length
            .checked_add(1)
            .and_then(|n| n.checked_mul(2))
            .map(|n| n.max(4))
            .context("Output path component capacity overflow")?;
        let stack = components
            .checked_mul(size_of::<std::path::Component<'_>>())
            .and_then(|n| n.checked_mul(2))
            .context("Output path component bytes overflow")?;
        let workspace = length
            .checked_mul(3)
            .and_then(|n| n.checked_add(stack))
            .context("Output path normalization workspace overflow")?;
        return Ok((text, workspace));
    }
    Ok((length, 0))
}

impl std::ops::Deref for BudgetPath {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.data
    }
}
impl AsRef<Path> for BudgetPath {
    fn as_ref(&self) -> &Path {
        &self.data
    }
}
impl std::fmt::Debug for BudgetPath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}
#[cfg(test)]
impl Clone for BudgetPath {
    fn clone(&self) -> Self {
        self.try_clone().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn owned_paths_match_standard_join_and_preserve_opaque_encoded_bytes() {
        let budget = budget(1 << 20);
        for base in [
            "",
            "/",
            "relative",
            "a/../b",
            "x/.",
            "λ🧬",
            "a//",
            r"\\?\C:\root\..",
        ] {
            for tail in [
                "",
                "file",
                "../other",
                "/absolute",
                "foo/bar",
                "δ",
                r"..\leaf",
            ] {
                let owned = BudgetPath::from_path(
                    Path::new(base),
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
                .unwrap();
                let joined = owned.try_join(tail).unwrap();
                assert_eq!(
                    joined.as_os_str().as_encoded_bytes(),
                    Path::new(base).join(tail).as_os_str().as_encoded_bytes()
                );
                assert_eq!(owned.as_os_str(), Path::new(base).as_os_str());
            }
        }
        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStrExt;
            let base = Path::new(std::ffi::OsStr::from_bytes(b"opaque\xff/root"));
            let tail = Path::new(std::ffi::OsStr::from_bytes(b"leaf\xfe"));
            let owned =
                BudgetPath::from_path(base, Some(&budget), MemoryClass::RetainedResults).unwrap();
            assert_eq!(
                owned.try_join(tail).unwrap().as_os_str().as_bytes(),
                base.join(tail).as_os_str().as_bytes()
            );
        }
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(
            budget.snapshot().peak_reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
    }

    #[test]
    fn owned_path_refusal_preserves_parent_and_retries_after_competing_owner_drops() {
        let (bound, workspace) = join_bounds(Path::new("parent"), Path::new("child"), 12).unwrap();
        let budget = budget((6 + bound + workspace) as u64);
        let parent = BudgetPath::from_path(
            Path::new("parent"),
            Some(&budget),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let other = budget.reserve(MemoryClass::RetainedResults, 1).unwrap();
        let error = parent.try_join("child").err().unwrap();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(parent.as_os_str(), "parent");
        assert_eq!(budget.snapshot().reserved_memory_bytes, 7);
        drop(other);
        let joined = parent.try_join("child").unwrap();
        assert_eq!(joined.as_ref(), Path::new("parent").join("child"));
        drop(joined);
        drop(parent);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn owned_path_reservations_cover_isolated_allocator_peaks_and_release_buffers() {
        const ENV: &str = "SHARKMER_PATH_ALLOCATION_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "kmer::budget::path::tests::owned_path_reservations_cover_isolated_allocator_peaks_and_release_buffers", "--nocapture", "--test-threads=1"])
                .env(ENV, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        for (case, base) in ["", "/absolute", "relative/../path", r"\\?\C:\root\.."]
            .into_iter()
            .enumerate()
        {
            for length in [0, 1, 31, 1024] {
                let tail = "δ".repeat(length);
                let expected = Path::new(base).join(&tail);
                let budget = budget(8 << 20);
                let _ = budget.snapshot();
                let baseline = crate::PEAK_ALLOC.current_usage();
                crate::PEAK_ALLOC.reset_peak_usage();
                let parent = BudgetPath::from_path(
                    Path::new(base),
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
                .unwrap();
                let joined = parent.try_join(&tail).unwrap();
                let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
                let reserved = budget.snapshot().peak_reserved_memory_bytes;
                assert_eq!(joined.as_os_str(), expected.as_os_str());
                assert!(
                    peak as u64 <= reserved,
                    "path heap {peak} exceeded {reserved}"
                );
                drop(joined);
                drop(parent);
                assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                println!(
                    "path_allocator case={case} tail_units={length} peak_bytes={peak} reserved_bytes={reserved}"
                );
            }
        }
    }
}
