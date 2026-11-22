//! Shared logical memory/disk reservations for v4 stages.
//!
//! Reservations are ownership guards, not allocator hooks or RSS measurements.
//! Production limits must include measured process/allocator overhead.

use anyhow::{Context, Result, ensure};
use std::ops::{Deref, DerefMut};
use std::sync::{Arc, Mutex, MutexGuard};

pub(crate) mod collections;
pub(crate) mod path;
mod process_files;
pub(crate) mod text;
pub(crate) use process_files::ProcessFileLimits;

/// Independent memory owners charged against the same run limit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
pub enum MemoryClass {
    CountSlots,
    InputBuffer,
    WorkerQueues,
    SortBuffer,
    LookupCache,
    TargetGraphs,
    ReadEvidence,
    ReplayBuffer,
    RetainedResults,
}

pub(crate) const MEMORY_CLASSES: usize = 9;
const CLASSES: usize = MEMORY_CLASSES;

/// The configured resource whose live reservation could not be extended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResourceKind {
    StageMemory(MemoryClass),
    RunMemory,
    TemporaryDisk,
}

impl std::fmt::Display for ResourceKind {
    fn fmt(&self, output: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::StageMemory(MemoryClass::CountSlots) => "Count-table memory",
            Self::StageMemory(MemoryClass::InputBuffer) => "Input-buffer memory",
            Self::StageMemory(MemoryClass::WorkerQueues) => "Worker-queue memory",
            Self::StageMemory(MemoryClass::SortBuffer) => "Sort-buffer memory",
            Self::StageMemory(MemoryClass::LookupCache) => "Lookup-cache memory",
            Self::StageMemory(MemoryClass::TargetGraphs) => "Target-graph memory",
            Self::StageMemory(MemoryClass::ReadEvidence) => "Read-evidence memory",
            Self::StageMemory(MemoryClass::ReplayBuffer) => "Replay-buffer memory",
            Self::StageMemory(MemoryClass::RetainedResults) => "Retained-result memory",
            Self::RunMemory => "Combined run memory",
            Self::TemporaryDisk => "Temporary disk space",
        };
        output.write_str(name)
    }
}

/// A recoverable configured-limit failure, distinct from invalid input,
/// arithmetic errors, allocator failures and broken accounting invariants.
///
/// Callers may downcast through `anyhow` context to decide whether to spill
/// counts or reduce concurrency. Required bytes use u128 so even a request
/// beyond u64 can be reported exactly without wrapping reservation accounting.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ResourceExhausted {
    pub resource: ResourceKind,
    pub required_bytes: u128,
    pub limit_bytes: u64,
}

impl std::fmt::Display for ResourceExhausted {
    fn fmt(&self, output: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            output,
            "{} reservation requires {} bytes; configured budget is {} bytes",
            self.resource, self.required_bytes, self.limit_bytes
        )?;
        match self.resource {
            ResourceKind::TemporaryDisk => {
                output.write_str(". Raise --temp-disk-budget if more storage is available")
            }
            ResourceKind::StageMemory(MemoryClass::TargetGraphs) | ResourceKind::RunMemory => {
                output.write_str(". Raise --memory-budget if more memory is available")
            }
            _ => Ok(()),
        }
    }
}

impl std::error::Error for ResourceExhausted {}

/// A temporary-file descriptor limit, distinct from byte reservations and I/O.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TemporaryFilesExhausted {
    pub required_files: u128,
    pub limit_files: u64,
}

impl std::fmt::Display for TemporaryFilesExhausted {
    fn fmt(&self, output: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            output,
            "Count/replay files require {} open descriptors; configured budget is {}",
            self.required_files, self.limit_files
        )
    }
}

impl std::error::Error for TemporaryFilesExhausted {}

/// Configured ceilings; simultaneous live reservations must also fit total RAM.
#[derive(Clone, Debug)]
pub struct ResourceLimits {
    pub memory_bytes: u64,
    pub process_overhead_bytes: u64,
    pub class_bytes: [u64; CLASSES],
    pub temporary_disk_bytes: u64,
    pub temporary_open_files: u64,
}

/// Observable logical usage/peaks, including the initial process reservation.
#[derive(Clone, Debug)]
pub struct ResourceSnapshot {
    pub reserved_memory_bytes: u64,
    pub peak_reserved_memory_bytes: u64,
    pub reserved_by_class: [u64; CLASSES],
    pub peak_reserved_by_class: [u64; CLASSES],
    pub reserved_temporary_disk_bytes: u64,
    pub peak_reserved_temporary_disk_bytes: u64,
    pub reserved_temporary_open_files: u64,
    pub peak_reserved_temporary_open_files: u64,
    pub valid: bool,
}

/// A shared run budget. Clones share accounting across stages and workers.
pub struct ResourceBudget {
    limits: ResourceLimits,
    state: Mutex<BudgetState>,
    process_file_limits: Option<ProcessFileLimits>,
    parent_memory: Option<MemoryReservation>,
    scope_finished: Option<Box<dyn FnOnce() + Send + Sync + std::panic::RefUnwindSafe>>,
}

struct BudgetState {
    snapshot: ResourceSnapshot,
    memory_peak: Option<u64>,
    disk_peak: Option<u64>,
}
impl Deref for BudgetState {
    type Target = ResourceSnapshot;
    fn deref(&self) -> &Self::Target {
        &self.snapshot
    }
}
impl DerefMut for BudgetState {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.snapshot
    }
}

/// Observe one non-nested phase without resetting the run's lifetime peak.
pub(crate) struct MemoryPeak {
    budget: Option<Arc<ResourceBudget>>,
}
impl MemoryPeak {
    /// End the phase after its last allocation and report simultaneous usage.
    pub(crate) fn finish(mut self) -> Result<u64> {
        let budget = self
            .budget
            .take()
            .context("Memory phase already finished")?;
        let mut state = budget.state();
        let peak = state
            .memory_peak
            .take()
            .context("Memory phase is not active")?;
        ensure!(state.valid, "Resource accounting is invalid");
        Ok(peak)
    }
}
impl Drop for MemoryPeak {
    fn drop(&mut self) {
        if let Some(budget) = &self.budget {
            budget.state().memory_peak = None;
        }
    }
}

/// Observe a temporary-disk phase while preserving the whole-run ledger peak.
pub(crate) struct DiskPeak {
    budget: Option<Arc<ResourceBudget>>,
}

impl DiskPeak {
    /// End observation after the phase's final file mutation.
    pub(crate) fn finish(mut self) -> Result<u64> {
        let budget = self.budget.take().context("Disk phase already finished")?;
        let mut state = budget.state();
        let peak = state.disk_peak.take().context("Disk phase is not active")?;
        ensure!(state.valid, "Resource accounting is invalid");
        Ok(peak)
    }
}

impl Drop for DiskPeak {
    fn drop(&mut self) {
        if let Some(budget) = &self.budget {
            budget.state().disk_peak = None;
        }
    }
}

impl ResourceBudget {
    /// Validate limits before any stage starts reserving memory or disk.
    pub fn new(limits: ResourceLimits) -> Result<Arc<Self>> {
        Self::build(limits, None, None, None)
    }

    /// Apply actual process observations before any count/replay file opens.
    pub(crate) fn for_counting(
        mut limits: ResourceLimits,
        output_lease_allowance: u64,
    ) -> Result<Arc<Self>> {
        let file_limits =
            ProcessFileLimits::detect(limits.temporary_open_files, output_lease_allowance)?;
        limits.temporary_open_files = file_limits.temporary_file_limit;
        Self::build(limits, Some(file_limits), None, None)
    }

    /// A child ceiling backed by a complete parent reservation. The parent
    /// charge and completion callback survive until the last child owner drops.
    pub(crate) fn scoped(
        limits: ResourceLimits,
        parent_memory: MemoryReservation,
        scope_finished: Box<dyn FnOnce() + Send + Sync + std::panic::RefUnwindSafe>,
    ) -> Result<Arc<Self>> {
        ensure!(
            limits.memory_bytes <= parent_memory.bytes,
            "Child memory limit exceeds its parent reservation"
        );
        Self::build(limits, None, Some(parent_memory), Some(scope_finished))
    }

    /// Report performed startup observations, never a fabricated OS census.
    pub(crate) fn process_file_limits(&self) -> Option<ProcessFileLimits> {
        self.process_file_limits
    }

    fn build(
        limits: ResourceLimits,
        process_file_limits: Option<ProcessFileLimits>,
        parent_memory: Option<MemoryReservation>,
        scope_finished: Option<Box<dyn FnOnce() + Send + Sync + std::panic::RefUnwindSafe>>,
    ) -> Result<Arc<Self>> {
        ensure!(
            limits.memory_bytes > 0,
            "Run memory budget must be positive"
        );
        ensure!(
            limits.process_overhead_bytes <= limits.memory_bytes,
            "Process overhead exceeds run memory budget"
        );
        ensure!(
            limits
                .class_bytes
                .iter()
                .all(|limit| *limit <= limits.memory_bytes),
            "Stage memory ceiling exceeds run budget"
        );
        let memory_bytes = limits.process_overhead_bytes;
        Ok(Arc::new(Self {
            limits,
            process_file_limits,
            parent_memory,
            scope_finished,
            state: Mutex::new(BudgetState {
                snapshot: ResourceSnapshot {
                    reserved_memory_bytes: memory_bytes,
                    peak_reserved_memory_bytes: memory_bytes,
                    reserved_by_class: [0; CLASSES],
                    peak_reserved_by_class: [0; CLASSES],
                    reserved_temporary_disk_bytes: 0,
                    peak_reserved_temporary_disk_bytes: 0,
                    reserved_temporary_open_files: 0,
                    peak_reserved_temporary_open_files: 0,
                    valid: true,
                },
                memory_peak: None,
                disk_peak: None,
            }),
        }))
    }

    // No user callback runs while this lock is held, and validation precedes
    // mutation. Recovering after unwinding lets guards return reservations.
    fn state(&self) -> MutexGuard<'_, BudgetState> {
        self.state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Return one owner's ceiling (not necessarily currently available bytes).
    pub fn class_limit(&self, class: MemoryClass) -> u64 {
        self.limits.class_bytes[class as usize]
    }

    /// Return a consistent usage snapshot.
    pub fn snapshot(&self) -> ResourceSnapshot {
        self.state().snapshot.clone()
    }

    /// Include live baseline owners while excluding completed earlier phases.
    pub(crate) fn begin_memory_peak(self: &Arc<Self>) -> Result<MemoryPeak> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        ensure!(
            state.memory_peak.is_none(),
            "A memory phase is already active"
        );
        state.memory_peak = Some(state.reserved_memory_bytes);
        Ok(MemoryPeak {
            budget: Some(Arc::clone(self)),
        })
    }

    /// Include live baseline files, excluding completed startup/publication.
    pub(crate) fn begin_disk_peak(self: &Arc<Self>) -> Result<DiskPeak> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        ensure!(state.disk_peak.is_none(), "A disk phase is already active");
        state.disk_peak = Some(state.reserved_temporary_disk_bytes);
        Ok(DiskPeak {
            budget: Some(Arc::clone(self)),
        })
    }

    fn resize_memory(&self, class: MemoryClass, old: u64, new: u64) -> Result<()> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        let index = class as usize;
        let class_bytes = u128::from(
            state.reserved_by_class[index]
                .checked_sub(old)
                .context("Stage reservation release exceeds live charge")?,
        ) + u128::from(new);
        let memory_bytes = u128::from(
            state
                .reserved_memory_bytes
                .checked_sub(old)
                .context("Run reservation release exceeds live charge")?,
        ) + u128::from(new);
        if class_bytes > u128::from(self.class_limit(class)) {
            return Err(ResourceExhausted {
                resource: ResourceKind::StageMemory(class),
                required_bytes: class_bytes,
                limit_bytes: self.class_limit(class),
            }
            .into());
        }
        if memory_bytes > u128::from(self.limits.memory_bytes) {
            return Err(ResourceExhausted {
                resource: ResourceKind::RunMemory,
                required_bytes: memory_bytes,
                limit_bytes: self.limits.memory_bytes,
            }
            .into());
        }
        // Both projections are bounded by checked u64 limits before casting.
        let class_bytes = class_bytes as u64;
        let memory_bytes = memory_bytes as u64;
        state.reserved_by_class[index] = class_bytes;
        state.peak_reserved_by_class[index] = state.peak_reserved_by_class[index].max(class_bytes);
        state.reserved_memory_bytes = memory_bytes;
        state.peak_reserved_memory_bytes = state.peak_reserved_memory_bytes.max(memory_bytes);
        if let Some(peak) = state.memory_peak.as_mut() {
            *peak = (*peak).max(memory_bytes);
        }
        Ok(())
    }

    /// Reserve before allocating. The guard releases its charge on drop.
    pub fn reserve(self: &Arc<Self>, class: MemoryClass, bytes: u64) -> Result<MemoryReservation> {
        self.resize_memory(class, 0, bytes)?;
        Ok(MemoryReservation {
            budget: Arc::clone(self),
            class,
            bytes,
        })
    }

    fn resize_disk(&self, old: u64, new: u64) -> Result<()> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        let bytes = u128::from(
            state
                .reserved_temporary_disk_bytes
                .checked_sub(old)
                .context("Temporary disk release exceeds live charge")?,
        ) + u128::from(new);
        if bytes > u128::from(self.limits.temporary_disk_bytes) {
            return Err(ResourceExhausted {
                resource: ResourceKind::TemporaryDisk,
                required_bytes: bytes,
                limit_bytes: self.limits.temporary_disk_bytes,
            }
            .into());
        }
        let bytes = bytes as u64;
        state.reserved_temporary_disk_bytes = bytes;
        state.peak_reserved_temporary_disk_bytes =
            state.peak_reserved_temporary_disk_bytes.max(bytes);
        if let Some(peak) = state.disk_peak.as_mut() {
            *peak = (*peak).max(bytes);
        }
        Ok(())
    }

    /// Reserve space before a temporary-file write; retained files keep guards.
    pub fn reserve_disk(self: &Arc<Self>, bytes: u64) -> Result<DiskReservation> {
        self.resize_disk(0, bytes)?;
        Ok(DiskReservation {
            budget: Arc::clone(self),
            bytes,
        })
    }

    fn resize_temporary_files(&self, old: u64, new: u64) -> Result<()> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        let files = u128::from(
            state
                .reserved_temporary_open_files
                .checked_sub(old)
                .context("Temporary descriptor release exceeds live charge")?,
        ) + u128::from(new);
        if files > u128::from(self.limits.temporary_open_files) {
            return Err(TemporaryFilesExhausted {
                required_files: files,
                limit_files: self.limits.temporary_open_files,
            }
            .into());
        }
        state.reserved_temporary_open_files = files as u64;
        state.peak_reserved_temporary_open_files =
            state.peak_reserved_temporary_open_files.max(files as u64);
        Ok(())
    }

    /// Reserve one descriptor before creating an owned count/replay file.
    pub fn reserve_temporary_file(self: &Arc<Self>) -> Result<TemporaryFileReservation> {
        self.resize_temporary_files(0, 1)?;
        Ok(TemporaryFileReservation {
            budget: Arc::clone(self),
        })
    }
}

impl Drop for ResourceBudget {
    fn drop(&mut self) {
        // Child guards hold an Arc to this budget. Their allocations and
        // reservations have closed before the parent allowance is released.
        drop(self.parent_memory.take());
        if let Some(finished) = self.scope_finished.take() {
            finished();
        }
    }
}

/// Non-cloneable ownership of one live temporary-file descriptor.
pub struct TemporaryFileReservation {
    budget: Arc<ResourceBudget>,
}

impl Drop for TemporaryFileReservation {
    fn drop(&mut self) {
        if self.budget.resize_temporary_files(1, 0).is_err() {
            self.budget.state().valid = false;
            log::error!(
                "Failed to release temporary descriptor reservation; resource accounting disabled"
            );
        }
    }
}

/// Non-cloneable ownership of a live memory reservation.
pub struct MemoryReservation {
    budget: Arc<ResourceBudget>,
    class: MemoryClass,
    bytes: u64,
}

impl MemoryReservation {
    /// Adjust the same live reservation; failure preserves its previous charge.
    pub fn resize(&mut self, bytes: u64) -> Result<()> {
        self.budget.resize_memory(self.class, self.bytes, bytes)?;
        self.bytes = bytes;
        Ok(())
    }
}

impl Drop for MemoryReservation {
    fn drop(&mut self) {
        if self
            .budget
            .resize_memory(self.class, self.bytes, 0)
            .is_err()
        {
            self.budget.state().valid = false;
            log::error!("Failed to release memory reservation; resource accounting disabled");
        }
    }
}

/// Non-cloneable ownership of temporary disk space until file deletion.
pub struct DiskReservation {
    budget: Arc<ResourceBudget>,
    bytes: u64,
}

impl DiskReservation {
    /// Bytes still owned by this reservation, including retained allowances.
    pub(crate) fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Reserve a file's new total size before extending it.
    pub fn resize(&mut self, bytes: u64) -> Result<()> {
        self.budget.resize_disk(self.bytes, bytes)?;
        self.bytes = bytes;
        Ok(())
    }
}

impl Drop for DiskReservation {
    fn drop(&mut self) {
        if self.budget.resize_disk(self.bytes, 0).is_err() {
            self.budget.state().valid = false;
            log::error!(
                "Failed to release temporary disk reservation; resource accounting disabled"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        Barrier,
        atomic::{AtomicUsize, Ordering},
    };

    fn budget() -> Arc<ResourceBudget> {
        let mut classes = [0; CLASSES];
        classes[MemoryClass::CountSlots as usize] = 800;
        classes[MemoryClass::TargetGraphs as usize] = 500;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 1000,
            process_overhead_bytes: 100,
            class_bytes: classes,
            temporary_disk_bytes: 500,
            temporary_open_files: 128,
        })
        .unwrap()
    }

    #[test]
    fn simultaneous_stage_reservations_and_disk_have_independent_limits() {
        let budget = budget();
        let mut counter = budget.reserve(MemoryClass::CountSlots, 600).unwrap();
        assert!(budget.reserve(MemoryClass::TargetGraphs, 301).is_err());
        let graph = budget.reserve(MemoryClass::TargetGraphs, 300).unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 1000);
        assert!(counter.resize(601).is_err());
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::CountSlots as usize],
            600
        );
        let mut disk = budget.reserve_disk(400).unwrap();
        assert!(disk.resize(501).is_err());
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 400);
        assert!(budget.reserve(MemoryClass::ReplayBuffer, 1).is_err());
        drop(graph);
        counter.resize(800).unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 900);
        drop(counter);
        drop(disk);
        let snapshot = budget.snapshot();
        assert_eq!(snapshot.reserved_memory_bytes, 100);
        assert_eq!(snapshot.reserved_temporary_disk_bytes, 0);
        assert_eq!(snapshot.peak_reserved_memory_bytes, 1000);
        assert_eq!(snapshot.peak_reserved_temporary_disk_bytes, 400);
        assert!(snapshot.valid);
    }

    #[test]
    fn typed_exhaustion_survives_context_and_preserves_live_reservations() {
        let budget = budget();
        let mut counter = budget.reserve(MemoryClass::CountSlots, 600).unwrap();
        let stage_error = budget
            .reserve(MemoryClass::CountSlots, 201)
            .err()
            .unwrap()
            .context("adding another count owner");
        assert_eq!(
            stage_error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::StageMemory(MemoryClass::CountSlots),
                required_bytes: 801,
                limit_bytes: 800,
            })
        );
        let run_error = budget
            .reserve(MemoryClass::TargetGraphs, 301)
            .err()
            .unwrap();
        assert_eq!(
            run_error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::RunMemory,
                required_bytes: 1001,
                limit_bytes: 1000,
            })
        );
        let graph = budget.reserve(MemoryClass::TargetGraphs, 300).unwrap();
        let resize_error = counter.resize(601).unwrap_err();
        assert_eq!(
            resize_error.downcast_ref::<ResourceExhausted>(),
            run_error.downcast_ref::<ResourceExhausted>()
        );
        let mut disk = budget.reserve_disk(400).unwrap();
        let disk_error = disk.resize(501).unwrap_err();
        assert_eq!(
            disk_error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::TemporaryDisk,
                required_bytes: 501,
                limit_bytes: 500,
            })
        );
        let snapshot = budget.snapshot();
        assert_eq!(snapshot.reserved_memory_bytes, 1000);
        assert_eq!(snapshot.peak_reserved_memory_bytes, 1000);
        assert_eq!(
            snapshot.reserved_by_class[MemoryClass::CountSlots as usize],
            600
        );
        assert_eq!(
            snapshot.reserved_by_class[MemoryClass::TargetGraphs as usize],
            300
        );
        assert_eq!(snapshot.reserved_temporary_disk_bytes, 400);
        assert_eq!(snapshot.peak_reserved_temporary_disk_bytes, 400);
        drop((counter, graph, disk));
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn invalid_accounting_is_not_recoverable_resource_exhaustion() {
        let budget = budget();
        budget.state().valid = false;
        let error = budget.reserve(MemoryClass::CountSlots, 1).err().unwrap();
        assert!(error.downcast_ref::<ResourceExhausted>().is_none());
        assert!(error.to_string().contains("accounting is invalid"));
        assert!(
            budget
                .reserve_disk(1)
                .err()
                .unwrap()
                .downcast_ref::<ResourceExhausted>()
                .is_none()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
    }

    #[test]
    fn competing_workers_cannot_overcommit_and_guards_return_on_unwind() {
        let budget = budget();
        let barrier = Arc::new(Barrier::new(8));
        let successes = Arc::new(AtomicUsize::new(0));
        std::thread::scope(|scope| {
            for _ in 0..8 {
                let budget = Arc::clone(&budget);
                let barrier = Arc::clone(&barrier);
                let successes = Arc::clone(&successes);
                scope.spawn(move || {
                    let guard = budget.reserve(MemoryClass::TargetGraphs, 100).ok();
                    if guard.is_some() {
                        successes.fetch_add(1, Ordering::SeqCst);
                    }
                    barrier.wait();
                    drop(guard);
                });
            }
        });
        assert_eq!(successes.load(Ordering::SeqCst), 5);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
        let panicked = std::panic::catch_unwind(|| {
            let _guard = budget.reserve(MemoryClass::CountSlots, 400).unwrap();
            panic!("exercise owner unwind");
        });
        assert!(panicked.is_err());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn invalid_and_overflowing_limits_do_not_create_reservations() {
        assert!(
            ResourceBudget::new(ResourceLimits {
                memory_bytes: 0,
                process_overhead_bytes: 0,
                class_bytes: [0; CLASSES],
                temporary_disk_bytes: 0,
                temporary_open_files: 0,
            })
            .is_err()
        );
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: u64::MAX,
            process_overhead_bytes: 1,
            class_bytes: [u64::MAX; CLASSES],
            temporary_disk_bytes: u64::MAX,
            temporary_open_files: 128,
        })
        .unwrap();
        let error = budget
            .reserve(MemoryClass::InputBuffer, u64::MAX)
            .err()
            .unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::RunMemory,
                required_bytes: u128::from(u64::MAX) + 1,
                limit_bytes: u64::MAX,
            })
        );
        let disk = budget.reserve_disk(u64::MAX).unwrap();
        let error = budget.reserve_disk(1).err().unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>(),
            Some(&ResourceExhausted {
                resource: ResourceKind::TemporaryDisk,
                required_bytes: u128::from(u64::MAX) + 1,
                limit_bytes: u64::MAX,
            })
        );
        drop(disk);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 1);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }
    #[test]
    fn temporary_descriptors_are_shared_and_rejected_requests_preserve_live_usage() {
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1,
            process_overhead_bytes: 0,
            class_bytes: [0; CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 2,
        })
        .unwrap();
        let first = budget.reserve_temporary_file().unwrap();
        let second = budget.reserve_temporary_file().unwrap();
        let before = budget.snapshot();
        let error = budget
            .reserve_temporary_file()
            .err()
            .unwrap()
            .context("merge output");
        assert_eq!(
            error.downcast_ref::<TemporaryFilesExhausted>(),
            Some(&TemporaryFilesExhausted {
                required_files: 3,
                limit_files: 2,
            })
        );
        assert_eq!(
            budget.snapshot().reserved_temporary_open_files,
            before.reserved_temporary_open_files
        );
        assert_eq!(
            budget.snapshot().peak_reserved_temporary_open_files,
            before.peak_reserved_temporary_open_files
        );
        drop(first);
        let replacement = budget.reserve_temporary_file().unwrap();
        assert_eq!(budget.snapshot().reserved_temporary_open_files, 2);
        drop((replacement, second));
        let snapshot = budget.snapshot();
        assert_eq!(snapshot.reserved_temporary_open_files, 0);
        assert_eq!(snapshot.peak_reserved_temporary_open_files, 2);
        assert_eq!(snapshot.reserved_memory_bytes, 0);
        assert_eq!(snapshot.reserved_temporary_disk_bytes, 0);
        assert!(snapshot.valid);
    }

    #[test]
    fn concurrent_file_owners_obey_one_shared_limit() {
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 1,
            process_overhead_bytes: 0,
            class_bytes: [0; CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 3,
        })
        .unwrap();
        let barrier = Arc::new(Barrier::new(12));
        let threads = (0..12)
            .map(|_| {
                let budget = Arc::clone(&budget);
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    let reservation = budget.reserve_temporary_file();
                    if let Err(error) = &reservation {
                        assert!(error.downcast_ref::<TemporaryFilesExhausted>().is_some());
                    }
                    // Every successful owner stays live until all requests finish.
                    barrier.wait();
                    reservation.is_ok()
                })
            })
            .collect::<Vec<_>>();
        let accepted = threads
            .into_iter()
            .map(|thread| usize::from(thread.join().unwrap()))
            .sum::<usize>();
        assert_eq!(accepted, 3);
        let snapshot = budget.snapshot();
        assert_eq!(snapshot.reserved_temporary_open_files, 0);
        assert_eq!(snapshot.peak_reserved_temporary_open_files, 3);
        assert!(snapshot.valid);
    }

    #[test]
    fn disk_phase_peak_keeps_live_files_and_excludes_finished_metadata() {
        let budget = budget();
        drop(budget.reserve_disk(400).unwrap());
        let baseline = budget.reserve_disk(25).unwrap();
        let phase = budget.begin_disk_peak().unwrap();
        assert!(budget.begin_disk_peak().is_err());
        let barrier = Arc::new(Barrier::new(5));
        std::thread::scope(|scope| {
            for _ in 0..4 {
                let budget = Arc::clone(&budget);
                let barrier = Arc::clone(&barrier);
                scope.spawn(move || {
                    let _file = budget.reserve_disk(50).unwrap();
                    barrier.wait();
                    barrier.wait();
                });
            }
            barrier.wait();
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 225);
            barrier.wait();
        });
        assert_eq!(phase.finish().unwrap(), 225);
        assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 400);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 25);
        drop(baseline);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn disk_phase_drop_cleans_on_refusal_and_unwind_without_resetting_lifetime_peak() {
        let budget = budget();
        let failed = (|| -> Result<()> {
            let _phase = budget.begin_disk_peak()?;
            let _file = budget.reserve_disk(501)?;
            Ok(())
        })();
        assert!(failed.is_err());
        let panicked = std::panic::catch_unwind(|| {
            let _phase = budget.begin_disk_peak().unwrap();
            let _file = budget.reserve_disk(450).unwrap();
            panic!("test disk-phase unwind");
        });
        assert!(panicked.is_err());
        assert_eq!(budget.begin_disk_peak().unwrap().finish().unwrap(), 0);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
        assert_eq!(budget.snapshot().peak_reserved_temporary_disk_bytes, 450);
        assert!(budget.snapshot().valid);
        assert_eq!(Arc::strong_count(&budget), 1);
    }

    #[test]
    fn memory_phase_peak_keeps_live_baseline_and_excludes_finished_startup() {
        let budget = budget();
        drop(budget.reserve(MemoryClass::CountSlots, 600).unwrap());
        let baseline = budget.reserve(MemoryClass::TargetGraphs, 50).unwrap();
        let phase = budget.begin_memory_peak().unwrap();
        assert!(budget.begin_memory_peak().is_err());
        let barrier = Arc::new(Barrier::new(5));
        std::thread::scope(|scope| {
            for _ in 0..4 {
                let budget = Arc::clone(&budget);
                let barrier = Arc::clone(&barrier);
                scope.spawn(move || {
                    let _owner = budget.reserve(MemoryClass::CountSlots, 100).unwrap();
                    barrier.wait();
                    barrier.wait();
                });
            }
            barrier.wait();
            assert_eq!(budget.snapshot().reserved_memory_bytes, 550);
            barrier.wait();
        });
        assert_eq!(phase.finish().unwrap(), 550);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 700);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 150);
        drop(baseline);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
    }

    #[test]
    fn memory_phase_drop_cleans_on_refusal_and_unwind_without_resetting_lifetime_peak() {
        let budget = budget();
        let failed = (|| -> Result<()> {
            let _phase = budget.begin_memory_peak()?;
            let _owner = budget.reserve(MemoryClass::CountSlots, u64::MAX)?;
            Ok(())
        })();
        assert!(failed.is_err());
        let panicked = std::panic::catch_unwind(|| {
            let _phase = budget.begin_memory_peak().unwrap();
            let _owner = budget.reserve(MemoryClass::CountSlots, 600).unwrap();
            panic!("test phase unwind");
        });
        assert!(panicked.is_err());
        assert_eq!(budget.begin_memory_peak().unwrap().finish().unwrap(), 100);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 100);
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 700);
        assert!(budget.snapshot().valid);
        let phase = budget.begin_memory_peak().unwrap();
        budget.state().valid = false;
        let err = phase.finish().unwrap_err();
        assert!(err.downcast_ref::<ResourceExhausted>().is_none());
        assert!(err.to_string().contains("accounting is invalid"));
        budget.state().valid = true;
        assert_eq!(budget.begin_memory_peak().unwrap().finish().unwrap(), 100);
    }
}
