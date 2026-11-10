//! Shared logical memory/disk reservations for v4 stages.
//!
//! Reservations are ownership guards, not allocator hooks or RSS measurements.
//! Production limits must include measured process/allocator overhead.

use anyhow::{Context, Result, ensure};
use std::sync::{Arc, Mutex, MutexGuard};

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
}

const CLASSES: usize = 8;

/// Configured ceilings; simultaneous live reservations must also fit total RAM.
#[derive(Clone, Debug)]
pub struct ResourceLimits {
    pub memory_bytes: u64,
    pub process_overhead_bytes: u64,
    pub class_bytes: [u64; CLASSES],
    pub temporary_disk_bytes: u64,
}

/// Observable logical usage/peaks, including the initial process reservation.
#[derive(Clone, Debug)]
pub struct ResourceSnapshot {
    pub reserved_memory_bytes: u64,
    pub peak_reserved_memory_bytes: u64,
    pub reserved_by_class: [u64; CLASSES],
    pub reserved_temporary_disk_bytes: u64,
    pub peak_reserved_temporary_disk_bytes: u64,
    pub valid: bool,
}

/// A shared run budget. Clones share accounting across stages and workers.
pub struct ResourceBudget {
    limits: ResourceLimits,
    state: Mutex<ResourceSnapshot>,
}

impl ResourceBudget {
    /// Validate limits before any stage starts reserving memory or disk.
    pub fn new(limits: ResourceLimits) -> Result<Arc<Self>> {
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
            state: Mutex::new(ResourceSnapshot {
                reserved_memory_bytes: memory_bytes,
                peak_reserved_memory_bytes: memory_bytes,
                reserved_by_class: [0; CLASSES],
                reserved_temporary_disk_bytes: 0,
                peak_reserved_temporary_disk_bytes: 0,
                valid: true,
            }),
        }))
    }

    // No user callback runs while this lock is held, and validation precedes
    // mutation. Recovering after unwinding lets guards return reservations.
    fn state(&self) -> MutexGuard<'_, ResourceSnapshot> {
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
        self.state().clone()
    }

    fn resize_memory(&self, class: MemoryClass, old: u64, new: u64) -> Result<()> {
        let mut state = self.state();
        ensure!(state.valid, "Resource accounting is invalid");
        let index = class as usize;
        let class_bytes = state.reserved_by_class[index]
            .checked_sub(old)
            .and_then(|n| n.checked_add(new))
            .context("Stage reservation arithmetic overflow")?;
        let memory_bytes = state
            .reserved_memory_bytes
            .checked_sub(old)
            .and_then(|n| n.checked_add(new))
            .context("Run reservation arithmetic overflow")?;
        ensure!(
            class_bytes <= self.class_limit(class),
            "{class:?} memory reservation exceeds its ceiling"
        );
        ensure!(
            memory_bytes <= self.limits.memory_bytes,
            "Combined live memory reservations exceed the run limit"
        );
        state.reserved_by_class[index] = class_bytes;
        state.reserved_memory_bytes = memory_bytes;
        state.peak_reserved_memory_bytes = state.peak_reserved_memory_bytes.max(memory_bytes);
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
        let bytes = state
            .reserved_temporary_disk_bytes
            .checked_sub(old)
            .and_then(|n| n.checked_add(new))
            .context("Temporary disk reservation arithmetic overflow")?;
        ensure!(
            bytes <= self.limits.temporary_disk_bytes,
            "Temporary disk reservation exceeds its limit"
        );
        state.reserved_temporary_disk_bytes = bytes;
        state.peak_reserved_temporary_disk_bytes =
            state.peak_reserved_temporary_disk_bytes.max(bytes);
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
                temporary_disk_bytes: 0
            })
            .is_err()
        );
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: u64::MAX,
            process_overhead_bytes: 1,
            class_bytes: [u64::MAX; CLASSES],
            temporary_disk_bytes: u64::MAX,
        })
        .unwrap();
        assert!(budget.reserve(MemoryClass::InputBuffer, u64::MAX).is_err());
        let disk = budget.reserve_disk(u64::MAX).unwrap();
        assert!(budget.reserve_disk(1).is_err());
        drop(disk);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 1);
        assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
    }
}
