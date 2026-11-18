//! Bounded target admission with per-target graph-buffer ownership.

use crate::kmer::budget::{
    MemoryClass, ResourceBudget, ResourceExhausted, ResourceKind, ResourceLimits,
};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::sync::{Arc, Condvar, Mutex};

/// Performed graph-buffer measurements for one target.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct TargetGraphBuffers {
    pub(crate) limit_bytes: u64,
    pub(crate) peak_reserved_buffer_bytes: u64,
}

/// Admission limits and observed concurrent target ownership.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct GraphBufferSummary {
    pub(crate) shared_limit_bytes: u64,
    pub(crate) per_target_limit_bytes: u64,
    pub(crate) concurrency_limit: usize,
    pub(crate) peak_active_targets: usize,
    pub(crate) peak_reserved_target_allowances_bytes: u64,
}

struct Usage {
    active: usize,
    peak: usize,
    valid: bool,
}

struct Inner {
    budget: Arc<ResourceBudget>,
    target_bytes: u64,
    slots: usize,
    usage: Mutex<Usage>,
    available: Condvar,
}

/// Share the same parent ledger as immutable counts; each admitted target
/// gets an independent child ceiling within its complete parent reservation.
pub(crate) struct GraphBufferPool {
    inner: Arc<Inner>,
}

impl GraphBufferPool {
    pub(crate) fn new(budget: Arc<ResourceBudget>, target_bytes: u64) -> Result<Self> {
        ensure!(
            target_bytes > 0,
            "Per-target graph-buffer limit must be positive"
        );
        let shared = budget.class_limit(MemoryClass::TargetGraphs);
        ensure!(
            target_bytes <= shared,
            "Per-target graph-buffer limit exceeds shared graph-buffer budget"
        );
        let slots = usize::try_from(shared / target_bytes)
            .context("Graph-buffer concurrency exceeds the supported range")?;
        Ok(Self {
            inner: Arc::new(Inner {
                budget,
                target_bytes,
                slots,
                usage: Mutex::new(Usage {
                    active: 0,
                    peak: 0,
                    valid: true,
                }),
                available: Condvar::new(),
            }),
        })
    }

    /// Borrow the parent for results that must outlive target admission.
    pub(crate) fn result_budget(&self) -> &Arc<ResourceBudget> {
        &self.inner.budget
    }

    pub(crate) fn acquire(&self) -> Result<TargetGraphLease> {
        let mut usage = self
            .inner
            .usage
            .lock()
            .map_err(|_| anyhow::anyhow!("Graph admission lock poisoned"))?;
        loop {
            ensure!(usage.valid, "Graph admission accounting is invalid");
            if usage.active < self.inner.slots {
                match self
                    .inner
                    .budget
                    .reserve(MemoryClass::TargetGraphs, self.inner.target_bytes)
                {
                    Ok(reservation) => {
                        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
                        classes[MemoryClass::TargetGraphs as usize] = self.inner.target_bytes;
                        let finished_pool = Arc::clone(&self.inner);
                        let budget = ResourceBudget::scoped(
                            ResourceLimits {
                                memory_bytes: self.inner.target_bytes,
                                process_overhead_bytes: 0,
                                class_bytes: classes,
                                temporary_disk_bytes: 0,
                                temporary_open_files: 0,
                            },
                            reservation,
                            Box::new(move || finished_pool.release()),
                        )?;
                        usage.active += 1;
                        usage.peak = usage.peak.max(usage.active);
                        return Ok(TargetGraphLease { budget });
                    }
                    Err(error) => {
                        let competing =
                            error
                                .downcast_ref::<ResourceExhausted>()
                                .is_some_and(|error| {
                                    matches!(
                                        error.resource,
                                        ResourceKind::RunMemory
                                            | ResourceKind::StageMemory(MemoryClass::TargetGraphs)
                                    )
                                });
                        // If another admitted target can release enough space,
                        // reduce concurrency instead of rejecting a target that
                        // fits by itself. Never wait on an impossible empty pool.
                        if !competing || usage.active == 0 {
                            return Err(error);
                        }
                    }
                }
            }
            usage = self
                .inner
                .available
                .wait(usage)
                .map_err(|_| anyhow::anyhow!("Graph admission lock poisoned"))?;
        }
    }

    pub(crate) fn summary(&self) -> Result<GraphBufferSummary> {
        let usage = self
            .inner
            .usage
            .lock()
            .map_err(|_| anyhow::anyhow!("Graph admission lock poisoned"))?;
        ensure!(
            usage.valid && usage.active == 0,
            "Graph admission is incomplete or invalid"
        );
        let snapshot = self.inner.budget.snapshot();
        ensure!(snapshot.valid, "Run resource accounting is invalid");
        ensure!(
            snapshot.reserved_by_class[MemoryClass::TargetGraphs as usize] == 0,
            "Graph buffer owners still hold target allowances"
        );
        Ok(GraphBufferSummary {
            shared_limit_bytes: self.inner.budget.class_limit(MemoryClass::TargetGraphs),
            per_target_limit_bytes: self.inner.target_bytes,
            concurrency_limit: self.inner.slots,
            peak_active_targets: usage.peak,
            peak_reserved_target_allowances_bytes: snapshot.peak_reserved_by_class
                [MemoryClass::TargetGraphs as usize],
        })
    }
}

/// Its parent allowance remains held through all graph work in this target.
pub(crate) struct TargetGraphLease {
    budget: Arc<ResourceBudget>,
}

impl TargetGraphLease {
    pub(crate) fn budget(&self) -> &Arc<ResourceBudget> {
        &self.budget
    }
    pub(crate) fn summary(&self) -> TargetGraphBuffers {
        TargetGraphBuffers {
            limit_bytes: self.budget.class_limit(MemoryClass::TargetGraphs),
            peak_reserved_buffer_bytes: self.budget.snapshot().peak_reserved_by_class
                [MemoryClass::TargetGraphs as usize],
        }
    }
}

impl Inner {
    fn release(&self) {
        // The scoped child has released its parent bytes before this callback.
        let mut usage = self
            .usage
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if usage.active == 0 {
            usage.valid = false;
            log::error!("Graph admission release exceeds active targets");
        } else {
            usage.active -= 1;
        }
        self.available.notify_all();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;
    use std::time::Duration;

    fn resources(total: u64, graph: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::TargetGraphs as usize] = graph;
        classes[MemoryClass::CountSlots as usize] = total;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: total,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn count_ownership_reduces_admission_and_dropping_a_target_wakes_the_next() {
        let budget = resources(600, 512);
        let counts = budget.reserve(MemoryClass::CountSlots, 200).unwrap();
        let pool = Arc::new(GraphBufferPool::new(Arc::clone(&budget), 256).unwrap());
        let first = pool.acquire().unwrap();
        assert!(pool.summary().is_err());
        let (started, ready) = mpsc::channel();
        let (entered, received) = mpsc::channel();
        let other = Arc::clone(&pool);
        let worker = std::thread::spawn(move || {
            started.send(()).unwrap();
            let lease = other.acquire().unwrap();
            entered.send(lease.summary()).unwrap();
        });
        ready.recv_timeout(Duration::from_secs(5)).unwrap();
        assert!(received.recv_timeout(Duration::from_millis(25)).is_err());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 456);
        drop(first);
        assert_eq!(
            received
                .recv_timeout(Duration::from_secs(5))
                .unwrap()
                .limit_bytes,
            256
        );
        worker.join().unwrap();
        assert_eq!(pool.summary().unwrap().peak_active_targets, 1);
        drop(counts);
        let first = pool.acquire().unwrap();
        let second = pool.acquire().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 512);
        drop(first);
        drop(second);
        assert_eq!(pool.summary().unwrap().peak_active_targets, 2);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn impossible_admission_fails_instead_of_waiting_and_child_growth_is_bounded() {
        let budget = resources(600, 512);
        let counts = budget.reserve(MemoryClass::CountSlots, 400).unwrap();
        let pool = GraphBufferPool::new(Arc::clone(&budget), 256).unwrap();
        let error = pool.acquire().err().unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::RunMemory
        );
        assert_eq!(pool.summary().unwrap().peak_active_targets, 0);
        drop(counts);
        let lease = pool.acquire().unwrap();
        let first = lease
            .budget()
            .reserve(MemoryClass::TargetGraphs, 200)
            .unwrap();
        assert!(
            lease
                .budget()
                .reserve(MemoryClass::TargetGraphs, 57)
                .is_err()
        );
        assert_eq!(lease.summary().peak_reserved_buffer_bytes, 200);
        drop(first);
        drop(lease);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn surviving_graph_owners_keep_parent_bytes_and_admission_until_their_drop() {
        let budget = resources(512, 256);
        let pool = Arc::new(GraphBufferPool::new(Arc::clone(&budget), 256).unwrap());
        let lease = pool.acquire().unwrap();
        let graph = lease
            .budget()
            .reserve(MemoryClass::TargetGraphs, 128)
            .unwrap();
        drop(lease);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 256);
        assert!(pool.summary().is_err());
        let (started, ready) = mpsc::channel();
        let (entered, received) = mpsc::channel();
        let other = Arc::clone(&pool);
        let worker = std::thread::spawn(move || {
            started.send(()).unwrap();
            let _lease = other.acquire().unwrap();
            entered.send(()).unwrap();
        });
        ready.recv_timeout(Duration::from_secs(5)).unwrap();
        assert!(received.recv_timeout(Duration::from_millis(25)).is_err());
        drop(graph);
        received.recv_timeout(Duration::from_secs(5)).unwrap();
        worker.join().unwrap();
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert_eq!(pool.summary().unwrap().peak_active_targets, 1);
    }
}
