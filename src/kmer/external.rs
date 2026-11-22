//! Exact external sorting from a bounded reusable buffer and bounded run slots.
//!
//! Online binary compaction keeps at most 64 run slots for a u64-sized stream.
//! Each immutable run owns one charged page and one temporary file. Merge
//! outputs remain separate from their inputs until completely written/verified.

use super::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use super::identity::{CanonicalKey, KmerSize};
use super::sorted_file::{DiskCountStore, PAGE_BYTES, SortedCountWriter};
use super::store::CountStore;
use anyhow::{Context, Result, ensure};
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};
use std::sync::Arc;

const RUN_SLOTS: usize = 64;

/// Actual producer work, distinct from emitted biological sequences.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ExternalMetrics {
    pub(crate) peak_table_bytes: usize,
    pub(crate) sorted_runs: u64,
    pub(crate) merges: u64,
}

/// Streaming exact producer; never collects the global distinct-key table.
pub(crate) struct ExternalCountBuilder {
    k: KmerSize,
    directory: PathBuf,
    budget: Arc<ResourceBudget>,
    pairs: Vec<(u64, u64)>,
    pair_limit: usize,
    runs: Vec<Option<DiskCountStore>>,
    occurrences: u64,
    metrics: ExternalMetrics,
    max_lookup_batch: usize,
    failed: bool,
    sort_reservation: MemoryReservation,
    _run_reservation: MemoryReservation,
}

impl ExternalCountBuilder {
    /// Fixed upper bound for the slots themselves; page buffers are separate.
    pub(crate) fn run_slot_bytes() -> Result<u64> {
        u64::try_from(RUN_SLOTS * size_of::<Option<DiskCountStore>>())
            .context("Run-slot byte size overflow")
    }

    /// Bound all retained run pages plus a carry and its merge output.
    pub(crate) fn cache_limit_bytes() -> u64 {
        (RUN_SLOTS as u64 + 2) * PAGE_BYTES as u64
    }

    /// Allocate the sort buffer and fixed run slots under the shared ledger.
    pub(crate) fn new(
        directory: &Path,
        k: KmerSize,
        budget: Arc<ResourceBudget>,
        sort_bytes: usize,
        max_lookup_batch: usize,
    ) -> Result<Self> {
        let pair_limit = sort_bytes / size_of::<(u64, u64)>();
        ensure!(
            pair_limit > 0,
            "External sort buffer must hold a key/count pair"
        );
        ensure!(
            max_lookup_batch > 0,
            "Count lookup batch reservation must be positive"
        );
        let mut sort_reservation = budget.reserve(
            MemoryClass::SortBuffer,
            u64::try_from(pair_limit * size_of::<(u64, u64)>())?,
        )?;
        let mut pairs = Vec::new();
        pairs
            .try_reserve_exact(pair_limit)
            .context("Unable to reserve external sort buffer")?;
        sort_reservation.resize(u64::try_from(pairs.capacity() * size_of::<(u64, u64)>())?)?;
        let mut run_reservation =
            budget.reserve(MemoryClass::WorkerQueues, Self::run_slot_bytes()?)?;
        let mut runs = Vec::new();
        runs.try_reserve_exact(RUN_SLOTS)
            .context("Unable to reserve external run slots")?;
        run_reservation.resize(u64::try_from(
            runs.capacity() * size_of::<Option<DiskCountStore>>(),
        )?)?;
        runs.resize_with(RUN_SLOTS, || None);
        Ok(Self {
            k,
            directory: directory.to_path_buf(),
            budget,
            pairs,
            pair_limit,
            runs,
            occurrences: 0,
            metrics: ExternalMetrics {
                peak_table_bytes: 0,
                sorted_runs: 0,
                merges: 0,
            },
            max_lookup_batch,
            failed: false,
            sort_reservation,
            _run_reservation: run_reservation,
        })
    }

    /// Return the checked total of accepted windows, including spilled runs.
    pub(crate) fn occurrences(&self) -> u64 {
        self.occurrences
    }

    /// Invalid identities/overflow preserve the accepted prefix. A spill I/O
    /// or resource failure invalidates finalization instead of publishing a subset.
    pub(crate) fn add(&mut self, key: CanonicalKey, count: NonZeroU64) -> Result<()> {
        ensure!(!self.failed, "External count producer previously failed");
        ensure!(
            key.size() == self.k,
            "Count key k differs from external producer k"
        );
        let total = self
            .occurrences
            .checked_add(count.get())
            .context("Exact count occurrence total overflow")?;
        if self.pairs.len() == self.pair_limit
            && let Err(error) = self.flush_run()
        {
            self.failed = true;
            return Err(error);
        }
        self.pairs.push((key.bits(), count.get()));
        self.occurrences = total;
        Ok(())
    }

    fn writer(&self) -> Result<SortedCountWriter> {
        SortedCountWriter::new(&self.directory, self.k, &self.budget, self.max_lookup_batch)
    }

    fn merge(&mut self, left: DiskCountStore, right: DiskCountStore) -> Result<DiskCountStore> {
        ensure!(
            left.metadata().k == self.k && right.metadata().k == self.k,
            "External run identity mismatch"
        );
        let next_merge = self
            .metrics
            .merges
            .checked_add(1)
            .context("External merge count overflow")?;
        let mut writer = self.writer()?;
        let (mut li, mut ri) = (0, 0);
        let next = |store: &DiskCountStore,
                    index: &mut u64|
         -> Result<Option<(CanonicalKey, NonZeroU64)>> {
            if *index == store.metadata().distinct_keys {
                return Ok(None);
            }
            let item = store.entry_at(*index)?;
            *index += 1;
            Ok(Some(item))
        };
        let mut a = next(&left, &mut li)?;
        let mut b = next(&right, &mut ri)?;
        loop {
            match (a, b) {
                (None, None) => break,
                (Some((key, count)), None) => {
                    writer.push(key, count)?;
                    a = next(&left, &mut li)?;
                }
                (None, Some((key, count))) => {
                    writer.push(key, count)?;
                    b = next(&right, &mut ri)?;
                }
                (Some((ka, ca)), Some((kb, cb))) => match ka.bits().cmp(&kb.bits()) {
                    std::cmp::Ordering::Less => {
                        writer.push(ka, ca)?;
                        a = next(&left, &mut li)?;
                    }
                    std::cmp::Ordering::Greater => {
                        writer.push(kb, cb)?;
                        b = next(&right, &mut ri)?;
                    }
                    std::cmp::Ordering::Equal => {
                        let sum = ca
                            .get()
                            .checked_add(cb.get())
                            .context("Exact count addition overflow during merge")?;
                        writer.push(ka, NonZeroU64::new(sum).context("Merged count is zero")?)?;
                        a = next(&left, &mut li)?;
                        b = next(&right, &mut ri)?;
                    }
                },
            }
        }
        let store = writer.finish([0; 32])?;
        self.metrics.merges = next_merge;
        // Inputs remain live until the complete output passes format validation.
        left.close()?;
        right.close()?;
        Ok(store)
    }

    fn insert_run(&mut self, mut run: DiskCountStore) -> Result<()> {
        for index in 0..RUN_SLOTS {
            if let Some(previous) = self.runs[index].take() {
                run = self.merge(previous, run)?;
            } else {
                self.runs[index] = Some(run);
                return Ok(());
            }
        }
        anyhow::bail!("External run count exceeds the u64 stream bound")
    }

    fn flush_run(&mut self) -> Result<()> {
        if self.pairs.is_empty() {
            return Ok(());
        }
        let next_run = self
            .metrics
            .sorted_runs
            .checked_add(1)
            .context("External sorted-run count overflow")?;
        self.pairs.sort_unstable_by_key(|&(key, _)| key);
        let mut writer = self.writer()?;
        let mut index = 0;
        while index < self.pairs.len() {
            let (bits, mut total) = self.pairs[index];
            index += 1;
            while index < self.pairs.len() && self.pairs[index].0 == bits {
                total = total
                    .checked_add(self.pairs[index].1)
                    .context("Exact count addition overflow during reduction")?;
                index += 1;
            }
            writer.push(
                CanonicalKey::new(bits, self.k)?,
                NonZeroU64::new(total).context("Reduced count is zero")?,
            )?;
        }
        let run = writer.finish([0; 32])?;
        self.insert_run(run)?;
        self.pairs.clear();
        self.metrics.sorted_runs = next_run;
        Ok(())
    }

    /// Reduce all runs, bind the complete selected-input digest and move the
    /// final disk store into PCR without materializing a global memory table.
    pub(crate) fn finish(
        mut self,
        selected_input_sha256: [u8; 32],
    ) -> Result<(DiskCountStore, ExternalMetrics)> {
        ensure!(!self.failed, "External count producer previously failed");
        self.flush_run()?;
        self.pairs = Vec::new();
        self.sort_reservation.resize(0)?;
        let mut reduced = None;
        for index in 0..RUN_SLOTS {
            if let Some(run) = self.runs[index].take() {
                reduced = Some(match reduced {
                    None => run,
                    Some(previous) => self.merge(previous, run)?,
                });
            }
        }
        // Intermediate files bind complete run contents, not the whole input.
        // Re-streaming the final reduction binds the actual selected-input digest.
        let mut writer = self.writer()?;
        if let Some(run) = reduced.as_ref() {
            ensure!(
                run.metadata().occurrences == self.occurrences,
                "External reduced total differs from selected windows"
            );
            run.visit(&mut |key, count| writer.push(key, count))?;
        }
        let store = writer.finish(selected_input_sha256)?;
        ensure!(
            store.metadata().occurrences == self.occurrences,
            "Final external total differs from selected windows"
        );
        if let Some(run) = reduced {
            run.close()?;
        }
        Ok((store, self.metrics))
    }
}

#[cfg(test)]
mod tests;
