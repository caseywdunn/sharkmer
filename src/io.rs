use crate::input_identity::{IdentityWorkerSummary, IdentityWriter, SelectedDigest};
use crate::kmer::adaptive::{AdaptiveCountBuilder, AdaptiveStore};
use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{
    DiskPeak, MemoryClass, MemoryPeak, MemoryReservation, ResourceBudget, ResourceLimits,
};
use crate::kmer::external::{ExternalCountBuilder, ExternalMetrics};
use crate::kmer::identity::{CanonicalKey, KmerSize, visit_canonical_windows};
use crate::kmer::packed::{PackedCountBuilder, PackedCountStore};
use crate::kmer::parallel::{
    AdaptiveParallelCountBuilder, AdaptiveShardedCountStore, ParallelCountBuilder, ParallelMetrics,
    ShardedCountStore,
};
use crate::kmer::sorted_file::DiskCountStore;
use crate::kmer::store::{CountStore, MemoryCountBuilder, MemoryCountStore};
use anyhow::{Context, Result, bail, ensure};
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use sha2::{Digest, Sha256};
use std::io::BufRead;
use std::io::IsTerminal;
use std::io::Write;
use std::num::NonZeroU64;
use std::path::PathBuf;
use std::sync::Arc;

use crate::cli::Args;
use crate::format::{format_bytes, format_count, format_duration};
use crate::kmer;
use crate::kmer::{Chunk, KmerCounts};

pub(crate) const FASTA_LINE_WIDTH: usize = 80;
pub(crate) const N_READS_PER_BATCH: u64 = 1000;
const KMER_BATCH_CAPACITY: usize = 1024;

mod ena;
pub(crate) use ena::preview_budget as metadata_preview_budget;
pub(crate) use ena::sample_name as resolve_input_sample;
pub(crate) use ena::{EnaResult, get_ena_fastq_urls};

mod resources;
#[cfg(test)]
use resources::{CounterLayout, CounterPlan};
pub(crate) use resources::{PreparedIngestion, prepare_ingestion};

/// Describes how to re-acquire read data for Pass 2 (read threading).
pub(crate) enum ReadSourcePlan {
    /// Local files that can be reopened
    LocalFiles(Vec<PathBuf>),
    /// Remote files cached locally (paths to cached files)
    CachedRemote(Vec<PathBuf>),
    /// Remote files that must be re-downloaded (URLs)
    UncachedRemote(BudgetVec<BudgetString>),
    /// Cannot re-read (stdin or replay otherwise unavailable)
    Unavailable,
}

/// Original-pass replay for exact ingestion, or a legacy source-reopen plan.
pub(crate) struct ReadPlan {
    pub(crate) source: ReadSourcePlan,
    pub(crate) paired: bool,
    pub(crate) max_reads: u64,
    pub(crate) replay: Option<ReplayInput>,
}

/// Replay ownership keeps its directory until the file has been removed.
pub(crate) struct ReplayInput {
    pub(crate) store: crate::replay::ReplayStore,
    directory: tempfile::TempDir,
}

impl ReplayInput {
    /// Close the replay file before removing its containing directory.
    pub(crate) fn close(self) -> Result<()> {
        self.store.close()?;
        self.directory
            .close()
            .context("Unable to remove replay directory")
    }
}

/// Production consumers borrow a replay source; legacy ingestion retains its adapter.
pub(crate) enum ThreadingInput {
    Memory(Vec<ReadRecord>),
    Replay(Box<ReplayInput>),
}

/// A copyable borrowed view shared across target workers.
#[derive(Clone, Copy)]
pub(crate) enum ThreadingSource<'a> {
    Memory(&'a [ReadRecord]),
    Replay(&'a crate::replay::ReplayStore),
}

impl ThreadingSource<'_> {
    pub(crate) fn is_empty(self) -> bool {
        match self {
            Self::Memory(reads) => reads.is_empty(),
            Self::Replay(store) => store.record_count() == 0,
        }
    }
}

impl ThreadingInput {
    pub(crate) fn summary(
        &self,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Option<crate::replay::ReplaySummary>> {
        match self {
            Self::Memory(_) => Ok(None),
            Self::Replay(input) => input.store.summary(budget).map(Some),
        }
    }
    pub(crate) fn source(&self) -> ThreadingSource<'_> {
        match self {
            Self::Memory(reads) => ThreadingSource::Memory(reads),
            Self::Replay(input) => ThreadingSource::Replay(&input.store),
        }
    }
    pub(crate) fn close(self) -> Result<()> {
        if let Self::Replay(input) = self {
            input.close()?;
        }
        Ok(())
    }
}

/// A retained read for Pass 2 threading.
pub(crate) struct ReadRecord {
    /// Raw sequence (ASCII ACGT)
    pub(crate) sequence: String,
    /// Global read index (0-based, in order of ingestion)
    pub(crate) index: u64,
    /// Mate designation for paired-end reads
    pub(crate) mate: Mate,
}

/// Mate designation for paired-end reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mate {
    Unpaired,
    R1,
    R2,
}

/// Warn if an output file already exists (it will be overwritten).
pub(crate) fn warn_if_exists(path: &str) {
    if std::path::Path::new(path).exists() {
        warn!("Overwriting existing file {}", path);
    }
}

/// Write a FASTA record with sequence lines wrapped at FASTA_LINE_WIDTH characters.
pub(crate) fn write_fasta_fields(
    writer: &mut impl Write,
    id: &str,
    description: Option<&str>,
    sequence: &[u8],
) -> std::io::Result<()> {
    write!(writer, ">{}", id)?;
    if let Some(desc) = description {
        write!(writer, " {}", desc)?;
    }
    writeln!(writer)?;
    for chunk in sequence.chunks(FASTA_LINE_WIDTH) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Validate a FASTQ record (4 lines). Returns an error if the record is malformed.
fn validate_fastq_record(
    header: &str,
    _sequence: &str,
    separator: &str,
    quality: &str,
    sequence_len: usize,
    record_num: u64,
) -> Result<()> {
    if header.starts_with('>') {
        bail!(
            "Input appears to be FASTA format, not FASTQ (record {} starts with '>'). \
             sharkmer requires FASTQ input with quality scores.",
            record_num + 1,
        );
    }
    ensure!(
        header.starts_with('@'),
        "FASTQ record {} has invalid header (expected '@', got '{}'): {}",
        record_num + 1,
        header.chars().next().unwrap_or(' '),
        header,
    );
    ensure!(
        separator.starts_with('+'),
        "FASTQ record {} has invalid separator line (expected '+', got '{}'): {}",
        record_num + 1,
        separator.chars().next().unwrap_or(' '),
        separator,
    );
    ensure!(
        quality.len() == sequence_len,
        "FASTQ record {} has mismatched sequence ({}) and quality ({}) lengths",
        record_num + 1,
        sequence_len,
        quality.len(),
    );
    Ok(())
}

/// Mutable state for FASTQ reading, shared across multiple input sources.
pub(crate) struct FastqReadState {
    pub(crate) chunks: Vec<kmer::Chunk>,
    pub(crate) chunk_index: usize,
    pub(crate) seqs: Vec<String>,
    pub(crate) n_reads_read: u64,
    pub(crate) n_bases_read: u64,
    direct: Option<DirectCounter>,
}

enum DirectBuilder {
    Failed,
    Memory(MemoryCountBuilder),
    Adaptive(AdaptiveCountBuilder),
    Packed(PackedCountBuilder),
    External(ExternalCountBuilder),
    Parallel(ParallelCountBuilder),
    AdaptiveParallel(AdaptiveParallelCountBuilder),
}

impl DirectBuilder {
    fn occurrences(&self) -> u64 {
        match self {
            Self::Failed => 0,
            Self::Memory(builder) => builder.occurrences(),
            Self::Adaptive(builder) => builder.occurrences(),
            Self::Packed(builder) => builder.occurrences(),
            Self::External(builder) => builder.occurrences(),
            Self::Parallel(builder) => builder.occurrences(),
            Self::AdaptiveParallel(builder) => builder.occurrences(),
        }
    }

    fn finalize(self, digest: [u8; 32]) -> Result<ExactStore> {
        Ok(match self {
            Self::Failed => anyhow::bail!("Count producer failed during owner transfer"),
            Self::Memory(builder) => ExactStore::Memory(Box::new(builder.finalize(digest))),
            Self::Adaptive(builder) => match builder.finish(digest)? {
                AdaptiveStore::Memory(store) => ExactStore::Memory(Box::new(store)),
                AdaptiveStore::Disk(store, metrics) => ExactStore::Disk {
                    store: Box::new(store),
                    metrics,
                    directory: None,
                },
            },
            Self::Packed(builder) => ExactStore::Packed(Box::new(builder.finalize(digest))),
            Self::External(builder) => {
                let (store, metrics) = builder.finish(digest)?;
                ExactStore::Disk {
                    store: Box::new(store),
                    metrics,
                    directory: None,
                }
            }
            Self::AdaptiveParallel(builder) => {
                let (store, metrics) = builder.finish_adaptive(digest)?;
                ExactStore::AdaptiveSharded {
                    store: Box::new(store),
                    metrics,
                    directory: None,
                }
            }
            Self::Parallel(builder) => {
                let (store, metrics) = builder.finish(digest)?;
                ExactStore::Sharded {
                    store: Box::new(store),
                    metrics,
                }
            }
        })
    }
}

pub(crate) enum ExactStore {
    Memory(Box<MemoryCountStore>),
    Packed(Box<PackedCountStore>),
    Sharded {
        store: Box<ShardedCountStore>,
        metrics: ParallelMetrics,
    },
    AdaptiveSharded {
        store: Box<AdaptiveShardedCountStore>,
        metrics: ParallelMetrics,
        directory: Option<tempfile::TempDir>,
    },
    Disk {
        store: Box<DiskCountStore>,
        metrics: ExternalMetrics,
        directory: Option<tempfile::TempDir>,
    },
}

impl ExactStore {
    fn as_store(&self) -> &dyn CountStore {
        match self {
            Self::Memory(store) => store.as_ref(),
            Self::Packed(store) => store.as_ref(),
            Self::Sharded { store, .. } => store.as_ref(),
            Self::Disk { store, .. } => store.as_ref(),
            Self::AdaptiveSharded { store, .. } => store.as_ref(),
        }
    }
}

/// Actual phase boundary for a completed worker-role transfer.
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub(crate) struct CountOwnerTransferSummary {
    pub(crate) selected_records_before_transfer: u64,
    pub(crate) occurrences_before_transfer: u64,
    pub(crate) current_table_bytes: usize,
    pub(crate) threshold_bytes: usize,
}

struct DirectCounter {
    builder: DirectBuilder,
    budget: Arc<ResourceBudget>,
    record_bytes: usize,
    k: usize,
    valid_bases: u64,
    occurrences: u64,
    selected_digest: SelectedDigest,
    identity_writer: Option<IdentityWriter>,
    completed_input_digest: Option<[u8; 32]>,
    identity_summary: Option<IdentityWorkerSummary>,
    identity_failed: bool,
    owner_transfer_threshold: Option<usize>,
    owner_transfer_summary: Option<CountOwnerTransferSummary>,
    key_batch: Vec<CanonicalKey>,
    next_key: usize,
    // Declared after the buffer so its charge outlives the owned allocation.
    _key_batch_reservation: MemoryReservation,
    replay: Option<crate::replay::ReplayWriter>,
    owned_replay_directory: Option<tempfile::TempDir>,
    // Last: file owners close before the containing directory is removed.
    owned_count_directory: Option<tempfile::TempDir>,
    memory_peak: Option<MemoryPeak>,
    disk_peak: Option<DiskPeak>,
}

impl DirectCounter {
    #[cfg(test)]
    fn resource_budget(
        classes: [u64; crate::kmer::budget::MEMORY_CLASSES],
        disk_bytes: u64,
        graph_bytes: u64,
    ) -> Result<Arc<ResourceBudget>> {
        ResourceBudget::for_counting(
            resources::resource_limits(classes, disk_bytes, graph_bytes)?,
            0,
        )
    }

    #[cfg(test)]
    fn from_plan(plan: CounterPlan, directory: Option<&std::path::Path>) -> Result<Self> {
        let budget = ResourceBudget::for_counting(plan.limits()?, 0)?;
        plan.start(budget, directory)
    }

    #[cfg(test)]
    fn new(
        k: usize,
        table_bytes: u64,
        record_bytes: usize,
        disk_bytes: u64,
        replay_bytes: u64,
        graph_buffer_bytes: u64,
    ) -> Result<Self> {
        Self::with_table_layout(
            k,
            table_bytes,
            record_bytes,
            disk_bytes,
            replay_bytes,
            false,
            graph_buffer_bytes,
        )
    }

    #[cfg(test)]
    fn with_table_layout(
        k: usize,
        table_bytes: u64,
        record_bytes: usize,
        disk_bytes: u64,
        replay_bytes: u64,
        packed: bool,
        graph_buffer_bytes: u64,
    ) -> Result<Self> {
        Self::from_plan(
            CounterPlan {
                k,
                record_bytes,
                disk_bytes,
                replay_bytes,
                graph_bytes: graph_buffer_bytes,
                layout: CounterLayout::Table {
                    table_bytes,
                    packed,
                },
            },
            None,
        )
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn spillable_with_workers(
        k: usize,
        table_bytes: u64,
        sort_bytes: u64,
        record_bytes: usize,
        disk_bytes: u64,
        directory: &std::path::Path,
        replay_bytes: u64,
        workers: usize,
        identity_worker: bool,
        graph_buffer_bytes: u64,
    ) -> Result<Self> {
        Self::spillable_with_policy(
            k,
            table_bytes,
            sort_bytes,
            record_bytes,
            disk_bytes,
            directory,
            replay_bytes,
            workers,
            identity_worker,
            None,
            graph_buffer_bytes,
        )
    }

    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    fn spillable_with_policy(
        k: usize,
        table_bytes: u64,
        sort_bytes: u64,
        record_bytes: usize,
        disk_bytes: u64,
        directory: &std::path::Path,
        replay_bytes: u64,
        workers: usize,
        identity_worker: bool,
        owner_transfer_threshold: Option<usize>,
        graph_buffer_bytes: u64,
    ) -> Result<Self> {
        Self::from_plan(
            CounterPlan {
                k,
                record_bytes,
                disk_bytes,
                replay_bytes,
                graph_bytes: graph_buffer_bytes,
                layout: CounterLayout::Spill {
                    table_bytes,
                    sort_bytes,
                    workers,
                    identity_worker,
                    owner_transfer_threshold,
                },
            },
            Some(directory),
        )
    }

    fn with_builder(
        k: usize,
        record_bytes: usize,
        budget: Arc<ResourceBudget>,
        builder: DirectBuilder,
        memory_peak: MemoryPeak,
        disk_peak: DiskPeak,
    ) -> Result<Self> {
        let key_bytes = u64::try_from(KMER_BATCH_CAPACITY * size_of::<CanonicalKey>())?;
        let mut reservation = budget.reserve(MemoryClass::WorkerQueues, key_bytes)?;
        let mut key_batch = Vec::new();
        key_batch
            .try_reserve_exact(KMER_BATCH_CAPACITY)
            .context("Unable to reserve k-mer batch buffer")?;
        let allocated = key_batch
            .capacity()
            .checked_mul(size_of::<CanonicalKey>())
            .context("K-mer batch byte size overflow")?;
        reservation.resize(u64::try_from(allocated)?)?;
        let mut selected_digest = SelectedDigest::new();
        selected_digest.update(b"sharkmer-selected-fastq-v1\0")?;
        Ok(Self {
            builder,
            budget,
            record_bytes,
            k,
            valid_bases: 0,
            occurrences: 0,
            selected_digest,
            identity_writer: None,
            completed_input_digest: None,
            identity_summary: None,
            identity_failed: false,
            owner_transfer_threshold: None,
            owner_transfer_summary: None,
            key_batch,
            next_key: 0,
            _key_batch_reservation: reservation,
            replay: None,
            owned_replay_directory: None,
            owned_count_directory: None,
            memory_peak: Some(memory_peak),
            disk_peak: Some(disk_peak),
        })
    }

    // Keep encoding and table updates in separate tight loops. The fixed
    // buffer is reused across records; its cursor preserves accepted updates
    // if allocation pressure interrupts a flush before the next key.
    #[inline(never)]
    fn flush_keys(&mut self) -> Result<()> {
        let one = NonZeroU64::new(1).context("Positive count increment is zero")?;
        match &mut self.builder {
            DirectBuilder::Failed => anyhow::bail!("Count producer failed during owner transfer"),
            DirectBuilder::Memory(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key, one)
                })
            }
            DirectBuilder::Adaptive(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key, one)
                })
            }
            DirectBuilder::Packed(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key, one)
                })
            }
            DirectBuilder::External(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key, one)
                })
            }
            DirectBuilder::AdaptiveParallel(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key)
                })
            }
            DirectBuilder::Parallel(builder) => {
                Self::flush_to(&mut self.key_batch, &mut self.next_key, |key| {
                    builder.add(key)
                })
            }
        }
    }

    #[inline]
    fn flush_to(
        keys: &mut Vec<CanonicalKey>,
        next: &mut usize,
        mut add: impl FnMut(CanonicalKey) -> Result<()>,
    ) -> Result<()> {
        while *next < keys.len() {
            add(keys[*next])?;
            *next += 1;
        }
        keys.clear();
        *next = 0;
        Ok(())
    }

    fn input_digest(&mut self) -> Result<[u8; 32]> {
        ensure!(
            !self.identity_failed,
            "Input identity producer previously failed"
        );
        if let Some(digest) = self.completed_input_digest {
            return Ok(digest);
        }
        if let Some(writer) = self.identity_writer.take() {
            let (digest, summary) = match writer.finish() {
                Ok(completed) => completed,
                Err(error) => {
                    self.identity_failed = true;
                    return Err(error);
                }
            };
            self.completed_input_digest = Some(digest);
            self.identity_summary = Some(summary);
            return Ok(digest);
        }
        let digest = self.selected_digest.clone().finalize()?;
        self.completed_input_digest = Some(digest);
        Ok(digest)
    }

    fn maybe_transfer_owner(&mut self) -> Result<()> {
        let Some(threshold) = self.owner_transfer_threshold else {
            return Ok(());
        };
        let bytes = match &self.builder {
            DirectBuilder::Adaptive(builder) => builder.memory_table_bytes(),
            _ => None,
        };
        let Some(bytes) = bytes else {
            self.owner_transfer_threshold = None;
            return Ok(());
        };
        if bytes < threshold {
            return Ok(());
        }
        self.flush_keys()?;
        let writer = self
            .identity_writer
            .take()
            .context("Identity worker missing during owner transfer")?;
        let (state, summary) = writer.finish_state()?;
        self.selected_digest = state;
        self.identity_summary = Some(summary);
        // The failed placeholder makes every later path reject publication if
        // allocation or worker startup fails after moving the accepted prefix.
        let DirectBuilder::Adaptive(prefix) =
            std::mem::replace(&mut self.builder, DirectBuilder::Failed)
        else {
            anyhow::bail!("Owner transfer requires an adaptive memory prefix");
        };
        let owner = AdaptiveParallelCountBuilder::from_prefix(prefix)?;
        ensure!(
            owner.occurrences() == self.occurrences,
            "Transferred owner prefix total differs"
        );
        self.builder = DirectBuilder::AdaptiveParallel(owner);
        self.owner_transfer_summary = Some(CountOwnerTransferSummary {
            selected_records_before_transfer: summary.selected_records,
            occurrences_before_transfer: self.occurrences,
            current_table_bytes: bytes,
            threshold_bytes: threshold,
        });
        self.owner_transfer_threshold = None;
        Ok(())
    }

    fn ingest(&mut self, fields: [&[u8]; 4], source: &str, index: u64, mate: Mate) -> Result<()> {
        ensure!(
            !self.identity_failed && self.completed_input_digest.is_none(),
            "Input identity is already finished or failed"
        );
        ensure!(
            self.key_batch.is_empty(),
            "Previous k-mer batch is unfinished"
        );
        let valid_bases = visit_canonical_windows(fields[1], self.k, |key| {
            self.key_batch.push(key);
            if self.key_batch.len() == KMER_BATCH_CAPACITY {
                self.flush_keys()?;
            }
            Ok(())
        })?;
        self.flush_keys()?;
        // The builder already checks the global total on every accepted add.
        self.occurrences = self.builder.occurrences();
        self.valid_bases = self
            .valid_bases
            .checked_add(u64::try_from(valid_bases)?)
            .context("Input valid-base count overflow")?;
        if let Some(identity) = &mut self.identity_writer {
            if let Err(error) = identity.record(index, mate, source, fields) {
                self.identity_failed = true;
                return Err(error);
            }
        } else {
            self.selected_digest.update(index.to_le_bytes())?;
            self.selected_digest.update([match mate {
                Mate::Unpaired => 0,
                Mate::R1 => 1,
                Mate::R2 => 2,
            }])?;
            // Length framing prevents ambiguous concatenations; delimiters are not
            // part of logical identity, so LF and CRLF describe the same records.
            for value in std::iter::once(source.as_bytes()).chain(fields) {
                self.selected_digest
                    .update(u64::try_from(value.len())?.to_le_bytes())?;
                self.selected_digest.update(value)?;
            }
        }
        if let Some(replay) = &mut self.replay {
            replay.add(crate::replay::ReplayRecord {
                index,
                fragment: if mate == Mate::Unpaired {
                    index
                } else {
                    index / 2
                },
                mate,
                source,
                fields,
            })?;
        }
        if let Err(error) = self.maybe_transfer_owner() {
            self.identity_failed = true;
            return Err(error);
        }
        Ok(())
    }
}

/// Finalized global counts; the legacy path retains its released precision.
pub(crate) enum FinalCounts {
    Legacy(KmerCounts),
    Exact {
        store: ExactStore,
        run_budget: Arc<ResourceBudget>,
        peak_record_buffer_bytes: u64,
        peak_kmer_buffer_bytes: u64,
        peak_count_input_reserved_bytes: u64,
        peak_sort_buffer_bytes: u64,
        peak_temporary_disk_bytes: u64,
        input_identity_worker: Option<IdentityWorkerSummary>,
        count_owner_transfer: Option<CountOwnerTransferSummary>,
        process_file_limits: Option<Box<crate::kmer::budget::ProcessFileLimits>>,
    },
}

impl FinalCounts {
    /// Shared exact-run ledger, omitted by the legacy adapter.
    pub(crate) fn resource_budget(&self) -> Option<&Arc<ResourceBudget>> {
        match self {
            Self::Exact { run_budget, .. } => Some(run_budget),
            Self::Legacy(_) => None,
        }
    }

    /// Borrow an assembly view without copying global counts.
    pub(crate) fn filtered_view(&self, minimum: u32) -> kmer::FilteredKmerCounts<'_> {
        match self {
            Self::Legacy(counts) => counts.filtered_view(minimum),
            Self::Exact {
                store: ExactStore::Memory(store),
                ..
            } => kmer::FilteredKmerCounts::from_memory(store.as_ref(), u64::from(minimum)),
            Self::Exact {
                store: ExactStore::Packed(store),
                ..
            } => kmer::FilteredKmerCounts::from_packed(store.as_ref(), u64::from(minimum)),
            Self::Exact {
                store: ExactStore::Sharded { store, .. },
                ..
            } => kmer::FilteredKmerCounts::from_sharded(store.as_ref(), u64::from(minimum)),
            Self::Exact {
                store: ExactStore::Disk { store, .. },
                ..
            } => kmer::FilteredKmerCounts::from_store(store.as_ref(), u64::from(minimum)),
            Self::Exact {
                store: ExactStore::AdaptiveSharded { store, .. },
                ..
            } => kmer::FilteredKmerCounts::from_store(store.as_ref(), u64::from(minimum)),
        }
    }

    /// Report usable entry capacity, not a distinct-key count.
    pub(crate) fn capacity(&self) -> usize {
        match self {
            Self::Legacy(counts) => counts.capacity(),
            Self::Exact {
                store: ExactStore::Memory(store),
                ..
            } => store.capacity(),
            Self::Exact {
                store: ExactStore::Packed(store),
                ..
            } => store.capacity(),
            Self::Exact {
                store: ExactStore::Sharded { store, .. },
                ..
            } => store.capacity(),
            Self::Exact {
                store: ExactStore::AdaptiveSharded { store, .. },
                ..
            } => store.capacity(),
            Self::Exact {
                store: ExactStore::Disk { .. },
                ..
            } => 0,
        }
    }

    /// Confirm disk cleanup before publishing the final output transaction.
    pub(crate) fn close(self) -> Result<()> {
        let directory = match self {
            Self::Exact {
                store: ExactStore::Disk {
                    store, directory, ..
                },
                ..
            } => {
                store.close()?;
                directory
            }
            Self::Exact {
                store:
                    ExactStore::AdaptiveSharded {
                        store, directory, ..
                    },
                ..
            } => {
                store.close()?;
                directory
            }
            _ => None,
        };
        if let Some(directory) = directory {
            directory
                .close()
                .context("Unable to remove external-count run directory")?;
        }
        Ok(())
    }

    /// Describe actual exact counting when that backend was used.
    pub(crate) fn summary(
        &self,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Option<crate::stats::CountStoreSummary>> {
        let Self::Exact {
            store,
            peak_record_buffer_bytes,
            peak_kmer_buffer_bytes,
            peak_count_input_reserved_bytes,
            peak_sort_buffer_bytes,
            peak_temporary_disk_bytes,
            input_identity_worker,
            count_owner_transfer,
            process_file_limits,
            ..
        } = self
        else {
            return Ok(None);
        };
        let metadata = store.as_store().metadata();
        let (owned_table_bytes, peak_table_bytes, external) = match store {
            ExactStore::Memory(store) => (store.table_bytes(), store.peak_table_bytes(), None),
            ExactStore::Packed(store) => (store.table_bytes(), store.peak_table_bytes(), None),
            ExactStore::Sharded { store, .. } => {
                (store.table_bytes(), store.peak_table_bytes(), None)
            }
            ExactStore::AdaptiveSharded { store, .. } => {
                (store.table_bytes(), store.peak_table_bytes(), None)
            }
            ExactStore::Disk { store, metrics, .. } => (
                0,
                metrics.peak_table_bytes,
                Some(crate::stats::ExternalCountSummary {
                    count_file_bytes: store.file_bytes(),
                    lookup_cache_bytes: store.cache_bytes(),
                    peak_count_page_bytes: store.peak_page_bytes(),
                    peak_sort_buffer_bytes: *peak_sort_buffer_bytes,
                    peak_temporary_disk_bytes: *peak_temporary_disk_bytes,
                    sorted_run_count: metrics.sorted_runs,
                    merge_count: metrics.merges,
                }),
            ),
        };
        Ok(Some(crate::stats::CountStoreSummary {
            process_file_limits: process_file_limits.as_deref().copied(),
            format_version: metadata.format_version,
            backend: metadata.backend,
            hash_algorithm: metadata.hash_algorithm,
            count_width_bits: metadata.count_width_bits,
            overflow_policy: metadata.overflow_policy,
            preprocessing: metadata.preprocessing,
            canonicalization: metadata.canonicalization,
            selected_input_sha256: crate::kmer::budget::text::BudgetString::hex(
                &metadata.selected_input_sha256,
                budget,
                MemoryClass::RetainedResults,
            )?,
            distinct_keys: metadata.distinct_keys,
            occurrences: metadata.occurrences,
            owned_table_bytes,
            peak_table_bytes,
            peak_record_buffer_bytes: *peak_record_buffer_bytes,
            peak_kmer_buffer_bytes: *peak_kmer_buffer_bytes,
            peak_count_input_reserved_bytes: *peak_count_input_reserved_bytes,
            external,
            input_identity_worker: *input_identity_worker,
            count_owner_transfer: *count_owner_transfer,
            adaptive_parallel: match store {
                ExactStore::AdaptiveSharded { store, .. } => {
                    Some(crate::stats::AdaptiveParallelSummary {
                        shards: store.metrics(),
                        peak_sort_buffer_bytes: *peak_sort_buffer_bytes,
                        peak_temporary_disk_bytes: *peak_temporary_disk_bytes,
                    })
                }
                _ => None,
            },
            packed: match store {
                ExactStore::Packed(store) => Some(store.metrics()),
                _ => None,
            },
            parallel: match store {
                ExactStore::Sharded { metrics, .. }
                | ExactStore::AdaptiveSharded { metrics, .. } => Some(*metrics),
                _ => None,
            },
        }))
    }

    /// Order-independent fingerprint of all full key/count pairs. This opt-in
    /// verification scan is excluded from ordinary counting/timing runs.
    pub(crate) fn fingerprint(
        &self,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<crate::stats::CountFingerprint> {
        let view = self.filtered_view(0);
        let mut xor = [0u8; 32];
        let mut distinct = 0u64;
        let mut occurrences = 0u64;
        view.visit(&mut |bits, count| {
            let mut hash = Sha256::new();
            hash.update(b"sharkmer-count-pair-v1\0");
            hash.update(u64::try_from(view.get_k())?.to_le_bytes());
            hash.update(bits.to_le_bytes());
            hash.update(count.to_le_bytes());
            for (destination, byte) in xor.iter_mut().zip(hash.finalize()) {
                *destination ^= byte;
            }
            distinct = distinct
                .checked_add(1)
                .context("Fingerprint cardinality overflow")?;
            occurrences = occurrences
                .checked_add(count)
                .context("Fingerprint occurrence overflow")?;
            Ok(())
        })?;
        Ok(crate::stats::CountFingerprint {
            algorithm: "sha256_pair_xor_v1",
            distinct_keys: distinct,
            occurrences,
            pair_digest: crate::kmer::budget::text::BudgetString::hex(
                &xor,
                budget,
                MemoryClass::RetainedResults,
            )?,
        })
    }
}

/// Build a context-rich error message for an I/O error that occurred mid-stream
/// during FASTQ reading. Distinguishes abrupt connection drops (the common ENA
/// failure mode) from actual FASTQ parse issues, and suggests the cache as a
/// workaround when streaming fails.
pub(crate) fn stream_io_error(
    line_role: &str,
    err: &std::io::Error,
    source_name: &str,
    n_reads_read: u64,
) -> anyhow::Error {
    use std::io::ErrorKind;
    let is_stream_drop = matches!(
        err.kind(),
        ErrorKind::UnexpectedEof
            | ErrorKind::ConnectionReset
            | ErrorKind::ConnectionAborted
            | ErrorKind::BrokenPipe
            | ErrorKind::TimedOut
            | ErrorKind::Interrupted
    );
    let is_remote = source_name.starts_with("http://") || source_name.starts_with("https://");

    if is_stream_drop {
        if is_remote {
            anyhow::anyhow!(
                "Stream from {} dropped while reading {} line of record {} (I/O error: {} — kind {:?}).\n\
                 This is usually a transient network interruption, not a bad FASTQ file.\n\
                 Retry the run. If it repeats, use the read cache (the default for --ena); \
                 cached downloads are verified by SHA-256 and do not suffer mid-stream drops.",
                source_name,
                line_role,
                n_reads_read + 1,
                err,
                err.kind(),
            )
        } else {
            anyhow::anyhow!(
                "Local read stream ended unexpectedly while reading {} line of record {} in {} \
                 (I/O error: {} — kind {:?}). The file may be truncated or corrupted.",
                line_role,
                n_reads_read + 1,
                source_name,
                err,
                err.kind(),
            )
        }
    } else {
        anyhow::anyhow!(
            "Failed to read {} line of record {} in {}: {} (kind {:?})",
            line_role,
            n_reads_read + 1,
            source_name,
            err,
            err.kind(),
        )
    }
}

/// Read FASTQ records from a buffered reader, ingesting sequences into chunks.
///
/// Reads 4 lines at a time (header, sequence, separator, quality) and validates
/// the format. Returns true if max_reads was reached.
fn read_fastq<R: BufRead>(
    reader: R,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source_name: &str,
    progress: &ProgressBar,
) -> Result<bool> {
    if state.direct.is_some() {
        return read_fastq_direct(
            reader,
            state,
            max_reads,
            validate_every,
            source_name,
            progress,
        );
    }
    let mut lines = reader.lines();
    let n_chunks = state.chunks.len();

    while let Some(header_result) = lines.next() {
        // Read 4 lines for one FASTQ record
        let header = header_result
            .map_err(|e| stream_io_error("header", &e, source_name, state.n_reads_read))?;

        let sequence = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("sequence", &e, source_name, state.n_reads_read))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing sequence line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let separator = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("separator", &e, source_name, state.n_reads_read))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing separator line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let quality = match lines.next() {
            Some(line) => {
                line.map_err(|e| stream_io_error("quality", &e, source_name, state.n_reads_read))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing quality line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        // Validate the record
        let should_validate = state.n_reads_read == 0
            || (validate_every > 0 && state.n_reads_read % validate_every == 0);
        if should_validate {
            validate_fastq_record(
                &header,
                &sequence,
                &separator,
                &quality,
                sequence.len(),
                state.n_reads_read,
            )?;
        }

        // Process the sequence
        state.n_bases_read += sequence.len() as u64;
        state.seqs.push(sequence);
        state.n_reads_read += 1;

        // If we have read enough reads, ingest them into current chunk
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            drain_batch(state, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

fn direct_reader<R: BufRead>(
    reader: R,
    state: &FastqReadState,
) -> Result<crate::fastq::FastqReader<R>> {
    let counter = state.direct.as_ref().context("Missing direct counter")?;
    crate::fastq::FastqReader::new(reader, counter.record_bytes, Arc::clone(&counter.budget))
}

fn ingest_direct_record(
    fields: [&[u8]; 4],
    state: &mut FastqReadState,
    validate_every: u64,
    source: &str,
    mate: Mate,
) -> Result<()> {
    let should_validate =
        state.n_reads_read == 0 || (validate_every > 0 && state.n_reads_read % validate_every == 0);
    if should_validate {
        validate_fastq_record(
            std::str::from_utf8(fields[0])?,
            std::str::from_utf8(fields[1])?,
            std::str::from_utf8(fields[2])?,
            std::str::from_utf8(fields[3])?,
            fields[1].len(),
            state.n_reads_read,
        )?;
    }
    let records = state
        .n_reads_read
        .checked_add(1)
        .context("FASTQ record count overflow")?;
    let bases = state
        .n_bases_read
        .checked_add(u64::try_from(fields[1].len())?)
        .context("FASTQ base count overflow")?;
    state
        .direct
        .as_mut()
        .context("Missing direct counter")?
        .ingest(fields, source, state.n_reads_read, mate)?;
    state.n_reads_read = records;
    state.n_bases_read = bases;
    Ok(())
}

fn read_fastq_direct<R: BufRead>(
    reader: R,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source: &str,
    progress: &ProgressBar,
) -> Result<bool> {
    let mut reader = direct_reader(reader, state)?;
    while reader.next(source, state.n_reads_read)? {
        ingest_direct_record(
            reader.fields(),
            state,
            validate_every,
            source,
            Mate::Unpaired,
        )?;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments)]
fn read_fastq_paired_direct<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source1: &str,
    source2: &str,
    progress: &ProgressBar,
) -> Result<bool> {
    let mut reader1 = direct_reader(reader1, state)?;
    let mut reader2 = direct_reader(reader2, state)?;
    let mut counts = [0u64; 2];
    loop {
        if !reader1.next(source1, state.n_reads_read)? {
            // Preserve the released selection at uneven EOF: the first extra
            // record from the longer file is included before the warning.
            if reader2.next(source2, state.n_reads_read)? {
                ingest_direct_record(reader2.fields(), state, validate_every, source2, Mate::R2)?;
                counts[1] += 1;
            }
            break;
        }
        ingest_direct_record(reader1.fields(), state, validate_every, source1, Mate::R1)?;
        counts[0] += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            return Ok(true);
        }
        if !reader2.next(source2, state.n_reads_read)? {
            break;
        }
        ingest_direct_record(reader2.fields(), state, validate_every, source2, Mate::R2)?;
        counts[1] += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            return Ok(true);
        }
    }
    progress.set_position(state.n_reads_read);
    if counts[0] != counts[1] {
        warn!(
            "Paired-end input length mismatch: {} has {} reads, {} has {} reads. The first unmatched read was included; later unmatched reads were not selected.",
            source1, counts[0], source2, counts[1]
        );
    }
    Ok(false)
}

/// Drain accumulated sequences into the current chunk.
fn drain_batch(state: &mut FastqReadState, n_chunks: usize) -> Result<()> {
    if state.direct.is_some() {
        return Ok(());
    }
    for seq in state.seqs.drain(..) {
        state.chunks[state.chunk_index].ingest_seq(&seq)?;
    }
    state.chunk_index = (state.chunk_index + 1) % n_chunks;
    Ok(())
}

/// Ingest FASTQ reads from all input sources (ENA, files, or stdin).
/// Returns the read state with populated chunks, summary statistics, and a
/// `ReadPlan` describing how to re-read the same data in Pass 2.
pub(crate) fn ingest_reads(
    prepared: PreparedIngestion<'_>,
    mut cached_ena_result: Option<EnaResult>,
    cache_config: Option<&crate::cache::CacheConfig>,
    show_progress: bool,
) -> Result<(FastqReadState, u64, u64, u64, ReadPlan)> {
    let start = std::time::Instant::now();
    info!("Ingesting reads...");
    let (args, direct, replay_requested, workspace_root) = prepared.start()?;
    let k = args.k;
    // Legacy chunks retain their separate analysis path.
    let n_chunks = if args.count_backend() != "legacy" {
        0
    } else if args.chunks == 0 {
        1
    } else {
        args.chunks
    };
    let chunks = (0..n_chunks).map(|_| Chunk::new(&k)).collect();
    let mut state = FastqReadState {
        chunks,
        chunk_index: 0,
        seqs: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
        direct,
    };

    let mut max_reads = args.max_reads.unwrap_or(0);
    if replay_requested && let Some(counter) = &mut state.direct {
        let fallback;
        let parent = match workspace_root.as_deref().or(args.count_temp_dir.as_deref()) {
            Some(parent) => parent,
            None => {
                fallback = std::env::temp_dir();
                &fallback
            }
        };
        let directory = tempfile::Builder::new()
            .prefix(".sharkmer-replay-run-")
            .tempdir_in(parent)
            .context("Unable to create replay directory")?;
        let cap = if args.paired && max_reads % 2 != 0 {
            max_reads
                .checked_add(1)
                .context("Paired read cap overflow")?
        } else {
            max_reads
        };
        let limit = args
            .fastq_record_bytes
            .checked_add(65536 + 96)
            .context("Replay frame limit overflow")?;
        counter.replay = Some(crate::replay::ReplayWriter::new(
            directory.path(),
            &counter.budget,
            limit,
            args.paired,
            cap,
        )?);
        counter.owned_replay_directory = Some(directory);
    }

    // Create progress indicator: bar with ETA if max_reads is known, spinner otherwise.
    let progress = if show_progress && max_reads > 0 {
        let pb = ProgressBar::new(max_reads);
        pb.set_style(
            ProgressStyle::with_template(
                "Ingesting reads {bar:30} {human_pos}/{human_len} [{per_sec}] [ETA {eta}]",
            )
            .expect("valid progress template"),
        );
        pb
    } else if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("Ingesting reads {spinner} {human_pos} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    // Track the input source for Pass 2 re-reading
    let mut read_source_plan = ReadSourcePlan::Unavailable;

    if let Some(accession) = &args.ena {
        // Stream reads from ENA (use cached result if available)
        let ena_result = match cached_ena_result.take() {
            Some(r) => r,
            None => get_ena_fastq_urls(
                accession,
                state.direct.as_ref().map(|counter| &counter.budget),
            )?,
        };
        let mut cached_paths: Vec<PathBuf> = Vec::new();
        let mut is_cached = false;
        for url in &ena_result.urls {
            let reader: Box<dyn BufRead> = if let Some(cache) = cache_config {
                // Use cache: lookup or download
                let cached = if let Some(counter) = &state.direct {
                    cache.lookup_bounded(url, max_reads, &counter.budget)?
                } else {
                    cache.lookup(url, max_reads)?
                };
                let local_path = match cached {
                    Some(path) => {
                        info!("Cache hit for {}", url);
                        path
                    }
                    None => {
                        info!("Cache miss for {}, downloading...", url);
                        if let Some(counter) = &state.direct {
                            cache.download_to_cache_bounded(
                                url,
                                max_reads,
                                counter.record_bytes,
                                &counter.budget,
                            )?
                        } else {
                            cache.download_to_cache(url, max_reads)?
                        }
                    }
                };
                cached_paths.push(local_path.clone());
                is_cached = true;
                let file = std::fs::File::open(&local_path).with_context(|| {
                    format!("Failed to open cached file: {}", local_path.display())
                })?;
                Box::new(crate::codec::Decoder::new(
                    file,
                    state.direct.as_ref().map(|counter| &counter.budget),
                )?)
            } else {
                // No cache: stream directly
                info!("Streaming from {} (no cache)...", url);
                let response =
                    crate::http::get(url).with_context(|| format!("Failed to download {}", url))?;
                Box::new(crate::codec::Decoder::new(
                    response,
                    state.direct.as_ref().map(|counter| &counter.budget),
                )?)
            };
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                url,
                &progress,
            )?;
            if reached_max {
                break;
            }
        }
        read_source_plan = if is_cached {
            ReadSourcePlan::CachedRemote(cached_paths)
        } else {
            warn!("Read threading will require re-downloading reads from ENA (no cache)");
            ReadSourcePlan::UncachedRemote(ena_result.urls)
        };
    } else if let Some(input_files) = &args.input {
        if args.paired {
            // Paired-end: alternate reads from R1 and R2
            let budget = state.direct.as_ref().map(|counter| &counter.budget);
            let reader1 = open_fastq_reader_bounded(&input_files[0], budget)?;
            let reader2 = open_fastq_reader_bounded(&input_files[1], budget)?;
            let name1 = input_files[0].to_string_lossy().to_string();
            let name2 = input_files[1].to_string_lossy().to_string();
            // Round max_reads up to even so Pass 1 consumes balanced
            // pairs. The rounded value is also recorded in the ReadPlan
            // below so Pass 2 (re-read for threading) consumes the same
            // number of reads; leaving the unrounded value here would
            // cause Pass 2 to stop one R2 short of Pass 1 and drop the
            // final read pair silently.
            if max_reads > 0 && max_reads % 2 != 0 {
                max_reads += 1;
            }
            read_fastq_paired(
                reader1,
                reader2,
                &mut state,
                max_reads,
                args.validate_every,
                &name1,
                &name2,
                &progress,
            )?;
        } else {
            // Read from one or more files sequentially
            for file_path in input_files.iter() {
                let reader = open_fastq_reader_bounded(
                    file_path,
                    state.direct.as_ref().map(|counter| &counter.budget),
                )?;
                let file_name = file_path.to_string_lossy();
                let reached_max = read_fastq(
                    reader,
                    &mut state,
                    max_reads,
                    args.validate_every,
                    &file_name,
                    &progress,
                )?;
                if reached_max {
                    break;
                }
            }
        }
        read_source_plan = ReadSourcePlan::LocalFiles(input_files.clone());
    } else {
        // Read from stdin (cannot re-read for Pass 2)
        let stdin = std::io::stdin();
        ensure!(
            !stdin.is_terminal(),
            "No input files specified and stdin is a terminal.\n\
             Provide FASTQ files as arguments, use --ena, or pipe data via stdin.\n\
             Example: sharkmer -s sample -k 21 reads.fastq\n\
             Example: sharkmer -s sample --ena SRR5324768\n\
             Example: zcat reads.fastq.gz | sharkmer -s sample -k 21"
        );
        let handle = stdin.lock();

        read_fastq(
            handle,
            &mut state,
            max_reads,
            args.validate_every,
            "stdin",
            &progress,
        )?;
        // read_source_plan remains Unavailable for stdin
    }

    progress.finish_and_clear();

    // Ingest any remaining sequences
    let n_chunks = state.chunks.len();
    drain_batch(&mut state, n_chunks)?;

    let mut n_reads_ingested: u64 = 0;
    let mut n_bases_ingested: u64 = 0;
    let mut n_kmers_ingested: u64 = 0;
    for chunk in state.chunks.iter() {
        n_reads_ingested += chunk.get_n_reads();
        n_bases_ingested += chunk.get_n_bases();
        n_kmers_ingested += chunk.get_n_kmers();
    }
    if let Some(counter) = &state.direct {
        n_reads_ingested = state.n_reads_read;
        n_bases_ingested = counter.valid_bases;
        n_kmers_ingested = counter.occurrences;
    }

    info!(
        "Read {} reads, {} bases",
        format_count(state.n_reads_read),
        format_bytes(state.n_bases_read)
    );
    info!(
        "Ingested {} subreads, {} bases, {} kmers",
        format_count(n_reads_ingested),
        format_bytes(n_bases_ingested),
        format_count(n_kmers_ingested)
    );

    // Warn if very few reads for sPCR
    let has_pcr = !args.pcr_panel.is_empty()
        || !args.pcr_panel_file.is_empty()
        || !args.pcr_primers.is_empty();
    if has_pcr && state.n_reads_read < 10_000 {
        warn!(
            "Only {} FASTQ records read. sPCR typically needs many more reads to produce results.",
            state.n_reads_read
        );
    }
    info!("Time to ingest reads: {}", format_duration(start.elapsed()));

    if n_reads_ingested == 0 {
        bail!("No reads were ingested. Check that input files contain valid FASTQ records.");
    }

    let replay = if let Some(counter) = &mut state.direct {
        if let Some(writer) = counter.replay.take() {
            let digest = counter.input_digest()?;
            Some(ReplayInput {
                store: writer.finish(digest)?,
                directory: counter
                    .owned_replay_directory
                    .take()
                    .context("Missing replay directory")?,
            })
        } else {
            None
        }
    } else {
        None
    };
    let read_plan = ReadPlan {
        replay,
        source: read_source_plan,
        paired: args.paired,
        max_reads,
    };

    Ok((
        state,
        n_reads_ingested,
        n_bases_ingested,
        n_kmers_ingested,
        read_plan,
    ))
}

/// Open a FASTQ file as a buffered reader, with automatic gzip detection.
fn open_fastq_reader(file_path: &std::path::Path) -> Result<Box<dyn BufRead>> {
    open_fastq_reader_bounded(file_path, None)
}

fn open_fastq_reader_bounded(
    file_path: &std::path::Path,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Box<dyn BufRead>> {
    let file_name = file_path.to_string_lossy();
    let has_gz_ext = file_name.ends_with(".gz") || file_name.ends_with(".gzip");

    // Open the file once. If we need to peek at the gzip magic bytes, do so
    // through the BufReader: fill_buf() loads bytes into the internal buffer
    // without consuming them, so subsequent reads (including through a
    // MultiGzDecoder wrapper) see the stream from position 0. This avoids a
    // redundant second open that showed up as two syscalls per FASTQ file.
    let file = std::fs::File::open(file_path)
        .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
    let mut buf_reader = crate::codec::Buffered::new(file, budget)?;

    let use_gzip = if has_gz_ext {
        true
    } else {
        let magic = buf_reader.fill_buf().context("Failed to peek at file")?;
        magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b
    };

    if use_gzip {
        Ok(Box::new(crate::codec::Decoder::new(buf_reader, budget)?))
    } else {
        Ok(Box::new(buf_reader))
    }
}

/// Read FASTQ sequences from two files in alternating order (R1, R2, R1, R2, ...).
/// Returns true if max_reads was reached.
#[allow(clippy::too_many_arguments)]
fn read_fastq_paired<R1: BufRead, R2: BufRead>(
    reader1: R1,
    reader2: R2,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source1_name: &str,
    source2_name: &str,
    progress: &ProgressBar,
) -> Result<bool> {
    if state.direct.is_some() {
        return read_fastq_paired_direct(
            reader1,
            reader2,
            state,
            max_reads,
            validate_every,
            source1_name,
            source2_name,
            progress,
        );
    }
    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();
    let n_chunks = state.chunks.len();
    let mut r1_records: u64 = 0;
    let mut r2_records: u64 = 0;
    let mismatch_exit: bool;

    loop {
        // Read one record from R1
        let r1_done = read_one_fastq_record(&mut lines1, state, validate_every, source1_name)?;
        if r1_done {
            // R1 exhausted. Verify R2 is also at EOF so we don't silently drop
            // unmatched R2 records (truncated R1, mismatched inputs).
            let r2_extra = read_one_fastq_record(&mut lines2, state, validate_every, source2_name)?;
            mismatch_exit = !r2_extra;
            if !r2_extra {
                r2_records += 1;
            }
            break;
        }
        r1_records += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            drain_batch(state, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }

        // Read one record from R2
        let r2_done = read_one_fastq_record(&mut lines2, state, validate_every, source2_name)?;
        if r2_done {
            // R2 ended mid-pair: R1 already consumed one more record than R2.
            mismatch_exit = true;
            break;
        }
        r2_records += 1;
        if state.n_reads_read % N_READS_PER_BATCH == 0 {
            drain_batch(state, n_chunks)?;
            progress.set_position(state.n_reads_read);
        }
        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true);
        }
    }

    if mismatch_exit && r1_records != r2_records {
        warn!(
            "Paired-end input length mismatch: {} has {} reads, {} has {} reads. \
             Extra reads in the longer file were not processed.",
            source1_name, r1_records, source2_name, r2_records
        );
    }

    Ok(false)
}

/// Read a single FASTQ record (4 lines) from a lines iterator.
/// Returns true if EOF was reached (no more records).
fn read_one_fastq_record<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    state: &mut FastqReadState,
    validate_every: u64,
    source_name: &str,
) -> Result<bool> {
    let header = match lines.next() {
        Some(result) => {
            result.map_err(|e| stream_io_error("header", &e, source_name, state.n_reads_read))?
        }
        None => return Ok(true), // EOF
    };

    let sequence = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("sequence", &e, source_name, state.n_reads_read))?
        }
        None => bail!(
            "Truncated FASTQ record at record {} in {}: missing sequence line",
            state.n_reads_read + 1,
            source_name,
        ),
    };

    let separator = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("separator", &e, source_name, state.n_reads_read))?
        }
        None => bail!(
            "Truncated FASTQ record at record {} in {}: missing separator line",
            state.n_reads_read + 1,
            source_name,
        ),
    };

    let quality = match lines.next() {
        Some(line) => {
            line.map_err(|e| stream_io_error("quality", &e, source_name, state.n_reads_read))?
        }
        None => bail!(
            "Truncated FASTQ record at record {} in {}: missing quality line",
            state.n_reads_read + 1,
            source_name,
        ),
    };

    let should_validate =
        state.n_reads_read == 0 || (validate_every > 0 && state.n_reads_read % validate_every == 0);
    if should_validate {
        validate_fastq_record(
            &header,
            &sequence,
            &separator,
            &quality,
            sequence.len(),
            state.n_reads_read,
        )?;
    }

    state.n_bases_read += sequence.len() as u64;
    state.seqs.push(sequence);
    state.n_reads_read += 1;

    Ok(false) // not EOF
}

/// Re-read FASTQ sequences for Pass 2 (read threading).
/// Opens the same sources used in Pass 1 and collects sequences into `ReadRecord`s.
pub(crate) fn reread_sequences(plan: &ReadPlan, show_progress: bool) -> Result<Vec<ReadRecord>> {
    let start = std::time::Instant::now();
    info!("Pass 2: re-reading sequences for read threading...");

    // Temp files held here are auto-deleted when this vec drops at end of
    // function. Each NamedTempFile has a unique path (random suffix), so
    // concurrent sharkmer invocations and same-length URLs cannot collide.
    let mut _pass2_tempfiles: Vec<tempfile::NamedTempFile> = Vec::new();
    let files: Vec<PathBuf> = match &plan.source {
        ReadSourcePlan::LocalFiles(paths) => paths.clone(),
        ReadSourcePlan::CachedRemote(paths) => paths.clone(),
        ReadSourcePlan::UncachedRemote(urls) => {
            // Re-download to temporary cache paths
            warn!("Re-downloading reads for Pass 2 (use --cache-dir to avoid this)");
            let mut paths = Vec::new();
            for url in urls {
                info!("Downloading {} for Pass 2...", url);
                let response = ureq::get(url)
                    .call()
                    .with_context(|| format!("Failed to download {} for Pass 2", url))?;
                // Create a unique temp file (random suffix) in the system temp
                // directory. The NamedTempFile guard is kept alive in
                // `_pass2_tempfiles` below so the file persists through Pass 2
                // reading and is deleted on function return.
                let mut tmp = tempfile::Builder::new()
                    .prefix("sharkmer_pass2_")
                    .suffix(".fastq.gz")
                    .tempfile()
                    .context("Failed to create temp file for Pass 2")?;
                std::io::copy(&mut response.into_reader(), tmp.as_file_mut())
                    .context("Failed to write Pass 2 temp file")?;
                paths.push(tmp.path().to_path_buf());
                _pass2_tempfiles.push(tmp);
            }
            paths
        }
        ReadSourcePlan::Unavailable => {
            return Ok(Vec::new());
        }
    };

    let progress = if show_progress && plan.max_reads > 0 {
        let pb = ProgressBar::new(plan.max_reads);
        pb.set_style(
            ProgressStyle::with_template("Pass 2 {bar:30} {human_pos}/{human_len} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("Pass 2 {spinner} {human_pos} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    let mut records: Vec<ReadRecord> = Vec::new();
    let mut index: u64 = 0;

    if plan.paired && files.len() == 2 {
        // Paired-end: alternate R1 and R2. Track per-side record counts
        // so a truncated mate file can be surfaced to the user rather
        // than silently producing an unbalanced read set for Pass 2
        // threading. This mirrors the check in read_fastq_paired for
        // Pass 1.
        let reader1 = open_fastq_reader(&files[0])?;
        let reader2 = open_fastq_reader(&files[1])?;
        let name1 = files[0].to_string_lossy().to_string();
        let name2 = files[1].to_string_lossy().to_string();
        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();
        let mut r1_records: u64 = 0;
        let mut r2_records: u64 = 0;
        let mismatch_exit: bool;

        'outer: loop {
            // Read R1
            match read_sequence_only(&mut lines1)? {
                Some(seq) => {
                    records.push(ReadRecord {
                        sequence: seq,
                        index,
                        mate: Mate::R1,
                    });
                    index += 1;
                    r1_records += 1;
                }
                None => {
                    // R1 exhausted. Check whether R2 also reached EOF —
                    // if not, R2 has extra unmatched records.
                    let r2_extra = read_sequence_only(&mut lines2)?.is_some();
                    mismatch_exit = r2_extra;
                    if r2_extra {
                        r2_records += 1; // counted so the warning numbers are meaningful
                    }
                    break 'outer;
                }
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                mismatch_exit = false;
                break;
            }

            // Read R2
            match read_sequence_only(&mut lines2)? {
                Some(seq) => {
                    records.push(ReadRecord {
                        sequence: seq,
                        index,
                        mate: Mate::R2,
                    });
                    index += 1;
                    r2_records += 1;
                }
                None => {
                    // R2 ended mid-pair: R1 already consumed one more.
                    mismatch_exit = true;
                    break 'outer;
                }
            }
            progress.set_position(index);
            if plan.max_reads > 0 && index >= plan.max_reads {
                mismatch_exit = false;
                break;
            }
        }

        if mismatch_exit && r1_records != r2_records {
            warn!(
                "Pass 2 paired-end input length mismatch: {} has {} reads, {} has {} reads. \
                 Extra reads in the longer file were not re-read for threading.",
                name1, r1_records, name2, r2_records
            );
        }
    } else {
        // Unpaired: read files sequentially
        for file_path in &files {
            let reader = open_fastq_reader(file_path)?;
            let mut lines = reader.lines();

            while let Some(seq) = read_sequence_only(&mut lines)? {
                records.push(ReadRecord {
                    sequence: seq,
                    index,
                    mate: Mate::Unpaired,
                });
                index += 1;
                if index % 10_000 == 0 {
                    progress.set_position(index);
                }
                if plan.max_reads > 0 && index >= plan.max_reads {
                    break;
                }
            }
            if plan.max_reads > 0 && index >= plan.max_reads {
                break;
            }
        }
    }

    progress.finish_and_clear();
    info!(
        "Pass 2: collected {} reads for threading in {}",
        format_count(index),
        format_duration(start.elapsed())
    );

    Ok(records)
}

/// Read a single FASTQ record and return only the sequence line.
/// Returns None at EOF.
fn read_sequence_only<R: BufRead>(lines: &mut std::io::Lines<R>) -> Result<Option<String>> {
    // Header
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ header")?;
        }
        None => return Ok(None),
    }
    // Sequence
    let sequence = match lines.next() {
        Some(result) => result.context("Failed to read FASTQ sequence")?,
        None => bail!("Truncated FASTQ record: missing sequence line"),
    };
    // Separator
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ separator")?;
        }
        None => bail!("Truncated FASTQ record: missing separator line"),
    }
    // Quality
    match lines.next() {
        Some(result) => {
            result.context("Failed to read FASTQ quality")?;
        }
        None => bail!("Truncated FASTQ record: missing quality line"),
    }
    Ok(Some(sequence))
}

/// Consolidate chunks into a single KmerCounts table, optionally writing histograms.
/// Returns (kmer_counts, n_singleton_kmers).
pub(crate) fn consolidate_and_histogram(
    state: &mut FastqReadState,
    args: &Args,
    k: usize,
    sample: &str,
    directory: &str,
    n_kmers_ingested: u64,
    show_progress: bool,
) -> Result<(FinalCounts, Option<u64>)> {
    if let Some(mut counter) = state.direct.take() {
        ensure!(
            args.chunks == 0 && state.chunks.is_empty(),
            "Exact counter cannot finalize legacy chunks"
        );
        ensure!(
            counter.key_batch.is_empty() && counter.next_key == 0,
            "Exact counter has an unfinished k-mer batch"
        );
        let budget = Arc::clone(&counter.budget);
        let key_buffer_bytes =
            u64::try_from(counter.key_batch.capacity() * size_of::<CanonicalKey>())?;
        // Input arenas are gone; the key buffer is still charged until this
        // counter is consumed and its table moves into immutable storage.
        let digest = counter.input_digest()?;
        let input_identity_worker = counter.identity_summary;
        let count_owner_transfer = counter.owner_transfer_summary;
        let mut store = counter.builder.finalize(digest)?;
        if let ExactStore::Disk { directory, .. } | ExactStore::AdaptiveSharded { directory, .. } =
            &mut store
        {
            *directory = counter.owned_count_directory.take();
        } else if let Some(directory) = counter.owned_count_directory.take() {
            directory
                .close()
                .context("Unable to remove unused count run directory")?;
        }
        let peak_count_input_reserved_bytes = counter
            .memory_peak
            .take()
            .context("Exact counter has no memory measurement")?
            .finish()?;
        let snapshot = budget.snapshot();
        let peak_temporary_disk_bytes = counter
            .disk_peak
            .take()
            .context("Exact counter has no disk measurement")?
            .finish()?;
        let peak_record_buffer_bytes =
            snapshot.peak_reserved_by_class[MemoryClass::InputBuffer as usize];
        ensure!(
            store.as_store().metadata().occurrences == n_kmers_ingested,
            "Exact count total differs from selected windows"
        );
        info!(
            "Finalized {} exact unique kmers with a total count of {}",
            store.as_store().metadata().distinct_keys,
            n_kmers_ingested
        );
        return Ok((
            FinalCounts::Exact {
                store,
                run_budget: Arc::clone(&budget),
                peak_record_buffer_bytes,
                peak_kmer_buffer_bytes: key_buffer_bytes,
                peak_count_input_reserved_bytes,
                peak_sort_buffer_bytes: snapshot.peak_reserved_by_class
                    [MemoryClass::SortBuffer as usize],
                peak_temporary_disk_bytes,
                input_identity_worker,
                count_owner_transfer,
                process_file_limits: budget.process_file_limits().map(Box::new),
            },
            None,
        ));
    }
    let spinner_msg = if args.chunks > 0 {
        "Consolidating kmer counts..."
    } else {
        "Merging kmer counts..."
    };
    info!("{}", spinner_msg);
    let spinner_style =
        ProgressStyle::with_template("{spinner:.cyan} {msg}").expect("valid spinner template");
    let consolidate_spinner = if show_progress {
        let sp = ProgressBar::new_spinner();
        sp.set_style(spinner_style);
        sp.set_message(spinner_msg.to_string());
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        sp
    } else {
        ProgressBar::hidden()
    };
    let start = std::time::Instant::now();

    let estimated_capacity: usize = state.chunks.iter().map(|c| c.get_kmer_counts().len()).sum();
    let mut kmer_counts: KmerCounts = KmerCounts::new_with_capacity(&k, estimated_capacity);
    let mut n_singleton_kmers: Option<u64> = None;

    let histo_comment = format!(
        "# sharkmer {} k={} chunks={}",
        env!("CARGO_PKG_VERSION"),
        args.k,
        args.chunks
    );

    if args.chunks > 0 {
        // Incremental histogram mode: update histogram as chunks are merged.
        // Snapshot only the Vec<u64> at each step (not the full Histogram struct
        // with its `histo_large` AHashMap tail).
        let mut histo_vecs: Vec<Vec<u64>> = Vec::with_capacity(args.chunks);
        let mut running_histo = kmer::Histogram::new(&args.histo_max);

        for chunk in state.chunks.drain(..) {
            kmer_counts.extend_with_histogram(chunk.get_kmer_counts(), &mut running_histo)?;
            drop(chunk);

            histo_vecs.push(running_histo.get_vector()?);
        }
        consolidate_spinner.finish_and_clear();
        info!(
            "Chunks consolidated, time: {}",
            format_duration(start.elapsed())
        );

        let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
        info!(
            "{} unique kmers with a total count of {} were found",
            kmer_counts.get_n_unique_kmers(),
            n_hashed_kmers
        );

        ensure!(
            n_hashed_kmers == n_kmers_ingested,
            "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
            n_hashed_kmers,
            n_kmers_ingested,
        );

        // Write incremental histograms with header
        info!("Writing histograms to file...");
        let histo_path = format!("{}{}.histo", directory, sample);
        warn_if_exists(&histo_path);
        let mut file =
            std::fs::File::create(&histo_path).context("Failed to create histogram file")?;

        // Comment line with version and parameters
        writeln!(file, "{}", histo_comment).context("Failed to write histogram comment")?;

        // Header row
        let header: String = std::iter::once("count".to_string())
            .chain((1..=histo_vecs.len()).map(|i| format!("chunk_{}", i)))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(file, "{}", header).context("Failed to write histogram header")?;

        // Data rows — histo_vecs already contains the precomputed vectors
        for i in 1..args.histo_max as usize + 2 {
            let line: String = std::iter::once(i.to_string())
                .chain(histo_vecs.iter().map(|v| v[i].to_string()))
                .collect::<Vec<_>>()
                .join("\t");
            writeln!(file, "{}", line).context("Failed to write histogram data")?;
        }

        // Write the final histogram with header
        info!("Writing final histogram to file...");
        let last_histo_vec = histo_vecs.last().context("No histograms were produced")?;

        let final_histo_path = format!("{}{}.final.histo", directory, sample);
        warn_if_exists(&final_histo_path);
        let mut file = std::fs::File::create(&final_histo_path)
            .context("Failed to create final histogram file")?;

        writeln!(file, "{}", histo_comment).context("Failed to write final histogram comment")?;
        writeln!(file, "count\tfrequency").context("Failed to write final histogram header")?;

        for (i, value) in last_histo_vec
            .iter()
            .enumerate()
            .skip(1)
            .take(args.histo_max as usize + 1)
        {
            writeln!(file, "{}\t{}", i, value).context("Failed to write final histogram data")?;
        }

        let singleton_count = *last_histo_vec
            .get(1)
            .context("Histogram vector too short to contain singleton count")?;
        n_singleton_kmers = Some(singleton_count);

        // Warn if singleton rate is very high (>95% of unique kmers)
        let n_unique = running_histo.get_n_unique_kmers();
        if n_unique > 0 {
            let singleton_rate = singleton_count as f64 / n_unique as f64;
            if singleton_rate > 0.95 {
                warn!(
                    "Very high singleton rate ({:.1}%). This may indicate very low coverage \
                     or contamination. sPCR results may be unreliable.",
                    singleton_rate * 100.0
                );
            }
        }

        let n_unique_kmers_histo: u64 = running_histo.get_n_unique_kmers();
        let n_kmers_histo: u64 = running_histo.get_n_kmers();

        debug!("{} unique kmers in histogram", n_unique_kmers_histo);
        debug!("{} kmers in histogram", n_kmers_histo);

        ensure!(
            n_kmers_histo == n_kmers_ingested,
            "The total count of kmers in the histogram ({}) does not equal the total expected count of kmers ({})",
            n_kmers_histo,
            n_kmers_ingested,
        );

        ensure!(
            n_unique_kmers_histo == kmer_counts.get_n_unique_kmers(),
            "The total count of unique kmers in the histogram ({}) does not equal the total count of hashed kmers ({})",
            n_unique_kmers_histo,
            kmer_counts.get_n_unique_kmers(),
        );
    } else {
        // No histogram mode: merge all chunks into a single kmer count table
        for chunk in state.chunks.drain(..) {
            kmer_counts.extend(chunk.get_kmer_counts())?;
            drop(chunk);
        }
        consolidate_spinner.finish_and_clear();
        info!(
            "Kmer counts merged, time: {}",
            format_duration(start.elapsed())
        );

        let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
        info!(
            "{} unique kmers with a total count of {} were found",
            kmer_counts.get_n_unique_kmers(),
            n_hashed_kmers
        );

        ensure!(
            n_hashed_kmers == n_kmers_ingested,
            "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
            n_hashed_kmers,
            n_kmers_ingested,
        );
    }

    Ok((FinalCounts::Legacy(kmer_counts), n_singleton_kmers))
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECORD_ONE: &str = "@read-one\nACGT\n+\n!!!!\n";
    const RECORD_TWO: &str = "@read-two\nTGCA\n+\n####\n";
    const RECORD_THREE: &str = "@read-three\nGATC\n+\n$$$$\n";
    const RECORD_FOUR: &str = "@read-four\nCTAG\n+\n%%%%\n";

    fn gzip_member(contents: &str) -> Vec<u8> {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(contents.as_bytes()).unwrap();
        encoder.finish().unwrap()
    }

    fn concatenated_gzip(records: &[&str]) -> Vec<u8> {
        records
            .iter()
            .flat_map(|record| gzip_member(record))
            .collect()
    }

    fn fastq_state() -> FastqReadState {
        let kmer_length = 3;
        FastqReadState {
            chunks: vec![Chunk::new(&kmer_length)],
            chunk_index: 0,
            seqs: Vec::new(),
            n_reads_read: 0,
            n_bases_read: 0,
            direct: None,
        }
    }

    fn read_sequences(
        path: &std::path::Path,
        max_reads: u64,
    ) -> Result<(Vec<String>, u64, u64, bool)> {
        let reader = open_fastq_reader(path)?;
        let mut state = fastq_state();
        let reached_max = read_fastq(
            reader,
            &mut state,
            max_reads,
            0,
            &path.to_string_lossy(),
            &ProgressBar::hidden(),
        )?;
        let sequences = state.seqs.clone();
        let n_chunks = state.chunks.len();
        drain_batch(&mut state, n_chunks)?;
        let n_reads = state.chunks.iter().map(Chunk::get_n_reads).sum();
        let n_kmers = state.chunks.iter().map(Chunk::get_n_kmers).sum();
        Ok((sequences, n_reads, n_kmers, reached_max))
    }

    #[test]
    fn retained_result_ceiling_shares_the_existing_combined_budget() {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::CountSlots as usize] = 200;
        classes[MemoryClass::InputBuffer as usize] = 50;
        classes[MemoryClass::RetainedResults as usize] = 999;
        let budget = DirectCounter::resource_budget(classes, 0, 100).unwrap();
        assert_eq!(budget.class_limit(MemoryClass::RetainedResults), 350);
        let retained = budget.reserve(MemoryClass::RetainedResults, 350).unwrap();
        let error = budget.reserve(MemoryClass::CountSlots, 1).err().unwrap();
        assert_eq!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .unwrap()
                .resource,
            crate::kmer::budget::ResourceKind::RunMemory
        );
        drop(retained);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn gzip_magic_detection_matches_plain_fastq() {
        let directory = tempfile::tempdir().unwrap();
        let plain_path = directory.path().join("reads.fastq");
        let compressed_path = directory.path().join("reads-compressed.fastq");
        let fastq = format!("{}{}", RECORD_ONE, RECORD_TWO);
        std::fs::write(&plain_path, &fastq).unwrap();
        std::fs::write(
            &compressed_path,
            concatenated_gzip(&[RECORD_ONE, RECORD_TWO]),
        )
        .unwrap();

        let compressed = read_sequences(&compressed_path, 0).unwrap();
        let plain = read_sequences(&plain_path, 0).unwrap();

        assert_eq!(compressed.0, plain.0);
        assert_eq!(compressed.1, plain.1);
        assert_eq!(compressed.2, plain.2);
    }

    #[test]
    fn concatenated_gzip_members_ingest_all_records() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.fastq.gz");
        std::fs::write(&path, concatenated_gzip(&[RECORD_ONE, RECORD_TWO])).unwrap();

        assert_eq!(
            read_sequences(&path, 0).unwrap().0,
            ["ACGT".to_string(), "TGCA".to_string()]
        );
    }

    #[test]
    fn corrupt_second_gzip_member_returns_an_error() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.fastq.gz");
        let mut corrupt_second_member = gzip_member(RECORD_TWO);
        *corrupt_second_member.last_mut().unwrap() ^= 0xff;
        let mut compressed = gzip_member(RECORD_ONE);
        compressed.extend(corrupt_second_member);
        std::fs::write(&path, compressed).unwrap();

        let error = read_sequences(&path, 0).unwrap_err();
        assert!(error.to_string().contains("Failed to read header line"));
    }

    #[test]
    fn paired_readers_continue_across_concatenated_members() {
        let directory = tempfile::tempdir().unwrap();
        let r1_path = directory.path().join("r1.fastq.gz");
        let r2_path = directory.path().join("r2.fastq.gz");
        std::fs::write(&r1_path, concatenated_gzip(&[RECORD_ONE, RECORD_TWO])).unwrap();
        std::fs::write(&r2_path, concatenated_gzip(&[RECORD_THREE, RECORD_FOUR])).unwrap();

        let mut state = fastq_state();
        let reached_max = read_fastq_paired(
            open_fastq_reader(&r1_path).unwrap(),
            open_fastq_reader(&r2_path).unwrap(),
            &mut state,
            0,
            0,
            "r1",
            "r2",
            &ProgressBar::hidden(),
        )
        .unwrap();

        assert!(!reached_max);
        assert_eq!(
            state.seqs,
            [
                "ACGT".to_string(),
                "GATC".to_string(),
                "TGCA".to_string(),
                "CTAG".to_string(),
            ]
        );
    }

    #[test]
    fn max_reads_stops_within_a_concatenated_member() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.fastq.gz");
        std::fs::write(
            &path,
            concatenated_gzip(&[
                &format!("{}{}", RECORD_ONE, RECORD_TWO),
                &format!("{}{}", RECORD_THREE, RECORD_FOUR),
            ]),
        )
        .unwrap();

        let (sequences, _, _, reached_max) = read_sequences(&path, 3).unwrap();
        assert!(reached_max);
        assert_eq!(
            sequences,
            ["ACGT".to_string(), "TGCA".to_string(), "GATC".to_string()]
        );
    }

    #[test]
    fn replay_matches_initial_ingestion_for_concatenated_members() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("cached.fastq.gz");
        std::fs::write(
            &path,
            concatenated_gzip(&[RECORD_ONE, RECORD_TWO, RECORD_THREE]),
        )
        .unwrap();

        let ingested_sequences = read_sequences(&path, 0).unwrap().0;
        let replayed = reread_sequences(
            &ReadPlan {
                replay: None,
                source: ReadSourcePlan::CachedRemote(vec![path]),
                paired: false,
                max_reads: 0,
            },
            false,
        )
        .unwrap();

        assert_eq!(
            replayed
                .iter()
                .map(|record| record.sequence.clone())
                .collect::<Vec<_>>(),
            ingested_sequences
        );
        assert!(replayed.iter().all(|record| record.mate == Mate::Unpaired));
        assert_eq!(
            replayed
                .iter()
                .map(|record| record.index)
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
    }
}

#[cfg(test)]
mod direct_tests;
