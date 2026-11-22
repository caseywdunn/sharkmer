//! Ordered selected-input hashing with bounded, reusable byte buffers.
//!
//! One worker consumes the existing length-framed byte stream. It does not
//! hash records independently or change the digest's serialization contract.

use crate::io::Mate;
use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
mod digest;
pub(crate) use digest::SelectedDigest;
use std::sync::{Arc, mpsc};
use std::thread::JoinHandle;

const BUFFER_BYTES: usize = 64 << 10;
const BUFFERS: usize = 3;
const QUEUED_BUFFERS: usize = 1;
const CHANNEL_ALLOWANCE_BYTES: usize = 4096;
const DOMAIN: &[u8] = b"sharkmer-selected-fastq-v1\0";

/// Actual input-identity worker ownership, omitted when hashing is inline.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct IdentityWorkerSummary {
    pub(crate) selected_records: u64,
    pub(crate) framed_bytes: u64,
    pub(crate) buffer_bytes: usize,
    pub(crate) queued_buffers: usize,
    pub(crate) owned_buffers: usize,
    pub(crate) reserved_bytes: u64,
}

struct CompletedIdentity {
    hash: SelectedDigest,
    framed_bytes: u64,
}

/// Own one coordinator buffer and the channels for one hashing worker.
pub(crate) struct IdentityWriter {
    work: Option<mpsc::SyncSender<Vec<u8>>>,
    recycled: Option<mpsc::Receiver<Vec<u8>>>,
    current: Vec<u8>,
    handle: Option<JoinHandle<Result<CompletedIdentity>>>,
    records: u64,
    framed_bytes: u64,
    failed: bool,
    // Drop closes channels and joins before any ownership charge is released.
    _pool: MemoryReservation,
}

impl IdentityWriter {
    /// Fixed pool and conservative channel allowance; not process RSS.
    pub(crate) fn reservation_bytes() -> Result<u64> {
        u64::try_from(BUFFERS * BUFFER_BYTES + CHANNEL_ALLOWANCE_BYTES + size_of::<Self>())
            .context("Input identity reservation overflow")
    }

    /// Reserve the whole pool before allocating buffers or starting a worker.
    pub(crate) fn new(budget: &Arc<ResourceBudget>) -> Result<Self> {
        let reservation = budget.reserve(MemoryClass::WorkerQueues, Self::reservation_bytes()?)?;
        let allocate = || -> Result<Vec<u8>> {
            let mut buffer = Vec::new();
            buffer
                .try_reserve_exact(BUFFER_BYTES)
                .context("Unable to allocate input identity buffer")?;
            ensure!(
                buffer.capacity() == BUFFER_BYTES,
                "Input identity allocation differs from its reservation"
            );
            Ok(buffer)
        };
        let current = allocate()?;
        let (work_sender, work_receiver) = mpsc::sync_channel::<Vec<u8>>(QUEUED_BUFFERS);
        let (recycle_sender, recycle_receiver) = mpsc::sync_channel::<Vec<u8>>(BUFFERS - 1);
        for _ in 1..BUFFERS {
            recycle_sender
                .send(allocate()?)
                .context("Unable to initialize input identity buffers")?;
        }
        let handle = std::thread::Builder::new()
            .name("input-identity".to_owned())
            .spawn(move || {
                let mut hash = SelectedDigest::new();
                hash.update(DOMAIN)?;
                let mut framed_bytes = 0u64;
                while let Ok(mut buffer) = work_receiver.recv() {
                    framed_bytes = framed_bytes
                        .checked_add(u64::try_from(buffer.len())?)
                        .context("Input identity byte total overflow")?;
                    hash.update(&buffer)?;
                    buffer.clear();
                    // Finish/drop can close recycling while queued work still drains.
                    let _ = recycle_sender.send(buffer);
                }
                Ok(CompletedIdentity { hash, framed_bytes })
            })
            .context("Unable to start input identity worker")?;
        Ok(Self {
            work: Some(work_sender),
            recycled: Some(recycle_receiver),
            current,
            handle: Some(handle),
            records: 0,
            framed_bytes: 0,
            failed: false,
            _pool: reservation,
        })
    }

    fn dispatch(&mut self, recycle: bool) -> Result<()> {
        if self.current.is_empty() {
            return Ok(());
        }
        self.work
            .as_ref()
            .context("Input identity queue is closed")?
            .send(std::mem::take(&mut self.current))
            .context("Input identity worker stopped before accepting bytes")?;
        if recycle {
            self.current = self
                .recycled
                .as_ref()
                .context("Input identity buffer pool is closed")?
                .recv()
                .context("Input identity worker stopped before returning a buffer")?;
            ensure!(
                self.current.is_empty() && self.current.capacity() == BUFFER_BYTES,
                "Invalid recycled input identity buffer"
            );
        }
        Ok(())
    }

    fn append(&mut self, mut bytes: &[u8]) -> Result<()> {
        while !bytes.is_empty() {
            let count = (BUFFER_BYTES - self.current.len()).min(bytes.len());
            self.current.extend_from_slice(&bytes[..count]);
            bytes = &bytes[count..];
            if self.current.len() == BUFFER_BYTES {
                self.dispatch(true)?;
            }
        }
        Ok(())
    }

    /// Append exactly the existing selected-fastq-v1 record representation.
    /// Large fields stream through the same fixed pool without a record clone.
    pub(crate) fn record(
        &mut self,
        index: u64,
        mate: Mate,
        source: &str,
        fields: [&[u8]; 4],
    ) -> Result<()> {
        ensure!(!self.failed, "Input identity producer previously failed");
        ensure!(self.work.is_some(), "Input identity queue is closed");
        ensure!(
            index == self.records,
            "Input identity record index is out of order"
        );
        let next_records = self
            .records
            .checked_add(1)
            .context("Input identity record count overflow")?;
        let values = [
            source.as_bytes(),
            fields[0],
            fields[1],
            fields[2],
            fields[3],
        ];
        let mut lengths = [0u64; 5];
        let mut added = 9u64 + 5 * 8;
        for (length, value) in lengths.iter_mut().zip(values) {
            *length = u64::try_from(value.len())?;
            added = added
                .checked_add(*length)
                .context("Input identity frame length overflow")?;
        }
        let next_bytes = self
            .framed_bytes
            .checked_add(added)
            .context("Input identity byte total overflow")?;
        let result = (|| -> Result<()> {
            let mut prefix = [0u8; 9];
            prefix[..8].copy_from_slice(&index.to_le_bytes());
            prefix[8] = match mate {
                Mate::Unpaired => 0,
                Mate::R1 => 1,
                Mate::R2 => 2,
            };
            self.append(&prefix)?;
            for (length, value) in lengths.into_iter().zip(values) {
                self.append(&length.to_le_bytes())?;
                self.append(value)?;
            }
            Ok(())
        })();
        if let Err(error) = result {
            self.failed = true;
            return Err(self.worker_error(error));
        }
        self.records = next_records;
        self.framed_bytes = next_bytes;
        Ok(())
    }

    fn close_channels(&mut self) {
        drop(self.work.take());
        drop(self.recycled.take());
    }

    fn worker_error(&mut self, channel_error: anyhow::Error) -> anyhow::Error {
        self.close_channels();
        let Some(handle) = self.handle.take() else {
            return channel_error;
        };
        match handle.join() {
            Ok(Err(error)) => error.context("Input identity worker failed"),
            Err(_) => anyhow::anyhow!("Input identity worker panicked"),
            Ok(Ok(_)) => channel_error,
        }
    }

    /// Drain the final partial buffer, join, and verify the worker's byte census.
    pub(crate) fn finish(self) -> Result<([u8; 32], IdentityWorkerSummary)> {
        let (hash, summary) = self.finish_state()?;
        Ok((hash.finalize()?, summary))
    }

    /// Drain and join without finalizing SHA, allowing the coordinator to
    /// continue the identical stream after a change in worker ownership.
    pub(crate) fn finish_state(mut self) -> Result<(SelectedDigest, IdentityWorkerSummary)> {
        ensure!(!self.failed, "Input identity producer previously failed");
        ensure!(self.work.is_some(), "Input identity queue is closed");
        if let Err(error) = self.dispatch(false) {
            return Err(self.worker_error(error));
        }
        self.close_channels();
        let completed = self
            .handle
            .take()
            .context("Missing input identity worker")?
            .join()
            .map_err(|_| anyhow::anyhow!("Input identity worker panicked"))??;
        ensure!(
            completed.framed_bytes == self.framed_bytes,
            "Input identity worker byte total differs from selected frames"
        );
        Ok((
            completed.hash,
            IdentityWorkerSummary {
                selected_records: self.records,
                framed_bytes: self.framed_bytes,
                buffer_bytes: BUFFER_BYTES,
                queued_buffers: QUEUED_BUFFERS,
                owned_buffers: BUFFERS,
                reserved_bytes: Self::reservation_bytes()?,
            },
        ))
    }
}

impl Drop for IdentityWriter {
    fn drop(&mut self) {
        self.close_channels();
        if let Some(handle) = self.handle.take() {
            match handle.join() {
                Ok(Ok(_)) => (),
                Ok(Err(error)) => log::debug!("Input identity worker stopped: {error:#}"),
                Err(_) => log::error!("Input identity worker panicked during cleanup"),
            }
        }
    }
}

#[cfg(test)]
mod tests;
