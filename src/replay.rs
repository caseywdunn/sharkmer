//! Bounded, integrity-checked replay of the selected original FASTQ records.

use crate::io::Mate;
use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use crate::kmer::temporary::OwnedDiskFile;
use anyhow::{Context, Result, ensure};
use rand::RngCore;
use sha2::{Digest, Sha256};
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};

const HEADER: usize = 64;
const PREFIX: usize = 64;
const CHECKSUM: usize = 32;
const FOOTER: usize = 96;
const MAGIC: &[u8; 8] = b"SHKRPL01";
const COMPLETE: &[u8; 8] = b"SHKREN01";

/// Report only actual replay storage and its owned-buffer measurements.
#[cfg_attr(test, derive(Clone))]
#[derive(serde::Serialize)]
pub(crate) struct ReplaySummary {
    format_version: u32,
    selected_records: u64,
    sequence_bases: u64,
    selected_input_sha256: crate::kmer::budget::text::BudgetString,
    file_bytes: u64,
    max_frame_bytes: usize,
    peak_replay_buffer_bytes: u64,
}

fn integer(bytes: &[u8], offset: usize) -> u64 {
    let mut result = [0; 8];
    result.copy_from_slice(&bytes[offset..offset + 8]);
    u64::from_le_bytes(result)
}

fn selected_hash() -> Sha256 {
    let mut hash = Sha256::new();
    hash.update(b"sharkmer-selected-fastq-v1\0");
    hash
}

/// A borrowed original record; sequence positions and quality bytes are intact.
pub(crate) struct ReplayRecord<'a> {
    pub(crate) index: u64,
    pub(crate) fragment: u64,
    pub(crate) mate: Mate,
    pub(crate) source: &'a str,
    pub(crate) fields: [&'a [u8]; 4],
}

impl ReplayRecord<'_> {
    fn mate_code(&self) -> u8 {
        match self.mate {
            Mate::Unpaired => 0,
            Mate::R1 => 1,
            Mate::R2 => 2,
        }
    }

    fn update_identity(&self, hash: &mut Sha256) -> Result<()> {
        hash.update(self.index.to_le_bytes());
        hash.update([self.mate_code()]);
        for field in std::iter::once(self.source.as_bytes()).chain(self.fields) {
            hash.update(u64::try_from(field.len())?.to_le_bytes());
            hash.update(field);
        }
        Ok(())
    }

    fn validate(&self, paired: bool) -> Result<()> {
        ensure!(
            self.fragment == if paired { self.index / 2 } else { self.index },
            "Replay fragment identity differs from selected-record policy"
        );
        ensure!(
            (self.mate != Mate::Unpaired) == paired,
            "Replay mate differs from selected-record policy"
        );
        ensure!(
            self.fields[0].starts_with(b"@"),
            "Invalid replay FASTQ header"
        );
        ensure!(
            self.fields[2].starts_with(b"+"),
            "Invalid replay FASTQ separator"
        );
        ensure!(
            self.fields[1].len() == self.fields[3].len(),
            "Replay sequence/quality length mismatch"
        );
        for field in self.fields {
            std::str::from_utf8(field).context("Invalid UTF-8 in replay FASTQ")?;
        }
        Ok(())
    }
}

struct Buffer {
    bytes: Vec<u8>,
    budget: Arc<ResourceBudget>,
    limit: usize,
    reservation: MemoryReservation,
}

impl Buffer {
    fn new(budget: &Arc<ResourceBudget>, limit: usize) -> Result<Self> {
        Ok(Self {
            bytes: Vec::new(),
            budget: Arc::clone(budget),
            limit,
            reservation: budget.reserve(MemoryClass::ReplayBuffer, 0)?,
        })
    }

    fn resize(&mut self, length: usize) -> Result<()> {
        ensure!(
            length <= self.limit,
            "Replay frame exceeds configured byte limit"
        );
        if length > self.bytes.capacity() {
            let capacity = length
                .checked_next_power_of_two()
                .unwrap_or(self.limit)
                .min(self.limit);
            let mut reservation = self
                .budget
                .reserve(MemoryClass::ReplayBuffer, u64::try_from(capacity)?)?;
            let mut bytes = Vec::new();
            bytes
                .try_reserve_exact(capacity)
                .context("Unable to allocate replay buffer")?;
            ensure!(
                bytes.capacity() <= self.limit,
                "Replay allocation exceeds frame limit"
            );
            reservation.resize(u64::try_from(bytes.capacity())?)?;
            self.bytes = bytes;
            self.reservation = reservation;
        }
        self.bytes.resize(length, 0);
        Ok(())
    }
}

/// An incomplete writer cannot be used as replay until its footer is verified.
pub(crate) struct ReplayWriter {
    file: OwnedDiskFile,
    buffer: Buffer,
    header: [u8; HEADER],
    records: u64,
    bases: u64,
    selected: Sha256,
    stream: Sha256,
    failed: bool,
}

impl ReplayWriter {
    /// Store the effective selection cap and pairing policy in the header.
    pub(crate) fn new(
        directory: &Path,
        budget: &Arc<ResourceBudget>,
        max_frame_bytes: usize,
        paired: bool,
        max_reads: u64,
    ) -> Result<Self> {
        ensure!(
            max_frame_bytes >= PREFIX + CHECKSUM,
            "Replay frame limit is too small"
        );
        let mut header = [0; HEADER];
        header[..8].copy_from_slice(MAGIC);
        header[8..12].copy_from_slice(&1u32.to_le_bytes());
        header[12] = u8::from(paired);
        header[13] = 1; // original UTF-8 FASTQ fields, quality encoding unchanged
        header[16..24].copy_from_slice(&u64::try_from(max_frame_bytes)?.to_le_bytes());
        header[24..32].copy_from_slice(&max_reads.to_le_bytes());
        rand::rngs::OsRng
            .try_fill_bytes(&mut header[32..48])
            .context("Unable to initialize replay file identity")?;
        let mut file = OwnedDiskFile::new(directory, budget)?;
        file.append(&header)?;
        let mut stream = Sha256::new();
        stream.update(header);
        Ok(Self {
            file,
            buffer: Buffer::new(budget, max_frame_bytes)?,
            header,
            records: 0,
            bases: 0,
            selected: selected_hash(),
            stream,
            failed: false,
        })
    }

    /// Copy just one bounded frame into reusable storage before appending it.
    pub(crate) fn add(&mut self, record: ReplayRecord<'_>) -> Result<()> {
        ensure!(!self.failed, "Replay writer previously failed");
        let result = self.add_record(record);
        if result.is_err() {
            self.failed = true;
        }
        result
    }

    fn add_record(&mut self, record: ReplayRecord<'_>) -> Result<()> {
        record.validate(self.header[12] == 1)?;
        ensure!(record.index == self.records, "Replay record order changed");
        let cap = integer(&self.header, 24);
        ensure!(
            cap == 0 || self.records < cap,
            "Replay exceeds effective selected-read cap"
        );
        let next = self
            .records
            .checked_add(1)
            .context("Replay record count overflow")?;
        let bases = self
            .bases
            .checked_add(u64::try_from(record.fields[1].len())?)
            .context("Replay base count overflow")?;
        let fields = [
            record.source.as_bytes(),
            record.fields[0],
            record.fields[1],
            record.fields[2],
            record.fields[3],
        ];
        let length = fields.iter().try_fold(PREFIX + CHECKSUM, |sum, field| {
            sum.checked_add(field.len())
                .context("Replay frame length overflow")
        })?;
        self.buffer.resize(length)?;
        let bytes = &mut self.buffer.bytes;
        bytes[..PREFIX].fill(0);
        bytes[..8].copy_from_slice(&record.index.to_le_bytes());
        bytes[8..16].copy_from_slice(&record.fragment.to_le_bytes());
        bytes[16] = record.mate_code();
        let mut offset = PREFIX;
        for (index, field) in fields.iter().enumerate() {
            bytes[24 + index * 8..32 + index * 8]
                .copy_from_slice(&u64::try_from(field.len())?.to_le_bytes());
            bytes[offset..offset + field.len()].copy_from_slice(field);
            offset += field.len();
        }
        let checksum = frame_checksum(&self.header, &bytes[..offset]);
        bytes[offset..].copy_from_slice(&checksum);
        self.file.append(bytes)?;
        self.stream.update(bytes);
        record.update_identity(&mut self.selected)?;
        self.records = next;
        self.bases = bases;
        Ok(())
    }

    /// Bind replay to the count pass's actual selected-record digest.
    pub(crate) fn finish(mut self, expected_digest: [u8; 32]) -> Result<ReplayStore> {
        ensure!(!self.failed, "Replay writer previously failed");
        let selected: [u8; 32] = self.selected.finalize().into();
        ensure!(
            selected == expected_digest,
            "Replay identity differs from counted input"
        );
        let mut footer = [0; FOOTER];
        footer[..8].copy_from_slice(COMPLETE);
        footer[8..16].copy_from_slice(&self.records.to_le_bytes());
        footer[16..24].copy_from_slice(&self.bases.to_le_bytes());
        footer[24..32].copy_from_slice(&self.file.len()?.to_le_bytes());
        footer[32..64].copy_from_slice(&selected);
        self.stream.update(&footer[..64]);
        footer[64..].copy_from_slice(&self.stream.finalize());
        self.file.append(&footer)?;
        self.file.sync()?;
        let budget = Arc::clone(&self.buffer.budget);
        let limit = self.buffer.limit;
        drop(self.buffer);
        ReplayStore::open(self.file, &budget, limit, expected_digest)
    }
}

fn frame_checksum(header: &[u8; HEADER], bytes: &[u8]) -> [u8; 32] {
    let mut hash = Sha256::new();
    hash.update(b"sharkmer-replay-frame-v1\0");
    hash.update(header);
    hash.update(bytes);
    hash.finalize().into()
}

/// Sharing a replay file does not share or retain all records in memory.
pub(crate) struct ReplayStore {
    file: Mutex<OwnedDiskFile>,
    budget: Arc<ResourceBudget>,
    header: [u8; HEADER],
    footer: [u8; FOOTER],
    limit: usize,
}

impl ReplayStore {
    fn open(
        mut file: OwnedDiskFile,
        budget: &Arc<ResourceBudget>,
        limit: usize,
        expected_digest: [u8; 32],
    ) -> Result<Self> {
        file.verify_extent()?;
        ensure!(
            file.len()? >= (HEADER + FOOTER) as u64,
            "Truncated replay file"
        );
        let mut header = [0; HEADER];
        file.read_exact_at(0, &mut header)?;
        ensure!(
            &header[..8] == MAGIC && header[8..12] == 1u32.to_le_bytes(),
            "Incompatible replay format"
        );
        ensure!(
            header[12] <= 1
                && header[13] == 1
                && header[14..16] == [0; 2]
                && header[48..] == [0; 16],
            "Incompatible replay policy"
        );
        ensure!(
            integer(&header, 16) == u64::try_from(limit)?,
            "Replay frame limit changed"
        );
        let mut footer = [0; FOOTER];
        let end = file.len()? - FOOTER as u64;
        file.read_exact_at(end, &mut footer)?;
        ensure!(
            &footer[..8] == COMPLETE && integer(&footer, 24) == end,
            "Incomplete replay file"
        );
        ensure!(
            footer[32..64] == expected_digest,
            "Replay selected-input identity mismatch"
        );
        let result = Self {
            file: Mutex::new(file),
            budget: Arc::clone(budget),
            header,
            footer,
            limit,
        };
        result.visit(&mut |_| Ok(()))?;
        Ok(result)
    }

    /// The declared pairing policy comes from the original selection pass.
    pub(crate) fn paired(&self) -> bool {
        self.header[12] == 1
    }

    /// Number of verified selected records, rather than the requested cap.
    pub(crate) fn record_count(&self) -> u64 {
        integer(&self.footer, 8)
    }

    /// Measurements and identity of the verified original-pass replay.
    pub(crate) fn summary(&self, budget: Option<&Arc<ResourceBudget>>) -> Result<ReplaySummary> {
        Ok(ReplaySummary {
            format_version: 1,
            selected_records: self.record_count(),
            sequence_bases: integer(&self.footer, 16),
            selected_input_sha256: crate::kmer::budget::text::BudgetString::hex(
                &self.footer[32..64],
                budget,
                MemoryClass::RetainedResults,
            )?,
            file_bytes: integer(&self.footer, 24) + FOOTER as u64,
            max_frame_bytes: self.limit,
            peak_replay_buffer_bytes: self.budget.snapshot().peak_reserved_by_class
                [MemoryClass::ReplayBuffer as usize],
        })
    }

    /// Validate frames before callbacks and verify completion at the end.
    pub(crate) fn visit(
        &self,
        visitor: &mut dyn FnMut(ReplayRecord<'_>) -> Result<()>,
    ) -> Result<()> {
        let mut buffer = Buffer::new(&self.budget, self.limit)?;
        let mut stream = Sha256::new();
        stream.update(self.header);
        let mut selected = selected_hash();
        let (mut offset, mut records, mut bases) = (HEADER as u64, 0u64, 0u64);
        let end = integer(&self.footer, 24);
        {
            let mut file = self
                .file
                .lock()
                .map_err(|_| anyhow::anyhow!("Replay file lock poisoned"))?;
            file.verify_extent()?;
            ensure!(file.len()? == end + FOOTER as u64, "Replay extent changed");
            let mut header = [0; HEADER];
            let mut footer = [0; FOOTER];
            file.read_exact_at(0, &mut header)?;
            file.read_exact_at(end, &mut footer)?;
            ensure!(
                header == self.header && footer == self.footer,
                "Replay identity or completion changed"
            );
        }
        while offset < end {
            ensure!(
                end - offset >= (PREFIX + CHECKSUM) as u64,
                "Truncated replay frame"
            );
            let mut prefix = [0; PREFIX];
            let mut ranges: [Range<usize>; 5] = std::array::from_fn(|_| 0..0);
            {
                let mut file = self
                    .file
                    .lock()
                    .map_err(|_| anyhow::anyhow!("Replay file lock poisoned"))?;
                file.read_exact_at(offset, &mut prefix)?;
                let mut length = PREFIX;
                for (index, range) in ranges.iter_mut().enumerate() {
                    let size = usize::try_from(integer(&prefix, 24 + index * 8))
                        .context("Replay field exceeds address space")?;
                    let next = length
                        .checked_add(size)
                        .context("Replay frame length overflow")?;
                    *range = length..next;
                    length = next;
                }
                length = length
                    .checked_add(CHECKSUM)
                    .context("Replay checksum extent overflow")?;
                ensure!(
                    u64::try_from(length)? <= end - offset,
                    "Replay frame exceeds completed extent"
                );
                // Both the caller's limit and the completed extent are checked
                // before allocating storage for any declared field length.
                buffer.resize(length)?;
                buffer.bytes[..PREFIX].copy_from_slice(&prefix);
                file.read_exact_at(offset + PREFIX as u64, &mut buffer.bytes[PREFIX..])?;
            }
            let bytes = &buffer.bytes;
            let payload_end = bytes.len() - CHECKSUM;
            ensure!(
                frame_checksum(&self.header, &bytes[..payload_end]) == bytes[payload_end..],
                "Replay frame checksum mismatch"
            );
            ensure!(prefix[17..24] == [0; 7], "Invalid replay frame flags");
            let mate = match prefix[16] {
                0 => Mate::Unpaired,
                1 => Mate::R1,
                2 => Mate::R2,
                _ => anyhow::bail!("Invalid replay mate"),
            };
            let record = ReplayRecord {
                index: integer(&prefix, 0),
                fragment: integer(&prefix, 8),
                mate,
                source: std::str::from_utf8(&bytes[ranges[0].clone()])?,
                fields: std::array::from_fn(|index| &bytes[ranges[index + 1].clone()]),
            };
            ensure!(record.index == records, "Replay record order changed");
            record.validate(self.paired())?;
            record.update_identity(&mut selected)?;
            records = records
                .checked_add(1)
                .context("Replay record total overflow")?;
            bases = bases
                .checked_add(u64::try_from(record.fields[1].len())?)
                .context("Replay base total overflow")?;
            offset = offset
                .checked_add(u64::try_from(bytes.len())?)
                .context("Replay extent overflow")?;
            stream.update(bytes);
            visitor(record)?;
        }
        stream.update(&self.footer[..64]);
        let digest: [u8; 32] = selected.finalize().into();
        ensure!(
            records == integer(&self.footer, 8)
                && bases == integer(&self.footer, 16)
                && digest == self.footer[32..64],
            "Replay totals or input identity mismatch"
        );
        let hash: [u8; 32] = stream.finalize().into();
        ensure!(
            hash == self.footer[64..],
            "Replay completed-stream checksum mismatch"
        );
        let cap = integer(&self.header, 24);
        ensure!(
            cap == 0 || records <= cap,
            "Replay exceeds effective read cap"
        );
        Ok(())
    }

    /// Confirm removal before allowing final output publication.
    pub(crate) fn close(self) -> Result<()> {
        self.file
            .into_inner()
            .map_err(|_| anyhow::anyhow!("Replay file lock poisoned"))?
            .close()
    }
}

#[cfg(test)]
mod tests;
