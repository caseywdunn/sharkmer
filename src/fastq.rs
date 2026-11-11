//! Reusable, bounded four-line FASTQ parsing for the normal counter.

use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result, bail, ensure};
use std::io::BufRead;
use std::ops::Range;
use std::sync::Arc;

/// One reader owns one arena; returned fields borrow it until the next record.
pub(crate) struct FastqReader<R> {
    reader: R,
    bytes: Vec<u8>,
    fields: [Range<usize>; 4],
    limit: usize,
    budget: Arc<ResourceBudget>,
    reservation: MemoryReservation,
}

impl<R: BufRead> FastqReader<R> {
    /// Charge both buffers during growth, before allocating either replacement.
    pub(crate) fn new(reader: R, limit: usize, budget: Arc<ResourceBudget>) -> Result<Self> {
        ensure!(limit > 0, "FASTQ record byte limit must be positive");
        let reservation = budget.reserve(MemoryClass::InputBuffer, 0)?;
        Ok(Self {
            reader,
            bytes: Vec::new(),
            fields: std::array::from_fn(|_| 0..0),
            limit,
            budget,
            reservation,
        })
    }

    fn grow(&mut self, minimum: usize) -> Result<()> {
        if minimum <= self.bytes.capacity() {
            return Ok(());
        }
        let capacity = minimum
            .checked_next_power_of_two()
            .unwrap_or(self.limit)
            .min(self.limit);
        let mut reservation = self
            .budget
            .reserve(MemoryClass::InputBuffer, u64::try_from(capacity)?)?;
        let mut replacement = Vec::new();
        replacement
            .try_reserve_exact(capacity)
            .context("Unable to reserve FASTQ record buffer")?;
        ensure!(
            replacement.capacity() <= self.limit,
            "Allocated FASTQ arena exceeds its byte limit"
        );
        reservation.resize(u64::try_from(replacement.capacity())?)?;
        replacement.extend_from_slice(&self.bytes);
        self.bytes = replacement;
        self.reservation = reservation;
        Ok(())
    }

    /// Read one record, preserving original bytes except CRLF/LF delimiters.
    /// No line or frame can grow beyond the configured arena limit.
    pub(crate) fn next(&mut self, source: &str, index: u64) -> Result<bool> {
        self.bytes.clear();
        for (field, role) in ["header", "sequence", "separator", "quality"]
            .iter()
            .enumerate()
        {
            let start = self.bytes.len();
            let mut terminated = false;
            loop {
                let block = self
                    .reader
                    .fill_buf()
                    .map_err(|error| crate::io::stream_io_error(role, &error, source, index))?;
                if block.is_empty() {
                    break;
                }
                let n = block
                    .iter()
                    .position(|byte| *byte == b'\n')
                    .map_or(block.len(), |i| i + 1);
                let needed = self
                    .bytes
                    .len()
                    .checked_add(n)
                    .context("FASTQ record length overflow")?;
                ensure!(
                    needed <= self.limit,
                    "FASTQ record {} in {} exceeds its {} byte buffer limit",
                    index + 1,
                    source,
                    self.limit
                );
                // End the fill_buf borrow before arena growth.
                self.grow(needed)?;
                let block = self
                    .reader
                    .fill_buf()
                    .map_err(|error| crate::io::stream_io_error(role, &error, source, index))?;
                terminated = block[n - 1] == b'\n';
                self.bytes.extend_from_slice(&block[..n]);
                self.reader.consume(n);
                if terminated {
                    break;
                }
            }
            if self.bytes.len() == start {
                if field == 0 {
                    return Ok(false);
                }
                bail!(
                    "Truncated FASTQ record at record {} in {}: missing {} line",
                    index + 1,
                    source,
                    role
                );
            }
            let mut end = self.bytes.len();
            if terminated {
                end -= 1;
                if end > start && self.bytes[end - 1] == b'\r' {
                    end -= 1;
                }
            }
            // Preserve the released UTF-8 rejection policy even when periodic
            // structural validation is disabled.
            std::str::from_utf8(&self.bytes[start..end]).with_context(|| {
                format!(
                    "Invalid UTF-8 in FASTQ {role} line of record {} in {source}",
                    index + 1
                )
            })?;
            self.fields[field] = start..end;
        }
        Ok(true)
    }

    /// Return header, sequence, separator and quality without copying them.
    pub(crate) fn fields(&self) -> [&[u8]; 4] {
        std::array::from_fn(|i| &self.bytes[self.fields[i].clone()])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::ResourceLimits;
    use std::io::{BufReader, Cursor};

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; 8];
        classes[MemoryClass::InputBuffer as usize] = bytes;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
        })
        .unwrap()
    }

    #[test]
    fn arena_reuses_storage_and_preserves_fields_across_tiny_blocks_and_crlf() {
        let pool = budget(1024);
        let data = b"@first\r\nACNT\r\n+\r\n!@+#\r\n@b\nTA\n+\nII";
        let mut reader = FastqReader::new(
            BufReader::with_capacity(2, Cursor::new(data)),
            128,
            Arc::clone(&pool),
        )
        .unwrap();
        assert!(reader.next("source", 0).unwrap());
        assert_eq!(
            reader.fields(),
            [b"@first".as_slice(), b"ACNT", b"+", b"!@+#"]
        );
        let pointer = reader.bytes.as_ptr();
        let capacity = reader.bytes.capacity();
        assert!(reader.next("source", 1).unwrap());
        assert_eq!(reader.fields(), [b"@b".as_slice(), b"TA", b"+", b"II"]);
        assert_eq!(reader.bytes.as_ptr(), pointer);
        assert_eq!(reader.bytes.capacity(), capacity);
        assert!(!reader.next("source", 2).unwrap());
        assert!(
            pool.snapshot().peak_reserved_by_class[MemoryClass::InputBuffer as usize] <= 2 * 128
        );
        drop(reader);
        assert_eq!(pool.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn byte_and_shared_resize_limits_fail_before_unbounded_allocation() {
        let pool = budget(32);
        let data = b"@id\nACGTACGT\n+\nIIIIIIII\n";
        let mut reader = FastqReader::new(Cursor::new(data), 128, Arc::clone(&pool)).unwrap();
        assert!(reader.next("source", 0).is_err());
        assert!(pool.snapshot().peak_reserved_memory_bytes <= 32);
        drop(reader);
        assert_eq!(pool.snapshot().reserved_memory_bytes, 0);
        let pool = budget(1024);
        let mut reader = FastqReader::new(Cursor::new(data), 8, Arc::clone(&pool)).unwrap();
        assert!(
            reader
                .next("source", 0)
                .unwrap_err()
                .to_string()
                .contains("byte buffer limit")
        );
        assert!(reader.bytes.capacity() <= 8);
    }

    #[test]
    fn truncation_and_non_utf8_fail_with_record_and_source_context() {
        for (data, role) in [
            (b"@id\nAC\n+\n".as_slice(), "quality"),
            (b"@id\nAC\n", "separator"),
            (b"@id\n", "sequence"),
            (b"@id\nAC\n+\n\xff\xff\n", "quality"),
        ] {
            let mut reader = FastqReader::new(Cursor::new(data), 128, budget(1024)).unwrap();
            let error = reader.next("broken.fastq", 7).unwrap_err().to_string();
            assert!(error.contains(role), "{error}");
            assert!(error.contains("record 8"), "{error}");
            assert!(error.contains("broken.fastq"), "{error}");
        }
    }
}
