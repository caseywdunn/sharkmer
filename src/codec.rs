//! Preadmitted backing allocations for buffered FASTQ and pinned gzip codecs.
//!
//! The envelope includes realloc overlap and a possibly boxed wrapper. Library
//! objects always drop before their guard. Transport/TLS allocations and
//! allocator/OS overhead use process headroom, outside this codec envelope.

use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::sync::Arc;

const BUFFER_BYTES: usize = 8192;
const FLATE_BUFFER_BYTES: usize = 32 << 10;
// flate2 1.1.9: XLEN <=65535, filename/comment Vec capacities <=65536;
// only one header is retained across concatenated members. Include the old
// 32768-byte vector during its final growth, plus a header-CRC allocation.
const HEADER_BYTES: usize = 65535 + 2 * 65536 + 32768 + size_of::<flate2::Crc>();
// miniz_oxide 0.8.9: HuffmanOxide's three symbol tables, LocalBuf,
// HashBuffers.dict and its two u16 hash-chain tables. All other compressor
// backing data is inline in the public CompressorOxide size below.
const DEFLATE_CHILD_BYTES: usize = 3 * 288 * 5 + (65536 * 13) / 10 + (32768 + 258) + 2 * 32768 * 2;

fn reserve(
    bytes: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Option<MemoryReservation>> {
    budget
        .map(|budget| budget.reserve(MemoryClass::RetainedResults, u64::try_from(bytes)?))
        .transpose()
}

fn envelope<T>(backing: usize) -> Result<usize> {
    backing
        .checked_add(size_of::<T>())
        .context("Codec allocation envelope overflow")
}

/// Buffered source bytes and their backing/wrapper admission.
pub(crate) struct Buffered<R> {
    inner: BufReader<R>,
    _owner: Option<MemoryReservation>,
}

impl<R: Read> Buffered<R> {
    /// Admit before allocating the source buffer or its boxed wrapper.
    pub(crate) fn new(reader: R, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let owner = reserve(envelope::<Self>(BUFFER_BYTES)?, budget)?;
        Ok(Self {
            inner: BufReader::with_capacity(BUFFER_BYTES, reader),
            _owner: owner,
        })
    }
}

impl<R: Read> Read for Buffered<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        self.inner.read(bytes)
    }
}
impl<R: Read> BufRead for Buffered<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }
    fn consume(&mut self, bytes: usize) {
        self.inner.consume(bytes);
    }
}

/// Concatenated-gzip decoding, including every member's bounded header.
pub(crate) struct Decoder<R> {
    inner: BufReader<flate2::read::MultiGzDecoder<R>>,
    _owner: Option<MemoryReservation>,
}

impl<R: Read> Decoder<R> {
    fn allocation_bytes() -> Result<usize> {
        envelope::<Self>(
            size_of::<miniz_oxide::inflate::stream::InflateState>()
                + FLATE_BUFFER_BYTES
                + BUFFER_BYTES
                + HEADER_BYTES,
        )
    }

    /// Admission precedes decoder allocation and its eager first-header parse.
    pub(crate) fn new(reader: R, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let owner = reserve(Self::allocation_bytes()?, budget)?;
        Ok(Self {
            inner: BufReader::with_capacity(
                BUFFER_BYTES,
                flate2::read::MultiGzDecoder::new(reader),
            ),
            _owner: owner,
        })
    }
}

impl<R: Read> Read for Decoder<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        self.inner.read(bytes)
    }
}
impl<R: Read> BufRead for Decoder<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }
    fn consume(&mut self, bytes: usize) {
        self.inner.consume(bytes);
    }
}

/// Fast cache recompression whose reservation outlives gzip and buffer flushing.
pub(crate) struct Encoder<W: Write> {
    inner: flate2::write::GzEncoder<BufWriter<W>>,
    _owner: Option<MemoryReservation>,
}

impl<W: Write> Encoder<W> {
    fn allocation_bytes() -> Result<usize> {
        // GzBuilder::new().write emits a fixed 10-byte header, without optional
        // strings. zio's 32KiB vector never grows; run_vec uses spare capacity.
        envelope::<Self>(
            size_of::<miniz_oxide::deflate::core::CompressorOxide>()
                + DEFLATE_CHILD_BYTES
                + FLATE_BUFFER_BYTES
                + BUFFER_BYTES
                + 10,
        )
    }

    /// Admit all nested compressor/buffer allocations before constructing them.
    pub(crate) fn new(writer: W, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let owner = reserve(Self::allocation_bytes()?, budget)?;
        Ok(Self {
            inner: flate2::write::GzEncoder::new(
                BufWriter::with_capacity(BUFFER_BYTES, writer),
                flate2::Compression::fast(),
            ),
            _owner: owner,
        })
    }

    /// Keep the guard through the trailer, final buffer flush and error drops.
    pub(crate) fn finish(self) -> io::Result<W> {
        let Self { inner, _owner } = self;
        let writer = inner
            .finish()?
            .into_inner()
            .map_err(|error| error.into_error())?;
        drop(_owner);
        Ok(writer)
    }
}

impl<W: Write> Write for Encoder<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.inner.write(bytes)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests;
