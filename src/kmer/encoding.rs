// kmer/encoding.rs
//! 2-bit DNA encoding, kmer extraction, and sequence conversion utilities.

use anyhow::{Result, bail, ensure};

/// A structure to hold a read in two bit encoding of bases.
///
/// Retained for tests and as a reference implementation; the hot path now uses
/// `kmers_from_ascii` to skip the encode/decode round-trip.
#[allow(dead_code)]
/// * `00` represents `A`
/// * `01` represents `C`
/// * `10` represents `G`
/// * `11` represents `T`
///
/// If the length is not divisible by 4, the least significant bits of the last byte will be ignored.
#[derive(Debug)]
pub struct Read {
    pub sequence: Vec<u8>,
    pub length: usize, // Number of bases in the read
}

#[allow(dead_code)]
impl Read {
    pub fn new(sequence: Vec<u8>, length: usize) -> Read {
        Read { sequence, length }
    }

    /// This function encodes a DNA sequence into a series of 8-bit integers,
    /// where each integer represents 4 consecutive bases from the sequence.
    /// The encoding uses a 2-bit representation for each base:
    ///
    /// * `A` -> `00`
    /// * `C` -> `01`
    /// * `G` -> `10`
    /// * `T` -> `11`
    ///
    /// Any other character will cause the function to panic.
    pub fn from_str(seq: &str) -> Result<Read> {
        let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4 + 1);
        let mut frame: u8 = 0; // The current 8-bit integer being constructed
        let mut length: usize = 0; // length of the current subread
        for c in seq.chars() {
            length += 1;
            let base = match c {
                'A' => 0, // 00
                'C' => 1, // 01
                'G' => 2, // 10
                'T' => 3, // 11
                _ => bail!(
                    "Invalid character {} in sequence {}. Only ACGT allowed.",
                    c,
                    seq
                ),
            };

            frame = (frame << 2) | base;
            if length % 4 == 0 {
                ints.push(frame);
                frame = 0; // Reset frame after pushing to the vector
            }
        }

        // Check if there is anything left in the frame, and rotate it into place if so
        let modulo = length % 4;
        if modulo != 0 {
            let shift = 2 * (4 - modulo);
            frame <<= shift;
            ints.push(frame);
        }

        // Create and return the read
        Ok(Read::new(ints, length))
    }

    pub fn validate(&self) -> bool {
        let mut valid = true;
        if self.length % 4 == 0 {
            if self.sequence.len() != self.length / 4 {
                valid = false;
            }
        } else if self.sequence.len() != (self.length / 4 + 1) {
            valid = false;
        }

        valid
    }

    /// Get canonical kmers from a Read.
    ///
    /// Given a sequence that's been encoded as a series of 8-bit integers (where each integer
    /// represents 4 bases), this function extracts canonical kmers from it. A canonical kmer
    /// is the lexicographically smaller of a kmer and its reverse complement.
    ///
    /// The encoding uses a 2-bit representation for each base:
    ///
    /// * `A` -> `00`
    /// * `C` -> `01`
    /// * `G` -> `10`
    /// * `T` -> `11`
    ///
    /// # Arguments
    ///
    /// * `ints` - A reference to the vector of 8-bit integers which encode a DNA sequence.
    /// * `k` - A reference to the length of the desired kmers (i.e., the number of bases in each kmer).
    ///
    /// # Returns
    ///
    /// A vector of 64-bit integers, where each integer represents a canonical kmer extracted
    /// from the encoded sequence.
    pub fn get_kmers(&self, k: &usize) -> Result<Vec<u64>> {
        ensure!(*k > 0 && *k < 32, "k must be between 1 and 31, got {}", k);
        let ints = &self.sequence;
        let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k) + 1);
        let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
        let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
        let mut n_valid = 0; // number of valid bases in the frame
        let mask: u64 = (1 << (2 * k)) - 1;

        // No kmers if the read is shorter than k
        if self.length < *k {
            return Ok(kmers);
        }

        // Iterate over the bases
        for &int in ints.iter() {
            // Iterate over the bases in the integer
            for j in 0..4 {
                // Get the base from the left side of the integer,
                // move it to the least two significant bits and mask it
                let base = ((int >> ((3 - j) * 2)) & 3) as u64;

                frame = (frame << 2) | base;
                revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));
                n_valid += 1;
                if n_valid >= *k {
                    let forward = frame & mask;
                    let reverse = revframe & mask;

                    if forward < reverse {
                        kmers.push(forward);
                    } else {
                        kmers.push(reverse);
                    }
                }
            }
        }

        // If the read length is not divisible by 4, the last byte will have some extra bases
        // that result in extra kmers we need to remove
        let modulo = self.length % 4;
        if modulo != 0 {
            let n_extra_bases = 4 - modulo;
            kmers.truncate(kmers.len().saturating_sub(n_extra_bases));
        }

        let expected = self.length - *k + 1;
        if kmers.len() != expected {
            bail!(
                "Number of kmers ({}) does not match expected ({}) for read of length {}",
                kmers.len(),
                expected,
                self.length
            );
        }

        Ok(kmers)
    }
}

impl PartialEq for Read {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.length == other.length
    }
}

/// Returns the reverse complement of a specified kmer.
///
/// The function computes the reverse complement of a kmer represented as a 64-bit unsigned integer.
/// Each base in the kmer is encoded using 2 bits, with the following scheme:
///
/// * `00` represents `A`
/// * `01` represents `C`
/// * `10` represents `G`
/// * `11` represents `T`
///
/// # Arguments
///
/// * `kmer` - A reference to the 64-bit unsigned integer representation of the kmer.
/// * `k` - A reference to the length of the kmer (i.e., the number of bases).
///
/// # Returns
///
/// A 64-bit unsigned integer representing the reverse complement of the given kmer.
pub fn revcomp_kmer(kmer: &u64, k: &usize) -> u64 {
    let mut revcomp = 0;
    for i in 0..*k {
        let base = (kmer >> (2 * i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

/// Converts a DNA sequence into a vector of Reads (really subreads).
///
///
/// Any sequence containing the base `N` is split into multiple subsequences at that point,
/// and each sub-sequence is encoded separately. The result is a vector of Reads,
/// with each Read holding the integer representation of a sub-sequence.
///
/// # Arguments
///
/// * `seq` - A reference to the DNA sequence string to be encoded.
///
/// # Returns
///
/// A vector of Reads, where each inner vector represents an encoded
/// DNA sub-sequence.
///
/// # Panics
///
/// The function will exit prematurely if it encounters a base other than `A`, `C`, `G`, `T`, or `N`.
#[allow(dead_code)]
pub fn seq_to_reads(seq: &str) -> Result<Vec<Read>> {
    let mut reads: Vec<Read> = Vec::new();

    // Split seq on N
    for subseq in seq.split('N') {
        if !subseq.is_empty() {
            let read = Read::from_str(subseq)?;
            if !read.validate() {
                bail!("Invalid read: {:?} from subsequence {}", read, subseq);
            }
            reads.push(read);
        }
    }
    Ok(reads)
}

pub fn kmer_to_seq(kmer: &u64, k: &usize) -> String {
    let mut seq = String::new();
    for i in 0..*k {
        let base = (kmer >> (2 * (k - i - 1))) & 3;
        let base = match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!("2-bit encoding only produces values 0-3"),
        };
        seq.push(base);
    }
    seq
}

/// Extract canonical kmers directly from an ASCII sequence in a single pass.
///
/// Splits on `N` (equivalent to `seq_to_reads` + `get_kmers`), but avoids the
/// intermediate `Read` struct and 2-bit packing round-trip. Returns the same
/// canonical kmers in the same order.
pub fn kmers_from_ascii(seq: &str, k: usize) -> Result<Vec<u64>> {
    ensure!(k > 0 && k < 32, "k must be between 1 and 31, got {}", k);
    let mask: u64 = (1 << (2 * k)) - 1;
    let mut kmers: Vec<u64> = Vec::with_capacity(seq.len().saturating_sub(k - 1));
    let mut frame: u64 = 0;
    let mut revframe: u64 = 0;
    let mut n_valid: usize = 0;

    for &b in seq.as_bytes() {
        let base: u64 = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => {
                // Reset — equivalent to starting a new subread
                n_valid = 0;
                frame = 0;
                revframe = 0;
                continue;
            }
            _ => bail!(
                "Invalid character '{}' in sequence. Only ACGTN allowed.",
                b as char
            ),
        };

        frame = (frame << 2) | base;
        revframe = (revframe >> 2) | ((3 - base) << (2 * (k - 1)));
        n_valid += 1;

        if n_valid >= k {
            let forward = frame & mask;
            let reverse = revframe & mask;
            kmers.push(forward.min(reverse));
        }
    }

    Ok(kmers)
}

/// Count the number of non-N bases in a sequence.
pub fn count_valid_bases(seq: &str) -> u64 {
    seq.as_bytes().iter().filter(|&&b| b != b'N').count() as u64
}

#[allow(dead_code)]
pub fn seq_to_kmer(seq: &str) -> Result<u64> {
    let mut kmer = 0;
    for c in seq.chars() {
        let base: u64 = match c {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => bail!("Invalid base '{}' in sequence '{}'", c, seq),
        };
        kmer = (kmer << 2) | base;
    }
    Ok(kmer)
}
