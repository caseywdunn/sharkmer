//! This module provides kmer functions.

use rustc_hash::FxHashMap;

// A structure with a hashmap for kmer counts and a u64 for the number of singleton kmers
pub struct KmerSummary {
    pub kmer_counts: FxHashMap<u64, u64>,
    pub n_singletons: u64,
}

/// A structure to hold a read in two bit encoding of bases
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

impl Read {
    pub fn new(sequence: Vec<u8>, length: usize) -> Read {
        Read { sequence, length }
    }

    pub fn validate(&self) -> bool {
        let mut valid = true;
        if self.length % 4 == 0 {
            if self.sequence.len() != self.length / 4 {
                valid = false;
            }
        } else {
            if self.sequence.len() != (self.length / 4 + 1) {
                valid = false;
            }
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
    ///
    /// # Example
    ///
    /// ```rust
    /// let encoded_sequence: Vec<u8> = /* some encoded sequence */;
    /// let kmer_length: usize = /* desired kmer length */;
    /// let canonical_kmers = get_kmers(&encoded_sequence, &kmer_length);
    /// ```
    pub fn get_kmers(&self, k: &usize) -> Vec<u64> {
        let ints = &self.sequence; // TODO: handle case where length is not divisible by 4
        let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k) + 1);
        let mut frame: u64 = 0; // read the bits for each base into the least significant end of this integer
        let mut revframe: u64 = 0; // read the bits for complement into the least significant end of this integer
        let mut n_valid = 0; // number of valid bases in the frame
        let mask: u64 = (1 << (2 * k)) - 1;

        // Iterate over the bases
        for (_i, &int) in ints.iter().enumerate() {
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

                    // assert_eq!(forward, revcomp_kmer(reverse, k));

                    if forward < reverse {
                        kmers.push(forward);
                    } else {
                        kmers.push(reverse);
                    }
                }
            }
        }
        kmers
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
///
/// # Example
///
/// ```rust
/// let kmer: u64 = /* some kmer encoding */;
/// let k: usize = /* kmer length */;
/// let reverse_complement = revcomp_kmer(&kmer, &k);
/// ```
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
/// This function encodes a DNA sequence into a series of 8-bit integers,
/// where each integer represents 4 consecutive bases from the sequence.
/// The encoding uses a 2-bit representation for each base:
///
/// * `A` -> `00`
/// * `C` -> `01`
/// * `G` -> `10`
/// * `T` -> `11`
///
/// Any sequence containing the base `N` is split into multiple sub-sequences at that point,
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
///
/// # Example
///
/// ```rust
/// let sequence = "ACGTNAGCT";
/// let reads = seq_to_ints(&sequence);
/// ```
pub fn seq_to_reads(seq: &str) -> Vec<Read> {
    let mut reads: Vec<Read> = Vec::new();
    let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4 + 1);
    let mut frame: u8 = 0; // The current 8-bit integer being constructed
    let mut position: usize = 0; // position in the sequence. not including Ns
    let mut length: usize = 0; // length of the read
    for c in seq.chars() {
        length += 1;
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            'N' => {
                if !ints.is_empty() {
                    // Check if there is anything left in the frame,
                    // if so shift it left to the most significan bits and push it to the vector
                    let modulo = (position + 1) % 4;
                    if modulo != 0 {
                        let shift = 2 * (4 - modulo);
                        frame = frame << shift;
                        ints.push(frame);
                    }

                    // Create and push the read to the vector
                    let read = Read::new(ints, length);
                    reads.push(read);
                    length = 0;
                    ints = Vec::with_capacity(seq.len() / 4 + 1);
                    frame = 0; // Reset frame before starting a new subread
                    position = 0;
                }
                continue;
            }
            _ => 5,
        };
        if base > 3 {
            break;
        }
        frame = (frame << 2) | base;
        if ((position + 1) % 4 == 0) & (position > 0) {
            ints.push(frame);
            frame = 0; // Reset frame after pushing to the vector
        }
        position += 1;
    }
    if !ints.is_empty() || reads.is_empty() {
        // Check if there is anything left in the frame,
        // if so shift it left to the most significan bits and push it to the vector
        let modulo = (position + 1) % 4;
        if modulo != 0 {
            let shift = 2 * (4 - modulo);
            frame = frame << shift;
            ints.push(frame);
        }

        // Create and push the read to the vector
        let read = Read::new(ints, position);
        reads.push(read);
        length = 0;
    }
    reads
}



/// Generates a histogram from kmer counts.
///
/// This function produces a histogram where the indices represent the counts of a kmer,
/// and the values at those indices represent the number of kmers with that count.
/// If a kmer's count exceeds the specified `histo_max`, it gets placed in the final bucket.
///
/// # Arguments
///
/// * `kmer_counts` - A reference to a `FxHashMap` where keys are kmers (as 64-bit integers)
///   and values are the corresponding counts of each kmer.
/// * `histo_max` - A reference to the maximum count to be considered for individual bins
///   in the histogram. Kmer counts exceeding this value will be lumped into a single
///   overflow bin.
///
/// # Returns
///
/// A vector representing the histogram. The value at index `i` represents the number of
/// kmers that appeared `i` times. The last value in the vector represents the number of
/// kmers with counts greater than `histo_max`.
///
/// # Example
///
/// ```rust
/// use rustc_hash::FxHashMap;
///
/// let mut kmer_counts = FxHashMap::default();
/// kmer_counts.insert(12345, 3);
/// kmer_counts.insert(67890, 5);
/// // ... populate more kmers and counts ...
///
/// let histo_max = 10;
/// let histogram = count_histogram(&kmer_counts, &histo_max);
/// // histogram will have 12 entries: one for each count from 0 to 10,
/// // plus an additional one for counts greater than 10.
/// ```
pub fn count_histogram(kmer_counts: &FxHashMap<u64, u64>, histo_max: &u64) -> Vec<u64> {
    // Create a histogram of counts
    let mut histo: Vec<u64> = vec![0; *histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for count in kmer_counts.values() {
        if *count <= *histo_max {
            histo[*count as usize] += 1;
        } else {
            histo[*histo_max as usize + 1] += 1;
        }
    }
    histo
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
            _ => panic!("Invalid base"),
        };
        seq.push(base);
    }
    seq
}

// Returns the count of a specified kmer, or 0 if the kmer is not present.
// Uses canonical kmer
pub fn get_kmer_count(kmer_counts: &FxHashMap<u64, u64>, kmer: &u64, k: &usize) -> u64 {
    let revcomp = revcomp_kmer(kmer, k);
    let canonical = if *kmer < revcomp { *kmer } else { revcomp };
    match kmer_counts.get(&canonical) {
        Some(count) => *count,
        None => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tests() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_seq_to_reads() {
        // 'A' => 0, // 00
        // 'C' => 1, // 01
        // 'G' => 2, // 10
        // 'T' => 3, // 11
        let seq = "CGTAATGCGGCGA";
        let expected = vec![Read::new(vec![0b01101100, 0b00111001, 0b10100110, 0b00000000], 13)];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_seq_to_reads_n() {
        // 'A' => 0, // 00
        // 'C' => 1, // 01
        // 'G' => 2, // 10
        // 'T' => 3, // 11
        let seq = "CGTANATGCGGCGA";
        let expected = vec![Read::new(vec![0b01101100], 4), Read::new(vec![0b00111001, 0b10100110, 0b00000000], 9)];
        let actual = seq_to_reads(seq);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_revcomp_kmer() {
        let kmer = 0b0010_0110;
        let actual = revcomp_kmer(&kmer, &3usize);

        // Check that the reverse complement of the reverse complement is the original
        assert_eq!(kmer, revcomp_kmer(&actual, &3usize));

        // Check against hard coded expected value
        let expected = 0b0001_1001;
        assert_eq!(actual, expected);

        let kmer = 0b0110_1100_0011_1001_1010_0110;
        let actual = revcomp_kmer(&kmer, &12usize);
        assert_eq!(kmer, revcomp_kmer(&actual, &12usize));

        let expected = 0b0110_0101_1001_0011_1100_0110;

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_get_kmers() {
        let ints = vec![0b01101100, 0b00111001, 0b10100110];
        let read = Read::new(ints, 12);
        // original	                // reverse complement
        // 0b01101100_00111001_10   0b01_10010011_11000110  >
        // 0b101100_00111001_1010   0b0101_10010011_110001  >
        // 0b1100_00111001_101001   0b100101_10010011_1100  >
        // 0b00_00111001_10100110   0b01100101_10010011_11  <
        //
        let expected = vec![
            0b01_1001_0011_1100_0110,
            0b01_0110_0100_1111_0001,
            0b10_0101_1001_0011_1100,
            0b00_0011_1001_1010_0110, 
        ];
        let actual = read.get_kmers(&9usize);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_kmer_to_seq() {
        let kmer = 0b1001_1000;
        let seq = kmer_to_seq(&kmer, &4usize);
        assert_eq!(seq, "GCGA");

        let kmer = 0b1001_1000_1001_1000;
        let seq = crate::kmer::kmer_to_seq(&kmer, &8usize);
        assert_eq!(seq, "GCGAGCGA");
    }
}
