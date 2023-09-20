//! This module provides kmer functions.
//!


use rustc_hash::FxHashMap;

// Create a structure with a hashmap for kmer counts and a u64 for the number of singleton kmers
pub struct KmerSummary {
    pub kmer_counts: FxHashMap<u64, u64>,
    pub n_singletons: u64,
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

/// Converts a DNA sequence into a vector of integer representations.
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
/// and each sub-sequence is encoded separately. The result is a vector of vectors,
/// with each inner vector holding the integer representation of a sub-sequence.
///
/// # Arguments
///
/// * `seq` - A reference to the DNA sequence string to be encoded.
///
/// # Returns
///
/// A vector of 8-bit integer vectors, where each inner vector represents an encoded
/// DNA sub-sequence. Sub-sequences are separated by any occurrence of the base `N`.
///
/// # Panics
///
/// The function will exit prematurely if it encounters a base other than `A`, `C`, `G`, `T`, or `N`.
///
/// # Example
///
/// ```rust
/// let sequence = "ACGTNAGCT";
/// let encoded = seq_to_ints(&sequence);
/// // Here, encoded would be a vector containing two inner vectors.
/// // The first inner vector would represent "ACGT" and the second would represent "AGCT".
/// ```
pub fn seq_to_ints(seq: &str) -> Vec<Vec<u8>> {
    let mut result: Vec<Vec<u8>> = Vec::new();
    let mut ints: Vec<u8> = Vec::with_capacity(seq.len() / 4);
    let mut frame: u8 = 0;
    let mut position: usize = 0; // position in the sequence. not including Ns
    for c in seq.chars() {
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            'N' => {
                if !ints.is_empty() {
                    result.push(ints);
                    ints = Vec::with_capacity(seq.len() / 4);
                    frame = 0; // Reset frame before starting a new subread
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
    if !ints.is_empty() || result.is_empty() {
        // Don't miss the last part
        result.push(ints);
    }
    result
}

/// Converts a vector of 8-bit integers into a vector of canonical kmers.
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
/// let canonical_kmers = ints_to_kmers(&encoded_sequence, &kmer_length);
/// ```
pub fn ints_to_kmers(ints: &Vec<u8>, k: &usize) -> Vec<u64> {
    let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k ) + 1);
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

#[cfg(test)]
mod tests {
    use super::*;
	
	#[test]
	fn test_tests() {
		assert_eq!(2 + 2, 4);
	}

    #[test]
	fn test_seq_to_ints() {
		// 'A' => 0, // 00
		// 'C' => 1, // 01
		// 'G' => 2, // 10
		// 'T' => 3, // 11
		let seq = "CGTAATGCGGCGA";
		let expected = vec![vec![0b01101100, 0b00111001, 0b10100110]];
		let actual = seq_to_ints(seq);
		assert_eq!(actual, expected);
	}

	#[test]
	fn test_seq_to_ints_n() {
		// 'A' => 0, // 00
		// 'C' => 1, // 01
		// 'G' => 2, // 10
		// 'T' => 3, // 11
		let seq = "CGTANATGCGGCGA";
		let expected = vec![vec![0b01101100], vec![0b00111001, 0b10100110]];
		let actual = seq_to_ints(seq);
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
	fn test_ints_to_kmers() {
		let ints = vec![0b01101100, 0b00111001, 0b10100110];
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
		let actual = ints_to_kmers(&ints, &9usize);
		assert_eq!(actual, expected);
	}

    #[test]
	fn test_kmer_to_seq(){
		let kmer = 0b1001_1000;
		let seq = kmer_to_seq(&kmer, &4usize);
		assert_eq!(seq, "GCGA");

		let kmer = 0b1001_1000_1001_1000;
		let seq = crate::kmer::kmer_to_seq(&kmer, &8usize);
		assert_eq!(seq, "GCGAGCGA");
	}

}