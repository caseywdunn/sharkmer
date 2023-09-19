use rustc_hash::FxHashMap;

// Create a structure with a hashmap for kmer counts and a u64 for the number of singleton kmers
pub struct KmerSummary {
    pub kmer_counts: FxHashMap<u64, u64>,
    pub n_singletons: u64,
}

// Get the kmer that is a reverse complement of a specified kmer
pub fn revcomp_kmer(kmer: u64, k: u8) -> u64 {
    let mut revcomp = 0;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 3;
        revcomp = (revcomp << 2) | (3 - base);
    }
    revcomp
}

// Convert read to integer encoded subreads, split on N in original sequence
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

pub fn ints_to_kmers(ints: &Vec<u8>, k: u8) -> Vec<u64> {
    let mut kmers: Vec<u64> = Vec::with_capacity((ints.len() * 4 / k as usize) + 1);
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
            if n_valid >= k {
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

pub fn count_histogram(kmer_counts: &FxHashMap<u64, u64>, histo_max: u64) -> Vec<u64> {
    // Create a histogram of counts
    let mut histo: Vec<u64> = vec![0; histo_max as usize + 2]; // +2 to allow for 0 and for >histo_max
    for count in kmer_counts.values() {
        if *count <= histo_max {
            histo[*count as usize] += 1;
        } else {
            histo[histo_max as usize + 1] += 1;
        }
    }
    histo
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
		let actual = revcomp_kmer(kmer, 3);

		// Check that the reverse complement of the reverse complement is the original
		assert_eq!(kmer, revcomp_kmer(actual, 3));

		// Check against hard coded expected value
		let expected = 0b0001_1001;
		assert_eq!(actual, expected);

		let kmer = 0b0110_1100_0011_1001_1010_0110;
		let actual = revcomp_kmer(kmer, 12);
		assert_eq!(kmer, revcomp_kmer(actual, 12));

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
		let actual = ints_to_kmers(&ints, 9);
		assert_eq!(actual, expected);
	}

}