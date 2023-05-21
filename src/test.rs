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
fn test_revcomp_kmer(){

	let kmer = 0b00_100110;
	let expected = 0b00_011001;
	let actual = revcomp_kmer(kmer, 3);
	assert_eq!(actual, expected);


	let kmer = 0b01101100_00111001_10100110;
	let expected = 0b01100101_10010011_11000110;
	let actual = revcomp_kmer(kmer, 9);
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
	let expected = vec![0b01_10010011_11000110, 0b0101_10010011_110001, 0b100101_10010011_1100, 0b00_00111001_10100110];
	let actual = ints_to_kmers(ints, 9);
	assert_eq!(actual, expected);
}