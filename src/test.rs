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