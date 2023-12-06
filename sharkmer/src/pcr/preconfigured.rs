use super::PCRParams; // Assuming PCRParams is defined in mod.rs

pub const CNIDARIA: &[PCRParams] = &[
	// --pcr "GACTGTTTACCAAAAACATA,GACTGTTTACCAAAAACATA,1000,16s" \
	PCRParams {
        forward_seq: "GACTGTTTACCAAAAACATA".to_string(),
        reverse_seq: "GACTGTTTACCAAAAACATA".to_string(),
        max_length: 1000,
        gene_name: "16s".to_string(),
        coverage: 3,
        mismatches: 2,
        trim: 15,
    },
  	// --pcr "TCATAAAGATATTGG,ATGCCCGAAAAACCA,2000,co1" \
	PCRParams {
        forward_seq: "TCATAAAGATATTGG".to_string(),
        reverse_seq: "ATGCCCGAAAAACCA".to_string(),
        max_length: 2000,
        gene_name: "co1".to_string(),
        coverage: 3,
        mismatches: 2,
        trim: 15,
    },	
  	// --pcr "AACCTGGTTGATCCTGCCAGT,TGATCCTTCTGCAGGTTCACCTAC,2500,18s" \
	PCRParams {
        forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),
        reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
        max_length: 2500,
        gene_name: "18s".to_string(),
        coverage: 3,
        mismatches: 2,
        trim: 15,
    },	
  	// --pcr "CCYYAGTAACGGCGAGT,SWACAGATGGTAGCTTCG,4000,28s"  \
	PCRParams {
        forward_seq: "CCYYAGTAACGGCGAGT".to_string(),
        reverse_seq: "SWACAGATGGTAGCTTCG".to_string(),
        max_length: 4000,
        gene_name: "28s".to_string(),
        coverage: 3,
        mismatches: 2,
        trim: 15,
    },	
  	// --pcr "TACACACCGCCCGTCGCTACTA,ACTCGCCGTTACTRRGG,1000,ITSfull" \
	PCRParams {
        forward_seq: "TACACACCGCCCGTCGCTACTA".to_string(),
        reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
        max_length: 1000,
        gene_name: "ITS".to_string(),
        coverage: 3,
        mismatches: 2,
        trim: 15,
    },	
];

pub fn get_panel(panel: &str) -> Option<&[PCRParams]> {
	match panel {
		"cnidaria" => Some(CNIDARIA),
		_ => None,
	}
}
