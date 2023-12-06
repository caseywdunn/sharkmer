use super::PCRParams;

#[allow(dead_code)]
struct PCRPanel {
    name: String,
    description: String,
    params: Vec<PCRParams>,
}

fn get_preconfigured_panels() -> Vec<PCRPanel> {
  vec![
    PCRPanel {
        name: "cnidaria".to_string(),
        description: "18S, 28S, ITS, 16S, and CO1".to_string(),
        params: get_cnidaria(),
    },
  ]
}

fn get_cnidaria() -> Vec<PCRParams> {
  vec![
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
  ]
}

pub fn get_panel(panel_name: &str) -> Result<Vec<PCRParams>, String> {
  let panels = get_preconfigured_panels();
  panels
      .iter()
      .find(|panel| panel.name == panel_name)
      .map(|panel| panel.params.clone())
      .ok_or_else(|| format!("Invalid preconfigured PCR panel name '{}'", panel_name))
}
