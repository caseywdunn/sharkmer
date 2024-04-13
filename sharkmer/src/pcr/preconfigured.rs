use super::PCRParams;

#[allow(dead_code)]
struct PCRPanel {
    name: String,
    description: String,
    params: Vec<PCRParams>,
}

fn get_preconfigured_panels() -> Vec<PCRPanel> {
  let mut panels = vec![
    PCRPanel {
        name: "cnidaria".to_string(),
        description: "18S, 28S, ITS, 16S, and CO1".to_string(),
        params: get_cnidaria(),
    },
    PCRPanel {
      name: "teleostei".to_string(),
      description: "12S and CO1".to_string(),
      params: get_teleostei(),
  },
  ];

  // Loop over the panels and prepend the PCRPanel name to the PCRParams genename, using _ as a deimiter
  for panel in panels.iter_mut() {
    for param in panel.params.iter_mut() {
      param.gene_name = format!("{}_{}", panel.name, param.gene_name);
      }
  }

  panels
}

fn get_cnidaria() -> Vec<PCRParams> {
  vec![
    // --pcr "GACTGTTTACCAAAAACATA,AATTCAACATCGAGG,1000,16s" \
    PCRParams {
          forward_seq: "GACTGTTTACCAAAAACATA".to_string(),
          reverse_seq: "AATTCAACATCGAGG".to_string(),
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

fn get_teleostei() -> Vec<PCRParams> {
  vec![
    // https://doi.org/10.1098/rstb.2005.1716
    // expected product is 655bp
    PCRParams {
          forward_seq: "TCAACCAACCACAAAGACATTGGCAC".to_string(),
          reverse_seq: "TAGACTTCTGGGTGGCCAAAGAATCA".to_string(),
          max_length: 1000,
          gene_name: "co1".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
      },
    // https://doi.org/10.1371/journal.pone.0266720 evaluates 
    // multiple primer pairs below
    
    // https://doi.org/10.1093/nar/gkr732
    // expected product is 106bp
    PCRParams {
          forward_seq: "ACTGGGATTAGATACCCC".to_string(),
          reverse_seq: "TAGAACAGGCTCCTCTAG".to_string(),
          max_length: 200,
          gene_name: "12s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
      },
      // https://doi.org/10.1002/ece3.3123
      // expected product is 219bp
      PCRParams {
          forward_seq: "GACCCTATGGAGCTTTAGAC".to_string(),
          reverse_seq: "CGCTGTTATCCCTADRGTAACT".to_string(),
          max_length: 300,
          gene_name: "16s".to_string(),
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
