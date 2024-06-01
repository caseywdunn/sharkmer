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

    PCRParams {
          forward_seq: "GRCTGTTTACCAAAAACATA".to_string(),
          reverse_seq: "AATTCAACATMGAGG".to_string(),
          max_length: 1000,
          gene_name: "16s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "Modified from Cunningham and Buss 1993 https://doi.org/10.1016/0305-1978(93)90009-G ".to_string(),
          notes: "Amplifies portions of domains IV and V".to_string(),
      },

    PCRParams {
          forward_seq: "TCATAARGATATHGG".to_string(),
          reverse_seq: "RTGNCCAAAAAACCA".to_string(),
          max_length: 2000,
          gene_name: "co1".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "Modified from Folmer et al. 1994".to_string(),
          notes: "".to_string(),
      },	

    PCRParams {
          forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),
          reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
          max_length: 2500,
          gene_name: "18s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "Medlin et al. 1988 https://doi.org/10.1016/0378-1119(88)90066-2".to_string(),
          notes: "".to_string(),
      },	

    PCRParams {
          forward_seq: "CCYYAGTAACGGCGAGT".to_string(),
          reverse_seq: "SWACAGATGGTAGCTTCG".to_string(),
          max_length: 4000,
          gene_name: "28s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "Evans et al. 2008. https://doi.org/10.1186/1471-2148-8-139".to_string(),
          notes: "Forward primer is F97, reverse is R3238".to_string(),
      },	

    PCRParams {
          forward_seq: "TACACACCGCCCGTCGCTACTA".to_string(), // CAS18sF1
          // forward_seq: "GTAGGTGAACCTGCAGAAGGATCA".to_string(), // Reverse complement of 18s reverse primer
          reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
          max_length: 1000,
          gene_name: "ITS".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
          notes: "Forward primer is CAS18sF1, position 1843, from Ji et al 2003. The reverse primer is the reverse complement of the 28s forward primer.".to_string(),
      },	
  ]
}

fn get_teleostei() -> Vec<PCRParams> {
  vec![

    PCRParams {
          forward_seq: "TCAACCAACCACAAAGACATTGGCAC".to_string(),
          reverse_seq: "TAGACTTCTGGGTGGCCAAAGAATCA".to_string(),
          max_length: 1000,
          gene_name: "co1".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "https://doi.org/10.1098/rstb.2005.1716".to_string(),
          notes: "expected product is 655bp".to_string(),
      },

    // https://doi.org/10.1371/journal.pone.0266720 evaluates 
    // multiple primer pairs below
    PCRParams {
          forward_seq: "ACTGGGATTAGATACCCC".to_string(),
          reverse_seq: "TAGAACAGGCTCCTCTAG".to_string(),
          max_length: 200,
          gene_name: "12s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "https://doi.org/10.1093/nar/gkr732".to_string(),
          notes: "expected product is 106bp".to_string(),
      },

      PCRParams {
          forward_seq: "GACCCTATGGAGCTTTAGAC".to_string(),
          reverse_seq: "CGCTGTTATCCCTADRGTAACT".to_string(),
          max_length: 300,
          gene_name: "16s".to_string(),
          coverage: 3,
          mismatches: 2,
          trim: 15,
          citation: "https://doi.org/10.1002/ece3.3123".to_string(),
          notes: "expected product is 219bp".to_string(),
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
