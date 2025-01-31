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
          min_length: 500,
          max_length: 700,
          gene_name: "16s".to_string(),
          min_coverage: 2,
          mismatches: 2,
          trim: 15,
          citation: "Modified from Cunningham and Buss 1993 https://doi.org/10.1016/0305-1978(93)90009-G ".to_string(),
          notes: "Amplifies portions of domains IV and V. Product is about 580bp".to_string(),
    },

    PCRParams {
          forward_seq: "TCATAARGATATHGG".to_string(),
          reverse_seq: "RTGNCCAAAAAACCA".to_string(),
          min_length: 600,
          max_length: 800,
          gene_name: "co1".to_string(),
          min_coverage: 2,
          mismatches: 2,
          trim: 15,
          citation: "Modified from Folmer et al. 1994".to_string(),
          notes: "Product is about 700bp long".to_string(),
    },	

    PCRParams {
          forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),
          reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
          min_length: 1600,
          max_length: 2000,
          gene_name: "18s".to_string(),
          min_coverage: 2,
          mismatches: 2,
          trim: 15,
          citation: "Medlin et al. 1988 https://doi.org/10.1016/0378-1119(88)90066-2".to_string(),
          notes: "Product is about 1790bp".to_string(),
    },	

    PCRParams {
          forward_seq: "CCYYAGTAACGGCGAGT".to_string(),
          reverse_seq: "SWACAGATGGTAGCTTCG".to_string(),
          min_length: 2900,
          max_length: 3500,
          gene_name: "28s".to_string(),
          min_coverage: 2,
          mismatches: 2,
          trim: 15,
          citation: "Evans et al. 2008. https://doi.org/10.1186/1471-2148-8-139".to_string(),
          notes: "Forward primer is F97, reverse is R3238. Product is about 3230bp.".to_string(),
    },	

    PCRParams {
          forward_seq: "TACACACCGCCCGTCGCTACTA".to_string(), // CAS18sF1
          reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
          min_length: 600,
          max_length: 1000,
          gene_name: "ITS".to_string(),
          min_coverage: 2,
          mismatches: 2,
          trim: 15,
          citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
          notes: "Forward primer is CAS18sF1, position 1843, from Ji et al 2003. The reverse primer is the reverse complement of the 28s forward primer. Product is about 770-880bp.".to_string(),
    },	
      
    PCRParams {
        forward_seq: "GTAGGTGAACCTGCAGAAGGATCA".to_string(), // Reverse complement of 18s reverse primer
        reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
        min_length: 600,
        max_length: 1000,
        gene_name: "ITSalt".to_string(),
        min_coverage: 2,
        mismatches: 2,
        trim: 15,
        citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
        notes: "Forward primer is reverse complement of 18s reverse primer. The reverse primer is the reverse complement of the 28s forward primer. Product is about 620-730bp.".to_string(),
    },

    PCRParams {
        // forward_seq: "ACGTGGTATGGTTGCCTCTG".to_string(),
        // reverse_seq: "CTTGATAACGCCAACGGCWAC".to_string(),
        // Custom primers based on alignment of a few cnidarian sequences
        forward_seq: "AMGWGGHATGGTDGCTGGTG".to_string(),
        reverse_seq: "YTTRATNAYDCCAACAGCWAC".to_string(),
        min_length: 200,
        max_length: 3000, // Exprected 350, but could have introns
        gene_name: "ef1a".to_string(),
        min_coverage: 2,
        mismatches: 2,
        trim: 15,
        citation: "".to_string(),
        notes: "".to_string(),
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
