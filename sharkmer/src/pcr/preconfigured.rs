use super::PCRParams;

#[allow(dead_code)]
struct PCRPanel {
    name: String,
    description: String,
    params: Vec<PCRParams>,
}

fn get_preconfigured_panels() -> Vec<PCRPanel> {
<<<<<<< Updated upstream
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
=======
    let mut panels = vec![
        PCRPanel {
            name: "cnidaria".to_string(),
            description: "18S, 28S, ITS, 16S, and CO1".to_string(),
            params: get_cnidaria(),
        },
        PCRPanel {
            name: "insecta".to_string(),
            description: "Three primer sets".to_string(),
            params: get_insecta(),
        },
    ];

  // Loop over the panels and prepend the PCRPanel name to the PCRParams genename, using _ as a deimiter
    for panel in panels.iter_mut() {
        for param in panel.params.iter_mut() {
        param.gene_name = format!("{}_{}", panel.name, param.gene_name);
        }
    }
>>>>>>> Stashed changes

  panels
}

<<<<<<< Updated upstream
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
=======

fn get_insecta() -> Vec<PCRParams> {
    vec![

    // Universal insect mitochondrial primers from doi: 10.1111/1755-0998.12942
        PCRParams {
            forward_seq: "ACTWTGTTACGACTTDTY".to_string(),
            reverse_seq: "AGGATTAGATACCCTDBT".to_string(),
            min_length: 300,
            max_length: 500,
            gene_name: "12s".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "Hex12SF2–Hex12SR2. Product is about 391bp".to_string(),
        },

        PCRParams {
            forward_seq: "TARTYCAACATCGRGGTC".to_string(),
            reverse_seq: "CYGTRCDAAGGTAGCATA".to_string(),
            min_length: 300,
            max_length: 500,
            gene_name: "16s".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "Chiar16SF–Chiar16SR. Product is about 348bp".to_string(),
        },

        PCRParams {
            forward_seq: "HCCHGAYATRGCHTTYCC".to_string(),
            reverse_seq: "TATDGTRATDGCHCCNGC".to_string(),
            min_length: 300,
            max_length: 500,
            gene_name: "CO1".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Modified from Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "HexCOIF4–HexCOIR4. Product is about 322bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "GGNCRHCARTGRTAYTGA".to_string(),
            reverse_seq: "RATYTCDGARCAYTGNCC".to_string(),
            min_length: 200,
            max_length: 400,
            gene_name: "CO2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "	HexCOX2F3–HexCOX2R3. Product is about 260".to_string(),
        },
        
        PCRParams {
            forward_seq: "NCAAATRTCNTTHTGRGG".to_string(),
            reverse_seq: "YCAYTCDGGYTKRATRTG".to_string(),
            min_length: 300,
            max_length: 500,
            gene_name: "CytB".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "HexCytBF3–HexCytBR3. Product is about 373".to_string(),
        },
        
        PCRParams {
            forward_seq: "ATHARYTTATCRTANCGR".to_string(),
            reverse_seq: "NTTYGAYTTTKCDGARGG".to_string(),
            min_length: 150,
            max_length: 300,
            gene_name: "ND1".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "HexND1F4–HexND1R4. Product is about 210bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "HGGDGCYTCNACATGDGC".to_string(),
            reverse_seq: "RGGNTAYCARCCDGARCG".to_string(),
            min_length: 150,
            max_length: 300,
            gene_name: "ND4".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "HexND4F4–HexND4R4. Product is about 211bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "RTCYYTNGARTAAAAHCC".to_string(),
            reverse_seq: "NGCHAAYTWTGARTWTGA".to_string(),
            min_length: 200,
            max_length: 400,
            gene_name: "ND5".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Marquina et al doi: 10.1111/1755-0998.12942".to_string(),
            notes: "HexND5F3–HexND5R3. Product is about 271bp".to_string(),
        },
    // 
    
    
    // Drosophila-specific mitochondrial primers from https://doi.org/10.1016/j.ympev.2010.11.022
    // Provide longer fragments than universal primers above
        PCRParams {
            forward_seq: "CCGGTTTGAACTCAGATCACGT".to_string(),
            reverse_seq: "CGCCTGTTTAACAAAAACAT".to_string(),
            min_length: 400,
            max_length: 600,
            gene_name: "16s_v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "O'Grady et al https://doi.org/10.1016/j.ympev.2010.11.022 ".to_string(),
            notes: "HexND5F3–HexND5R3. Product is about 515bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "CAACATTTATTTTGATTTTTTGG".to_string(),
            reverse_seq: "TYCATTGCACTAATCTGCCATATTAG".to_string(),
            min_length: 700,
            max_length: 950,
            gene_name: "CO1_v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "O'Grady et al https://doi.org/10.1016/j.ympev.2010.11.022 ".to_string(),
            notes: "Product is about 830bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "ATGGCAGATTAGTGCAATGG".to_string(),
            reverse_seq: "GTTTAAGAGACCAGTACTTG".to_string(),
            min_length: 600,
            max_length: 850,
            gene_name: "CO2_v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "O'Grady et al https://doi.org/10.1016/j.ympev.2010.11.022 ".to_string(),
            notes: "Product is about 765bp".to_string(),
        },

        PCRParams {
            forward_seq: "AGCTATTGGGTTCAGACCCC".to_string(),
            reverse_seq: "GAAGTTTGGTTTAAACCTCC".to_string(),
            min_length: 400,
            max_length: 600,
            gene_name: "NADH".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "O'Grady et al https://doi.org/10.1016/j.ympev.2010.11.022 ".to_string(),
            notes: "Product is about 520bp".to_string(),
        },
    //

    // Drosophila-specific nuclear primers from https://doi.org/10.1016/j.ympev.2010.11.022
        PCRParams {
            forward_seq: "GCTTWTGAGACCGCTGATGG".to_string(),
            reverse_seq: "ATCTTRTCGAGACGCTGGAA".to_string(),
            min_length: 700,
            max_length: 1000,
            gene_name: "EF1g".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Magnacca et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Product is about 856bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "GCGTCTTTCTATTGCGCTACTAT".to_string(),
            reverse_seq: "GCTTGTACGGACTGCTGATTATT".to_string(),
            min_length: 800,
            max_length: 1100,
            gene_name: "Fz4".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Magnacca et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Product is about 943bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "CCCGACCTGGTTGAGGCTGCCAAGAATGC".to_string(),
            reverse_seq: "ACATATGCTCAGGGTGATTGCGTATGCA".to_string(),
            min_length: 900,
            max_length: 1200,
            gene_name: "Gpdh".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Magnacca et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Product is about 1071bp".to_string(),
        },

        PCRParams {
            forward_seq: "GCCATGTTCTSYGGMCAGCAYAT".to_string(),
            reverse_seq: "TAACGACCTCCNACCCARTCCCA".to_string(),
            min_length: 500,
            max_length: 800,
            gene_name: "Pgi".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Magnacca et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Product is about 633bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "CAGCAGCGTTACAATCTCCAGCC".to_string(),
            reverse_seq: "CCGAAGGGGCTCTTGGAGTTCAC".to_string(),
            min_length: 600,
            max_length: 900,
            gene_name: "Yp2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Magnacca et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Product is about 757bp".to_string(),
        },
    //    

    // Universal metazoan nuclear ribosomal primers from https://doi.org/10.1371/journal.pone.0134314
        PCRParams {
            forward_seq: "CTGGTGCCAGCAGCCGCGGYAA".to_string(),
            reverse_seq: "TCCGTCAATTYCTTTAAGTT".to_string(),
            min_length: 600,
            max_length: 900,
            gene_name: "18s".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Machida et al https://doi.org/10.1371/journal.pone.0134314".to_string(),
            notes: "Product is about 757bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "GGGAAAGAAGACCCTGTTGAG".to_string(),
            reverse_seq: "GCTTGGCBGCCACAAGCCAGTTA".to_string(),
            min_length: 600,
            max_length: 900,
            gene_name: "28s".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Machida et al https://doi.org/10.1371/journal.pone.0134314".to_string(),
            notes: "Product is about 757bp".to_string(),
        },
        
        PCRParams {
            forward_seq: "AACTTAAAGRAATTGACGGA".to_string(),
            reverse_seq: "CTCAACAGGGTCTTCTTTCCC".to_string(),
            min_length: 3500,
            max_length: 5500,
            gene_name: "ITS".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Machida et al https://doi.org/10.1016/j.ympev.2015.06.014".to_string(),
            notes: "Forward primer is reverse complement of 18s reverse primer. The reverse primer is the reverse complement of the 28s forward primer. Product is about 4800bp.".to_string(),
        },
    //
    
    // Ribosomal primers from Cnidarian panel, may work here as well  
        PCRParams {
            forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),
            reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
            min_length: 1600,
            max_length: 2000,
            gene_name: "18s_v2".to_string(),
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
            gene_name: "28s_v2".to_string(),
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
            gene_name: "ITS_v2".to_string(),
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
            gene_name: "ITS_v3".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
            notes: "Forward primer is reverse complement of 18s reverse primer. The reverse primer is the reverse complement of the 28s forward primer. Product is about 620-730bp.".to_string(),
        },
    ]
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
>>>>>>> Stashed changes
  ]
}

pub fn get_panel(panel_name: &str) -> Result<Vec<PCRParams>, String> {
<<<<<<< Updated upstream
  let panels = get_preconfigured_panels();
  panels
      .iter()
      .find(|panel| panel.name == panel_name)
      .map(|panel| panel.params.clone())
      .ok_or_else(|| format!("Invalid preconfigured PCR panel name '{}'", panel_name))
=======
    let panels = get_preconfigured_panels();
    panels
        .iter()
        .find(|panel| panel.name == panel_name)
        .map(|panel| panel.params.clone())
        .ok_or_else(|| format!("Invalid preconfigured PCR panel name '{}'", panel_name))
>>>>>>> Stashed changes
}
