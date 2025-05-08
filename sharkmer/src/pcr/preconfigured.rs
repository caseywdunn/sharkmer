use super::PCRParams;
use super::pcrparams_string;

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
            description: "Universal and custom primers targeting mitochondrial and nuclear products.".to_string(),
            params: get_cnidaria(),
        },

        PCRPanel {
            name: "human".to_string(),
            description: "Human mitochondrial primers from Ramos et al. 2009. doi: 10.1002/elps.200800601".to_string(),
            params: get_human(),
        },

        PCRPanel {
            name: "teleostei".to_string(),
            description: "Universal fish primers for mitochondrial and nuclear products".to_string(),
            params: get_teleostei(),
        },

        PCRPanel {
            name: "angiospermae".to_string(),
            description: "Universal plant primers for chloroplast and ITS products".to_string(),
            params: get_angiospermae(),
        },

        PCRPanel {
            name: "insecta".to_string(),
            description: "Multiple primer sets, including universal and Drosophila specific, as well as for broad metazoan markers".to_string(),
            params: get_insecta(),
        },

        PCRPanel {
            name: "bacteria".to_string(),
            description: "Bacterial 16S metabarcoding primers".to_string(),
            params: get_bacteria(),
        },

        PCRPanel {
            name: "metazoa".to_string(),
            description: "Metazoa 18S metabarcoding primers".to_string(),
            params: get_metazoa(),
        },

    ];

    // Loop over the panels and prepend the PCRPanel name to the PCRParams genename, using _ as a deimiter
    for panel in panels.iter_mut() {
        for param in panel.params.iter_mut() {
        param.gene_name = format!("{}_{}", panel.name, param.gene_name);
        }
    }

    panels.sort_by_key(|panel| panel.name.clone());

    panels
}

fn get_cnidaria() -> Vec<PCRParams> {
    vec![

        PCRParams {
            forward_seq: "GRCTGTTTACCAAAAACATA".to_string(),
            reverse_seq: "AATTCAACATMGAGG".to_string(),
            min_length: 500,
            max_length: 700,
            gene_name: "16S".to_string(),
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
            gene_name: "CO1".to_string(),
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
            gene_name: "18S".to_string(),
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
            gene_name: "28S".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Evans et al. 2008. https://doi.org/10.1186/1471-2148-8-139".to_string(),
            notes: "Forward primer is F97, reverse is R3238. Product is about 3230bp.".to_string(),
        },	

        PCRParams {
            forward_seq: "TACACACCGCCCGTCGCTACTA".to_string(), // CAS18SF1
            reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
            min_length: 600,
            max_length: 1000,
            gene_name: "ITS".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
            notes: "Forward primer is CAS18SF1, position 1843, from Ji et al 2003. The reverse primer is the reverse complement of the 28S forward primer. Product is about 770-880bp.".to_string(),
        },	

        PCRParams {
            forward_seq: "GTAGGTGAACCTGCAGAAGGATCA".to_string(), // Reverse complement of 18S reverse primer
            reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
            min_length: 600,
            max_length: 1000,
            gene_name: "ITS-v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
            notes: "Forward primer is reverse complement of 18S reverse primer. The reverse primer is the reverse complement of the 28S forward primer. Product is about 620-730bp.".to_string(),
        },

        PCRParams {
            // forward_seq: "ACGTGGTATGGTTGCCTCTG".to_string(),
            // reverse_seq: "CTTGATAACGCCAACGGCWAC".to_string(),
            // Custom primers based on alignment of a few cnidarian sequences
            forward_seq: "AMGWGGHATGGTDGCTGGTG".to_string(),
            reverse_seq: "YTTRATNAYDCCAACAGCWAC".to_string(),
            min_length: 200,
            max_length: 3000, // Exprected 350, but could have introns
            gene_name: "EF1A".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "".to_string(),
            notes: "".to_string(),
        },	

        PCRParams {
            forward_seq: "CGTGAAACCGYTRRAAGGG".to_string(),
            reverse_seq: "TTGGTCCGTGTTTCAAGACG".to_string(),
            min_length: 300,
            max_length: 700, // Exprected 350, but could have introns
            gene_name: "28S-v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "McCartin et al 2024. https://doi.org/10.7717/peerj.18607".to_string(),
            notes: "Anth-28S-eDNA; Product is about 540bp".to_string(),
        },	

    ]
}

fn get_human() -> Vec<PCRParams> {
    vec![

    // Human mt primers
        PCRParams {
            gene_name: "mt1404-3947".to_string(),
            forward_seq: "ACTTAAGGGTCGAAGGTGGATT".to_string(),
            reverse_seq: "TCGATGTTGAAGCCTGAGACTA".to_string(),
            min_length: 2400,
            max_length: 2600,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt3734-6739".to_string(),
            forward_seq: "AAGTCACCCTAGCCATCATTCTA".to_string(),
            reverse_seq: "GATATCATAGCTCAGACCATACC".to_string(),
            min_length: 2800,
            max_length: 3200,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt6511-9220".to_string(),
            forward_seq: "CTGCTGGCATCACTATACTACTA".to_string(),
            reverse_seq: "GATTGGTGGGTCATTATGTGTTG".to_string(),
            min_length: 2500,
            max_length: 2900,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt8910-10648".to_string(),
            forward_seq: "CTTACCACAAGGCACACCTACA".to_string(),
            reverse_seq: "GGCACAATATTGGCTAAGAGGG".to_string(),
            min_length: 1600,
            max_length: 1900,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt10360-12226".to_string(),
            forward_seq: "GTCTGGCCTATGAGTGACTACA".to_string(),
            reverse_seq: "CAGTTCTTGTGAGCTTTCTCGG".to_string(),
            min_length: 1700,
            max_length: 2000,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt11977-13830".to_string(),
            forward_seq: "CTCCCTCTACATATTTACCACAAC".to_string(),
            reverse_seq: "AAGTCCTAGGAAAGTGACAGCGA".to_string(),
            min_length: 1700,
            max_length: 2000,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt13477-15349".to_string(),
            forward_seq: "GCAGGAATACCTTTCCTCACAG".to_string(),
            reverse_seq: "GTGCAAGAATAGGAGGTGGAGT".to_string(),
            min_length: 1600,
            max_length: 2000,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt14898-151".to_string(),
            forward_seq: "TAGCCATGCACTACTCACCAGA".to_string(),
            reverse_seq: "GGATGAGGCAGGAATCAAAGAC".to_string(),
            min_length: 1700,
            max_length: 2000,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

        PCRParams {
            gene_name: "mt16488-1677".to_string(),
            forward_seq: "CTGTATCCGACATCTGGTTCCT".to_string(),
            reverse_seq: "GTTTAGCTCAGAGCGGTCAAGT".to_string(),
            min_length: 1600,
            max_length: 1900,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ramos et al 2009. doi: 10.1002/elps.200800601".to_string(),
            notes: "mt region".to_string(),
        },

    ]
}

fn get_teleostei() -> Vec<PCRParams> {
    vec![

    // Universal fish primers
        PCRParams {
            gene_name: "18S".to_string(),
            forward_seq: "TAACATATGCTTGTCTCAAAG".to_string(),
            reverse_seq: "CCTGTATTGTTATTTTTCGTCAC".to_string(),
            min_length: 300,
            max_length: 700,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Karabanov et al 2014. https://doi.org/10.3390/w14030437".to_string(),
            notes: "if18S1".to_string(),
        },

        PCRParams {
            gene_name: "16S".to_string(),
            forward_seq: "CGCCTGTTTATCAAAAACAT".to_string(),
            reverse_seq: "CCGGTCTGAACTCAGATCACGT".to_string(),
            min_length: 400,
            max_length: 700,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ivanova et al 2007. doi: 10.1111/j.1471-8286.2007.01748.x".to_string(),
            notes: "16Sar".to_string(),
        },

        PCRParams {
            gene_name: "CO1".to_string(),
            forward_seq: "TCAACCAACCACAAAGACATTGGCAC".to_string(),
            reverse_seq: "TAGACTTCTGGGTGGCCAAAGAATCA".to_string(),
            min_length: 500,
            max_length: 800,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Naz et al 2023. doi: 10.2478/aoas-2022-0073".to_string(),
            notes: "Fish-R1".to_string(),
        },

        PCRParams {
            gene_name: "CytB".to_string(),
            forward_seq: "AAAGCTTCCATCCAACATCTCAGCATGATGAAA".to_string(),
            reverse_seq: "AAACTGCAGCCCCTCAGAATGATATTTGTCCTCA".to_string(),
            min_length: 300,
            max_length: 600,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Nugroho et al 2007. https://www.bioflux.com.ro/docs/2019.1074-1079.pdf".to_string(),
            notes: "L14841; H15149".to_string(),
        },

        PCRParams {
            gene_name: "12S".to_string(),
            forward_seq: "ACTGGGATTAGATACCCCACTATG".to_string(),
            reverse_seq: "GAGAGTGACGGGCGGTGT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Zhang et al 2020. doi: 10.1111/2041-210X.13485 ".to_string(),
            notes: "Ac12S".to_string(),
        },

    ]
}


fn get_angiospermae() -> Vec<PCRParams> {
    vec![

    // Universal flowering plant primers
        PCRParams {
            gene_name: "psbA-trnH".to_string(),
            forward_seq: "GTTATGCATGAACGTAATGCTC".to_string(),
            reverse_seq: "CGCGCATGGTGGATTCACAATCC".to_string(),
            min_length: 100,
            max_length: 700,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "psbA3'f; trnHf".to_string(),
        },

        PCRParams {
            gene_name: "rpl36-infA-rps8".to_string(),
            forward_seq: "CACAAATTTTACGAACGAAG".to_string(),
            reverse_seq: "TAATGACAGAYCGAGARGCTCGAC".to_string(),
            min_length: 400,
            max_length: 600,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "rpl36f; rps8r".to_string(),
        },

        PCRParams {
            gene_name: "trnK-rps16".to_string(),
            forward_seq: "TACTCTACCRTTGAGTTAGCAAC".to_string(),
            reverse_seq: "AAAGGKGCTCAACCTACARGAAC".to_string(),
            min_length: 500,
            max_length: 900,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "trnK5'r; rps16-4546mod".to_string(),
        },

        PCRParams {
            gene_name: "trnV-atpE".to_string(),
            forward_seq: "GTGTAAACGAGTTGCTCTACCA".to_string(),
            reverse_seq: "CGACATTTGCACATTTAGATGCTAC".to_string(),
            min_length: 400,
            max_length: 1000,
             min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "trnV5f; S1022".to_string(),
        },

        PCRParams {
            gene_name: "trnC-ycf6".to_string(),
            forward_seq: "CCAGTTCAAATCTGGGTGTC".to_string(),
            reverse_seq: "CCCAAGCAAGACTTACTATATCC".to_string(),
            min_length: 100,
            max_length: 1200,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "trnC; petN1r".to_string(),
        },

        PCRParams {
            gene_name: "ycf6-psbM".to_string(),
            forward_seq: "GGATATAGTAAGTCTTGCTTGGG".to_string(),
            reverse_seq: "TTCTTGCATTTATTGCTACTGC".to_string(),
            min_length: 300,
            max_length: 1600,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "petN1; psbM2r".to_string(),
        },

        PCRParams {
            gene_name: "psbM-trnD".to_string(),
            forward_seq: "GCGGTAGGAACTAGAATAAATAG".to_string(),
            reverse_seq: "GGGATTGTAGTTCAATTGGT".to_string(),
            min_length: 500,
            max_length: 1400,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "psbMA1; trnD".to_string(),
        },

        PCRParams {
            gene_name: "atpB-rbcL".to_string(),
            forward_seq: "AGAAGTAGTAGGATTGATTCTCATA".to_string(),
            reverse_seq: "GAATCCAACACTTGCTTTAGTCTCT".to_string(),
            min_length: 500,
            max_length: 1000,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "S2r; RBCL1".to_string(),
        },

        PCRParams {
            gene_name: "trnL-F".to_string(),
            forward_seq: "GGTTCAAGTCCCTCTATCCC".to_string(),
            reverse_seq: "ATTTGAACTGGTGACACGAG".to_string(),
            min_length: 200,
            max_length: 600,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "e; f".to_string(),
        },

        PCRParams {
            gene_name: "ITS".to_string(),
            forward_seq: "CCTTATCATTTAGAGGAAGGAG".to_string(),
            reverse_seq: "TCCTCCGCTTATTGATATGC".to_string(),
            min_length: 400,
            max_length: 900,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Kress et al 2005. https://doi.org/10.1073/pnas.050312310".to_string(),
            notes: "ITS5a; ITS4".to_string(),
        },

    ]
}


fn get_insecta() -> Vec<PCRParams> {
    vec![

    // Universal insect mitochondrial primers from doi: 10.1111/1755-0998.12942
        PCRParams {
            forward_seq: "ACTWTGTTACGACTTDTY".to_string(),
            reverse_seq: "AGGATTAGATACCCTDBT".to_string(),
            min_length: 300,
            max_length: 500,
            gene_name: "12S".to_string(),
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
            gene_name: "16S".to_string(),
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


    // Drosophila-specific mitochondrial primers doi: https://doi.org/10.1016/j.ympev.2010.11.022
    // Provide longer fragments than universal primers above
        PCRParams {
            forward_seq: "CCGGTTTGAACTCAGATCACGT".to_string(),
            reverse_seq: "CGCCTGTTTAACAAAAACAT".to_string(),
            min_length: 400,
            max_length: 600,
            gene_name: "16S-v2".to_string(),
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
            gene_name: "CO1-v2".to_string(),
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
            gene_name: "CO2-v2".to_string(),
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

    // Drosophila-specific nuclear primers doi: https://doi.org/10.1016/j.ympev.2010.11.022
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

    // Universal metazoan nuclear ribosomal primers doi: https://doi.org/10.1371/journal.pone.0134314
        PCRParams {
            forward_seq: "CTGGTGCCAGCAGCCGCGGYAA".to_string(),
            reverse_seq: "TCCGTCAATTYCTTTAAGTT".to_string(),
            min_length: 600,
            max_length: 900,
            gene_name: "18S".to_string(),
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
            gene_name: "28S".to_string(),
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
            notes: "Forward primer is reverse complement of 18S reverse primer. The reverse primer is the reverse complement of the 28S forward primer. Product is about 4800bp.".to_string(),
        },
    //

    // Nuclear ribosomal primers from cnidarian panel, to test if they work for insects as well
        PCRParams {
            forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),
            reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
            min_length: 1600,
            max_length: 2000,
            gene_name: "18S-v2".to_string(),
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
            gene_name: "28S-v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Evans et al. 2008. https://doi.org/10.1186/1471-2148-8-139".to_string(),
            notes: "Forward primer is F97, reverse is R3238. Product is about 3230bp.".to_string(),
        },	

        PCRParams {
            forward_seq: "TACACACCGCCCGTCGCTACTA".to_string(), // CAS18SF1
            reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
            min_length: 600,
            max_length: 1000,
            gene_name: "ITS-v2".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
            notes: "Forward primer is CAS18SF1, position 1843, from Ji et al 2003. The reverse primer is the reverse complement of the 28S forward primer. Product is about 770-880bp.".to_string(),
        },	

        PCRParams {
            forward_seq: "GTAGGTGAACCTGCAGAAGGATCA".to_string(), // Reverse complement of 18S reverse primer
            reverse_seq: "ACTCGCCGTTACTRRGG".to_string(),
            min_length: 600,
            max_length: 1000,
            gene_name: "ITS_v3".to_string(),
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Ji et al 2003. https://doi.org/10.1046/j.1471-8286.2003.00519.x".to_string(),
            notes: "Forward primer is reverse complement of 18S reverse primer. The reverse primer is the reverse complement of the 28S forward primer. Product is about 620-730bp.".to_string(),
        },
    ]
}


fn get_metazoa() -> Vec<PCRParams> {
    vec![

    // Metabarcoding 18S primers
        PCRParams {
            gene_name: "18S-V3".to_string(),
            forward_seq: "AACGGCTACCACATCCAAGG".to_string(),
            reverse_seq: "CACCAGACTTGCCCTCCAAT".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Within V3".to_string(),
        },

        PCRParams {
            gene_name: "18S-V5-V7S".to_string(),
            forward_seq: "TGACGGAAGGGCACCACCAG".to_string(),
            reverse_seq: "TCCACCAACTAAGAACGGCC".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Between V5 and beginning of V7 (short amplicon)".to_string(),
        },

        PCRParams {
            gene_name: "18S-V5-V7L".to_string(),
            forward_seq: "AAACGATGCCGACTAGCGAT".to_string(),
            reverse_seq: "TCCACCAACTAAGAACGGCC".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Between V5 and beginning of V7 (long amplicon)".to_string(),
        },
        PCRParams {
            gene_name: "18S-V7".to_string(),
            forward_seq: "GGCCGTTCTTAGTTGGTGGA".to_string(),
            reverse_seq: "TGCGGCCCAGAACATCTAAG".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Within V7".to_string(),
        },

        PCRParams {
            gene_name: "18S-V7p-V8".to_string(),
            forward_seq: "AACAGGTCTGTGATGCCCTT".to_string(),
            reverse_seq: "TGTGTACAAAGGGCAGGGAC".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Part of V7 and most of V8".to_string(),
        },

        PCRParams {
            gene_name: "18S-V9".to_string(),
            forward_seq: "CTTTGTACACACCGCCCGTC".to_string(),
            reverse_seq: "CCTTGTTACGACTTTTACTTCCTCT".to_string(),
            min_length: 50,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Damian-Serrano et al. 2022 https://doi.org/10.1371/journal.pone.0267761".to_string(),
            notes: "Within V9".to_string(),
        },

        PCRParams {
            gene_name: "18S-V9-v2".to_string(),
            forward_seq: "CCCTGCCHTTTGTACACAC".to_string(),
            reverse_seq: "CCTTCYGCAGGTTCACCTAC".to_string(),
            min_length: 1,
            max_length: 300,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "From Amaral-Zettler 2009, https://doi.org/10.1371/journal.pone.0006372".to_string(),
            notes: "1380F;1510R".to_string(),
        },
    ]
}

fn get_bacteria() -> Vec<PCRParams> {
    vec![

    // Bacterial 16S variable region primers
        PCRParams {
            gene_name: "16S-27F-338R".to_string(),
            forward_seq: "AGAGTTTGATCCTGGCTCAG".to_string(),
            reverse_seq: "TGCTGCCTCCCGTAGGAGT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V1-V2".to_string(),
        },

        PCRParams {
            gene_name: "16S-V2f-V3r".to_string(),
            forward_seq: "AGTGGCGGACGGGTGAGTAA".to_string(),
            reverse_seq: "CCGCGGCTGCTGGCAC".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V2-V3".to_string(),
        },

        PCRParams {
            gene_name: "16S-341F-785R".to_string(),
            forward_seq: "CCTACGGGNGGCWGCAG".to_string(),
            reverse_seq: "GACTACHVGGGTATCTAATCC".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V3-V4".to_string(),
        },

        PCRParams {
            gene_name: "16S-PRK341F-PRK806R".to_string(),
            forward_seq: "CCTACGGGRBGCASCAG".to_string(),
            reverse_seq: "GGACTACYVGGGTATCTAAT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V3-V4-v2".to_string(),
        },

        PCRParams {
            gene_name: "16S-515F-806R".to_string(),
            forward_seq: "GTGCCAGCMGCCGCGGTAA".to_string(),
            reverse_seq: "GGACTACHVGGGTWTCTAAT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V4".to_string(),
        },

        PCRParams {
            gene_name: "16S-515F-806RB".to_string(),
            forward_seq: "GTGCCAGCMGCCGCGGTAA".to_string(),
            reverse_seq: "GGACTACNVGGGTWTCTAAT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V4".to_string(),
        },

        PCRParams {
            gene_name: "16S-515F-Y-926R".to_string(),
            forward_seq: "GTGYCAGCMGCCGCGGTAA".to_string(),
            reverse_seq: "CCGYCAATTYMTTTRAGTTT".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V4-V5".to_string(),
        },

        PCRParams {
            gene_name: "16S-B969F-BA1406R".to_string(),
            forward_seq: "ACGCGHNRAACCTTACC".to_string(),
            reverse_seq: "ACGGGCRGTGWGTRCAA".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al 2023. doi: 10.3389/fmars.2023.1199116".to_string(),
            notes: "V6-V8".to_string(),
        },

        PCRParams {
            gene_name: "16S-799F-1391R".to_string(),
            forward_seq: "AACMGGATTAGATACCCKG".to_string(),
            reverse_seq: "GACGGGCGGTGWGTRCA".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Chelius and Triplett, 2001, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V5-V6-V7".to_string(),
        },

        PCRParams {
            gene_name: "16S-967F-1391R".to_string(),
            forward_seq: "CAACGCGAAGAACCTTACC".to_string(),
            reverse_seq: "GACGGGCGGTGWGTRCA".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Sogin et al., 2006, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V6-V7".to_string(),
        },

        PCRParams {
            gene_name: "16S-799F-1193R".to_string(),
            forward_seq: "AACMGGATTAGATACCCKG".to_string(),
            reverse_seq: "ACGTCATCCCCACCTTCC".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Bodenhausen et al., 2013, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V5-V6-V7".to_string(),
        },

        PCRParams {
            gene_name: "16S-68F-783Rabc".to_string(),
            forward_seq: "TNANACATGCAAGTCGRRCG".to_string(),
            reverse_seq: "CTACCAGGGTATCTAATCCTG".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "McAllister et al., 2011, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V1-V4".to_string(),
        },

        PCRParams {
            gene_name: "16S-68F-518R".to_string(),
            forward_seq: "TNANACATGCAAGTCGRRCG".to_string(),
            reverse_seq: "WTTACCGCGGCTGCTGG".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Lee et al., 2010, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V1-V3".to_string(),
        },

        PCRParams {
            gene_name: "16S-341F-783Rabc".to_string(),
            forward_seq: "CCTACGGGNGGCWGCAG".to_string(),
            reverse_seq: "CTACCAGGGTATCTAATCCTG".to_string(),
            min_length: 200,
            max_length: 500,
            min_coverage: 2,
            mismatches: 2,
            trim: 15,
            citation: "Sakai et al., 2004, doi: http://dx.doi.org/10.3389/fmicb.2016.00650".to_string(),
            notes: "V3-V4".to_string(),
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

pub fn print_pcr_panels(){
    let panels = get_preconfigured_panels();
    for panel in panels {
        println!("{}: {}", panel.name, panel.description);
        for param in panel.params {
            println!("  {}", param.gene_name);
            println!("    forward: {}", param.forward_seq);
            println!("    reverse: {}", param.reverse_seq);
            println!("    min-length: {}", param.min_length);
            println!("    max-length: {}", param.max_length);
            println!("    min-coverage: {}", param.min_coverage);
            println!("    mismatches: {}", param.mismatches);
            println!("    trim: {}", param.trim);
            println!("    citation: {}", param.citation);
            println!("    notes: {}", param.notes);
            println!("    arguments: --pcr \"{}\"", pcrparams_string(&param));
            println!();
        }
    }
}
