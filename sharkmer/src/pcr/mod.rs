use crate::kmer::*;
use bio::io::fasta;
use std::io::Write;
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use petgraph::Graph;
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::algo::all_simple_paths;

// Create a structure to hold a kmer representing an oligo up to 32 nucleotides long in the 
// length*2 least significant bits
struct Oligo {
    length: usize,
    kmer: u64,
}

#[derive(PartialEq)]
enum PrimerDirection {
    Forward,
    Reverse,
}


// De Bruijn graph node
struct DBNode {
    sub_kmer: u64, // k-1 mer that contains overlap between kmers
    is_start: bool,  // Contains the forward primer
    is_end: bool,  // Contains the reverse complement of the reverse primer
    is_terminal: bool, // Is a terminal node
    visited: bool, // Has been visited during graph traversal
}

struct DBEdge {
    kmer: u64, // kmer that contains overlap between sub_kmers
    count: u64, // Number of times this kmer was observed
}

// Given an oligo as a String, return a Oligo struct representing it
fn string_to_oligo (seq: &str) -> Oligo {
    let mut kmer: u64 = 0;
    let mut length: usize = 0;
    for c in seq.chars() {
        let base = match c {
            'A' => 0, // 00
            'C' => 1, // 01
            'G' => 2, // 10
            'T' => 3, // 11
            _ => panic!("Invalid nucleotide {} in {}", c, seq),
        };

        kmer = (kmer << 2) | base as u64;
        length += 1;
    }
    Oligo {
        length,
        kmer,
    }
}

// Given a primer that may include ambiguous nucleotides, return a vector 
// of sequences that include all possible resolutions of the ambiguity. If 
// there are no ambiguous nucleotides, the vector contains only the original
// sequence.
fn resolve_primer(primer: String) -> Vec<String> {
    // Initial set of sequences to be expanded
    let mut sequences: Vec<String> = vec![primer.clone()];

    // For each nucleotide in the primer
    for (idx, nuc) in primer.chars().enumerate() {
        let possible_nucs = match nuc {
            // IUPAC nucleotide codes for ambiguous bases
            'R' => vec!["A".to_string(), "G".to_string()],
            'Y' => vec!["C".to_string(), "T".to_string()],
            'S' => vec!["G".to_string(), "C".to_string()],
            'W' => vec!["A".to_string(), "T".to_string()],
            'K' => vec!["G".to_string(), "T".to_string()],
            'M' => vec!["A".to_string(), "C".to_string()],
            'B' => vec!["C".to_string(), "G".to_string(), "T".to_string()],
            'D' => vec!["A".to_string(), "G".to_string(), "T".to_string()],
            'H' => vec!["A".to_string(), "C".to_string(), "T".to_string()],
            'V' => vec!["A".to_string(), "C".to_string(), "G".to_string()],
            'N' => vec!["A".to_string(), "C".to_string(), "G".to_string(), "T".to_string()],
            // Return the same nucleotide if it's not ambiguous
            _ => vec![nuc.to_string()],
        };

        let mut new_sequences = Vec::new();
        for seq in &sequences {
            for possible_nuc in &possible_nucs {
                let mut new_seq: Vec<char> = seq.chars().collect();
                new_seq[idx] = possible_nuc.chars().next().unwrap();
                new_sequences.push(new_seq.into_iter().collect());
            }
        }
        sequences = new_sequences;
    }

    sequences
}


fn permute_sequences(sequences: Vec<String>) -> Vec<String> {
    let mut unique_sequences = HashSet::new();

    for seq in sequences {
        unique_sequences.insert(seq.clone()); // Add original sequence

        // Do not permute the last nucleotide
        let last_index = seq.len() - 1;

        for (idx, nuc) in seq.chars().enumerate() {
            // Skip permutation for the last nucleotide
            if idx == last_index {
                continue;
            }

            for &replacement in ["A", "T", "C", "G"].iter() {
                if replacement != nuc.to_string() { // ensure we aren't replacing with the same nucleotide
                    let mut new_seq: Vec<char> = seq.chars().collect();
                    new_seq[idx] = replacement.chars().next().unwrap();
                    unique_sequences.insert(new_seq.into_iter().collect());
                }
            }
        }
    }

    unique_sequences.into_iter().collect()
}


fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        'Y' => 'R', // C or T
        'R' => 'Y', // A or G
        'S' => 'S', // C or G
        'W' => 'W', // A or T
        'K' => 'M', // G or T
        'M' => 'K', // A or C
        'B' => 'V', // C or G or T
        'V' => 'B', // A or C or G
        'D' => 'H', // A or G or T
        'H' => 'D', // A or C or T
        'N' => 'N', // A or C or G or T
        _ => panic!("Invalid nucleotide: {}", c),
    }).collect()
}


fn find_oligos_in_kmers (oligos: &[Oligo], kmers: &HashSet<u64>, k: &usize, dir: PrimerDirection) -> HashSet<u64>{
    // Find the kmers that contain the oligos.
    // If direction is forward, match the oligo at the start of the kmer.
    // If direction is reverse, match the oligo at the end of the kmer.

    // Assume all oligos have the same length
    let oligo_length = oligos[0].length;

    // Create hash set of oligos
    let oligo_set: HashSet<u64> = match dir {
        PrimerDirection::Forward => {
            // Rotate the oligo kmers to the start of the kmer
            oligos.iter().map(|oligo| oligo.kmer << (2*(*k-oligo_length))).collect()
        },
        PrimerDirection::Reverse => {
            // Just add the oligo kmers
            oligos.iter().map(|oligo| oligo.kmer).collect()
        },
    };    

    // Create mask for kmers so they can be compared to the oligo subsequence
    // If dir is forward, mask is set to 1 starting at the first position of the kmer for 2*oligo_length
    // If dir is reverse, mask is set to 1 for the 2*oligo_length least significant bits
    let mut mask: u64 = 0;
    match dir {
        PrimerDirection::Forward => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
            mask <<= 2 * *k - 2 * oligo_length;
        },
        PrimerDirection::Reverse => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
        },
    }

    // Create a mutable copy of kmers
    let mut kmers_match = kmers.clone();

    // Retain only the elements of kmers that intersect with oligo_set
    kmers_match.retain(|kmer| oligo_set.contains(&(kmer & mask)));

    kmers_match
}


fn n_nonterminal_nodes_in_graph(graph: &Graph<DBNode, DBEdge>) -> usize {
    let mut n_nonterminal_nodes = 0;
    for node in graph.node_indices() {
        if !graph[node].is_terminal {
            n_nonterminal_nodes += 1;
        }
    }
    n_nonterminal_nodes
}

fn n_unvisited_nodes_in_graph(graph: &Graph<DBNode, DBEdge>) -> usize {
    let mut n = 0;
    for node in graph.node_indices() {
        if !graph[node].visited {
            n += 1;
        }
    }
    n
}

fn get_path_length(graph: &Graph<DBNode, DBEdge>, new_node: NodeIndex ) -> usize {
    // Get the length of the path from the start node to the new node
    let mut path_length = 0;
    let mut current_node = new_node;
    loop {

        //current_node = NodeIndex::new(current_node_index).unwrap();
        let current_node_data = graph.node_weight(current_node).unwrap();
        if current_node_data.is_start {
            break;
        }

        path_length += 1;
        if let Some(next_node) = graph.neighbors_directed(current_node, Direction::Incoming).next() {
            current_node = next_node;
        } else {
            // Handle the case where there's no valid path to the start node from the given node
            // This could be a panic, a default value, or an error return, depending on your needs
            panic!("No path to the start node from the given node.");
        }
    }
    path_length
}

fn get_dbedge(kmer: &u64, kmer_counts: &FxHashMap<u64, u64>, k: &usize) -> DBEdge {
    
    let mut canonical_kmer = *kmer;

    let reverse_kmer = revcomp_kmer(kmer, k);
    if reverse_kmer < canonical_kmer {
        canonical_kmer = reverse_kmer;
    }
    
    DBEdge{
        kmer: *kmer,
        count: *kmer_counts.get(&canonical_kmer).unwrap(),
    }
}

pub fn do_pcr(
    kmer_counts: &FxHashMap<u64, u64>, 
    k: &usize, 
    max_length: &usize, 
    forward_seq: &str, 
    reverse_seq: &str, 
    run_name: &str,
    verbosity: usize,
) -> Vec<bio::io::fasta::Record> {

    // Create a vector to hold the records
    let mut records: Vec<fasta::Record> = Vec::new();

    // Preprocess the primers
    let mut forward = forward_seq.to_string();
    let mut reverse = reverse_seq.to_string();

    // Check if either is longer than k, if so retain only the last k nucleotides
    if forward.len() > *k {
        forward = forward[forward.len() - *k..].to_string();
        println!("Truncated the forward primer to {}", forward);
    }

    if reverse.len() > *k {
        reverse = reverse[reverse.len() - *k..].to_string();
        println!("Truncated the reverse primer to {}", reverse);
    }

    // Expand ambigous nucleotides
    let mut forward_variants = resolve_primer(forward);
    let mut reverse_variants = resolve_primer(reverse);

    // Get all possible variants of the primers
    forward_variants = permute_sequences(forward_variants);
    reverse_variants = permute_sequences(reverse_variants);

    // Replace the reverse variants with their reverse complements
    let mut reverse_variants_revcomp = Vec::new();
    for variant in reverse_variants.iter() {
        reverse_variants_revcomp.push(reverse_complement(variant));
    }
    reverse_variants = reverse_variants_revcomp;

    println!("There are {} variants of the forward primer", forward_variants.len());
    println!("There are {} variants of the reverse primer", reverse_variants.len());

    // Get the Oligos from the primer variants
    let mut forward_oligos: Vec<Oligo> = Vec::new();
    let mut reverse_oligos: Vec<Oligo> = Vec::new();
    for variant in forward_variants.iter() {
        forward_oligos.push(string_to_oligo(variant));
    }
    for variant in reverse_variants.iter() {
        reverse_oligos.push(string_to_oligo(variant));
    }
    

    // Create a hash set of the keys of kmer_counts_filtered
    print!("Creating hash set of kmers for assembly...");
    std::io::stdout().flush().unwrap();
    let mut kmers: std::collections::HashSet<u64> = kmer_counts.keys().copied().collect();
    
    // Add the reverse complement of each key to the hash set with revcomp_kmer()
    for kmer in kmer_counts.keys() {
        kmers.insert(revcomp_kmer(kmer, k));
    }

    println!(" done");

    // Find the kmers that contain the forward and reverse primers
    let start = std::time::Instant::now();
    print!("Finding kmers that contain the forward primer...");
    std::io::stdout().flush().unwrap();
    let forward_matches = find_oligos_in_kmers (&forward_oligos, &kmers, k, PrimerDirection::Forward);

    println!(" done, time: {:?}", start.elapsed());

    println!("  There are {} forward matches", forward_matches.len());
    for f in &forward_matches{
        let count = crate::kmer::get_kmer_count(kmer_counts, f, k);
        println!("  {}, count {}", crate::kmer::kmer_to_seq(f, k), count);
    }

    let start = std::time::Instant::now();
    print!("Finding kmers that contain the reverse primer...");
    std::io::stdout().flush().unwrap();
    let reverse_matches = find_oligos_in_kmers (&reverse_oligos, &kmers, k, PrimerDirection::Reverse);

    println!(" done, time: {:?}", start.elapsed());

    println!("  There are {} reverse matches", reverse_matches.len());
    for f in &reverse_matches{
        let count = crate::kmer::get_kmer_count(kmer_counts, f, k);
        println!("  {}, count {}", crate::kmer::kmer_to_seq(f, k), count);
    }

    // If the forward_matches or the reverse_matches are empty, exit
    if forward_matches.is_empty() | reverse_matches.is_empty() {
        println!("Binding sites were not found for both primers. Not searching for products.");
        return records;
    }

    println!("Creating graph, seeding with nodes that contain primer matches...");
    // Construct the graph
    let mut graph: Graph<DBNode, DBEdge> = Graph::new();

    // Create a suffix mask with 1 in the (2 * (k - 1)) least significant bits
    let suffix_mask: u64 = (1 << (2 * (*k - 1))) - 1;

    // Add the forward matches to the graph
    let mut start_nodes: Vec<NodeIndex> = Vec::new();
    for &kmer in &forward_matches {
        let prefix = kmer >> 2;

        // If the node with sub_kmer == suffix already exists, update the node so that is_start = true
        // Otherwise, create a new node

        let mut node_exists = false;
        for node in graph.node_indices() {
            if graph[node].sub_kmer == prefix {
                graph[node].is_start = true;
                node_exists = true;
                break;
            }
        }

        if !node_exists {
            let new_node = graph.add_node(DBNode {
                sub_kmer: prefix,
                is_start: true,
                is_end: false,
                is_terminal: false,
                visited: false, // Will not be extended
            });
            start_nodes.push(new_node);
        }
    }
    
    // Add the reverse matches to the graph
    let mut end_nodes: Vec<NodeIndex> = Vec::new();
    for &kmer in &reverse_matches {
        let suffix = kmer & suffix_mask;

        // If the node with sub_kmer == suffix already exists, update the node so that is_end = true
        // Otherwise, create a new node

        let mut node_exists = false;
        for node in graph.node_indices() {
            if graph[node].sub_kmer == suffix {
                graph[node].is_end = true;
                graph[node].is_terminal = true;
                node_exists = true;
                break;
            }
        }

        if !node_exists {
            let new_node = graph.add_node(DBNode {
                sub_kmer: suffix,
                is_start: false,
                is_end: true,
                is_terminal: true,
                visited: true,
            });
            end_nodes.push(new_node);
        }
    }

    

    // Loop over the nodes and see if any of the nodes are start and end nodes
    // If so, print a warning
    for node in graph.node_indices() {
        if graph[node].is_start && graph[node].is_end {
            println!("Warning: node {} is both a start and end node", node.index());
        }
    }

    // Print a summary of how many start and end nodes there are
    let mut n_start_nodes = 0;
    let mut n_end_nodes = 0;
    for node in graph.node_indices() {
        if graph[node].is_start {
            n_start_nodes += 1;
        }
        if graph[node].is_end {
            n_end_nodes += 1;
        }
    }
    println!("There are {} start nodes", n_start_nodes);
    println!("There are {} end nodes", n_end_nodes);

    // Print the information for each node
    for node in graph.node_indices() {
        println!("Node {}:", node.index());
        println!("  sub_kmer: {}", crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k-1)));
        println!("  is_start: {}", graph[node].is_start);
        println!("  is_end: {}", graph[node].is_end);
        println!("  is_terminal: {}", graph[node].is_terminal);
    }

    let start = std::time::Instant::now();
    println!("Extending the assembly graph...");
    
    
    
    // Extend graph by adding edges and new nodes to non-terminal nodes.
    // Terminal nodes have any of the following properties:
    // - is_end = true
    // - kmers that extend the node have been searched for but not found
    // - The node has a path length longer than max_length-k+1 from a start node

    // Extension continues until all nodes have been visited

    // Extension proceeds by:
    // - For each non-terminal node, find the kmers that extend the node
    // - Get the suffix of each kmer that extends the node
    // - If a node with sub_kmer == suffix already exists, add an edge to the existing node
    // - Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node
    
    // where:
    // - The prefix of the kmer of the node is the sub_kmer of the parent node in the graph
    // - The suffix of the kmer of the node is the sub_kmer of the new node in the graph
    // - If a node with the sub_kmer already exists, add a new edge to the existing node
    
    while n_unvisited_nodes_in_graph(&graph) > 0 {
        // Iterate over the nodes
        for node in graph.node_indices() {
            if !(graph[node].visited) {
                // Get the suffix of the kmer of the node
                let sub_kmer = graph[node].sub_kmer;

                if verbosity > 1 {
                    print!("  {} sub_kmer being extended for node {}. ", crate::kmer::kmer_to_seq(&sub_kmer, &(*k-1)), node.index());
                    std::io::stdout().flush().unwrap();
                }
                
                // Get the kmers that could extend the node
                let mut candidate_kmers: HashSet<u64> = HashSet::new();

                for base in 0..4 {
                    let kmer = (sub_kmer << 2) | base;
                    candidate_kmers.insert(kmer);
                }

                // Retain only the candidate kmers that are in the set of kmers
                candidate_kmers.retain(|kmer| kmers.contains(kmer));

                if verbosity > 1 {
                    print!("There are {} candidate kmers for extension. ", candidate_kmers.len());
                    std::io::stdout().flush().unwrap();
                }

                // If there are no candidate kmers, the node is terminal
                if candidate_kmers.is_empty() {
                    if verbosity > 1 {
                        print!("Marking node as terminal because there are no candidates for extension. ");
                        std::io::stdout().flush().unwrap();
                    }
                    graph[node].is_terminal = true;
                    graph[node].visited = true;
                    continue;
                }

                // Add new nodes if needed, and new edges
                for kmer in candidate_kmers.iter() {
                    let suffix = kmer & suffix_mask;

                    // Check if the node extends by itself and mark it as terminal if it does
                    if suffix == sub_kmer {
                        graph[node].is_terminal = true;
                        graph[node].visited = true;
                        if verbosity > 1 {
                            print!("Node {} extends itself. Marking as terminal. ", node.index());
                            std::io::stdout().flush().unwrap();
                        }
                        break;
                    }

                    // If the node with sub_kmer == suffix already exists, add an edge to the existing node
                    // Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node

                    let mut node_exists = false;
                    for existing_node in graph.node_indices() {
                        if graph[existing_node].sub_kmer == suffix {
                            let edge = get_dbedge(kmer, &kmer_counts, k);
                            graph.add_edge(node, existing_node, edge);
                            node_exists = true;
                            break;
                        }
                    }

                    if !node_exists {
                        let new_node = graph.add_node(DBNode {
                            sub_kmer: suffix,
                            is_start: false,
                            is_end: false,
                            is_terminal: false,
                            visited: false,
                        });
                        let edge = get_dbedge(kmer, &kmer_counts, k);
                        graph.add_edge(node, new_node, edge);

                        if verbosity > 1 {
                            print!("Added sub_kmer {} for new node {}. ", crate::kmer::kmer_to_seq(&suffix, &(*k-1)), new_node.index());
                            std::io::stdout().flush().unwrap();
                        }

                        // Check if the new node is max_length-k+1 from a start node
                        // If so, mark the new node as terminal
                        let path_length = get_path_length(&graph, new_node);
                        if path_length >= *max_length - (*k) + 1 {
                            graph[new_node].is_terminal = true;
                            if verbosity > 1 {
                                print!("Marking new node {} as terminal because it exceeds max_length from start. ", new_node.index());
                                std::io::stdout().flush().unwrap();
                            }
                        }
                        
                    }
                }
                graph[node].visited = true;

                if verbosity > 1 {
                    println!("There are now {} unvisited and {} non-terminal nodes in the graph. ", n_unvisited_nodes_in_graph(&graph), n_nonterminal_nodes_in_graph(&graph));
                    std::io::stdout().flush().unwrap();
                }

                
            }
        }
    }

    println!("done.  Time to extend graph: {:?}", start.elapsed());

    let start = std::time::Instant::now();
    println!("Traversing the assembly graph...");

    // Get all paths from start nodes to terminal nodes
    let mut all_paths = Vec::new();

    for start in &start_nodes {
        for end in &end_nodes {
            let paths_for_this_pair = all_simple_paths::<Vec<NodeIndex>, &Graph<DBNode, DBEdge>>(
                &graph,
                *start,
                *end,
                1,
                Some(*max_length - (*k) + 1),
            );
            
            all_paths.extend(paths_for_this_pair);
        }
    }

    // For each path, get the sequence of the path
    let mut i = 0;
    for path in all_paths {
        let mut sequence = String::new();
        // The first time through the loop add the whole sequence, after that just add the last base
        for node in path.iter() {
            let node_data = graph.node_weight(*node).unwrap();
            let subread = crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(*k-1));
            if sequence.is_empty() {
                sequence = subread;
            } else {
                sequence = format!("{}{}", sequence, subread.chars().last().unwrap(), );
            }
        }
        println!("{}", sequence);
        let id = format!("{} product {} length {}", run_name, i, sequence.len());
        let record = fasta::Record::with_attrs(&id, None, &(sequence.as_bytes()));
        records.push(record);
        i += 1;
    }

    println!("done.  Time to traverse graph: {:?}", start.elapsed());

    records
}

#[cfg(test)]
mod tests {
    use super::*;

	#[test]
	fn test_resolve_primers() {
		// Check with no ambiguous nucleotides
		let seq1 = "CGTAATGCGGCGA".to_string();
		let mut expected1 = vec![seq1.clone()];
		let mut result1 = resolve_primer(seq1);
		result1.sort();
		expected1.sort();
		assert_eq!(result1, expected1);

		// Check with one ambiguous nucleotide
		let seq2 = "CGTAATGCGGCGN".to_string();
		let mut expected2 = vec![
			"CGTAATGCGGCGA".to_string(),
			"CGTAATGCGGCGC".to_string(),
			"CGTAATGCGGCGG".to_string(),
			"CGTAATGCGGCGT".to_string(),
		];
		let mut result2 = resolve_primer(seq2);
		result2.sort();
		expected2.sort();
		assert_eq!(result2, expected2);

		// Check with one ambiguous nucleotide
		let seq3 = "CGTAATRCGGCGA".to_string();
		let mut expected3 = vec![
			"CGTAATACGGCGA".to_string(),
			"CGTAATGCGGCGA".to_string(),
		];
		let mut result3 = resolve_primer(seq3);
		result3.sort();
		expected3.sort();
		assert_eq!(result3, expected3);

		// Check with two ambiguous nucleotides
		let seq4 = "CGTAATRCGGCGY".to_string();
		let mut expected4 = vec![
			"CGTAATACGGCGC".to_string(),
			"CGTAATGCGGCGC".to_string(),
			"CGTAATACGGCGT".to_string(),
			"CGTAATGCGGCGT".to_string(),
		];
		let mut result4 = resolve_primer(seq4);
		result4.sort();
		expected4.sort();
		assert_eq!(result4, expected4);


	}

	#[test]
	fn test_permute_sequences(){
		let seq1 = vec!["CGA".to_string()];
		let mut expected1 = vec![
			"CAA".to_string(),
			"CCA".to_string(),
			"CGA".to_string(),
			"CTA".to_string(),
			"AGA".to_string(),
			"GGA".to_string(),
			"TGA".to_string(),
		];
		let mut result1 = permute_sequences(seq1);
		result1.sort();
		println!("Permutations: {}", result1.join(", "));
		expected1.sort();
		assert_eq!(result1, expected1);

	}

	#[test]
	fn test_string_to_oligo(){
		let oligo = string_to_oligo("GCGA");
		assert_eq!(oligo.kmer, 0b1001_1000);
		assert_eq!(oligo.length, 4);
	}

}