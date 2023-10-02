use crate::kmer::*;
use bio::io::fasta;
use petgraph::algo::all_simple_paths;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use petgraph::Graph;
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use std::io::Write;

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
    sub_kmer: u64,     // k-1 mer that contains overlap between kmers
    is_start: bool,    // Contains the forward primer
    is_end: bool,      // Contains the reverse complement of the reverse primer
    is_terminal: bool, // Is a terminal node
    visited: bool,     // Has been visited during graph traversal
}

struct DBEdge {
    _kmer: u64, // kmer that contains overlap between sub_kmers
    count: u64, // Number of times this kmer was observed
}

// Given an oligo as a String, return a Oligo struct representing it
fn string_to_oligo(seq: &str) -> Oligo {
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
    Oligo { length, kmer }
}

// Given a primer that may include ambiguous nucleotides, return a set
// of sequences that include all possible resolutions of the ambiguity. If
// there are no ambiguous nucleotides, the set contains only the original
// sequence.
fn resolve_primer(primer: String) -> HashSet<String> {
    // Add the original sequence to the set
    let mut sequences: HashSet<String> = HashSet::new();

    // For each nucleotide in the primer
    for nuc in primer.chars() {
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
            'N' => vec![
                "A".to_string(),
                "C".to_string(),
                "G".to_string(),
                "T".to_string(),
            ],
            // Return the same nucleotide if it's not ambiguous
            _ => vec![nuc.to_string()],
        };

        // If sequences is empty, add each of the possible nucleotides as its own sequence
        if sequences.is_empty() {
            for possible_nuc in &possible_nucs {
                sequences.insert(possible_nuc.to_string());
            }
        } else {
            // Otherwise, for each sequence in sequences, add a new sequence for each possible nucleotide
            let mut new_sequences: HashSet<String> = HashSet::new();
            for seq in &sequences {
                for possible_nuc in &possible_nucs {
                    let mut new_seq: Vec<char> = seq.chars().collect();
                    new_seq.push(possible_nuc.chars().next().unwrap());
                    new_sequences.insert(new_seq.into_iter().collect());
                }
            }
            // Replace the old shorter sequences with the new extended sequences
            sequences = new_sequences;
        }
    }

    sequences
}

/// Given a set of sequences, return a set of all sequences that differ
/// from each original sequence at up to r positions. Includes the original
/// sequences.
fn permute_sequences(sequences: HashSet<String>, r: &usize) -> HashSet<String> {
    let mut permutations = HashSet::new();

    for seq in &sequences {
        for positions in combinations(seq.len(), *r) {
            generate_recursive_permutations(seq, &positions, 0, &mut permutations);
        }
    }

    permutations
}

fn combinations(n: usize, r: usize) -> Vec<Vec<usize>> {
    if r == 0 {
        return vec![Vec::new()];
    }

    if n == r {
        return vec![(0..r).collect()];
    }

    let without_last = combinations(n - 1, r);
    let mut with_last = combinations(n - 1, r - 1);
    for item in &mut with_last {
        item.push(n - 1);
    }

    without_last.into_iter().chain(with_last).collect()
}

fn generate_recursive_permutations(
    seq: &str,
    positions: &Vec<usize>,
    current: usize,
    unique_sequences: &mut HashSet<String>,
) {
    let nucleotides = ["A", "T", "C", "G"];

    if current == positions.len() {
        unique_sequences.insert(seq.to_string());
        return;
    }

    let pos = positions[current];
    for &nucleotide in nucleotides.iter() {
        let mut new_seq = seq.chars().collect::<Vec<_>>();
        new_seq[pos] = nucleotide.chars().next().unwrap();
        generate_recursive_permutations(
            &new_seq.iter().collect::<String>(),
            positions,
            current + 1,
            unique_sequences,
        );
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
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
        })
        .collect()
}

fn find_oligos_in_kmers(
    oligos: &[Oligo],
    kmers: &HashSet<u64>,
    k: &usize,
    dir: PrimerDirection,
) -> HashSet<u64> {
    // Find the kmers that contain the oligos.
    // If direction is forward, match the oligo at the start of the kmer.
    // If direction is reverse, match the oligo at the end of the kmer.

    // Assume all oligos have the same length
    let oligo_length = oligos[0].length;

    // Create hash set of oligos
    let oligo_set: HashSet<u64> = match dir {
        PrimerDirection::Forward => {
            // Rotate the oligo kmers to the start of the kmer
            oligos
                .iter()
                .map(|oligo| oligo.kmer << (2 * (*k - oligo_length)))
                .collect()
        }
        PrimerDirection::Reverse => {
            // Just add the oligo kmers
            oligos.iter().map(|oligo| oligo.kmer).collect()
        }
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
        }
        PrimerDirection::Reverse => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
        }
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

fn get_path_length(graph: &Graph<DBNode, DBEdge>, new_node: NodeIndex) -> Option<usize> {
    // Get the length of the path from the start node to the new node
    let mut path_length = 0;
    let mut current_node = new_node;

    // Create an empty set to contain the visited nodes
    let mut visited_nodes: HashSet<NodeIndex> = HashSet::new();

    loop {
        // If the current node has already been visited, break
        if visited_nodes.contains(&current_node) {
            return None;
        }

        // Add the current node to the visited nodes
        visited_nodes.insert(current_node);

        //current_node = NodeIndex::new(current_node_index).unwrap();
        let current_node_data = graph.node_weight(current_node).unwrap();
        if current_node_data.is_start {
            break;
        }

        path_length += 1;
        if let Some(next_node) = graph
            .neighbors_directed(current_node, Direction::Incoming)
            .next()
        {
            current_node = next_node;
        } else {
            // Handle the case where there's no valid path to the start node from the given node
            // This could be a panic, a default value, or an error return, depending on your needs
            panic!("No path to the start node from the given node.");
        }
    }
    Some(path_length)
}

fn get_dbedge(kmer: &u64, kmer_counts: &FxHashMap<u64, u64>, k: &usize) -> DBEdge {
    DBEdge {
        _kmer: *kmer,
        count: crate::kmer::get_kmer_count(kmer_counts, kmer, k),
    }
}

pub struct PCRParams {
    pub forward_seq: String,
    pub reverse_seq: String,
    pub max_length: usize,
    pub gene_name: String,
    pub coverage: u64,
    pub mismatches: usize,
    pub trim: usize,
}

pub fn do_pcr(
    kmer_counts: &FxHashMap<u64, u64>,
    k: &usize,
    verbosity: usize,
    params: &PCRParams,
) -> Vec<bio::io::fasta::Record> {
    // Create a vector to hold the fasta records
    let mut records: Vec<fasta::Record> = Vec::new();

    // Preprocess the primers
    let mut forward = params.forward_seq.clone();
    let mut reverse = params.reverse_seq.clone();

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
    forward_variants = permute_sequences(forward_variants, &params.mismatches);
    reverse_variants = permute_sequences(reverse_variants, &params.mismatches);

    // Replace the reverse variants with their reverse complements
    let mut reverse_variants_revcomp = HashSet::new();
    for variant in reverse_variants.iter() {
        reverse_variants_revcomp.insert(reverse_complement(variant));
    }
    reverse_variants = reverse_variants_revcomp;

    println!(
        "There are {} variants of the forward primer",
        forward_variants.len()
    );
    println!(
        "There are {} variants of the reverse primer",
        reverse_variants.len()
    );

    // Get the Oligos from the primer variants
    let mut forward_oligos: Vec<Oligo> = Vec::new();
    let mut reverse_oligos: Vec<Oligo> = Vec::new();
    for variant in forward_variants.iter() {
        forward_oligos.push(string_to_oligo(variant));
    }
    for variant in reverse_variants.iter() {
        reverse_oligos.push(string_to_oligo(variant));
    }

    // Create a hash set of the keys of kmer_counts
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
    let forward_matches =
        find_oligos_in_kmers(&forward_oligos, &kmers, k, PrimerDirection::Forward);

    println!(" done, time: {:?}", start.elapsed());

    let mut max_forward_count: u64 = 0;
    println!("  There are {} forward matches", forward_matches.len());
    for f in &forward_matches {
        let count = crate::kmer::get_kmer_count(kmer_counts, f, k);
        if count > max_forward_count {
            max_forward_count = count;
        }
        println!("  {}, count {}", crate::kmer::kmer_to_seq(f, k), count);
    }

    let start = std::time::Instant::now();
    print!("Finding kmers that contain the reverse primer...");
    std::io::stdout().flush().unwrap();
    let reverse_matches =
        find_oligos_in_kmers(&reverse_oligos, &kmers, k, PrimerDirection::Reverse);

    println!(" done, time: {:?}", start.elapsed());

    let mut max_reverse_count: u64 = 0;
    println!("  There are {} reverse matches", reverse_matches.len());
    for f in &reverse_matches {
        let count = crate::kmer::get_kmer_count(kmer_counts, f, k);
        if count > max_reverse_count {
            max_reverse_count = count;
        }
        println!("  {}, count {}", crate::kmer::kmer_to_seq(f, k), count);
    }

    // If the forward_matches or the reverse_matches are empty, exit
    if forward_matches.is_empty() | reverse_matches.is_empty() {
        println!("Binding sites were not found for both primers. Not searching for products.");
        return records;
    }

    // If the count of kmers containing primers is significantly higher than min_count, apply a higher min coverage
    let mut min_count = max_reverse_count;
    if max_forward_count < max_reverse_count {
        min_count = max_forward_count;
    }

    let coverage_multiplier = 5;
    let new_coverage = min_count / coverage_multiplier;
    if min_count > coverage_multiplier * params.coverage {
        println!("The count of kmers containing primers have high coverage {} relative to the coverage threshold of {}.  Increasing min coverage to {}.", min_count, params.coverage, new_coverage);

        // Create a hash set of the keys of kmer_counts
        println!("  Updating hash set of kmers to include only those that exceed updated coverage threshold.");
        // Remove all members of kmers
        kmers.clear();

        // Add each kmer and its reverse complement to kmers if the kmer count is >= new_coverage
        for kmer in kmer_counts.keys() {
            if kmer_counts[kmer] >= new_coverage {
                kmers.insert(*kmer);
                kmers.insert(revcomp_kmer(kmer, k));
            }
        }
    }

    // Construct the graph
    println!("Creating graph, seeding with nodes that contain primer matches...");
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
            println!(
                "Warning: node {} is both a start and end node",
                node.index()
            );
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

    if verbosity > 1 {
        // Print the information for each node
        for node in graph.node_indices() {
            println!("Node {}:", node.index());
            println!(
                "  sub_kmer: {}",
                crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k - 1))
            );
            println!("  is_start: {}", graph[node].is_start);
            println!("  is_end: {}", graph[node].is_end);
            println!("  is_terminal: {}", graph[node].is_terminal);
        }
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
                    print!(
                        "  {} sub_kmer being extended for node {}. ",
                        crate::kmer::kmer_to_seq(&sub_kmer, &(*k - 1)),
                        node.index()
                    );
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
                    print!(
                        "There are {} candidate kmers for extension. ",
                        candidate_kmers.len()
                    );
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
                            print!(
                                "Node {} extends itself. Marking as terminal. ",
                                node.index()
                            );
                            std::io::stdout().flush().unwrap();
                        }
                        break;
                    }

                    // If the node with sub_kmer == suffix already exists, add an edge to the existing node
                    // Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node

                    let mut node_exists = false;
                    for existing_node in graph.node_indices() {
                        if graph[existing_node].sub_kmer == suffix {
                            let edge = get_dbedge(kmer, kmer_counts, k);
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
                        let edge = get_dbedge(kmer, kmer_counts, k);
                        let edge_count = edge.count;
                        graph.add_edge(node, new_node, edge);

                        if verbosity > 1 {
                            print!(
                                "Added sub_kmer {} for new node {} with edge kmer count {}. ",
                                crate::kmer::kmer_to_seq(&suffix, &(*k - 1)),
                                new_node.index(),
                                edge_count
                            );
                            std::io::stdout().flush().unwrap();
                        }

                        // Check if the new node is max_length-k+1 from a start node
                        // If so, mark the new node as terminal
                        let path_length = get_path_length(&graph, new_node);

                        // If the path length is None, the node is part of a cycle and is marked terminal.
                        // If the path length is Some, is marked terminal if the path length is >= max_length-k+1
                        if let Some(path_length) = path_length {
                            if verbosity > 1 {
                                print!("Path length is {}. ", path_length);
                                std::io::stdout().flush().unwrap();
                            }

                            if path_length > params.max_length - (*k) {
                                graph[new_node].is_terminal = true;
                                if verbosity > 1 {
                                    print!("Marking new node {} as terminal because it exceeds max_length from start. ", new_node.index());
                                    std::io::stdout().flush().unwrap();
                                }
                            }
                        } else {
                            graph[new_node].is_terminal = true;
                            if verbosity > 1 {
                                print!("Marking new node {} as terminal because it is part of a cycle. ", new_node.index());
                                std::io::stdout().flush().unwrap();
                            }
                        }
                    }
                }
                graph[node].visited = true;

                if verbosity > 1 {
                    println!(
                        "There are now {} unvisited and {} non-terminal nodes in the graph. ",
                        n_unvisited_nodes_in_graph(&graph),
                        n_nonterminal_nodes_in_graph(&graph)
                    );
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
                Some(params.max_length - (*k) + 1),
            );

            all_paths.extend(paths_for_this_pair);
        }
    }

    // For each path, get the sequence of the path
    for (i, path) in all_paths.into_iter().enumerate() {
        let mut sequence = String::new();
        // The first time through the loop add the whole sequence, after that just add the last base
        for node in path.iter() {
            let node_data = graph.node_weight(*node).unwrap();
            let subread = crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(*k - 1));
            if sequence.is_empty() {
                sequence = subread;
            } else {
                sequence = format!("{}{}", sequence, subread.chars().last().unwrap(),);
            }
        }
        println!("{}", sequence);
        let id = format!(
            "{} product {} length {}",
            params.gene_name,
            i,
            sequence.len()
        );
        let record = fasta::Record::with_attrs(&id, None, sequence.as_bytes());
        records.push(record);
    }

    println!("done.  Time to traverse graph: {:?}", start.elapsed());

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    // These functions are used in the tests
    fn factorial(n: usize) -> usize {
        let mut result = 1;
        for i in 2..=n {
            result *= i;
        }
        result
    }

    fn n_combinations(n: usize, r: usize) -> usize {
        // Given a sequence of length n and r sites that can be permuted,
        // there are (n! / (r!(n-r)!)) combinations of r sites in the sequence.
        factorial(n) / (factorial(r) * factorial(n - r))
    }

    fn expected_permutations(n: usize, r: usize) -> usize {
        // Given a sequence of length n and r sites that can be permuted,
        // there are (n! / (r!(n-r)!)) combinations of r sites in the sequence and
        // 4^r permutations of the r sites.
        // So there are (n! / (r!(n-r)!)) * 4^r permutations.
        // The original sequence is included in each combination, but we only
        // want to include it once. Se we subtract the number of combinations
        // and add 1.

        n_combinations(n, r) * 4_usize.pow(r as u32) - n_combinations(n, r) + 1
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(2), 2);
        assert_eq!(factorial(3), 6);
        assert_eq!(factorial(8), 40320);
    }

    #[test]
    fn test_n_combinations() {
        assert_eq!(n_combinations(3, 2), 3);
        assert_eq!(n_combinations(5, 5), 1);
    }

    #[test]
    fn test_combinations() {
        // Make sure there are the correct numbers of combinations
        assert_eq!(combinations(3, 2).len(), n_combinations(3, 2));
        assert_eq!(combinations(8, 5).len(), n_combinations(8, 5));

        // Make sure the first combination has r elements
        assert_eq!(combinations(8, 5)[0].len(), 5);

        // Make sure that the combinations are unique
        let combos = combinations(8, 5);
        let combos_set: HashSet<Vec<usize>> = HashSet::from_iter(combos.clone());
        assert_eq!(combos.len(), combos_set.len());
    }

    #[test]
    fn test_expected_permutations() {
        assert_eq!(expected_permutations(2, 1), 7);
    }

    #[test]
    fn test_resolve_primers() {
        // Check with no ambiguous nucleotides
        let seq1 = "CGTAATGCGGCGA".to_string();
        let mut expected1: HashSet<String> = HashSet::new();
        expected1.insert(seq1.clone());
        let result1 = resolve_primer(seq1);
        assert_eq!(result1, expected1);

        // Check with one ambiguous nucleotide
        let seq2 = "CGTAATGCGGCGN".to_string();
        let mut expected2: HashSet<String> = HashSet::new();
        expected2.insert("CGTAATGCGGCGA".to_string());
        expected2.insert("CGTAATGCGGCGC".to_string());
        expected2.insert("CGTAATGCGGCGG".to_string());
        expected2.insert("CGTAATGCGGCGT".to_string());

        let result2 = resolve_primer(seq2);
        assert_eq!(result2, expected2);

        // Check with one ambiguous nucleotide
        let seq3 = "CGTAATRCGGCGA".to_string();
        let mut expected3: HashSet<String> = HashSet::new();
        expected3.insert("CGTAATACGGCGA".to_string());
        expected3.insert("CGTAATGCGGCGA".to_string());
        let result3 = resolve_primer(seq3);
        assert_eq!(result3, expected3);

        // Check with two ambiguous nucleotides
        let seq4 = "CGTAATRCGGCGY".to_string();
        let mut expected4: HashSet<String> = HashSet::new();
        expected4.insert("CGTAATACGGCGC".to_string());
        expected4.insert("CGTAATGCGGCGC".to_string());
        expected4.insert("CGTAATACGGCGT".to_string());
        expected4.insert("CGTAATGCGGCGT".to_string());

        let result4 = resolve_primer(seq4);
        assert_eq!(result4, expected4);

        // Check when the first nucleotide is ambiguous
        let seq5 = "RCGTAATCGGCGA".to_string();
        let mut expected5: HashSet<String> = HashSet::new();
        expected5.insert("ACGTAATCGGCGA".to_string());
        expected5.insert("GCGTAATCGGCGA".to_string());
        let result5 = resolve_primer(seq5);
        assert_eq!(result5, expected5);
    }

    #[test]
    fn test_permute_sequences() {
        // Check specific permutations for a tiny example
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CG".to_string());

            let mut expected: HashSet<String> = HashSet::new();
            expected.insert("CA".to_string());
            expected.insert("CC".to_string());
            expected.insert("CG".to_string());
            expected.insert("CT".to_string());
            expected.insert("AG".to_string());
            expected.insert("GG".to_string());
            expected.insert("TG".to_string());

            let result = permute_sequences(seq, &1_usize);
            println!(
                "Permutations: {}",
                result.iter().cloned().collect::<Vec<_>>().join(", ")
            );
            assert_eq!(result, expected);
        }

        // Check sequence length 3 and r=3
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());

            let result = permute_sequences(seq, &3_usize);

            // All sites permuted, so there should be 4^n sequences
            assert_eq!(result.len(), 64);
        }

        // Construct all the permutations procedurally
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());
            let r: usize = 2;
            let n: usize = 3;

            let bases = vec!['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("C{}{}", bases[i], bases[j]));
                    expected.insert(format!("{}G{}", bases[i], bases[j]));
                    expected.insert(format!("{}{}T", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }

        // Procedural test where n=4 and r=2
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("ACGT".to_string());
            let r: usize = 2;
            let n: usize = 4;

            let bases = vec!['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("{}{}GT", bases[i], bases[j]));
                    expected.insert(format!("{}C{}T", bases[i], bases[j]));
                    expected.insert(format!("{}CG{}", bases[i], bases[j]));
                    expected.insert(format!("A{}{}T", bases[i], bases[j]));
                    expected.insert(format!("A{}G{}", bases[i], bases[j]));
                    expected.insert(format!("AC{}{}", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }
    }

    #[test]
    fn test_string_to_oligo() {
        let oligo = string_to_oligo("GCGA");
        assert_eq!(oligo.kmer, 0b1001_1000);
        assert_eq!(oligo.length, 4);
    }
}
