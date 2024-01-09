// pcr/mod.rs

use bio::alignment::distance::simd::*;
use bio::io::fasta;
use colored::*;
use petgraph::algo::{all_simple_paths, is_cyclic_directed};
use petgraph::dot::{Config, Dot};
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::{Bfs, EdgeRef};
use petgraph::Direction;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::fs::File;
use std::io::Write;
use textwrap::{fill, indent};

use crate::{COLOR_FAIL, kmer};
use crate::COLOR_NOTE;
use crate::COLOR_SUCCESS;
use crate::COLOR_WARNING;
use crate::kmer::KmerCounts;

pub mod preconfigured;

// Constants that may require tuning

/// The multiplier for establishing when a kmer is considered to have high coverage,
/// relative to the coverage threshold. It is then used to also adjust the threshold.
const COVERAGE_MULTIPLIER: u64 = 2;

/// A multiplier for adjusting the threshold as it is applied.
const COVERAGE_STEPS: u64 = 4;

/// The maximum number of kmers containing the forward or reverse primers to maintain,
/// with only those with the highest count being retained
const MAX_NUM_PRIMER_KMERS: usize = 10;

/// How frequently to check extension graph for ballooning growth
const EXTENSION_EVALUATION_FREQUENCY: usize = 1_000;

/// The graph depth over which to evaluate ballooning growth
const EXTENSION_EVALUATION_DEPTH: usize = 4;

/// Sets the threshold for detecting ballooning, where threshold is
/// 4^(EXTENSION_EVALUATION_DEPTH-EXTENSION_EVALUATION_DIFF)
const EXTENSION_EVALUATION_DIFF: usize = 1;

/// Edit threshold for https://docs.rs/bio/latest/bio/alignment/distance/simd/fn.bounded_levenshtein.html
const DISTANCE_EDIT_THRESHOLD: u32 = 10;

/// If an edges count has more than BALLOONING_COUNT_THRESHOLD_MULTIPLIER * median_edge_count, it is
/// likely going to balloon and is not added to the graph
const BALLOONING_COUNT_THRESHOLD_MULTIPLIER: f64 = 10.0;

/// Give up if the graph gets too large
const MAX_NUM_NODES: usize = 50_000;

struct AssemblyRecord {
    fasta_record: fasta::Record,
    kmer_min_count: u64,
}

// Create a structure to hold a kmer representing an oligo up to 32 nucleotides long in the
// length*2 least significant bits
pub struct Oligo {
    pub length: usize,
    pub kmer: u64,
}

#[derive(PartialEq)]
enum PrimerDirection {
    Forward,
    Reverse,
}

// De Bruijn graph node
#[derive(Debug, Clone)]
pub struct DBNode {
    pub sub_kmer: u64,     // k-1 mer that contains overlap between kmers
    pub is_start: bool,    // Contains the forward primer
    pub is_end: bool,      // Contains the reverse complement of the reverse primer
    pub is_terminal: bool, // Is a terminal node
    pub visited: bool,     // Has been visited during graph traversal
}

#[derive(Debug, Clone)]
pub struct DBEdge {
    pub _kmer: u64, // kmer that contains overlap between sub_kmers
    pub count: u64, // Number of times this kmer was observed
}

// Iteratively remove nodes that do not have outgoing edges and are not end nodes
// These are terminal side branches.
pub fn remove_side_branches(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    let mut removed_nodes = 1;
    while removed_nodes > 0 {
        removed_nodes = 0;

        let mut nodes_to_remove: Vec<_> = graph
            .node_indices()
            .filter(|&node| {
                if graph[node].is_end {
                    false
                } else {
                    graph.neighbors_directed(node, Direction::Outgoing).count() == 0
                }
            })
            .collect();

        nodes_to_remove.sort_by(|a, b| b.cmp(a));

        for node in nodes_to_remove {
            graph.remove_node(node);
            removed_nodes += 1;
        }
    }
}

// Get a vector of edge counts by traversing the graph backwards from the focal node
#[allow(dead_code)]
fn get_backward_edge_counts(
    graph: &StableDiGraph<DBNode, DBEdge>,
    focal_node: NodeIndex,
    depth: usize,
) -> Vec<u64> {
    let mut edge_counts = Vec::new();
    let mut current_node = focal_node;
    let mut current_depth = 0;

    while current_depth < depth {
        let incoming_edges: Vec<_> = graph
            .edges_directed(current_node, Direction::Incoming)
            .collect();

        if incoming_edges.is_empty() {
            // Stop if the current node has no incoming edges
            break;
        }

        // Assuming there's only one incoming edge per node in this context
        if let Some(edge) = incoming_edges.first() {
            edge_counts.push(edge.weight().count);
            current_node = edge.source();
        }

        current_depth += 1;
    }

    edge_counts
}

// Get a vector of node degrees by traversing the graph backwards from the focal node
fn get_backward_node_degrees(
    graph: &StableDiGraph<DBNode, DBEdge>,
    focal_node: NodeIndex,
    depth: usize,
) -> Vec<usize> {
    let mut node_degrees = Vec::new();
    let mut current_node = focal_node;
    let mut current_depth = 0;

    while current_depth < depth {
        // Get the number of outgoing edges from the current node
        let degree = graph
            .edges_directed(current_node, Direction::Outgoing)
            .count();
        node_degrees.push(degree);

        let incoming_edges: Vec<_> = graph
            .edges_directed(current_node, Direction::Incoming)
            .collect();

        // Move to the next node (the source of the first incoming edge)
        if let Some(edge) = incoming_edges.first() {
            current_node = edge.source();
        } else {
            // Break if there are no incoming edges to follow
            break;
        }

        current_depth += 1;
    }

    node_degrees
}

pub fn compute_mean(numbers: &[u64]) -> f64 {
    let sum: u64 = numbers.iter().sum();
    sum as f64 / numbers.len() as f64
}

pub fn compute_median(numbers: &[u64]) -> f64 {
    let mut sorted = numbers.to_vec();
    sorted.sort();

    let mid = sorted.len() / 2;

    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[mid] as f64
    }
}

// Given an oligo sequence, return a Oligo struct representing it
pub fn string_to_oligo(seq: &str) -> Oligo {
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
fn resolve_primer(primer: &str) -> HashSet<String> {
    let mut sequences: HashSet<String> = HashSet::new();

    for nuc in primer.chars() {
        let possible_nucs = match nuc {
            'R' => vec!['A', 'G'],
            'Y' => vec!['C', 'T'],
            'S' => vec!['G', 'C'],
            'W' => vec!['A', 'T'],
            'K' => vec!['G', 'T'],
            'M' => vec!['A', 'C'],
            'B' => vec!['C', 'G', 'T'],
            'D' => vec!['A', 'G', 'T'],
            'H' => vec!['A', 'C', 'T'],
            'V' => vec!['A', 'C', 'G'],
            'N' => vec!['A', 'C', 'G', 'T'],
            _ => vec![nuc], // Use the nucleotide directly
        };

        if sequences.is_empty() {
            for possible_nuc in possible_nucs {
                sequences.insert(possible_nuc.to_string());
            }
        } else {
            let mut new_sequences = HashSet::new();
            for seq in sequences.iter() {
                for &possible_nuc in &possible_nucs {
                    let new_seq = format!("{}{}", seq, possible_nuc);
                    new_sequences.insert(new_seq);
                }
            }
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

// Find the kmers that contain the oligos.
// If direction is forward, match the oligo at the start of the kmer.
// If direction is reverse, match the oligo at the end of the kmer.
fn find_oligos_in_kmers(oligos: &[Oligo], kmers: &KmerCounts, dir: &PrimerDirection, min_count: &u64) -> KmerCounts {
    // Assume all oligos have the same length
    let oligo_length = oligos[0].length;

    // Create hash set of oligos
    let oligo_set: HashSet<u64> = match dir {
        PrimerDirection::Forward => {
            // Rotate the oligo kmers to the start of the kmer
            oligos
                .iter()
                .map(|oligo| oligo.kmer << (2 * (kmers.get_k() - oligo_length)))
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
            mask <<= 2 * kmers.get_k() - 2 * oligo_length;
        }
        PrimerDirection::Reverse => {
            for _i in 0..(2 * oligo_length) {
                mask = (mask << 1) | 1;
            }
        }
    }

    let mut kmers_match: KmerCounts = KmerCounts::new(&kmers.get_k());

    for (kmer, count) in kmers.iter() {
        if oligo_set.contains(&(kmer & mask)) && count >= min_count {
            kmers_match.insert(kmer, count);
        }
    }

    kmers_match
}

fn get_suffix_mask(k: &usize) -> u64 {
    let suffix_mask: u64 = (1 << (2 * *k - 1)) - 1;
    suffix_mask
}

pub fn n_nonterminal_nodes_in_graph(graph: &StableDiGraph<DBNode, DBEdge>) -> usize {
    let mut n_nonterminal_nodes = 0;
    for node in graph.node_indices() {
        if !graph[node].is_terminal {
            n_nonterminal_nodes += 1;
        }
    }
    n_nonterminal_nodes
}

pub fn n_unvisited_nodes_in_graph(graph: &StableDiGraph<DBNode, DBEdge>) -> usize {
    let mut n = 0;
    for node in graph.node_indices() {
        if !graph[node].visited {
            n += 1;
        }
    }
    n
}

pub fn get_path_length(
    graph: &StableDiGraph<DBNode, DBEdge>,
    new_node: NodeIndex,
) -> Option<usize> {
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

pub fn get_dbedge(kmer: &u64, kmer_counts: &kmer::KmerCounts) -> DBEdge {
    DBEdge {
        _kmer: *kmer,
        count: kmer_counts.get_canonical(kmer),
    }
}

pub fn would_form_cycle(
    graph: &StableDiGraph<DBNode, DBEdge>,
    parent: NodeIndex,
    child: NodeIndex,
) -> bool {
    // If there's a path from the `child` to `parent`, adding an edge from `parent` to `child` would create a cycle
    let mut bfs = Bfs::new(graph, child);
    while let Some(node) = bfs.next(graph) {
        if node == parent {
            return true;
        }
    }
    false
}

pub fn get_start_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    let mut nodes: Vec<NodeIndex> = Vec::new();
    for node in graph.node_indices() {
        if graph[node].is_start {
            nodes.push(node)
        }
    }
    nodes
}

pub fn get_end_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    let mut nodes: Vec<NodeIndex> = Vec::new();
    for node in graph.node_indices() {
        if graph[node].is_end {
            nodes.push(node)
        }
    }
    nodes
}

// Find the descendants of a node, each descendent no more than `depth` edges away
fn descendants(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    depth: usize,
) -> HashSet<NodeIndex> {
    let mut visited: HashSet<NodeIndex> = HashSet::new();
    let mut queue = VecDeque::new();

    queue.push_back((node, 0));
    visited.insert(node);

    let mut descendants: HashSet<NodeIndex> = HashSet::new();

    while let Some((current_node, current_depth)) = queue.pop_front() {
        if current_depth >= depth {
            continue;
        }

        for neighbor in graph.neighbors(current_node) {
            if visited.insert(neighbor) {
                // insert returns false if the item was already in the set
                descendants.insert(neighbor);
                queue.push_back((neighbor, current_depth + 1));
            }
        }
    }

    descendants
}

// Get all descendants of a node in a directed graph
fn get_descendants(graph: &StableDiGraph<DBNode, DBEdge>, node: NodeIndex) -> Vec<NodeIndex> {
    let mut visited = HashSet::new();
    let mut stack = vec![node];
    visited.insert(node); // mark the original node as visited

    while let Some(node) = stack.pop() {
        for neighbor in graph.neighbors(node) {
            if !visited.contains(&neighbor) {
                visited.insert(neighbor); // mark neighbor as visited
                stack.push(neighbor);
            }
        }
    }

    // Remove the original node from the set before converting it to a Vec
    visited.remove(&node);

    visited.into_iter().collect()
}

fn pop_balloons(graph: &mut StableDiGraph<DBNode, DBEdge>, k: &usize, verbosity: usize) {
    let mut to_clip: Vec<NodeIndex> = Vec::new();
    for node in graph.node_indices() {
        let d = descendants(graph, node, EXTENSION_EVALUATION_DEPTH);
        let n = d.len();
        // The maximum number of descendants would be 4^EXTENSION_EVALUATION_DEPTH
        if n > 4_usize.pow((EXTENSION_EVALUATION_DEPTH) as u32) {
            let seq = crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k - 1));
            println!("{}", format!("WARNING: Node {} with sequence {} has {} descendants at a depth of {}. This exceed the maximum of 4^{}={} that is expected", node.index(), seq, n, EXTENSION_EVALUATION_DEPTH, EXTENSION_EVALUATION_DEPTH, 4_usize.pow((EXTENSION_EVALUATION_DEPTH) as u32)).color(COLOR_WARNING));
            // println!("  Descendants: {:?}", d);

            // Get a vector of sequences of the descendants
            let mut seqs: Vec<String> = Vec::new();
            for descendant in d {
                seqs.push(crate::kmer::kmer_to_seq(
                    &graph[descendant].sub_kmer,
                    &(*k - 1),
                ));
            }

            println!("  Sequences: {:?}", seqs);
        }

        if n > 4_usize.pow((EXTENSION_EVALUATION_DEPTH - EXTENSION_EVALUATION_DIFF) as u32) {
            to_clip.push(node);
            if verbosity > 1 {
                println!(
                    "  Node {} with sequence {} has {} descendants at a depth of {}, descendents will be clipped",
                    node.index(),
                    crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k - 1)),
                    n,
                    EXTENSION_EVALUATION_DEPTH
                );
            }
        }
    }

    // Mark the clipped nodes as terminal
    for node in &to_clip {
        graph[*node].is_terminal = true;
    }

    // Vector to hold nodes to be pruned
    let mut to_prune: Vec<NodeIndex> = Vec::new();
    for node in &to_clip {
        // Node may have been removed already, so check if it is in graph
        if graph.node_weight(*node).is_some() {
            // prune away all the descendants of the node, but keep the node
            let mut descendants = get_descendants(graph, *node);
            to_prune.append(&mut descendants);
        }
    }

    if !to_prune.is_empty() {
        println!(
            "    Removing {} nodes descended from {} nodes with ballooning graph extension",
            to_prune.len(),
            to_clip.len()
        );
    }

    // Sort in descending order. This is because node indices following pruned node are decremented,
    // so the highest ones need to be pruned first or the remaining indices are no longer valid
    to_prune.sort_by(|a, b| b.cmp(a));

    // Remove the nodes in to_prune
    for node in to_prune {
        graph.remove_node(node);
    }
}

pub fn summarize_extension(graph: &StableDiGraph<DBNode, DBEdge>, pad: &str) {
    // Print the number of nodes and edges in the graph
    println!("{}There are {} nodes in the graph", pad, graph.node_count());
    println!("{}There are {} edges in the graph", pad, graph.edge_count());

    let has_cycles = is_cyclic_directed(graph);
    if has_cycles {
        println!("{}The graph has cycles", pad);
    } else {
        println!("{}The graph does not have cycles", pad);
    }

    // Print the mean, median, and max degree of all nodes
    let mut degrees: Vec<usize> = Vec::new();
    for node in graph.node_indices() {
        degrees.push(graph.neighbors(node).count());
    }

    let max_degree = degrees.iter().max().unwrap();
    let degrees_u64: Vec<u64> = degrees.iter().map(|&x| x as u64).collect();
    let mean_degree = compute_mean(&degrees_u64);
    let median_degree = compute_median(&degrees_u64);

    println!(
        "{}Max node degree {}, mean {:.2} median {:.1}",
        pad, max_degree, mean_degree, median_degree
    );

    // Print the mean, median, and max count of edges
    let mut counts: Vec<u64> = Vec::new();
    for edge in graph.edge_indices() {
        counts.push(graph[edge].count);
    }

    let max_count = counts.iter().max().unwrap();
    let mean_count = compute_mean(&counts);
    let median_count = compute_median(&counts);

    println!(
        "{}Max edge count {}, mean {:.2} median {:.1}",
        pad, max_count, mean_count, median_count
    );

    // Create a vector of n_descendants for each node given a depth of EXTENSION_EVALUATION_DEPTH
    let mut n_descendants_vec: Vec<u64> = Vec::new();
    for node in graph.node_indices() {
        let n = descendants(graph, node, EXTENSION_EVALUATION_DEPTH).len();
        n_descendants_vec.push(n as u64);
    }

    let max_n_descendants = n_descendants_vec.iter().max().unwrap();
    let mean_n_descendants = compute_mean(&n_descendants_vec);
    let median_n_descendants = compute_median(&n_descendants_vec);

    println!(
        "{}Number of descendants to a depth of {}, max {} mean {:.2} median {:.1}",
        pad,
        EXTENSION_EVALUATION_DEPTH,
        max_n_descendants,
        mean_n_descendants,
        median_n_descendants
    );
}

fn pairwise_sequence_distances(records: &Vec<fasta::Record>) -> Vec<Vec<Option<u32>>> {
    // https://docs.rs/bio/latest/bio/alignment/distance/simd/fn.bounded_levenshtein.html

    let n = records.len();
    let mut matrix: Vec<Vec<Option<u32>>> = vec![vec![Some(0); n]; n];

    for i in 0..n {
        for j in i + 1..n {
            let dist =
                bounded_levenshtein(records[i].seq(), records[j].seq(), DISTANCE_EDIT_THRESHOLD);
            matrix[i][j] = dist;
            matrix[j][i] = dist; // Symmetric matrix
        }
    }

    matrix
}

// Individual steps in the PCR process

/// Ingest pcr parameters, and return a set of all possible variants of the
/// forward or reverse primer. The variants are generated by resolving
/// ambiguous nucleotides, and then permuting the sequences to include
/// up to `mismatches` mismatches.
fn preprocess_primer(
    params: &PCRParams,
    dir: PrimerDirection,
    k: &usize,
    _verbosity: &usize,
) -> HashSet<String> {
    let mut primer = params.forward_seq.clone();
    if dir == PrimerDirection::Reverse {
        primer = params.reverse_seq.clone();
    }

    let mut trim = params.trim;
    if trim > (k - 1) {
        println!("  Primer must not be longer than k-1. Trim length is {}, k is {}, so adjusting trim length to {}", trim, k, k-1);
        trim = k - 1;
    }

    // Check if either is longer than trim, if so retain only the last trim nucleotides
    if primer.len() > trim {
        primer = primer[primer.len() - trim..].to_string();
        println!(
            "  Trimming the primer to {} so that it is within the trim length of {}.",
            primer, trim
        );
    }

    // Expand ambiguous nucleotides
    let mut primer_variants = resolve_primer(&primer);

    // Get all possible variants of the primers
    primer_variants = permute_sequences(primer_variants, &params.mismatches);

    if dir == PrimerDirection::Reverse {
        // Replace the reverse variants with their reverse complements
        let mut primer_variants_revcomp = HashSet::new();
        for variant in primer_variants.iter() {
            primer_variants_revcomp.insert(reverse_complement(variant));
        }
        primer_variants = primer_variants_revcomp;
    }

    println!(
        "  There are {} variants of the primer",
        primer_variants.len()
    );

    primer_variants
}

/// Given a set of primer variants, return a set of kmers from the data that contain the primers
fn get_kmers_from_primers(
    primer_variants: &HashSet<String>,
    kmer_counts: &KmerCounts,
    dir: PrimerDirection,
    min_count: &u64,
) -> KmerCounts {
    // Get the kmers that contain the primers
    let mut oligos: Vec<Oligo> = Vec::new();
    for variant in primer_variants.iter() {
        oligos.push(string_to_oligo(variant));
    }

    find_oligos_in_kmers(&oligos, kmer_counts, &dir, min_count)
}

/// Given a set of kmers that contain primers, filter them to retain only those with the highest counts
/// If there are more than MAX_NUM_PRIMER_KMERS, retain only the MAX_NUM_PRIMER_KMERS with the highest counts
/// If there are less than MAX_NUM_PRIMER_KMERS, retain all of them
/// Returns a hash map of the kmers and their counts
fn filter_primer_kmers(matches: KmerCounts) -> KmerCounts {
    if matches.is_empty() {
        return matches;
    }

    let mut counts: Vec<u64> = matches.counts();

    counts.sort();
    counts.reverse();

    // set equal to last element
    let mut top_count_cutoff = counts.last().unwrap();

    if counts.len() > MAX_NUM_PRIMER_KMERS {
        top_count_cutoff = &counts[MAX_NUM_PRIMER_KMERS - 1];
    }

    // If there are less than MAX_NUM_PRIMER_KMERS matches, this is the lowest count
    // If there are more than MAX_NUM_PRIMER_KMERS matches, get the count of the
    // MAX_NUM_PRIMER_KMERS highest count matches and use that as the cutoff

    let mut matches_keep: KmerCounts = KmerCounts::new(&matches.get_k());
    for (kmer, count) in matches.iter() {
        let mut keep: bool = false;
        if count >= top_count_cutoff {
            // Add the kmer to the hash map
            matches_keep.insert(kmer, count);
            keep = true;
        }
        println!(
            "  {}, count {}, keep {}",
            crate::kmer::kmer_to_seq(kmer, &matches.get_k()),
            count,
            keep,
        );
    }

    // Replace matches with matches_keep
    matches_keep
}

fn get_median_edge_count(graph: &StableDiGraph<DBNode, DBEdge>) -> Option<f64> {
    let mut counts: Vec<u64> = Vec::new();
    for edge in graph.edge_indices() {
        counts.push(graph[edge].count);
    }

    if counts.is_empty() {
        return None;
    }

    let median_edge_count = compute_median(&counts);
    Some(median_edge_count)
}

fn get_primer_kmers(params: &PCRParams, kmer_counts: &KmerCounts, verbosity: &usize) -> (KmerCounts, KmerCounts) {
    // Preprocess the primers to get all variants to be considered
    let forward_variants = preprocess_primer(
        params,
        PrimerDirection::Forward,
        &kmer_counts.get_k(),
        verbosity,
    );
    let reverse_variants = preprocess_primer(
        params,
        PrimerDirection::Reverse,
        &kmer_counts.get_k(),
        verbosity,
    );

    // Get the kmers that contain the primers
    let mut forward_primer_kmers =
        get_kmers_from_primers(&forward_variants, &kmer_counts, PrimerDirection::Forward, &params.coverage);
    forward_primer_kmers = filter_primer_kmers(forward_primer_kmers);

    let mut reverse_primer_kmers =
        get_kmers_from_primers(&reverse_variants, &kmer_counts, PrimerDirection::Reverse, &params.coverage);
    reverse_primer_kmers = filter_primer_kmers(reverse_primer_kmers);

    (forward_primer_kmers, reverse_primer_kmers)
}


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
fn extend_graph(seed_graph: &StableDiGraph<DBNode, DBEdge>, kmer_counts: &KmerCounts, min_count: &u64, params: &PCRParams, verbosity: &usize) -> StableDiGraph<DBNode, DBEdge> {

    let suffix_mask: u64 = get_suffix_mask( &kmer_counts.get_k() );

    let mut graph = seed_graph.clone();

    let mut last_check: usize = 0;
    while n_unvisited_nodes_in_graph(&graph) > 0 {
        let n_nodes = graph.node_count();

        if n_nodes > MAX_NUM_NODES {
            println!("{}",
                format!("WARNING: There are {} nodes in the graph. This exceeds the maximum of {}, abandoning search.", n_nodes, MAX_NUM_NODES).color(COLOR_WARNING)
            );
            break;
        }

        // Get the median edge count, or min_count if there are no edges
        let edge_count_summary: f64 = match get_median_edge_count(&graph) {
            Some(count) => count,
            None => *min_count as f64,
        };

        // Periodically evaluate extension
        // After pruning, n_nodes can be less than last_check and there will be an overflow on subtracting
        if (n_nodes > last_check) && ((n_nodes - last_check) > EXTENSION_EVALUATION_FREQUENCY) {
            last_check = n_nodes - (n_nodes % EXTENSION_EVALUATION_FREQUENCY);

            println!("  Evaluating extension:");
            summarize_extension(&graph, "    ");

            // Some graphs balloon in size and get to hundreds of thousands of nodes while extension gets
            // slower and slower because there are so many growing tips. This may be due to a sequencing
            // adapter becoming integrated into the graph, for example. So periodically check for a region
            // of high degree and prune it if found
            pop_balloons(&mut graph, &kmer_counts.get_k(), *verbosity);
        }

        // Iterate over the nodes
        let node_indices: Vec<_> = graph.node_indices().collect();
        for node in node_indices {
            if !(graph[node].visited) {
                // Get the suffix of the kmer of the node
                let sub_kmer = graph[node].sub_kmer;

                if *verbosity > 1 {
                    print!(
                        "  {} sub_kmer being extended for node {}. ",
                        crate::kmer::kmer_to_seq(&sub_kmer, &(kmer_counts.get_k() - 1)),
                        node.index()
                    );
                    std::io::stdout().flush().unwrap();
                }

                // Get the kmers that could extend the node
                let mut candidate_kmers: HashSet<u64> = HashSet::new();

                for base in 0..4 {
                    let kmer = (sub_kmer << 2) | base;

                    // If the kmer is in the kmer_counts hash and has count >= min_count, 
                    // add it to the candidate kmers
                    if kmer_counts.contains(&kmer) {
                        if kmer_counts.get_count(&kmer) >= *min_count {
                            candidate_kmers.insert(kmer);
                        }
                    }
                    
                }

                if *verbosity > 1 {
                    print!(
                        "There are {} candidate kmers for extension. ",
                        candidate_kmers.len()
                    );
                    std::io::stdout().flush().unwrap();
                }

                // If there are no candidate kmers, the node is terminal
                if candidate_kmers.is_empty() {
                    if *verbosity > 1 {
                        println!("Marking node as terminal because there are no candidates for extension. ");
                        std::io::stdout().flush().unwrap();
                    }
                    graph[node].is_terminal = true;
                    graph[node].visited = true;
                    continue;
                }

                // Get the degrees of ancestor nodes, skipping degree of this node
                let node_degrees = get_backward_node_degrees(&graph, node, 20);
                let node_degrees_slice = &node_degrees[1..];

                // If the node is in a rapidly ballooning region of the graph, don't add it
                if candidate_kmers.len() > 2 {
                    if node_degrees_slice.len() >= 2 {
                        if (node_degrees_slice[0] > 2) && (node_degrees_slice[1] > 2) {
                            if *verbosity > 1 {
                                println!("Marking node as terminal because it and immediate ancestors have high degree. ");
                                std::io::stdout().flush().unwrap();
                            }
                            graph[node].is_terminal = true;
                            graph[node].visited = true;
                            continue;
                        }
                    }
                }

                // Check for a lower growth rate over a longer path, and terminate if found
                if node_degrees_slice.len() >= 15 {
                    // Get the number of elements of node_degrees_slice that are greater than 1
                    let n_high_degree = node_degrees_slice.iter().filter(|&x| *x > 1).count();
                    if n_high_degree >= 3 {
                        if *verbosity > 1 {
                            println!("Marking node as terminal because its recent ancestors have moderately elevated degree. ");
                            std::io::stdout().flush().unwrap();
                        }
                        graph[node].is_terminal = true;
                        graph[node].visited = true;
                        continue;
                    }
                }

                // Add new nodes if needed, and new edges
                for kmer in candidate_kmers.iter() {
                    let suffix = kmer & suffix_mask;

                    // Check if the node extends by itself and mark it as terminal if it does
                    if suffix == sub_kmer {
                        graph[node].is_terminal = true;
                        graph[node].visited = true;
                        if *verbosity > 1 {
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
                            if !would_form_cycle(&graph, node, existing_node) {
                                let edge = get_dbedge(kmer, &kmer_counts);
                                graph.add_edge(node, existing_node, edge);
                                if graph[existing_node].is_end {
                                    println!("End node incorporated into graph, complete PCR product found.");
                                }
                                let outgoing =
                                    graph.neighbors_directed(node, Direction::Outgoing).count();
                                if outgoing > 4 {
                                    println!("{}",
                                        format!("WARNING: Node {} has {} outgoing edges. This exceed the maximum of 4 that is expected", node.index(), outgoing).color(COLOR_WARNING)
                                    );
                                }
                            } else {
                                graph[node].is_terminal = true;

                                if *verbosity > 1 {
                                    print!(
                                        "Adding edge to node {} would form cycle. Not adding edge, and marking current node as terminal. ",
                                        node.index()
                                    );
                                    std::io::stdout().flush().unwrap();
                                }
                            }

                            node_exists = true;
                            break;
                        }
                    }

                    if !node_exists {
                        let edge = get_dbedge(kmer, &kmer_counts);
                        let edge_count = edge.count;

                        // Don't add node and edge if the edge count is very high
                        if (edge_count as f64)
                            > (edge_count_summary * BALLOONING_COUNT_THRESHOLD_MULTIPLIER)
                        {
                            if *verbosity > 1 {
                                print!(
                                    "Edge count {} exceeds {} * median edge count {}. Not adding edge with kmer {} or node with sub_kmer {}. ",
                                    edge_count,
                                    BALLOONING_COUNT_THRESHOLD_MULTIPLIER,
                                    edge_count_summary,
                                    crate::kmer::kmer_to_seq(kmer, &kmer_counts.get_k()),
                                    crate::kmer::kmer_to_seq(&suffix, &(kmer_counts.get_k() - 1))
                                );
                                std::io::stdout().flush().unwrap();
                            }
                            continue;
                        }

                        let new_node = graph.add_node(DBNode {
                            sub_kmer: suffix,
                            is_start: false,
                            is_end: false,
                            is_terminal: false,
                            visited: false,
                        });

                        graph.add_edge(node, new_node, edge);
                        let outgoing = graph.neighbors_directed(node, Direction::Outgoing).count();
                        if outgoing > 4 {
                            println!("{}",
                                format!("WARNING: Node {} has {} outgoing edges when adding new node. This exceed the maximum of 4 that is expected", node.index(), outgoing).color(COLOR_WARNING)
                            );
                        }

                        if *verbosity > 1 {
                            print!(
                                "Added sub_kmer {} for new node {} with edge kmer count {}. ",
                                crate::kmer::kmer_to_seq(&suffix, &(kmer_counts.get_k() - 1)),
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
                            if *verbosity > 1 {
                                print!("Path length is {}. ", path_length);
                                std::io::stdout().flush().unwrap();
                            }

                            if path_length > params.max_length - (kmer_counts.get_k()) {
                                graph[new_node].is_terminal = true;
                                graph[new_node].visited = true;
                                if *verbosity > 1 {
                                    print!("Marking new node {} as terminal because it exceeds max_length from start. ", new_node.index());
                                    std::io::stdout().flush().unwrap();
                                }
                            }
                        } else {
                            graph[new_node].is_terminal = true;
                            if *verbosity > 1 {
                                print!("Marking new node {} as terminal because it is part of a cycle. ", new_node.index());
                                std::io::stdout().flush().unwrap();
                            }
                        }
                    }
                }
                graph[node].visited = true;

                if *verbosity > 1 {
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

    graph

}

#[derive(Clone)]
pub struct PCRParams {
    pub forward_seq: String,
    pub reverse_seq: String,
    pub max_length: usize,
    pub gene_name: String,
    pub coverage: u64,
    pub mismatches: usize,
    pub trim: usize,
}

// The primary function for PCR
pub fn do_pcr(
    kmer_counts: &KmerCounts,
    sample_name: &str,
    verbosity: usize,
    params: &PCRParams,
) -> Vec<bio::io::fasta::Record> {
    println!(
        "{}",
        format!("Running PCR on gene {}", params.gene_name).color(COLOR_NOTE)
    );

    println!("Preprocessing primers");
    let (forward_primer_kmers, reverse_primer_kmers) = get_primer_kmers(params, kmer_counts, &verbosity);
    
    if forward_primer_kmers.is_empty() {
        println!(
            "{}",
            format!(
                "For gene {}, binding sites were not found for the forward primer. Abandoning PCR.",
                params.gene_name
            )
            .color(COLOR_FAIL)
        );
        println!(
            "{}",
            "  Suggested action: optimize primer sequence."
                .to_string()
                .color(COLOR_FAIL)
        );
        let records: Vec<fasta::Record> = Vec::new();
        return records;
    }

    if reverse_primer_kmers.is_empty() {
        println!(
            "{}",
            format!(
                "For gene {}, binding sites were not found for the reverse primer. Abandoning PCR.",
                params.gene_name
            )
            .color(COLOR_FAIL)
        );
        println!(
            "{}",
            "  Suggested action: optimize primer sequence."
                .to_string()
                .color(COLOR_FAIL)
        );
        let records: Vec<fasta::Record> = Vec::new();
        return records;
    }

    // Construct the graph
    println!("Creating graph, seeding with nodes that contain primer matches...");
    let mut seed_graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();

    // Create a suffix mask with 1 in the (2 * (k - 1)) least significant bits
    let suffix_mask: u64 = get_suffix_mask( &kmer_counts.get_k() );

    // Add the forward primer matches to the graph
    for kmer in forward_primer_kmers.kmers() {
        let prefix = kmer >> 2;

        // If the node with sub_kmer == suffix already exists, update the node so that is_start = true
        // Otherwise, create a new node

        let mut node_exists = false;
        for node in seed_graph.node_indices() {
            if seed_graph[node].sub_kmer == prefix {
                seed_graph[node].is_start = true;
                node_exists = true;
                break;
            }
        }

        if !node_exists {
            seed_graph.add_node(DBNode {
                sub_kmer: prefix,
                is_start: true,
                is_end: false,
                is_terminal: false,
                visited: false, // Will not be extended
            });
        }
    }

    // Add the reverse primer matches to the graph
    for kmer in reverse_primer_kmers.kmers() {
        let suffix = kmer & suffix_mask;

        // If the node with sub_kmer == suffix already exists, update the node so that is_end = true
        // Otherwise, create a new node

        let mut node_exists = false;
        for node in seed_graph.node_indices() {
            if seed_graph[node].sub_kmer == suffix {
                seed_graph[node].is_end = true;
                seed_graph[node].is_terminal = true;
                node_exists = true;
                break;
            }
        }

        if !node_exists {
            seed_graph.add_node(DBNode {
                sub_kmer: suffix,
                is_start: false,
                is_end: true,
                is_terminal: true,
                visited: true,
            });
        }
    }

    // Loop over the nodes and see if any of the nodes are start and end nodes
    // If so, print a warning
    for node in seed_graph.node_indices() {
        if seed_graph[node].is_start && seed_graph[node].is_end {
            println!(
                "{}",
                format!(
                    "WARNING: node {} is both a start and end node",
                    node.index()
                )
                .color(COLOR_WARNING)
            );
        }
    }

    println!("There are {} start nodes", get_start_nodes(&seed_graph).len());
    println!("There are {} end nodes", get_end_nodes(&seed_graph).len());

    if verbosity > 1 {
        // Print the information for each node
        for node in seed_graph.node_indices() {
            println!("Node {}:", node.index());
            println!(
                "  sub_kmer: {}",
                crate::kmer::kmer_to_seq(&seed_graph[node].sub_kmer, &(kmer_counts.get_k() - 1))
            );
            println!("  is_start: {}", seed_graph[node].is_start);
            println!("  is_end: {}", seed_graph[node].is_end);
            println!("  is_terminal: {}", seed_graph[node].is_terminal);
        }
    }

    let start = std::time::Instant::now();
    println!("Extending the assembly graph...");

    // Create a vector to hold the fasta records
    let mut assembly_records: Vec<AssemblyRecord> = Vec::new();

    // Get the minimum of (max_forward_count, max_reverse_count), consider this as the observed count of the primers
    // and a preliminary expectation for the coverage of the PCR product
    let max_forward_count = forward_primer_kmers.get_max_count();
    let max_reverse_count = reverse_primer_kmers.get_max_count();

    let mut primer_count = max_reverse_count;
    if max_forward_count < max_reverse_count {
        primer_count = max_forward_count;
    }
    
    println!("Observed primer coverage is {}, user specified minimum coverage is {}", primer_count, params.coverage);

    // Creates a vector of coverage thresholds, starting with coverage_high_threshold and decreasing to params.coverage
    // in COVERAGE_STEPS steps
    let coverage_high_threshold = primer_count / COVERAGE_MULTIPLIER;
    let mut coverage_thresholds: Vec<u64> = Vec::new();
    let step_size = (coverage_high_threshold - params.coverage) / (COVERAGE_STEPS-1);

    for i in 0..COVERAGE_STEPS {
        coverage_thresholds.push(coverage_high_threshold - (i * step_size));
    }

    // Make sure the last element is params.coverage, even if there are rounding errors
    coverage_thresholds[COVERAGE_STEPS as usize -1] = params.coverage;

    for min_count in coverage_thresholds.iter() {
        println!("  Extending graph with minimum kmer count {}", min_count);
        let mut graph = extend_graph(&seed_graph, &kmer_counts, &min_count, &params, &verbosity);

        // Make hashmaps of the start and end nodes, where the key is the node index and the value is the number of edges
        let mut start_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();
        let mut end_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();

        for node in graph.node_indices() {
            if graph[node].is_start {
                start_nodes_map.insert(
                    node,
                    graph.neighbors_directed(node, Direction::Outgoing).count(),
                );
            }
            if graph[node].is_end {
                end_nodes_map.insert(
                    node,
                    graph.neighbors_directed(node, Direction::Incoming).count(),
                );
            }
        }

        if verbosity > 0 {
            // Print the number of edges for each start node
            for (node, edge_count) in &start_nodes_map {
                // Print the node subkmer and number of edges
                println!(
                    "  Start node {} with subkmer {} has {} edges",
                    node.index(),
                    crate::kmer::kmer_to_seq(&graph[*node].sub_kmer, &(kmer_counts.get_k() - 1)),
                    edge_count
                );
            }

            // Print the number of edges for each end node
            for (node, edge_count) in &end_nodes_map {
                // Print the node subkmer and number of edges
                println!(
                    "  End node {} with subkmer {} has {} edges",
                    node.index(),
                    crate::kmer::kmer_to_seq(&graph[*node].sub_kmer, &(kmer_counts.get_k() - 1)),
                    edge_count
                );
            }
        }

        // Drop entries that have no edges
        start_nodes_map.retain(|_node, edge_count| *edge_count > 0);
        end_nodes_map.retain(|_node, edge_count| *edge_count > 0);

        println!("Final extension statistics:");
        summarize_extension(&graph, "    ");
        println!(
            "  There are {} start nodes with edges",
            start_nodes_map.len()
        );
        println!("  There are {} end nodes with edges", end_nodes_map.len());

        println!("done.  Time to extend graph: {:?}", start.elapsed());

        if verbosity > 5 {
            let dot_format = format!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));

            // Write the DOT format to a file
            let file_name = format!("{}_{}.dot", sample_name, params.gene_name); // Concatenating the file extension
            println!("Writing dot file {}", file_name);
            let mut file = File::create(&file_name).expect("Unable to create file");
            file.write_all(dot_format.as_bytes())
                .expect("Unable to write data");
        }

        // Simplify the graph
        // all_simple_paths() hangs if the input graph is too complex
        // Also want to regularize some graph features

        let start = std::time::Instant::now();
        println!("Pruning the assembly graph...");

        remove_side_branches(&mut graph);

        // Remove end nodes without edges
        let mut nodes_to_remove: Vec<NodeIndex> = graph
            .node_indices()
            .filter(|&node| {
                if graph[node].is_end {
                    graph.neighbors_directed(node, Direction::Incoming).count() == 0
                } else {
                    false
                }
            })
            .collect();

        nodes_to_remove.sort_by(|a, b| b.cmp(a));

        for node in nodes_to_remove {
            graph.remove_node(node);
        }

        println!("  There are {} nodes in the graph", graph.node_count());
        println!("  There are {} edges in the graph", graph.edge_count());

        println!("  There are {} start nodes", get_start_nodes(&graph).len());
        println!("  There are {} end nodes", get_end_nodes(&graph).len());

        println!("done.  Time to prune graph: {:?}", start.elapsed());

        // Get all paths from start nodes to terminal nodes
        let start = std::time::Instant::now();
        println!("Traversing the assembly graph to find paths from forward to reverse primers...");
        let mut all_paths = Vec::new();

        for start in get_start_nodes(&graph) {
            for end in get_end_nodes(&graph) {
                let paths_for_this_pair =
                    all_simple_paths::<Vec<NodeIndex>, &StableDiGraph<DBNode, DBEdge>>(
                        &graph,
                        start,
                        end,
                        1,
                        Some(params.max_length - (kmer_counts.get_k()) + 1),
                    );

                all_paths.extend(paths_for_this_pair);
            }
        }

        println!(
            "  There are {} paths from forward to reverse primers in the graph",
            all_paths.len()
        );
        println!("done.  Time to traverse graph: {:?}", start.elapsed());

        if all_paths.is_empty() {
            println!("{}", format!("For gene {}, no path was found from a forward primer binding site to a reverse binding \nsite. Abandoning PCR.", params.gene_name).color(COLOR_FAIL));
            println!("{}", format!("  Suggested actions:").color(COLOR_FAIL));
            println!("{}", format!("    - The max_length for the PCR product of {} my be too short. Consider increasing it.", params.max_length).color(COLOR_FAIL));
            println!("{}", format!("    - The primers may have non-specific binding and are not close enough to generate a \n      product. Consider increasing the primer TRIM length from the default to \n      create a more specific primer.").color(COLOR_FAIL));
            println!("{}", format!("      - The maximum count of a forward kmer is {} and of a reverse kmer is {}. Large \n        differences in value can indicate non-specific binding of one of the primers.", max_forward_count, max_reverse_count).color(COLOR_FAIL));

            let count_threshold = 5;
            if (max_forward_count < count_threshold) | (max_reverse_count < count_threshold) {
                println!("{}", format!("    - The maximum count of a forward kmer is {} and of a reverse kmer is {}. A \n      low value, in this case less than {}, for either can indicate that read \n      coverage for this gene is too low to traverse from a forward to reverse \n      primer. Consider increasing coverage.", max_forward_count, max_reverse_count, count_threshold).color(COLOR_FAIL));
            }

            let records: Vec<fasta::Record> = Vec::new();
            return records;
        }

        println!("Generating sequences from paths...");

        // For each path, get the sequence of the path
        for (i, path) in all_paths.into_iter().enumerate() {
            let mut sequence = String::new();
            let mut edge_counts: Vec<u64> = Vec::new();
            let mut parent_node: NodeIndex = NodeIndex::new(0);
            // The first time through the loop add the whole sequence, after that just add the last base
            for node in path.iter() {
                let node_data = graph.node_weight(*node).unwrap();
                let subread = crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(kmer_counts.get_k() - 1));
                if sequence.is_empty() {
                    sequence = subread;
                    parent_node = *node;
                } else {
                    sequence = format!("{}{}", sequence, subread.chars().last().unwrap(),);

                    // Get the edge count for the edge from the parent node to this node
                    let edge = graph.find_edge(parent_node, *node).unwrap();
                    let edge_data = graph.edge_weight(edge).unwrap();
                    edge_counts.push(edge_data.count);
                    parent_node = *node;
                }
            }

            // Get some stats on the path counts
            let count_mean = compute_mean(&edge_counts);
            let count_median = compute_median(&edge_counts);
            let count_min = edge_counts.iter().min().unwrap();
            let count_max = edge_counts.iter().max().unwrap();

            let id = format!(
                "{} {} product {} length {} kmer count stats mean {:.2} median {} min {} max {}",
                sample_name,
                params.gene_name,
                i,
                sequence.len(),
                count_mean,
                count_median,
                count_min,
                count_max
            );
            println!(">{}", id);
            println!("{}", sequence);
            let record = fasta::Record::with_attrs(&id, None, sequence.as_bytes());
            // Create fasta record and add to vector
            let assembly_record: AssemblyRecord = AssemblyRecord {
                fasta_record: record,
                kmer_min_count: *count_min,
            };
            assembly_records.push(assembly_record);
        }
        if assembly_records.is_empty() {
            print!("Did not obtain PCR product.");
        } else {
            print!("Obtained PCR product.");
            break;
        }
    }

    println!("done.");

    // Order the records by descending kmer_min_count
    assembly_records.sort_by(|a, b| b.kmer_min_count.cmp(&a.kmer_min_count));

    let mut records: Vec<fasta::Record> = Vec::new();
    for fasta_record in assembly_records {
        records.push(fasta_record.fasta_record);
    }

    let num_records_all = records.len();

    // Thin the records to remove nearly-duplicate sequences that are highly similar. These can be abundant when coverage is high.
    // Calculate pairwise distances between all sequences
    let distances = pairwise_sequence_distances(&records);

    // Get the first column of the matrix
    let first_column: Vec<Option<u32>> = distances
        .iter()
        .filter_map(|row| row.get(0).cloned())
        .collect();

    // Loop through the indices and values of first_column in reverse order. If the value is Some, drop the corresponding
    // element of records because it is so similar to the first.
    for (index, &value) in first_column.iter().enumerate().rev() {
        if index > 0 && value.is_some() {
            records.remove(index);
        }
    }

    if num_records_all == records.len() {
        println!(
            "{}",
            format!(
                "For gene {}, {} PCR products were generated and retained.",
                params.gene_name, num_records_all
            )
            .color(COLOR_SUCCESS)
        );
    } else {
        println!("{}", format!("For gene {}, {} PCR products were generated and {} were retained (the others were minor variants of the first).", params.gene_name, num_records_all, records.len()).color(COLOR_SUCCESS));
    }

    // Return the records
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

    // Create a test graph of known structure
    fn create_test_graph() -> (
        StableDiGraph<DBNode, DBEdge>,
        HashMap<&'static str, NodeIndex>,
    ) {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut nodes = HashMap::new();

        // Add nodes and store their indices in the HashMap
        nodes.insert(
            "a",
            graph.add_node(DBNode {
                sub_kmer: 0,
                is_start: true,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "b",
            graph.add_node(DBNode {
                sub_kmer: 1,
                is_start: false,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "c",
            graph.add_node(DBNode {
                sub_kmer: 2,
                is_start: false,
                is_end: false,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "d",
            graph.add_node(DBNode {
                sub_kmer: 3,
                is_start: false,
                is_end: true,
                is_terminal: false,
                visited: false,
            }),
        );
        nodes.insert(
            "e",
            graph.add_node(DBNode {
                sub_kmer: 4,
                is_start: false,
                is_end: true,
                is_terminal: false,
                visited: false,
            }),
        );

        // Add directed edges
        graph.add_edge(
            nodes["a"],
            nodes["b"],
            DBEdge {
                _kmer: 10,
                count: 5,
            },
        );
        graph.add_edge(
            nodes["b"],
            nodes["c"],
            DBEdge {
                _kmer: 11,
                count: 10,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["d"],
            DBEdge {
                _kmer: 12,
                count: 4,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["e"],
            DBEdge {
                _kmer: 13,
                count: 1,
            },
        );

        (graph, nodes)
    }

    #[test]
    fn test_get_descendants() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(get_descendants(&graph, nodes["a"]).len(), 4);
        assert_eq!(get_descendants(&graph, nodes["b"]).len(), 3);
        assert_eq!(get_descendants(&graph, nodes["c"]).len(), 2);
        assert_eq!(get_descendants(&graph, nodes["d"]).len(), 0);
        assert_eq!(get_descendants(&graph, nodes["e"]).len(), 0);
    }

    #[test]
    fn test_n_descendants() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(descendants(&graph, nodes["a"], 1).len(), 1); // Direct successor
        assert_eq!(descendants(&graph, nodes["a"], 2).len(), 2); // Including b's successors
        assert_eq!(descendants(&graph, nodes["a"], 3).len(), 4); // All nodes reachable from a within 3 steps
        assert_eq!(descendants(&graph, nodes["a"], 4).len(), 4); // All nodes reachable from a within 4 steps
        assert_eq!(descendants(&graph, nodes["b"], 2).len(), 3);
    }

    #[test]
    fn test_get_backward_edge_counts() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(get_backward_edge_counts(&graph, nodes["a"], 3).len(), 0);
        assert_eq!(get_backward_edge_counts(&graph, nodes["b"], 3), [5]);
        assert_eq!(get_backward_edge_counts(&graph, nodes["c"], 3), [10, 5]);
        assert_eq!(get_backward_edge_counts(&graph, nodes["d"], 3), [4, 10, 5]);
        assert_eq!(get_backward_edge_counts(&graph, nodes["e"], 3), [1, 10, 5]);
    }

    #[test]
    fn test_get_backward_node_degrees() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(get_backward_node_degrees(&graph, nodes["a"], 3), [1]);
        assert_eq!(get_backward_node_degrees(&graph, nodes["b"], 3), [1, 1]);
        assert_eq!(get_backward_node_degrees(&graph, nodes["c"], 3), [2, 1, 1]);
        assert_eq!(
            get_backward_node_degrees(&graph, nodes["d"], 4),
            [0, 2, 1, 1]
        );
        assert_eq!(
            get_backward_node_degrees(&graph, nodes["e"], 4),
            [0, 2, 1, 1]
        );
        assert_eq!(get_backward_node_degrees(&graph, nodes["d"], 3), [0, 2, 1]);
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
        let result1 = resolve_primer(&seq1);
        assert_eq!(result1, expected1);

        // Check with one ambiguous nucleotide
        let seq2 = "CGTAATGCGGCGN".to_string();
        let mut expected2: HashSet<String> = HashSet::new();
        expected2.insert("CGTAATGCGGCGA".to_string());
        expected2.insert("CGTAATGCGGCGC".to_string());
        expected2.insert("CGTAATGCGGCGG".to_string());
        expected2.insert("CGTAATGCGGCGT".to_string());

        let result2 = resolve_primer(&seq2);
        assert_eq!(result2, expected2);

        // Check with one ambiguous nucleotide
        let seq3 = "CGTAATRCGGCGA".to_string();
        let mut expected3: HashSet<String> = HashSet::new();
        expected3.insert("CGTAATACGGCGA".to_string());
        expected3.insert("CGTAATGCGGCGA".to_string());
        let result3 = resolve_primer(&seq3);
        assert_eq!(result3, expected3);

        // Check with two ambiguous nucleotides
        let seq4 = "CGTAATRCGGCGY".to_string();
        let mut expected4: HashSet<String> = HashSet::new();
        expected4.insert("CGTAATACGGCGC".to_string());
        expected4.insert("CGTAATGCGGCGC".to_string());
        expected4.insert("CGTAATACGGCGT".to_string());
        expected4.insert("CGTAATGCGGCGT".to_string());

        let result4 = resolve_primer(&seq4);
        assert_eq!(result4, expected4);

        // Check when the first nucleotide is ambiguous
        let seq5 = "RCGTAATCGGCGA".to_string();
        let mut expected5: HashSet<String> = HashSet::new();
        expected5.insert("ACGTAATCGGCGA".to_string());
        expected5.insert("GCGTAATCGGCGA".to_string());
        let result5 = resolve_primer(&seq5);
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

        // Check for inclusion of particular strings as a spot check
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("TGCAGGTTCACCTAC".to_string());
            let r: usize = 2;
            let result = permute_sequences(seq, &r);

            // Check for inclusion of original sequence
            assert!(result.contains("TGCAGGTTCACCTAC"));

            // Check for inclusion of "GGCAGGTTCACCTAC"
            assert!(result.contains("GGCAGGTTCACCTAC"));

            
        }
    }

    #[test]
    fn test_string_to_oligo() {
        let oligo = string_to_oligo("GCGA");
        assert_eq!(oligo.kmer, 0b1001_1000);
        assert_eq!(oligo.length, 4);
    }

    #[test]
    fn test_get_end_nodes() {
        let (graph, nodes) = create_test_graph();

        let end_nodes = get_end_nodes(&graph);
        assert_eq!(end_nodes.len(), 2);
        assert!(end_nodes.contains(&nodes["d"]));
        assert!(end_nodes.contains(&nodes["e"]));
    }

    #[test]
    fn test_get_start_nodes() {
        let (graph, nodes) = create_test_graph();

        let start_nodes = get_start_nodes(&graph);
        assert_eq!(start_nodes.len(), 1);
        assert!(start_nodes.contains(&nodes["a"]));
    }

    fn build_test_case() -> (String, usize, usize, KmerCounts, PCRParams) {
        // This is the 18s that was assembled from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR26955578, per the readme tutorial
        // Padded with C's at the ends
        let read_string = "CCCCCCCCCCCCGTTGATCCTGCCAGTATCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCACTTGTACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATTGTACTCTCTTACTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCGAAAAGTCCCGACTCTCGTGGAAGGGATGTATTTATTAGATTAAAAACCAATGCGGCTTAACGGCCGCTTACAAACTTGGTGATTCATAGTAACTGTTCGAATCGCATGGCCTTCTCTGTTCGTGCCGGCGATGTTTCATTCAAATTTCTGCCCTATCAACTGTCGATGGTAAGATAGTGGCTTACCATGGTCGCAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACGTGGGGAGGTAGTGACAAAAAATAACAATACAGGGCTTTTTTGTAGTCTTGTAATTGGAATGAGTACAATTTAAATCTCTTAACGAGGACCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATTGTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGGCTCGACGGCGACGGTCAGCCGCAAGGTATGTCACTGTCGACGTTGGCCTTCTTCGCGCAGACTTCGCGTGCTCTTAACTGAGTGTGCGTTGGATACGCGACGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTTGTGCTTGGATACATAAGCATGGAATAATGGAATAGGACTTTGGTTCTATTTTCCGTTGGTTTCTGGAACCGAAGTAATGATTAATAGGGACAGTTGGGGGCATTCGTATTTCGTTGTCAGAGGTGAAATTCTTGGATTTACGAAAGACGAACTAATGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGATCGAAGACGATCAGATACCGTCCTAGTTCTAACCATAAACGATGCCGACTAGGGATCAGCGAGTGTTATTTGATGACCTCGTTGGCACCTTATGGGAAACCAAAGTTTTTGGGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTCACCAGGTCCAGACATAGTAAGGATTGACAGATTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTGACCTGCTAAATAGTCAGACGATTCTCGAATCGCTCTCGACTTCTTAGAGGGACTGTTGCGTGTTTAACCAAAGTCAGGAAGGCAATAACAGGTCTGTGATGCCCTTAGATGTCCTGGGCCGCACGCGCGCTACACTGACGATGGCAACGAGTCGCTCCTTCACCGAAAGGTGTGGGTAATCTTGTGAATCATCGTCGTGCTGGGGATAGATCATTGTAATTCTTGATCTTGAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCGGATTGGCACTGTCAGATGGGCTTCGGTCCATCCGACGGACGTCAAAAAGTTGGTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGCCCCCCCCCCCC".to_string();

        let k:usize = 21;
        let replicates: usize = 10;
        let mut kmer_counts = KmerCounts::new(&k);
        for _i in 0..replicates {
            let reads = kmer::seq_to_reads(&read_string);
            kmer_counts.ingest_reads(&reads);
        }

        kmer_counts.add_reverse_complements();

        let params = PCRParams {
            forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),

            // AGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCG 3' end of 18S above
            // CGCAGGTTCACCTACGGAAACCTTGTTACGACTTTTACT revcomp of above

            reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
            max_length: 2500,
            gene_name: "18s".to_string(),
            coverage: 3,
            mismatches: 2,
            trim: 15,
        };

        (read_string, k, replicates, kmer_counts, params)
    }

    #[test]
    fn test_primer_preprocessing_steps() {
        let (_, _, _, kmer_counts, params) = build_test_case();
        
        let verbosity: usize = 0;
        let reverse_variants = preprocess_primer(
            &params,
            PrimerDirection::Reverse,
            &kmer_counts.get_k(),
            &verbosity,
        );

        // There should be 991 variants of the reverse primer when r=2
        assert_eq!(reverse_variants.len(), 991);

        println!("Reverse variants: {:?}", reverse_variants);

        // Check for inclusion of revcomp of original sequence TGCAGGTTCACCTAC
        let member = "GTAGGTGAACCTGCA".to_string();
        assert!(reverse_variants.contains(&member));

        // Check for inclusion of revcomp of off by one variant GGCAGGTTCACCTAC
        let member = "GTAGGTGAACCTGCC".to_string();
        assert!(reverse_variants.contains(&member));

        let mut reverse_primer_kmers =
        get_kmers_from_primers(&reverse_variants, &kmer_counts, PrimerDirection::Reverse, &params.coverage);
        
        // Check for kmer
        assert_eq!(reverse_primer_kmers.len(), 1);
        
        reverse_primer_kmers = filter_primer_kmers(reverse_primer_kmers);

        // Check for kmer after filtering
        assert_eq!(reverse_primer_kmers.len(), 1);
    }

    #[test]
    fn test_integration() {

        let (read_string, k, replicates, kmer_counts, params) = build_test_case();

        // Check the number of kmers
        // Times 2 on right since reverse complements are added
        assert_eq!(kmer_counts.len(), (read_string.len() - k + 1)*2);

        // Check the total count of kmers
        assert_eq!(kmer_counts.get_n_kmers(), ((read_string.len() - k + 1) * replicates * 2) as u64);

        let verbosity = 0;
        let (forward_primer_kmers, reverse_primer_kmers) = get_primer_kmers(&params, &kmer_counts, &verbosity);

        assert_eq!(forward_primer_kmers.len(), 1);
        assert_eq!(reverse_primer_kmers.len(), 1);


    }
}
