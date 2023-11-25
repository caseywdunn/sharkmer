use bio::bio_types::sequence;
use bio::io::fasta;
use colored::*;
use petgraph::algo::all_simple_paths;
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::Direction;
use petgraph::visit::{Dfs, Walker};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::Write;

use crate::COLOR_FAIL;
use crate::COLOR_NOTE;
use crate::COLOR_SUCCESS;
use crate::COLOR_WARNING;

const MAX_NUM_NODES: usize = 5_000;

fn get_consensus(sequences: &Vec<String>) -> String {
    let mut consensus = String::new();
    let mut consensus_vec: Vec<char> = Vec::new();

    // Get the length of the shortest sequence
    let shortest_sequence_length = sequences
        .iter()
        .map(|sequence| sequence.len())
        .min()
        .unwrap_or(0);

    // Get the consensus sequence
    for i in 0..shortest_sequence_length {
        let mut base_counts: HashMap<char, usize> = HashMap::new();
        for sequence in sequences {
            let base = sequence.chars().nth(i).unwrap();
            let count = base_counts.entry(base).or_insert(0);
            *count += 1;
        }

        let mut max_count = 0;
        let mut max_base = 'N';
        for (base, count) in base_counts {
            if count > max_count {
                max_count = count;
                max_base = base;
            }
        }

        consensus_vec.push(max_base);
    }

    consensus = consensus_vec.into_iter().collect();

    consensus
}

fn mode(numbers: &Vec<usize>) -> usize {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for number in numbers {
        let count = counts.entry(*number).or_insert(0);
        *count += 1;
    }

    let mut max_count = 0;
    let mut mode = 0;
    for (number, count) in counts {
        if count > max_count {
            max_count = count;
            mode = number;
        }
    }

    mode
}

/// Removes all nodes from the graph that are not ancestors of the specified node.
///
/// # Arguments
///
/// * `graph` - The directed graph.
/// * `node` - The node whose ancestors will be kept.
fn remove_non_ancestors (graph: &mut StableDiGraph<crate::pcr::DBNode, crate::pcr::DBEdge>, node: &NodeIndex) {
    // First, find all the ancestors of the given node
    let mut ancestors = HashSet::new();
    let mut dfs = Dfs::new(&*graph, *node);
    while let Some(nx) = dfs.next(&*graph) {
        ancestors.insert(nx);
    }

    // Then, remove all nodes that are not ancestors of the given node
    let nodes_to_remove: Vec<_> = graph.node_indices()
                                      .filter(|&n| !ancestors.contains(&n))
                                      .collect();
    
    for n in nodes_to_remove {
        graph.remove_node(n);
    }
}

/// For a given starting kmer, generate fasta records
fn generate_fastas(
    params: &RADParams,
    sample_name: &str,
    kmer_counts: &FxHashMap<u64, u64>,
    kmers: &std::collections::HashSet<u64>,
    starting_kmer: &u64,
    k: &usize,
    verbosity: usize,
) -> Vec<fasta::Record> {
    let mut to_print = String::new();
    let mut records_vec: Vec<fasta::Record> = Vec::new();
    let prefix = starting_kmer >> 2;
    let starting_seq = crate::kmer::kmer_to_seq(starting_kmer, k);
    let suffix_mask: u64 = (1 << (2 * (*k - 1))) - 1;

    // Get the kmer and mask for cut1 at the end of the kmer
    let cut1_kmer_terminal = crate::pcr::string_to_oligo(&params.cut1).kmer;
    let mut cut1_mask_terminal: u64 = 0;
    for _i in 0..(2 * params.cut1.len()) {
        cut1_mask_terminal = (cut1_mask_terminal << 1) | 1;
    }

    // Get the kmer and mask for cut2 at the end of the kmer
    let cut2_kmer_terminal = crate::pcr::string_to_oligo(&params.cut2).kmer;
    let mut cut2_mask_terminal: u64 = 0;
    for _i in 0..(2 * params.cut2.len()) {
        cut2_mask_terminal = (cut2_mask_terminal << 1) | 1;
    }

    to_print = format!("{}{} Starting kmer analysis\n", to_print, starting_seq);

    // Create a graph with the starting kmer as the only node
    let mut graph: StableDiGraph<crate::pcr::DBNode, crate::pcr::DBEdge> = StableDiGraph::new();
    graph.add_node(crate::pcr::DBNode {
        sub_kmer: prefix,
        is_start: true,
        is_end: false,
        is_terminal: false,
        visited: false, // Will not be extended
    });

    // Extend the graph until all nodes have been visited
    while crate::pcr::n_unvisited_nodes_in_graph(&graph) > 0 {
        let n_nodes = graph.node_count();

        if n_nodes > MAX_NUM_NODES {
            to_print = format!("{}  There are {} nodes in the graph. This exceeds the maximum of {}, abandoning search.\n", to_print, n_nodes, MAX_NUM_NODES);
            if verbosity > 0 {
                print!("{}", to_print);
            }
            return records_vec;
        }

        let node_indices: Vec<_> = graph.node_indices().collect();
        'nodes: for node in node_indices {
            if !(graph[node].visited) {
                // Get the suffix of the kmer of the node
                let sub_kmer = graph[node].sub_kmer;

                if verbosity > 9 {
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

                // Retain only the candidate kmers that are in the kmers hash
                candidate_kmers.retain(|kmer| kmers.contains(kmer));

                if verbosity > 9 {
                    to_print = format!(
                        "{}  There are {} candidate kmers for extension.\n",
                        to_print,
                        candidate_kmers.len()
                    );
                }

                // If there are no candidate kmers, the node is terminal
                if candidate_kmers.is_empty() {
                    if verbosity > 9 {
                        to_print = format!("{}  There are no candidate kmers for extension. Marking node as terminal.\n", to_print);
                    }
                    graph[node].is_terminal = true;
                    graph[node].visited = true;
                    continue 'nodes;
                }

                // Add new nodes if needed, and new edges
                for kmer in candidate_kmers.iter() {
                    let suffix = kmer & suffix_mask;

                    // Check if the node extends by itself and mark it as terminal if it does
                    if suffix == sub_kmer {
                        graph[node].is_terminal = true;
                        graph[node].visited = true;
                        if verbosity > 9 {
                            to_print = format!(
                                "{}  Node {} extends itself. Marking as terminal.\n",
                                to_print,
                                node.index()
                            );
                        }
                        continue 'nodes;
                    }

                    // If the node with sub_kmer == suffix already exists, add an edge to the existing node
                    // Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node

                    let mut node_exists = false;
                    for existing_node in graph.node_indices() {
                        if graph[existing_node].sub_kmer == suffix {
                            if !crate::pcr::would_form_cycle(&graph, node, existing_node) {
                                let edge = crate::pcr::get_dbedge(kmer, kmer_counts, k);
                                graph.add_edge(node, existing_node, edge);

                                if graph[existing_node].is_end {
                                    if verbosity > 5 {
                                        to_print =
                                            format!("{}  New edge added to End node.\n", to_print);
                                    }
                                }
                                let outgoing =
                                    graph.neighbors_directed(node, Direction::Outgoing).count();
                                if outgoing > 4 {
                                    println!("{}",
										format!("WARNING: Node {} has {} outgoing edges. This exceed the maximum of 4 that is expected", node.index(), outgoing)
									);
                                }
                            } else {
                                graph[node].is_terminal = true;

                                if verbosity > 2 {
                                    to_print = format!("{}  Adding edge to node {} would form cycle. Not adding edge, and marking current node as terminal.\n", to_print, node.index());
                                }
                            }

                            node_exists = true;
                            break;
                        }
                    }

                    if !node_exists {

                        // Add the new node, marking it as an end node if it matches the cut2 kmer
                        let new_node: NodeIndex;
                        if (kmer & cut2_mask_terminal) == cut2_kmer_terminal {
                            // kmer ends with cut2, so it is terminal and end
                            new_node = graph.add_node(crate::pcr::DBNode {
                                sub_kmer: suffix,
                                is_start: false,
                                is_end: true,
                                is_terminal: true,
                                visited: true,
                            });
                        } else if (kmer & cut1_mask_terminal) == cut1_kmer_terminal {
                            // kmer ends with cut1, so it is terminal but not end
                            new_node = graph.add_node(crate::pcr::DBNode {
                                sub_kmer: suffix,
                                is_start: false,
                                is_end: false,
                                is_terminal: true,
                                visited: true,
                            });
                        } else {
                            // kmer is neither start nor end
                            new_node = graph.add_node(crate::pcr::DBNode {
                                sub_kmer: suffix,
                                is_start: false,
                                is_end: false,
                                is_terminal: false,
                                visited: false,
                            });
                        }

                        let edge = crate::pcr::get_dbedge(kmer, kmer_counts, k);
                        let edge_count = edge.count;
                        graph.add_edge(node, new_node, edge);
                        let outgoing = graph.neighbors_directed(node, Direction::Outgoing).count();
                        if outgoing > 4 {
                            println!("{}",
								format!("WARNING: Node {} has {} outgoing edges when adding new node. This exceed the maximum of 4 that is expected", node.index(), outgoing).color(COLOR_WARNING)
							);
                        }

                        if verbosity > 9 {
                            
                            to_print = format!(
                                "{}  Added sub_kmer {} for new node {} with edge kmer count {}.\n",
                                to_print,
                                crate::kmer::kmer_to_seq(&suffix, &(*k - 1)),
                                new_node.index(),
                                edge_count
                            );
                        }

                        // Check if the new node is max_length-k+1 from a start node
                        // If so, mark the new node as terminal
                        let path_length = crate::pcr::get_path_length(&graph, new_node);
                        
                        
                        if verbosity > 1 {
                            let is_end = (kmer & cut2_mask_terminal) == cut2_kmer_terminal;
                            if is_end {
                                to_print = format!(
                                    "{}  New node {} is an end node, with path length {}.\n",
                                    to_print,
                                    new_node.index(),
                                    path_length.unwrap_or(0)
                                );
                            }
                        }

                        // If the path length is None, the node is part of a cycle and is marked terminal.
                        // If the path length is Some, is marked terminal if the path length is >= max_length-k+1
                        if let Some(path_length) = path_length {
                            if verbosity > 9 {
                                to_print =
                                    format!("{}  Path length is {}.\n", to_print, path_length);
                            }

                            if path_length > params.max_length - (*k) {
                                graph[new_node].is_terminal = true;
                                graph[new_node].visited = true;
                                if verbosity > 5 {
                                    to_print = format!("{}  Marking new node {} as terminal because it exceeds max_length from start.\n", to_print, new_node.index());
                                }
                            }
                        } else {
                            graph[new_node].is_terminal = true;
                            if verbosity > 5 {
                                to_print = format!("{}  Marking new node {} as terminal because it is part of a cycle.\n", to_print, new_node.index());
                            }
                        }
                    }
                }
                graph[node].visited = true;

                if verbosity > 2 {
                    to_print = format!(
                        "{}There are now {} unvisited and {} non-terminal nodes in the graph.\n",
                        to_print,
                        crate::pcr::n_unvisited_nodes_in_graph(&graph),
                        crate::pcr::n_nonterminal_nodes_in_graph(&graph)
                    );
                }
            }
        }
    }  // End graph extension

    // There should only be a single start node, the one that seeded the graph
    let start_nodes = crate::pcr::get_start_nodes(&graph);
    if start_nodes.len() != 1 {
        panic!("There should be a single start node,");
    }

    to_print = format!(
        "{}  There are {} nodes in the raw graph.\n",
        to_print,
        graph.node_count()
    );
    
    let end_nodes = crate::pcr::get_end_nodes(&graph);
    if end_nodes.is_empty() {
        to_print = format!(
            "{}  No end nodes found for cut1 kmer {}\n",
            to_print,
            crate::kmer::kmer_to_seq(starting_kmer, &(*k - 1))
        );
        if verbosity > 0 {
            print!("{}", to_print);
        }
        return records_vec;
    } else {
        to_print = format!(
            "{}  There are {} end nodes for cut1 kmer {}\n",
            to_print,
            end_nodes.len(),
            crate::kmer::kmer_to_seq(starting_kmer, &(*k - 1))
        );
    }

    // Simplify the graph by removing nodes that are not ancestors of any end node
    crate::pcr::remove_side_branches(&mut graph);

    if end_nodes.len() != crate::pcr::get_end_nodes(&graph).len() {
        panic!("The number of end nodes changed after removing terminal side branches");
    }

    to_print = format!(
        "{}  There are {} nodes in the simplified graph.\n",
        to_print,
        graph.node_count()
    );

    // Look for paths between each combination of start and end nodes
    for start in crate::pcr::get_start_nodes(&graph) {
        let start_seq = crate::kmer::kmer_to_seq(&graph[start].sub_kmer, &(*k - 1));
        'end_nodes: for end in crate::pcr::get_end_nodes(&graph) {

            let end_seq = crate::kmer::kmer_to_seq(&graph[end].sub_kmer, &(*k - 1));

			// Create a copy of the graph
			let mut sub_graph = graph.clone();
			remove_non_ancestors(&mut sub_graph, &end);

            let paths_for_this_pair =
                all_simple_paths::<
                    Vec<NodeIndex>,
                    &StableDiGraph<crate::pcr::DBNode, crate::pcr::DBEdge>,
                >(&sub_graph, start, end, 1, Some(params.max_length - (*k) + 1));

            // A vector of sequences for this pair of start and end nodes
            let mut sequences_for_this_pair: Vec<String> = Vec::new();
            // For each path, get the sequence of the path
            for path in paths_for_this_pair.into_iter() {
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

                if sequence.len() >= params.min_length && sequence.len() <= params.max_length {
                    sequences_for_this_pair.push(sequence);
                }
            }

            // Get a vector of sequence lengths
            let mut sequence_lengths: Vec<usize> = sequences_for_this_pair
                .iter()
                .map(|sequence| sequence.len())
                .collect();

            let sequence_lengths_set: HashSet<usize> = sequence_lengths.iter().cloned().collect();
            let sequence_lengths_mode = mode(&sequence_lengths);
            if sequence_lengths_set.len() == 0 {
                to_print = format!(
                    "{}  No viable sequences found between start node {} and end node {}\n",
                    to_print,
                    start.index(),
                    end.index()
                );
                continue 'end_nodes;
            } else if sequence_lengths_set.len() == 1 {
                to_print = format!(
                    "{}  All {} sequences between start node {} and end node {} are {} bases long\n",
                    to_print,
                    sequences_for_this_pair.len(),
                    start.index(),
                    end.index(),
                    sequence_lengths_mode
                );
            } else {
                to_print = format!(
                    "{}  The {} sequences between start node {} and end node {} have variable length\n    Subsetting to those with mode length {}",
                    to_print,
                    sequences_for_this_pair.len(),
                    start.index(),
                    end.index(),
                    sequence_lengths_mode
                );
                sequences_for_this_pair.retain(|sequence| sequence.len() == sequence_lengths_mode);
            }
            let consensus_sequence = get_consensus(&sequences_for_this_pair);

            // Create a fasta record for this sequence
            let record = fasta::Record::with_attrs(
                format!(
                    "{}_{}-{}",
                    sample_name,
                    start_seq,
                    end_seq
                ).as_str(),
                None,
                &consensus_sequence.as_bytes(),
            );
            records_vec.push(record);
        }
    }

    records_vec
}

pub struct RADParams {
    pub cut1: String,
    pub cut2: String,
    pub min_length: usize,
    pub max_length: usize,
    pub name: String,
    pub coverage: u64,
}

pub fn do_rad(
    kmer_counts_map: &FxHashMap<u64, u64>,
    k: &usize,
    sample_name: &str,
    verbosity: usize,
    params: &RADParams,
) -> Vec<bio::io::fasta::Record> {
    //let mut records: Vec<fasta::Record> = Vec::new();

    println!(
        "{}",
        format!("Running RAD analysis {}", params.name).color(COLOR_NOTE)
    );

    // Remove kmer_counts entries with less than coverage
    println!(
        "Removing kmers with coverage less than {}...",
        params.coverage
    );

    // The kmer_counts_map includes only canonical kmers. Add the reverse complements,
    // while also filtering for coverage
    let mut kmer_counts: FxHashMap<u64, u64> = FxHashMap::default();
    let mut n_unique_kmers: u64 = 0;
    for (kmer, count) in kmer_counts_map {
        if count >= &params.coverage {
            kmer_counts.insert(*kmer, *count);
            kmer_counts.insert(crate::kmer::revcomp_kmer(kmer, k), *count);
            n_unique_kmers += 1;
        }
    }

    // Create a hash set of the keys of kmer_counts
    let kmers: std::collections::HashSet<u64> = kmer_counts.keys().copied().collect();

    println!(
        "  The number of unique kmers went from {} to {}",
        kmer_counts_map.len(),
        n_unique_kmers
    );

    let mut cut1_mask: u64 = 0;
    for _i in 0..(2 * params.cut1.len()) {
        cut1_mask = (cut1_mask << 1) | 1;
    }
    cut1_mask <<= 2 * *k - 2 * params.cut1.len();

    // Get the cut1 kmer and rotate it to the start of the kmer
    let mut cut1_kmer = crate::pcr::string_to_oligo(&params.cut1).kmer;
    cut1_kmer <<= 2 * *k - 2 * params.cut1.len();

    // Get all the kmers that start with cut1
    let mut cut1_kmers: Vec<u64> = Vec::new();
    for kmer in kmer_counts.keys() {
        if kmer & cut1_mask == cut1_kmer {
            cut1_kmers.push(*kmer);
        }
    }

    println!("  Found {} kmers that start with cut1", cut1_kmers.len());

    // Loop over all the cut1 kmers and try to extend them as a de bruijn graph.
    // Check if each extension matches the cut2 kmer. If it does, and the length
    // is less than the min length, continue to next cut1 kmer. If it is longer
    // than min length, add it to the list of records. If you get to max length,
    // stop and go to next cut1 kmer.

    let records: Vec<fasta::Record> = cut1_kmers
        .par_iter()
        .map(|starting_kmer| {
            generate_fastas(
                params,
                sample_name,
                &kmer_counts,
                &kmers,
                starting_kmer,
                k,
                verbosity,
            )
        })
        .flatten() // This flattens Vec<fasta::Record> into an iterator of fasta::Record
        .collect(); // Collect all records into a Vec

    records
}
