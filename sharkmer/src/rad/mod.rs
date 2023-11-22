use rustc_hash::FxHashMap;
use bio::io::fasta;
use colored::*;
use petgraph::stable_graph::StableDiGraph;
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::algo::all_simple_paths;
use std::collections::HashSet;
use std::collections::HashMap;
use std::io::Write;
use rayon::prelude::*;



use crate::COLOR_FAIL;
use crate::COLOR_NOTE;
use crate::COLOR_SUCCESS;
use crate::COLOR_WARNING;

const MAX_NUM_NODES: usize = 5_000;

fn generate_fastas(
	params: &RADParams,
	sample_name: &str,
	kmer_counts: &FxHashMap<u64, u64>,
	kmers: &std::collections::HashSet<u64>, 
	starting_kmer: &u64, 
	cut2_kmer: &u64,
	k: &usize, 
	verbosity: usize
) -> Vec<fasta::Record> {
	let mut to_print = String::new();
	let mut records_vec: Vec<fasta::Record> = Vec::new();
	let prefix = starting_kmer >> 2;
	let starting_seq = crate::kmer::kmer_to_seq(&starting_kmer, &k);
	let suffix_mask: u64 = (1 << (2 * (*k - 1))) - 1;

	let mut cut2_mask: u64 = 0;
	for _i in 0..(2 * params.cut2.len()) {
		cut2_mask = (cut2_mask << 1) | 1;
	}

	if verbosity > 0 {
		println!("{} Starting kmer analysis", starting_seq);
	}

	let mut graph: StableDiGraph<crate::pcr::DBNode, crate::pcr::DBEdge> = StableDiGraph::new();
	graph.add_node(crate::pcr::DBNode {
		sub_kmer: prefix,
		is_start: true,
		is_end: false,
		is_terminal: false,
		visited: false, // Will not be extended
	});
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
		for node in node_indices {
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
					print!(
						"There are {} candidate kmers for extension. ",
						candidate_kmers.len()
					);
					std::io::stdout().flush().unwrap();
				}

				// If there are no candidate kmers, the node is terminal
				if candidate_kmers.is_empty() {
					if verbosity > 9 {
						println!("Marking node as terminal because there are no candidates for extension. ");
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
						if verbosity > 9 {
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
							if !crate::pcr::would_form_cycle(&graph, node, existing_node) {
								let edge = crate::pcr::get_dbedge(kmer, &kmer_counts, k);
								graph.add_edge(node, existing_node, edge);
								
								if graph[existing_node].is_end {
									if verbosity > 5 {
										println!("New edge added to End node.");
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
						let edge = crate::pcr::get_dbedge(kmer, &kmer_counts, k);
						let edge_count = edge.count;
						
						// Add the new node, marking it as an end node if it matches the cut2 kmer
						let is_end = kmer & cut2_mask == *cut2_kmer;
						let new_node = graph.add_node(crate::pcr::DBNode {
							sub_kmer: suffix,
							is_start: false,
							is_end: is_end,
							is_terminal: is_end,
							visited: is_end,
						});

						graph.add_edge(node, new_node, edge);
						let outgoing = graph.neighbors_directed(node, Direction::Outgoing).count();
						if outgoing > 4 {
							println!("{}",
								format!("WARNING: Node {} has {} outgoing edges when adding new node. This exceed the maximum of 4 that is expected", node.index(), outgoing).color(COLOR_WARNING)
							);
						}

						if verbosity > 9 {
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
						let path_length = crate::pcr::get_path_length(&graph, new_node);

						if is_end {
							if verbosity > 1 {
								println!("  New node {} is an end node, with path length {}. ", new_node.index(), path_length.unwrap_or(0));
								std::io::stdout().flush().unwrap();
							}
						}

						// If the path length is None, the node is part of a cycle and is marked terminal.
						// If the path length is Some, is marked terminal if the path length is >= max_length-k+1
						if let Some(path_length) = path_length {
							if verbosity > 9 {
								print!("Path length is {}. ", path_length);
								std::io::stdout().flush().unwrap();
							}

							if path_length > params.max_length - (*k) {
								graph[new_node].is_terminal = true;
								graph[new_node].visited = true;
								if verbosity > 5 {
									print!("Marking new node {} as terminal because it exceeds max_length from start. ", new_node.index());
									std::io::stdout().flush().unwrap();
								}
							}
						} else {
							graph[new_node].is_terminal = true;
							if verbosity > 5 {
								print!("Marking new node {} as terminal because it is part of a cycle. ", new_node.index());
								std::io::stdout().flush().unwrap();
							}
						}
					}
				}
				graph[node].visited = true;

				if verbosity > 2 {
					println!(
						"There are now {} unvisited and {} non-terminal nodes in the graph. ",
						crate::pcr::n_unvisited_nodes_in_graph(&graph),
						crate::pcr::n_nonterminal_nodes_in_graph(&graph)
					);
					std::io::stdout().flush().unwrap();
				}
			}
		}
	}
	// Make hashmaps of the start and end nodes, where the key is the node index and the value is the number of edges
	let mut start_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();
	let mut end_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();

	if end_nodes_map.is_empty() {
		to_print = format!("{}  No end nodes found for cut1 kmer {}\n", to_print, crate::kmer::kmer_to_seq(&starting_kmer, &(*k - 1)));
		if verbosity > 0 {
			print!("{}", to_print);
		}
		return records_vec;
	}

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

	// Drop entries that have no edges
	start_nodes_map.retain(|_node, edge_count| *edge_count > 0);
	end_nodes_map.retain(|_node, edge_count| *edge_count > 0);

	println!("Final extension statistics:");
	crate::pcr::summarize_extension(&graph, "    ");
	println!(
		"  There are {} start nodes with edges",
		start_nodes_map.len()
	);
	println!("  There are {} end nodes with edges", end_nodes_map.len());	

	// Simplify the graph
	// all_simple_paths() hangs if the input graph is too complex
	// Also want to regularize some graph features

	let start = std::time::Instant::now();
	println!("Pruning the assembly graph...");

	// Iteratively remove nodes that do not have outgoing edges and are not end nodes
	// These are terminal side branches.
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

	// Get all paths from start nodes to terminal nodes
	let start = std::time::Instant::now();
	println!("Traversing the assembly graph to find paths from forward to reverse primers...");
	let mut all_paths = Vec::new();
	

	for start in crate::pcr::get_start_nodes(&graph) {
		for end in crate::pcr::get_end_nodes(&graph) {
			let paths_for_this_pair = all_simple_paths::<
				Vec<NodeIndex>,
				&StableDiGraph<crate::pcr::DBNode, crate::pcr::DBEdge>,
			>(
				&graph, start, end, 1, Some(params.max_length - (*k) + 1)
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
		to_print = format!("{}  No path found\n", to_print);
		if verbosity > 0 {
			print!("{}", to_print);
		}
		return records_vec;
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
			let subread = crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(*k - 1));
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

		if sequence.len() < params.min_length {
			to_print = format!("{}  RAD product {} is too short ({} bases). Skipping.\n", to_print, i, sequence.len());
			continue;
		}

		// Get some stats on the path counts
		let count_mean = crate::pcr::compute_mean(&edge_counts);
		let count_median = crate::pcr::compute_median(&edge_counts);
		let count_min = edge_counts.iter().min().unwrap();
		let count_max = edge_counts.iter().max().unwrap();

		let id = format!(
			"{} {} start {} product {} length {} kmer count stats mean {:.2} median {} min {} max {}",
			sample_name,
			params.name,
			starting_seq,
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
		records_vec.push(record);
	}
	return records_vec;
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
	let suffix_mask: u64 = (1 << (2 * (*k - 1))) - 1;

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
	

	// Get the cut2 kmer
	let cut2_kmer = crate::pcr::string_to_oligo(&params.cut2).kmer;

	// Get all the kmers that start with cut1
	let mut cut1_kmers: Vec<u64> = Vec::new();
	for (kmer, _count) in &kmer_counts {
		if kmer & cut1_mask == cut1_kmer {
			cut1_kmers.push(*kmer);
		}
	}

	print!("  Found {} kmers that start with cut1", cut1_kmers.len());

	// Loop over all the cut1 kmers and try to extend them as a de bruijn graph.
	// Check if each extension matches the cut2 kmer. If it does, and the length
	// is less than the min length, continue to next cut1 kmer. If it is longer
	// than min length, add it to the list of records. If you get to max length,
	// stop and go to next cut1 kmer.

	let records: Vec<fasta::Record> = cut1_kmers.par_iter()
		.map(|starting_kmer| generate_fastas(params, sample_name, &kmer_counts, &kmers, starting_kmer, &cut2_kmer, k, verbosity) )
		.flatten() // This flattens Vec<fasta::Record> into an iterator of fasta::Record
		.collect();  // Collect all records into a Vec

	records
}