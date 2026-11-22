//! Read-aware bubble detection with owned target-local buffers.
use super::memory::buffers::BudgetSet;
use super::threading::{BranchLink, ThreadingAnnotations};
use super::{DBEdge, DBNode};
use crate::kmer::budget::ResourceBudget;
use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use anyhow::{Context, Result};
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;
use std::sync::Arc;
#[cfg(test)]
mod reference;
#[cfg(test)]
use ahash::AHashMap;

const MAX_BUBBLE_DEPTH: usize = 50;

struct BranchRanking {
    edges: BudgetVec<EdgeIndex>,
    support: u32,
    key: u64,
    ordinal: usize,
}

#[cfg(test)]
pub fn resolve_bubbles(
    graph: &StableDiGraph<DBNode, DBEdge>,
    annotations: &ThreadingAnnotations,
) -> AHashMap<EdgeIndex, f64> {
    resolve_bubbles_budgeted(graph, annotations, None)
        .unwrap()
        .into_iter()
        .collect()
}

/// Resolve the existing simple-bubble policy under the graph's target ceiling.
pub(super) fn resolve_bubbles_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    annotations: &ThreadingAnnotations,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetMap<EdgeIndex, f64>> {
    let mut preferences = BudgetMap::new(budget);
    let bubbles = detect_simple_bubbles_budgeted(graph, budget)?;
    for bubble in bubbles.iter() {
        let rankings = rank_branches(graph, bubble, annotations, budget)?;
        if rankings.len() >= 2 {
            let max_support = rankings.iter().map(|r| r.support).max().unwrap_or(0);
            for ranking in rankings.iter() {
                let preference = if max_support > 0 {
                    ranking.support as f64 / max_support as f64
                } else {
                    1.0
                };
                for &edge in ranking.edges.iter() {
                    preferences.insert(edge, preference)?;
                }
            }
        }
    }
    Ok(preferences)
}

struct Bubble {
    source: NodeIndex,
    #[allow(dead_code)]
    sink: NodeIndex,
    branches: BudgetVec<BudgetVec<EdgeIndex>>,
}

#[cfg(test)]
fn detect_simple_bubbles(graph: &StableDiGraph<DBNode, DBEdge>) -> BudgetVec<Bubble> {
    detect_simple_bubbles_budgeted(graph, None).unwrap()
}

fn detect_simple_bubbles_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<Bubble>> {
    let mut bubbles = BudgetVec::new(budget);
    for source in graph.node_indices() {
        let outgoing =
            BudgetVec::collect(graph.edges_directed(source, Direction::Outgoing), budget)?;
        if outgoing.len() < 2 {
            continue;
        }
        let mut endpoints: BudgetMap<NodeIndex, BudgetVec<BudgetVec<EdgeIndex>>> =
            BudgetMap::new(budget);
        for edge in outgoing.iter() {
            let mut current = edge.target();
            let mut branch_edges = BudgetVec::new(budget);
            branch_edges.push(edge.id())?;
            let mut visited = BudgetSet::new(budget);
            visited.insert(source)?;
            let mut depth = 0;
            let mut terminated_naturally = false;
            loop {
                if depth >= MAX_BUBBLE_DEPTH {
                    break;
                }
                depth += 1;
                if visited.contains(&current) {
                    break;
                }
                visited.insert(current)?;
                // Counting up to two edges avoids allocating a temporary list
                // for every step of each bounded linear branch.
                let mut next = graph.edges_directed(current, Direction::Outgoing);
                match (next.next(), next.next()) {
                    (Some(edge), None) => {
                        branch_edges.push(edge.id())?;
                        current = edge.target();
                    }
                    _ => {
                        terminated_naturally = true;
                        break;
                    }
                }
            }
            // Cycles and depth-limited walks are not convergence evidence.
            if terminated_naturally {
                if endpoints.get(&current).is_none() {
                    endpoints.insert(current, BudgetVec::new(budget))?;
                }
                endpoints
                    .get_mut(&current)
                    .context("Missing bubble endpoint after insertion")?
                    .push(branch_edges)?;
            }
        }
        for (sink, branches) in endpoints {
            if branches.len() >= 2 && sink != source {
                bubbles.push(Bubble {
                    source,
                    sink,
                    branches,
                })?;
            }
        }
    }
    Ok(bubbles)
}

fn rank_branches(
    graph: &StableDiGraph<DBNode, DBEdge>,
    bubble: &Bubble,
    annotations: &ThreadingAnnotations,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<BranchRanking>> {
    let mut rankings = BudgetVec::new(budget);
    for (ordinal, branch) in bubble.branches.iter().enumerate() {
        let mut total = 0u32;
        for edge in branch.iter() {
            total = total
                .checked_add(
                    annotations
                        .edge_support
                        .get(edge)
                        .map_or(0, |s| s.read_support_total),
                )
                .context("Bubble read support overflow")?;
        }
        let mut phasing = 0u32;
        if let Some(&first) = branch.first() {
            for incoming in graph.edges_directed(bubble.source, Direction::Incoming) {
                let link = BranchLink {
                    incoming_edge: incoming.id(),
                    outgoing_edge: first,
                };
                phasing = phasing
                    .checked_add(annotations.branch_links.get(&link).copied().unwrap_or(0))
                    .context("Bubble phasing support overflow")?;
            }
        }
        let support = total
            .checked_add(phasing)
            .context("Combined bubble support overflow")?;
        rankings.push(BranchRanking {
            edges: branch.try_clone()?,
            support,
            key: branch_sort_key(graph, branch)?,
            ordinal,
        })?;
    }
    // The ordinal preserves stable order even in constructed graphs with
    // duplicate directional keys. In-place sorting needs no hidden scratch.
    rankings.sort_unstable_by(|a, b| {
        b.support
            .cmp(&a.support)
            .then(a.key.cmp(&b.key))
            .then(a.ordinal.cmp(&b.ordinal))
    });
    Ok(rankings)
}

fn branch_sort_key(graph: &StableDiGraph<DBNode, DBEdge>, edges: &[EdgeIndex]) -> Result<u64> {
    let Some(&first) = edges.first() else {
        return Ok(0);
    };
    let (source, target) = graph
        .edge_endpoints(first)
        .context("Missing bubble branch edge")?;
    Ok((graph[source].sub_kmer << 2) | (graph[target].sub_kmer & 3))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bubble_graph() -> StableDiGraph<DBNode, DBEdge> {
        // Create a diamond bubble:
        //   n0 --e0--> n1 --e2--> n3
        //   n0 --e1--> n2 --e3--> n3
        let mut graph = StableDiGraph::new();
        let n0 = graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: false,
        });
        let n2 = graph.add_node(DBNode {
            sub_kmer: 2,
            is_start: false,
            is_end: false,
        });
        let n3 = graph.add_node(DBNode {
            sub_kmer: 3,
            is_start: false,
            is_end: true,
        });

        graph.add_edge(
            n0,
            n1,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n0,
            n2,
            DBEdge {
                count: 3,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n1,
            n3,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n2,
            n3,
            DBEdge {
                count: 3,
                coverage_ratio: 1.0,
            },
        );

        graph
    }

    #[test]
    fn test_detect_simple_bubble() {
        let graph = make_bubble_graph();
        let bubbles = detect_simple_bubbles(&graph);
        assert_eq!(bubbles.len(), 1);
        assert_eq!(bubbles[0].branches.len(), 2);
    }

    #[test]
    fn test_resolve_with_read_support() {
        let graph = make_bubble_graph();

        let mut annotations = ThreadingAnnotations::new(None);

        // Give branch 1 (through n1) more read support
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();
        let e2 = graph.edge_indices().nth(2).unwrap();
        let e3 = graph.edge_indices().nth(3).unwrap();

        annotations
            .edge_support
            .insert(
                e0,
                super::super::threading::EdgeReadSupport {
                    read_support_total: 10,
                    read_support_unambiguous: 8,
                },
            )
            .unwrap();
        annotations
            .edge_support
            .insert(
                e2,
                super::super::threading::EdgeReadSupport {
                    read_support_total: 10,
                    read_support_unambiguous: 8,
                },
            )
            .unwrap();
        annotations
            .edge_support
            .insert(
                e1,
                super::super::threading::EdgeReadSupport {
                    read_support_total: 2,
                    read_support_unambiguous: 1,
                },
            )
            .unwrap();
        annotations
            .edge_support
            .insert(
                e3,
                super::super::threading::EdgeReadSupport {
                    read_support_total: 2,
                    read_support_unambiguous: 1,
                },
            )
            .unwrap();

        let prefs = resolve_bubbles(&graph, &annotations);

        // Branch through n1 (e0, e2) should have higher preference
        assert!(prefs[&e0] > prefs[&e1]);
        assert!(prefs[&e2] > prefs[&e3]);
    }

    #[test]
    fn test_no_bubble_linear() {
        // Linear graph has no bubbles
        let mut graph = StableDiGraph::new();
        let n0 = graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: true,
        });
        graph.add_edge(
            n0,
            n1,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );

        let bubbles = detect_simple_bubbles(&graph);
        assert!(bubbles.is_empty());
    }

    /// Two branches that never converge within MAX_BUBBLE_DEPTH must not be
    /// recorded as a bubble, even if their depth-limited endpoints happen
    /// to share a node. The old code treated the depth-cutoff position as
    /// a sink, which could fabricate bubbles from unrelated long paths.
    #[test]
    fn test_depth_limited_branches_not_a_bubble() {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        // Two long linear branches from a shared source. Neither ever
        // reaches a real convergence point — they just run off into
        // the distance. MAX_BUBBLE_DEPTH is 50, so we make each branch
        // 60 nodes long.
        let source = graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
        });
        let mk_chain = |graph: &mut StableDiGraph<DBNode, DBEdge>, base: u64| -> NodeIndex {
            let mut prev = source;
            for i in 0..60u64 {
                let n = graph.add_node(DBNode {
                    sub_kmer: base + i + 1,
                    is_start: false,
                    is_end: false,
                });
                graph.add_edge(
                    prev,
                    n,
                    DBEdge {
                        count: 5,
                        coverage_ratio: 1.0,
                    },
                );
                prev = n;
            }
            prev
        };
        let _ = mk_chain(&mut graph, 1000);
        let _ = mk_chain(&mut graph, 2000);

        let bubbles = detect_simple_bubbles(&graph);
        assert!(
            bubbles.is_empty(),
            "depth-limited branches must not be classified as a bubble"
        );
    }

    /// Tiebreaking must produce deterministic ranking when two branches
    /// have equal read support. The winner must be independent of the
    /// bubble's internal HashMap iteration order across runs.
    #[test]
    fn test_bubble_tiebreak_is_deterministic() {
        let graph = make_bubble_graph();
        // Both branches get identical support — the tie is resolved by
        // the content-derived branch_sort_key.
        let mut annotations = ThreadingAnnotations::new(None);
        for e in graph.edge_indices() {
            annotations
                .edge_support
                .insert(
                    e,
                    super::super::threading::EdgeReadSupport {
                        read_support_total: 10,
                        read_support_unambiguous: 5,
                    },
                )
                .unwrap();
        }

        // Resolve repeatedly; the full preference map must be bit-identical
        // across calls. If the tiebreak leaked HashMap iteration order, the
        // scores assigned to individual edges would vary between calls.
        let reference: Vec<(EdgeIndex, u64)> = {
            let prefs = resolve_bubbles(&graph, &annotations);
            let mut pairs: Vec<(EdgeIndex, u64)> =
                prefs.into_iter().map(|(e, s)| (e, s.to_bits())).collect();
            pairs.sort_by_key(|(e, _)| e.index());
            pairs
        };
        for _ in 0..10 {
            let prefs = resolve_bubbles(&graph, &annotations);
            let mut pairs: Vec<(EdgeIndex, u64)> =
                prefs.into_iter().map(|(e, s)| (e, s.to_bits())).collect();
            pairs.sort_by_key(|(e, _)| e.index());
            assert_eq!(pairs, reference, "tiebreak must be deterministic");
        }
    }
    fn bubble_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn owned_bubble_preferences_match_frozen_oracle_with_holes_ties_and_extra_branches() {
        for case in 0..48u64 {
            let mut graph = make_bubble_graph();
            let source = NodeIndex::new(0);
            let sink = NodeIndex::new(3);
            let hole = graph.add_node(DBNode {
                sub_kmer: 99,
                is_start: false,
                is_end: false,
            });
            graph.remove_node(hole);
            for extra in 0..case % 5 {
                let node = graph.add_node(DBNode {
                    sub_kmer: extra + 4,
                    is_start: false,
                    is_end: false,
                });
                graph.add_edge(
                    source,
                    node,
                    DBEdge {
                        count: 3,
                        coverage_ratio: 1.0,
                    },
                );
                graph.add_edge(
                    node,
                    sink,
                    DBEdge {
                        count: 3,
                        coverage_ratio: 1.0,
                    },
                );
            }
            let prefix = graph.add_node(DBNode {
                sub_kmer: 17,
                is_start: false,
                is_end: false,
            });
            let incoming = graph.add_edge(
                prefix,
                source,
                DBEdge {
                    count: 1,
                    coverage_ratio: 1.0,
                },
            );
            if case % 3 == 0 {
                graph.add_edge(
                    sink,
                    source,
                    DBEdge {
                        count: 1,
                        coverage_ratio: 1.0,
                    },
                );
            }
            let mut annotations = ThreadingAnnotations::new(None);
            let mut old = crate::pcr::threading::reference::ThreadingAnnotations {
                edge_support: AHashMap::new(),
                branch_links: AHashMap::new(),
                paired_links: Vec::new(),
            };
            for edge in graph.edge_indices() {
                let total = ((edge.index() as u64 * 7 + case) % 9) as u32;
                annotations
                    .edge_support
                    .insert(
                        edge,
                        crate::pcr::threading::EdgeReadSupport {
                            read_support_total: total,
                            read_support_unambiguous: 0,
                        },
                    )
                    .unwrap();
                old.edge_support.insert(
                    edge,
                    crate::pcr::threading::reference::EdgeReadSupport {
                        read_support_total: total,
                        read_support_unambiguous: 0,
                    },
                );
            }
            for edge in graph.edges_directed(source, Direction::Outgoing) {
                let count = (case % 4) as u32;
                annotations
                    .branch_links
                    .insert(
                        BranchLink {
                            incoming_edge: incoming,
                            outgoing_edge: edge.id(),
                        },
                        count,
                    )
                    .unwrap();
                old.branch_links.insert(
                    crate::pcr::threading::reference::BranchLink {
                        incoming_edge: incoming,
                        outgoing_edge: edge.id(),
                    },
                    count,
                );
            }
            let budget = bubble_budget(1 << 20);
            let actual: AHashMap<_, _> =
                resolve_bubbles_budgeted(&graph, &annotations, Some(&budget))
                    .unwrap()
                    .into_iter()
                    .collect();
            assert_eq!(
                actual,
                reference::resolve_bubbles(&graph, &old),
                "case {case}"
            );
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn bubble_refusal_and_support_overflow_drop_workspace_without_mutating_the_graph() {
        let graph = make_bubble_graph();
        let mut annotations = ThreadingAnnotations::new(None);
        let (mut refused, mut succeeded) = (0, 0);
        for bytes in (1..=80).map(|n| n * 64).chain([1 << 20]) {
            let budget = bubble_budget(bytes);
            match resolve_bubbles_budgeted(&graph, &annotations, Some(&budget)) {
                Ok(preferences) => {
                    assert_eq!(preferences.len(), 4);
                    succeeded += 1;
                }
                Err(error) => {
                    assert!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some(),
                        "{error:#}"
                    );
                    refused += 1;
                }
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(graph.node_count(), 4);
            assert_eq!(graph.edge_count(), 4);
        }
        assert!(refused > 4 && succeeded > 0);
        for edge in graph.edge_indices() {
            annotations
                .edge_support
                .insert(
                    edge,
                    crate::pcr::threading::EdgeReadSupport {
                        read_support_total: u32::MAX,
                        read_support_unambiguous: 0,
                    },
                )
                .unwrap();
        }
        let budget = bubble_budget(1 << 20);
        let error = resolve_bubbles_budgeted(&graph, &annotations, Some(&budget))
            .err()
            .unwrap();
        assert!(error.to_string().contains("support overflow"));
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_none()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}
