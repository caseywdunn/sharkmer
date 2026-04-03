// pcr/components.rs — Seed component identification, prioritization, and budget allocation

use petgraph::graph::NodeIndex;

/// A connected group of seed nodes identified during seed evaluation.
#[allow(dead_code)] // Fields populated in Step 2
pub(super) struct SeedComponent {
    /// Forward (start) seed node indices in this component
    pub start_seeds: Vec<NodeIndex>,
    /// Reverse (end) seed node indices in this component
    pub end_seeds: Vec<NodeIndex>,
    /// Whether seed eval found a forward seed that reached a reverse seed
    pub has_connected_seeds: bool,
    /// Average branching ratio across seeds in this component (lower is better)
    pub avg_branching_ratio: f64,
    /// Sum of seed eval node counts across seeds in this component
    pub total_extension_nodes: usize,
    /// Allocated node budget for extension
    pub node_budget: usize,
    /// Priority score for ordering (higher = extend first)
    pub priority_score: f64,
}

impl SeedComponent {
    /// All seed node indices in this component (start + end).
    #[allow(dead_code)] // Used in Step 3
    pub fn all_seeds(&self) -> Vec<NodeIndex> {
        let mut seeds = self.start_seeds.clone();
        seeds.extend_from_slice(&self.end_seeds);
        seeds
    }
}
