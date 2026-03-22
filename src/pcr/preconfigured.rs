use anyhow::{bail, Context, Result};

use super::pcrparams_string;
use super::PCRParams;

#[derive(serde::Deserialize)]
struct PanelFile {
    name: String,
    description: String,
    primers: Vec<PCRParams>,
}

/// Parse a YAML string into a PanelFile.
fn parse_panel_yaml(yaml_str: &str) -> Result<PanelFile> {
    serde_yaml::from_str(yaml_str).context("Failed to parse panel YAML")
}

/// Load a panel from a user-supplied YAML file path.
pub fn load_panel_file(path: &str) -> Result<Vec<PCRParams>> {
    let yaml_str = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read panel file: {}", path))?;
    let mut panel = parse_panel_yaml(&yaml_str)
        .with_context(|| format!("Failed to parse panel file: {}", path))?;

    // Prepend panel name to gene names
    for param in panel.primers.iter_mut() {
        param.gene_name = format!("{}_{}", panel.name, param.gene_name);
    }

    Ok(panel.primers)
}

/// Returns (name, raw_yaml) pairs for all built-in panels.
fn get_builtin_panel_sources() -> Vec<(&'static str, &'static str)> {
    vec![
        (
            "angiospermae",
            include_str!("../../panels/angiospermae.yaml"),
        ),
        ("bacteria", include_str!("../../panels/bacteria.yaml")),
        ("cnidaria", include_str!("../../panels/cnidaria.yaml")),
        ("human", include_str!("../../panels/human.yaml")),
        ("insecta", include_str!("../../panels/insecta.yaml")),
        ("metazoa", include_str!("../../panels/metazoa.yaml")),
        ("teleostei", include_str!("../../panels/teleostei.yaml")),
    ]
}

fn get_builtin_panels() -> Vec<PanelFile> {
    get_builtin_panel_sources()
        .iter()
        .map(|(_, yaml)| {
            parse_panel_yaml(yaml).expect("Built-in panel YAML should always be valid")
        })
        .collect()
}

fn get_preconfigured_panels() -> Vec<PanelFile> {
    let mut panels = get_builtin_panels();

    // Prepend panel name to gene names
    for panel in panels.iter_mut() {
        for param in panel.primers.iter_mut() {
            param.gene_name = format!("{}_{}", panel.name, param.gene_name);
        }
    }

    panels.sort_by_key(|panel| panel.name.clone());

    panels
}

pub fn get_panel(panel_name: &str) -> Result<Vec<PCRParams>, String> {
    let panels = get_preconfigured_panels();

    panels
        .iter()
        .find(|panel| panel.name == panel_name)
        .map(|panel| panel.primers.clone())
        .ok_or_else(|| {
            let available: Vec<&str> = panels.iter().map(|p| p.name.as_str()).collect();
            format!(
                "Unknown panel '{}'. Available panels: {}",
                panel_name,
                available.join(", ")
            )
        })
}

/// Export a built-in panel as raw YAML to stdout.
pub fn export_panel_yaml(panel_name: &str) -> Result<String> {
    for (name, yaml) in get_builtin_panel_sources() {
        if name == panel_name {
            return Ok(yaml.to_string());
        }
    }

    let available: Vec<&str> = get_builtin_panel_sources()
        .iter()
        .map(|(name, _)| *name)
        .collect();
    bail!(
        "Unknown panel '{}'. Available panels: {}",
        panel_name,
        available.join(", ")
    )
}

pub fn print_pcr_panels() {
    let panels = get_preconfigured_panels();
    println!("Available preconfigured PCR panels:\n");
    for panel in panels {
        let prefix = format!("{}_", panel.name);
        let gene_names: Vec<&str> = panel
            .primers
            .iter()
            .map(|p| p.gene_name.strip_prefix(&prefix).unwrap_or(&p.gene_name))
            .collect();
        println!("  {:<16} {} [{}]", panel.name, panel.description, gene_names.join(", "));
    }
    println!("\nUse --export-panel <name> to see full details for a panel.");
}
