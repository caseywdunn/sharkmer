use anyhow::{bail, Context, Result};

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

/// Check whether a string looks like a URL.
pub fn is_url(source: &str) -> bool {
    source.starts_with("http://") || source.starts_with("https://")
}

/// Load a panel from a user-supplied YAML file path.
fn load_panel_file(path: &str) -> Result<Vec<PCRParams>> {
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

/// Load a panel from a URL.
fn load_panel_url(url: &str) -> Result<Vec<PCRParams>> {
    log::info!("Downloading primer panel from {}", url);
    let response = ureq::get(url)
        .call()
        .with_context(|| format!("Failed to download panel from URL: {}", url))?;

    let yaml_str = response
        .into_string()
        .with_context(|| format!("Failed to read panel response from URL: {}", url))?;

    let mut panel = parse_panel_yaml(&yaml_str)
        .with_context(|| format!("Failed to parse panel YAML from URL: {}", url))?;

    // Prepend panel name to gene names
    for param in panel.primers.iter_mut() {
        param.gene_name = format!("{}_{}", panel.name, param.gene_name);
    }

    Ok(panel.primers)
}

/// Load a panel from a file path or URL.
pub fn load_panel_source(source: &str) -> Result<Vec<PCRParams>> {
    if is_url(source) {
        load_panel_url(source)
    } else {
        load_panel_file(source)
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_url() {
        assert!(is_url("https://example.com/panel.yaml"));
        assert!(is_url("http://example.com/panel.yaml"));
        assert!(!is_url("/path/to/panel.yaml"));
        assert!(!is_url("relative/panel.yaml"));
        assert!(!is_url("panel.yaml"));
    }

    #[test]
    fn test_load_panel_file_from_fixture() {
        let params = load_panel_file("tests/fixtures/test_panel.yaml").unwrap();
        assert_eq!(params.len(), 1);
        assert_eq!(params[0].gene_name, "test_panel_18S");
    }

    #[test]
    fn test_load_panel_source_local_file() {
        let params = load_panel_source("tests/fixtures/test_panel.yaml").unwrap();
        assert_eq!(params.len(), 1);
        assert_eq!(params[0].gene_name, "test_panel_18S");
    }

    #[test]
    fn test_load_panel_source_bad_url() {
        let result = load_panel_source("https://localhost:1/nonexistent_panel.yaml");
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("Failed to download panel from URL"));
    }
}

pub fn print_pcr_panels() {
    let panels = get_preconfigured_panels();
    println!("Available PCR panels (use --export-panel <name> for details):\n");
    for panel in &panels {
        let n = panel.primers.len();
        let noun = if n == 1 { "primer" } else { "primers" };
        println!(
            "  {:<16} {} ({} {})",
            panel.name, panel.description, n, noun
        );
    }
}
