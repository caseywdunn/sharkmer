use anyhow::{Context, Result, bail};
use std::collections::BTreeMap;

use super::PCRParams;

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct PanelFile {
    name: String,
    /// Panel semver version. Optional during migration; populated for all in-tree panels.
    #[serde(default)]
    version: Option<String>,
    description: String,
    #[serde(default)]
    #[allow(dead_code)]
    maintainers: Vec<Maintainer>,
    #[serde(default)]
    #[allow(dead_code)]
    changelog: Vec<ChangelogEntry>,
    primers: Vec<PCRParams>,
    /// Gold-standard reference amplicon sequences for validation. Not used by
    /// the Rust pipeline; consumed by the Python validation tooling.
    #[serde(default)]
    #[allow(dead_code)]
    references: Option<Vec<ReferenceGene>>,
    #[serde(default)]
    #[allow(dead_code)]
    validation: Option<ValidationBlock>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct Maintainer {
    name: String,
    #[serde(default)]
    orcid: Option<String>,
    #[serde(default)]
    contact: Option<String>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ChangelogEntry {
    version: String,
    date: String,
    #[serde(default)]
    sharkmer_version: Option<String>,
    changes: String,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ReferenceGene {
    gene_name: String,
    sequences: Vec<ReferenceSequence>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ReferenceSequence {
    taxon: String,
    #[serde(default)]
    accession: Option<String>,
    sequence: String,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ValidationBlock {
    #[serde(default)]
    last_validated: Option<LastValidated>,
    #[serde(default)]
    samples: Vec<ValidationSample>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct LastValidated {
    sharkmer_version: String,
    panel_version: String,
    date: String,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ValidationSample {
    accession: String,
    #[serde(default)]
    taxon: Option<String>,
    #[serde(default)]
    taxonomy: Option<String>,
    max_reads: Vec<u64>,
    /// Deprecated: validation thresholds formerly written by `--write`.
    /// Kept for backward compatibility with existing panel YAMLs.
    #[serde(default)]
    expected: BTreeMap<String, ExpectedGene>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ExpectedGene {
    #[serde(default)]
    min_identity: Option<f64>,
    #[serde(default)]
    length: Option<usize>,
    #[serde(default)]
    length_tolerance: Option<usize>,
}

/// Parse a YAML string into a PanelFile.
fn parse_panel_yaml(yaml_str: &str) -> Result<PanelFile> {
    serde_yaml_ng::from_str(yaml_str).context("Failed to parse panel YAML")
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
        .with_context(|| format!("Failed to parse panel file '{}'. Check for YAML syntax errors and ensure all primer fields are valid.", path))?;

    log_panel_version(&panel, path);

    // Prepend panel name to gene names
    for param in panel.primers.iter_mut() {
        param.gene_name = format!("{}_{}", panel.name, param.gene_name);
    }

    Ok(panel.primers)
}

/// Log the loaded panel's name and version at info level. A missing
/// `version` field produces a warning because versioning is expected for
/// reproducibility and a versionless panel breaks the panel development
/// cycle documented in PCR.md. Panel versions are deliberately *not*
/// checked against sharkmer's own version — PCR.md specifies they are
/// independent — but the values are surfaced so users can cross-reference
/// a run against the panel's changelog.
fn log_panel_version(panel: &PanelFile, source: &str) {
    match panel.version.as_deref() {
        Some(v) => log::info!(
            "Loaded panel '{}' v{} from {} ({} primer pair(s))",
            panel.name,
            v,
            source,
            panel.primers.len()
        ),
        None => log::warn!(
            "Panel '{}' from {} has no `version` field. Versioning is \
             recommended for reproducibility; see PCR.md.",
            panel.name,
            source
        ),
    }
}

/// Maximum size of a panel YAML file fetched over the network. Real panels
/// are a few kilobytes at most; this cap protects against a misconfigured
/// or malicious URL that would otherwise stream unbounded data.
const MAX_PANEL_YAML_BYTES: u64 = 10 * 1024 * 1024; // 10 MB

/// Load a panel from a URL.
fn load_panel_url(url: &str) -> Result<Vec<PCRParams>> {
    use std::io::Read;
    use std::time::Duration;

    log::info!("Downloading primer panel from {}", url);

    // Use an agent with explicit timeouts so a hung or slow endpoint cannot
    // block sharkmer indefinitely on startup.
    let agent = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(10))
        .timeout_read(Duration::from_secs(30))
        .build();

    let response = agent.get(url).call().with_context(|| {
        format!(
            "Failed to download panel from URL: {} (network error, timeout, or HTTP failure)",
            url
        )
    })?;

    // Bound the response size to avoid OOM on pathological content.
    let mut yaml_str = String::new();
    response
        .into_reader()
        .take(MAX_PANEL_YAML_BYTES + 1)
        .read_to_string(&mut yaml_str)
        .with_context(|| format!("Failed to read panel response body from URL: {}", url))?;
    if yaml_str.len() as u64 > MAX_PANEL_YAML_BYTES {
        bail!(
            "Panel YAML at {} exceeds maximum size of {} bytes",
            url,
            MAX_PANEL_YAML_BYTES
        );
    }

    let mut panel = parse_panel_yaml(&yaml_str)
        .with_context(|| format!("Downloaded panel from {} but failed to parse as YAML", url))?;

    log_panel_version(&panel, url);

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
        ("c_elegans", include_str!("../../panels/c_elegans.yaml")),
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

    panels.sort_by(|a, b| a.name.cmp(&b.name));

    panels
}

pub fn get_panel(panel_name: &str) -> Result<Vec<PCRParams>, String> {
    let panels = get_preconfigured_panels();

    match panels.iter().find(|panel| panel.name == panel_name) {
        Some(panel) => {
            log_panel_version(panel, &format!("built-in panel '{}'", panel_name));
            Ok(panel.primers.clone())
        }
        None => {
            let available: Vec<&str> = panels.iter().map(|p| p.name.as_str()).collect();
            Err(format!(
                "Unknown panel '{}'. Available panels: {}",
                panel_name,
                available.join(", ")
            ))
        }
    }
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
    fn test_builtin_panels_load_and_are_versioned() {
        // Every in-tree panel must parse and declare a version.
        let panels = get_builtin_panels();
        assert!(!panels.is_empty());
        for panel in &panels {
            assert!(
                panel.version.is_some(),
                "Panel '{}' missing version field",
                panel.name
            );
        }
    }

    #[test]
    fn test_deny_unknown_panel_field() {
        // A typo at the panel level must be rejected by deny_unknown_fields.
        let yaml = r#"
name: typo_panel
versoin: 1.0.0
description: "typo in version field"
primers:
  - gene_name: "X"
    forward_seq: "A"
    reverse_seq: "T"
"#;
        let result = parse_panel_yaml(yaml);
        assert!(result.is_err(), "expected rejection of unknown field");
    }

    #[test]
    fn test_deny_unknown_primer_field() {
        // A typo at the primer level must also be rejected.
        let yaml = r#"
name: typo_panel
version: 1.0.0
description: "typo in primer field"
primers:
  - gene_name: "X"
    forward_seq: "A"
    reverse_seq: "T"
    forward_sqe: "oops"
"#;
        let result = parse_panel_yaml(yaml);
        assert!(
            result.is_err(),
            "expected rejection of unknown primer field"
        );
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
        let version = panel.version.as_deref().unwrap_or("unversioned");
        println!(
            "  {:<16} v{:<8} {} ({} {})",
            panel.name, version, panel.description, n, noun
        );
    }
}
