use anyhow::{Context, Result, bail};
use std::collections::BTreeMap;

use super::PCRParams;

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct PanelFile {
    name: String,
    /// Schema version. Absent = v1 (legacy). "2" = current schema.
    #[serde(default)]
    #[allow(dead_code)]
    schema_version: Option<String>,
    /// Panel version (semver). Required in schema v2.
    #[serde(default)]
    panel_version: Option<String>,
    description: String,
    /// NCBI-preferred taxon name for the clade this panel targets.
    #[serde(default)]
    #[allow(dead_code)]
    clade: Option<String>,
    /// NCBI Taxonomy ID for the clade.
    #[serde(default)]
    #[allow(dead_code)]
    taxon_id: Option<u32>,
    /// Optional override for the prefix prepended to gene names in output.
    /// Falls back to `name` when absent.
    #[serde(default)]
    gene_prefix: Option<String>,
    /// Lifecycle status: "experimental", "stable", or "deprecated".
    #[serde(default)]
    #[allow(dead_code)]
    status: Option<String>,
    /// URL where this panel file can be obtained (for externally distributed panels).
    #[serde(default)]
    #[allow(dead_code)]
    source_url: Option<String>,
    /// SPDX license identifier (e.g. "CC-BY-4.0").
    #[serde(default)]
    #[allow(dead_code)]
    license: Option<String>,
    /// Citation for the panel itself (distinct from per-primer citations).
    #[serde(default)]
    #[allow(dead_code)]
    citation: Option<String>,
    /// Free-text notes about the panel.
    #[serde(default)]
    #[allow(dead_code)]
    notes: Option<String>,
    /// JSON Schema URL for editor validation (e.g. VS Code YAML extension).
    #[serde(rename = "$schema")]
    #[serde(default)]
    #[allow(dead_code)]
    json_schema: Option<String>,
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
    #[serde(default)]
    notes: Option<String>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ChangelogEntry {
    panel_version: String,
    date: String,
    #[serde(default)]
    sharkmer_version: Option<String>,
    changes: String,
    #[serde(default)]
    notes: Option<String>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
#[allow(dead_code)]
struct ReferenceGene {
    gene: String,
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
    #[serde(default)]
    notes: Option<String>,
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

/// Derive the internal `gene_name` identifier from `gene`, optional `region`,
/// and optional `index`. Format: `{gene}[-{region}][_{index}]`.
fn derive_gene_name(gene: &str, region: Option<&str>, index: Option<u32>) -> String {
    match (region, index) {
        (Some(r), Some(i)) => format!("{}-{}_{}", gene, r, i),
        (Some(r), None) => format!("{}-{}", gene, r),
        (None, Some(i)) => format!("{}_{}", gene, i),
        (None, None) => gene.to_string(),
    }
}

/// Validate character constraints on `gene`.
///
/// `_` is always forbidden (index delimiter).
/// `-` is forbidden only when `region` is also set, because the derived name
/// `{gene}-{region}` becomes ambiguous if `gene` itself contains `-`.
/// When no `region` is set, dashes in `gene` are unambiguous (e.g. `psbA-trnH`).
fn validate_gene_chars(gene: &str, has_region: bool) -> Result<()> {
    if gene.contains('_') {
        bail!(
            "gene '{}' must not contain '_' (reserved as index delimiter in output names).",
            gene
        );
    }
    if has_region && gene.contains('-') {
        bail!(
            "gene '{}' must not contain '-' when a `region` is also set, because the derived \
             output name `{{gene}}-{{region}}` would be ambiguous. Use alphanumeric characters \
             only for `gene` when pairing it with a `region` (e.g. 'CytB' not 'Cyt-b').",
            gene
        );
    }
    Ok(())
}

/// Validate that a `region` value contains no `_` characters.
fn validate_region_chars(region: &str) -> Result<()> {
    if region.contains('_') {
        bail!(
            "region '{}' must not contain '_' (reserved as index delimiter in output names).",
            region
        );
    }
    Ok(())
}

/// Validate that all (gene, region, index) combinations are unique within a panel.
fn validate_primer_uniqueness(primers: &[super::PCRParams], panel_name: &str) -> Result<()> {
    use std::collections::HashMap;
    let mut seen: HashMap<(String, Option<String>, Option<u32>), usize> = HashMap::new();
    for (i, p) in primers.iter().enumerate() {
        if let Some(gene) = &p.gene {
            let key = (gene.clone(), p.region.clone(), p.index);
            if let Some(prev) = seen.insert(key.clone(), i) {
                bail!(
                    "Panel '{}': duplicate primer entries for (gene={:?}, region={:?}, index={:?}) \
                     at positions {} and {}. Add an `index:` field to distinguish them.",
                    panel_name,
                    key.0,
                    key.1,
                    key.2,
                    prev,
                    i
                );
            }
        }
    }
    Ok(())
}

/// Post-process primers loaded from a panel YAML:
/// - Validate character constraints on `gene` and `region`
/// - Validate uniqueness of (gene, region, index) within the panel
/// - Derive `gene_name` from structured fields (required in schema v2)
/// - Warn when index is absent on a primer that shares (gene, region) with another
fn resolve_primer_gene_names(primers: &mut [super::PCRParams], panel_name: &str) -> Result<()> {
    // Validate characters and count (gene, region) groups
    for p in primers.iter() {
        if let Some(gene) = &p.gene {
            validate_gene_chars(gene, p.region.is_some())?;
        }
        if let Some(region) = &p.region {
            validate_region_chars(region)?;
        }
    }

    // Validate uniqueness
    validate_primer_uniqueness(primers, panel_name)?;

    // Derive gene_name for each primer that uses the structured fields
    for p in primers.iter_mut() {
        if let Some(gene) = &p.gene.clone() {
            p.gene_name = derive_gene_name(gene, p.region.as_deref(), p.index);
        }
        // If gene is absent the entry was loaded via CLI path and gene_name
        // is already set directly — leave it untouched.
    }

    Ok(())
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
    require_clade_for_v2(&panel, path)?;

    resolve_primer_gene_names(&mut panel.primers, &panel.name)
        .with_context(|| format!("Invalid primer specification in panel file '{}'", path))?;

    let prefix = panel.gene_prefix.as_deref().unwrap_or(&panel.name);
    for param in panel.primers.iter_mut() {
        param.gene_name = format!("{}_{}", prefix, param.gene_name);
    }

    Ok(filter_deprecated_primers(panel.primers, &panel.name))
}

/// Log the loaded panel's name and version at info level. Warns on missing
/// `panel_version` and on `status: "deprecated"`.
fn log_panel_version(panel: &PanelFile, source: &str) {
    if panel.status.as_deref() == Some("deprecated") {
        log::warn!(
            "Panel '{}' from {} has status 'deprecated'. Consider switching to a newer panel.",
            panel.name,
            source
        );
    }
    match panel.panel_version.as_deref() {
        Some(v) => log::info!(
            "Loaded panel '{}' v{} from {} ({} primer pair(s))",
            panel.name,
            v,
            source,
            panel.primers.len()
        ),
        None => log::warn!(
            "Panel '{}' from {} has no `panel_version` field. Versioning is \
             recommended for reproducibility; see PCR.md.",
            panel.name,
            source
        ),
    }
}

/// For schema v2 panels, require that `clade` is set.
fn require_clade_for_v2(panel: &PanelFile, source: &str) -> Result<()> {
    if panel.schema_version.as_deref() == Some("2") && panel.clade.is_none() {
        bail!(
            "Panel '{}' from {} declares schema_version: \"2\" but is missing the required \
             `clade` field. Set `clade` to the NCBI-preferred taxon name for the target clade \
             (e.g. `clade: \"Cnidaria\"`).",
            panel.name,
            source
        );
    }
    Ok(())
}

/// Remove deprecated primers from the list, emitting a warning for each one
/// skipped. Called after gene names are fully resolved and prefixed so that
/// the warning message uses the final output name.
fn filter_deprecated_primers(primers: Vec<PCRParams>, panel_name: &str) -> Vec<PCRParams> {
    let mut active = Vec::with_capacity(primers.len());
    for p in primers {
        if p.deprecated {
            let mut msg = format!(
                "Panel '{}': skipping deprecated primer '{}'.",
                panel_name, p.gene_name
            );
            if let Some(ref by) = p.deprecated_by {
                msg.push_str(&format!(" Use '{}' instead.", by));
            }
            if let Some(ref reason) = p.deprecated_reason {
                msg.push_str(&format!(" Reason: {}", reason));
            }
            log::warn!("{}", msg);
        } else {
            active.push(p);
        }
    }
    active
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
    require_clade_for_v2(&panel, url)?;

    resolve_primer_gene_names(&mut panel.primers, &panel.name).with_context(|| {
        format!(
            "Invalid primer specification in panel downloaded from '{}'",
            url
        )
    })?;

    let prefix = panel.gene_prefix.as_deref().unwrap_or(&panel.name);
    for param in panel.primers.iter_mut() {
        param.gene_name = format!("{}_{}", prefix, param.gene_name);
    }

    Ok(filter_deprecated_primers(panel.primers, &panel.name))
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
        ("hydrozoa", include_str!("../../panels/hydrozoa.yaml")),
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

    for panel in panels.iter_mut() {
        require_clade_for_v2(panel, &format!("built-in panel '{}'", panel.name))
            .unwrap_or_else(|e| panic!("{}", e));
        resolve_primer_gene_names(&mut panel.primers, &panel.name).unwrap_or_else(|e| {
            panic!("Built-in panel '{}' has invalid primers: {}", panel.name, e)
        });
        let prefix = panel.gene_prefix.as_deref().unwrap_or(&panel.name);
        let prefix = prefix.to_string();
        for param in panel.primers.iter_mut() {
            param.gene_name = format!("{}_{}", prefix, param.gene_name);
        }
    }

    panels.sort_by(|a, b| a.name.cmp(&b.name));

    panels
}

pub fn get_panel(panel_name: &str) -> Result<Vec<PCRParams>, String> {
    let panels = get_preconfigured_panels();

    match panels.iter().find(|panel| panel.name == panel_name) {
        Some(panel) => {
            let source = format!("built-in panel '{}'", panel_name);
            log_panel_version(panel, &source);
            Ok(filter_deprecated_primers(
                panel.primers.clone(),
                &panel.name,
            ))
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

pub fn print_pcr_panels() {
    let panels = get_preconfigured_panels();
    println!("Available PCR panels (use --export-panel <name> for details):\n");
    for panel in &panels {
        let n = panel.primers.len();
        let noun = if n == 1 { "primer" } else { "primers" };
        let version = panel.panel_version.as_deref().unwrap_or("unversioned");
        println!(
            "  {:<16} v{:<8} {} ({} {})",
            panel.name, version, panel.description, n, noun
        );
    }
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
        // Every in-tree panel must parse, declare panel_version, and — if
        // schema_version is "2" — declare clade.
        let panels = get_builtin_panels();
        assert!(!panels.is_empty());
        for panel in &panels {
            assert!(
                panel.panel_version.is_some(),
                "Panel '{}' missing panel_version field",
                panel.name
            );
            if panel.schema_version.as_deref() == Some("2") {
                assert!(
                    panel.clade.is_some(),
                    "Panel '{}' has schema_version \"2\" but is missing required `clade` field",
                    panel.name
                );
            }
        }
    }

    #[test]
    fn test_v2_panel_missing_clade_is_rejected() {
        let yaml = r#"
name: no_clade_panel
schema_version: "2"
panel_version: "1.0.0"
description: "v2 panel without clade"
primers:
  - gene: "X"
    forward_seq: "AAAA"
    reverse_seq: "TTTT"
    min_length: 100
    max_length: 500
    min_count: 2
    mismatches: 2
    trim: 15
"#;
        let mut panel = parse_panel_yaml(yaml).unwrap();
        resolve_primer_gene_names(&mut panel.primers, &panel.name).unwrap();
        let result = require_clade_for_v2(&panel, "test");
        assert!(result.is_err(), "expected error for v2 panel missing clade");
        assert!(
            result.unwrap_err().to_string().contains("clade"),
            "error message should mention clade"
        );
    }

    #[test]
    fn test_deny_unknown_panel_field() {
        // A typo at the panel level must be rejected by deny_unknown_fields.
        let yaml = r#"
name: typo_panel
versoin: 1.0.0
description: "typo in version field"
primers:
  - gene: "X"
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
panel_version: 1.0.0
description: "typo in primer field"
primers:
  - gene: "X"
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
    fn test_derive_gene_name() {
        assert_eq!(derive_gene_name("CO1", None, None), "CO1");
        assert_eq!(derive_gene_name("18S", Some("V9"), None), "18S-V9");
        assert_eq!(derive_gene_name("CO1", None, Some(2)), "CO1_2");
        assert_eq!(
            derive_gene_name("18S", Some("V5-V7"), Some(1)),
            "18S-V5-V7_1"
        );
    }

    #[test]
    fn test_validate_gene_chars_rejects_dash_with_region() {
        // Dash in gene is only an error when region is also set
        assert!(validate_gene_chars("Cyt-b", true).is_err());
        assert!(validate_gene_chars("CO-1", true).is_err());
    }

    #[test]
    fn test_validate_gene_chars_allows_dash_without_region() {
        // Established names like psbA-trnH are fine when no region is set
        assert!(validate_gene_chars("psbA-trnH", false).is_ok());
        assert!(validate_gene_chars("trnL-F", false).is_ok());
        assert!(validate_gene_chars("Cyt-b", false).is_ok());
    }

    #[test]
    fn test_validate_gene_chars_rejects_underscore() {
        assert!(validate_gene_chars("18S_rRNA", false).is_err());
        assert!(validate_gene_chars("18S_rRNA", true).is_err());
    }

    #[test]
    fn test_validate_gene_chars_accepts_valid() {
        assert!(validate_gene_chars("CO1", false).is_ok());
        assert!(validate_gene_chars("18S", false).is_ok());
        assert!(validate_gene_chars("5.8S", false).is_ok());
    }

    #[test]
    fn test_validate_region_chars_rejects_underscore() {
        assert!(validate_region_chars("V5_V7").is_err());
    }

    #[test]
    fn test_validate_region_chars_accepts_dash() {
        assert!(validate_region_chars("V5-V7").is_ok());
        assert!(validate_region_chars("V9").is_ok());
    }

    #[test]
    fn test_primer_uniqueness_conflict() {
        let yaml = r#"
name: dup_panel
panel_version: "1.0.0"
description: "duplicate primer test"
primers:
  - gene: "CO1"
    forward_seq: "AAAA"
    reverse_seq: "TTTT"
  - gene: "CO1"
    forward_seq: "CCCC"
    reverse_seq: "GGGG"
"#;
        let mut panel = parse_panel_yaml(yaml).unwrap();
        let result = resolve_primer_gene_names(&mut panel.primers, &panel.name);
        assert!(result.is_err(), "expected uniqueness conflict error");
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("duplicate"), "error should mention duplicate");
    }

    #[test]
    fn test_load_panel_source_bad_url() {
        let result = load_panel_source("https://localhost:1/nonexistent_panel.yaml");
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("Failed to download panel from URL"));
    }
}
