use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};
use std::collections::BTreeMap;

fn budget(bytes: u64) -> Arc<ResourceBudget> {
    ResourceBudget::new(ResourceLimits {
        memory_bytes: bytes,
        process_overhead_bytes: 0,
        class_bytes: [bytes; MEMORY_CLASSES],
        temporary_disk_bytes: 0,
        temporary_open_files: 0,
    })
    .unwrap()
}

#[derive(Debug, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Row {
    text: String,
    #[serde(default)]
    values: Vec<String>,
    #[serde(default)]
    optional: Option<String>,
}
#[derive(Debug, PartialEq, serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Model {
    name: String,
    rows: Vec<Row>,
    map: BTreeMap<String, Row>,
}

#[test]
fn tracked_derives_preserve_defaults_aliases_struct_sequences_and_rejections() {
    for text in [
        "name: 123\nrows: []\nmap: {}\n",
        "name: Unicode Σßλ\nrows: [&row {text: &scalar test, values: [*scalar, *scalar], optional: more}, *row]\nmap: {one: *row, two: *row}\n",
        "[name, [[first, [a, b], null], [second, [], yes]], {key: [third, [], null]}]",
        "name: x\nrows: []\nmap: {}\nextra: ignored\n",
        "name: x\nname: y\nrows: []\nmap: {}\n",
        "name: x\nrows: [{text: x, values: {bad: shape}}]\nmap: {}\n",
        "name: x\nrows: []\nmap: {}\n---\nname: y\nrows: []\nmap: {}\n",
    ] {
        let ledger = budget(64 << 20);
        let old = serde_yaml_ng::from_str::<Model>(text);
        let new = yaml::<Model>(text.as_bytes(), Some(&ledger));
        match (old, new) {
            (Ok(old), Ok(new)) => {
                assert_eq!(old, *new, "{text}");
                assert!(ledger.snapshot().reserved_memory_bytes > 0);
            }
            (Err(_), Err(_)) => {}
            (old, new) => panic!(
                "configuration compatibility differs old_ok={} new_ok={} input={text}",
                old.is_ok(),
                new.is_ok()
            ),
        }
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}

#[test]
fn alias_model_refusal_preserves_typed_resource_error_after_loader_admission() {
    let text = format!(
        "name: x\nrows: [{{text: &text {}, values: [{}]}}]\nmap: {{}}\n",
        "A".repeat(4096),
        vec!["*text"; 256].join(",")
    );
    let probe = budget(64 << 20);
    let workspace_bytes = {
        let _workspace =
            crate::output::parsing::workspace::Workspace::new(text.as_bytes(), Some(&probe))
                .unwrap();
        probe.snapshot().reserved_memory_bytes
    };
    let ledger = budget(workspace_bytes + 8192);
    let error = yaml::<Model>(text.as_bytes(), Some(&ledger)).err().unwrap();
    assert!(
        error.downcast_ref::<ResourceExhausted>().is_some(),
        "{error:#}"
    );
    assert!(ledger.snapshot().peak_reserved_memory_bytes >= workspace_bytes);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    assert!(ledger.snapshot().valid);
}

#[test]
fn owned_configuration_covers_isolated_allocator_peaks_and_releases() {
    const ENV: &str = "SHARKMER_CONFIGURATION_ALLOCATOR_CHILD";
    if std::env::var_os(ENV).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact", "configuration::tests::owned_configuration_covers_isolated_allocator_peaks_and_releases", "--nocapture", "--test-threads=1"]).env(ENV, "1").status().unwrap();
        assert!(status.success());
        return;
    }
    for copies in [0, 32, 512] {
        let text = format!(
            "name: test\nrows: [{{text: &text {}, values: [{}]}}]\nmap: {{one: {{text: *text}}, two: {{text: *text}}}}\n",
            "Σßλ".repeat(128),
            vec!["*text"; copies].join(",")
        );
        let ledger = budget(128 << 20);
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let parsed = yaml::<Model>(text.as_bytes(), Some(&ledger)).unwrap();
        let observed = crate::PEAK_ALLOC.peak_usage().saturating_sub(baseline);
        let reserved = ledger.snapshot().peak_reserved_memory_bytes;
        assert!(
            observed as u64 <= reserved,
            "copies={copies} observed={observed} reserved={reserved}"
        );
        assert_eq!(parsed.rows[0].values.len(), copies);
        println!("configuration_allocator copies={copies} observed={observed} reserved={reserved}");
        drop(parsed);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}

#[test]
fn retained_yaml_error_keeps_loader_admission_until_disposal_without_changing_text() {
    let text = format!(
        "name: demo\nrows: []\nmap: {{}}\n{}: value\n",
        "unknown".repeat(8192)
    );
    let previous = serde_yaml_ng::from_str::<Model>(&text)
        .unwrap_err()
        .to_string();
    let ledger = budget(64 << 20);
    let error = yaml::<Model>(text.as_bytes(), Some(&ledger)).err().unwrap();
    assert_eq!(error.to_string(), previous);
    assert!(
        ledger.snapshot().reserved_memory_bytes >= previous.len() as u64,
        "retained YAML error text has no surviving admission"
    );
    drop(error);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    assert!(ledger.snapshot().valid);
}
