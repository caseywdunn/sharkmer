//! Panel input, parser and retained primer owners on the original run ledger.
use super::*;
use crate::configuration::{self, Owned};
use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use std::io::Read;
use std::sync::Arc;

fn read_body(
    mut reader: impl Read,
    capacity: usize,
    limit: Option<u64>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<u8>> {
    let mut bytes =
        BudgetVec::with_capacity_in_class(capacity, budget, MemoryClass::RetainedResults)?;
    let mut scratch = [0; 16_384];
    loop {
        let allowance = limit.map_or(scratch.len(), |limit| {
            usize::try_from((limit + 1).saturating_sub(bytes.len() as u64))
                .unwrap_or(scratch.len())
                .min(scratch.len())
        });
        let n = match reader.read(&mut scratch[..allowance]) {
            Ok(0) => break,
            Ok(n) => n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error.into()),
        };
        if limit.is_some_and(|limit| bytes.len() as u64 + n as u64 > limit) {
            bail!(
                "Panel YAML exceeds maximum size of {} bytes",
                limit.unwrap_or(0)
            );
        }
        bytes.extend_from_slice(&scratch[..n])?;
    }
    // Preserve read_to_string's UTF-8 requirement rather than admitting other
    // encodings that the YAML byte loader can otherwise recognize.
    std::str::from_utf8(&bytes)
        .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidData, error))?;
    Ok(bytes)
}

fn decode(
    input: &[u8],
    source: &str,
    kind: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Owned<Vec<PCRParams>>> {
    let (data, mut owner) = configuration::yaml::<PanelFile>(input, budget)?.into_parts();
    let mut panel = data;
    {
        let label = if kind == "built-in panel" {
            Some(owner.text(format_args!("built-in panel '{source}'"))?)
        } else {
            None
        };
        let context = label.as_deref().unwrap_or(source);
        log_panel_version(&panel, context);
        require_clade_for_v2(&panel, context)?;
    }
    // Borrow keys while checking duplicates: no cloned name/region allocations.
    let mut seen = BudgetMap::in_class(budget, MemoryClass::RetainedResults);
    for (index, primer) in panel.primers.iter().enumerate() {
        if let Some(gene) = primer.gene.as_deref() {
            validate_gene_chars(gene, primer.region.is_some())?;
            let key = (gene, primer.region.as_deref(), primer.index);
            if let Some(previous) = seen.get(&key) {
                bail!(
                    "Panel '{}': duplicate primer entries for (gene={:?}, region={:?}, index={:?}) at positions {} and {}. Add an `index:` field to distinguish them.",
                    panel.name,
                    gene,
                    primer.region,
                    primer.index,
                    previous,
                    index
                );
            }
            seen.insert(key, index)?;
        }
        if let Some(region) = primer.region.as_deref() {
            validate_region_chars(region)?;
        }
    }
    drop(seen);
    let prefix = panel.gene_prefix.as_deref().unwrap_or(&panel.name);
    for primer in &mut panel.primers {
        // Derive and prefix in one admitted allocation, preserving the existing
        // spelling even for legacy entries without structured target fields.
        let gene = primer.gene.as_deref().unwrap_or(&primer.gene_name);
        primer.gene_name = match (primer.region.as_deref(), primer.index) {
            (Some(region), Some(index)) if primer.gene.is_some() => {
                owner.text(format_args!("{prefix}_{gene}-{region}_{index}"))?
            }
            (Some(region), None) if primer.gene.is_some() => {
                owner.text(format_args!("{prefix}_{gene}-{region}"))?
            }
            (None, Some(index)) if primer.gene.is_some() => {
                owner.text(format_args!("{prefix}_{gene}_{index}"))?
            }
            _ => owner.text(format_args!("{prefix}_{gene}"))?,
        };
        primer.source = owner.text(format_args!("{kind} '{source}'"))?;
        if primer.deprecated {
            let mut message = owner.text(format_args!(
                "Panel '{}': skipping deprecated primer '{}'.",
                panel.name, primer.gene_name
            ))?;
            if let Some(by) = &primer.deprecated_by {
                message = owner.text(format_args!("{message} Use '{by}' instead."))?;
            }
            if let Some(reason) = &primer.deprecated_reason {
                message = owner.text(format_args!("{message} Reason: {reason}"))?;
            }
            log::warn!("{message}");
        }
    }
    panel.primers.retain(|primer| !primer.deprecated);
    let primers = std::mem::take(&mut panel.primers);
    drop(panel);
    let bytes = primers.iter().try_fold(
        crate::kmer::budget::collections::bytes::<PCRParams>(primers.capacity())?,
        |total, primer| {
            total
                .checked_add(primer.configuration_bytes()?)
                .context("Panel retained size overflow")
        },
    )?;
    owner.retain_bytes(bytes)?;
    Ok(owner.finish(primers))
}

/// Load only the requested built-in panel, retaining its parsed allocations.
pub(crate) fn builtin(
    name: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Owned<Vec<PCRParams>>> {
    let sources = get_builtin_panel_sources();
    if let Some((_, text)) = sources.iter().find(|(candidate, _)| *candidate == name) {
        return decode(text.as_bytes(), name, "built-in panel", budget);
    }
    let mut owner = configuration::Owner::new(budget)?;
    let mut available = String::new();
    for (name, _) in sources {
        available = owner.text(format_args!(
            "{available}{}{name}",
            if available.is_empty() { "" } else { ", " }
        ))?;
    }
    bail!("Unknown panel '{}'. Available panels: {}", name, available)
}

/// Preserve local/URL panel rules with owned input and retained model data.
pub(crate) fn source(
    source: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Owned<Vec<PCRParams>>> {
    if is_url(source) {
        log::info!("Downloading primer panel from {}", source);
        let response = crate::http::panel(source).with_context(|| {
            format!(
                "Failed to download panel from URL: {} (network error, timeout, or HTTP failure)",
                source
            )
        })?;
        let input = read_body(response, 0, Some(MAX_PANEL_YAML_BYTES), budget)
            .with_context(|| format!("Failed to read panel response body from URL: {source}"))?;
        decode(&input, source, "panel URL", budget)
            .with_context(|| format!("Downloaded panel from {source} but failed to parse as YAML"))
    } else {
        let file = std::fs::File::open(source)
            .with_context(|| format!("Failed to read panel file: {source}"))?;
        let capacity = usize::try_from(file.metadata()?.len())
            .context("Panel file exceeds supported platform size")?;
        let input = read_body(file, capacity, None, budget)
            .with_context(|| format!("Failed to read panel file: {source}"))?;
        decode(&input, source, "panel file", budget).with_context(|| format!("Failed to parse panel file '{source}'. Check for YAML syntax errors and ensure all primer fields are valid."))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};

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

    #[test]
    fn bounded_builtin_and_local_panels_preserve_every_primer_field_and_release_discarded_metadata()
    {
        for (name, _) in get_builtin_panel_sources() {
            let ledger = budget(512 << 20);
            let mut old = get_panel(name).unwrap();
            for primer in &mut old {
                primer.source = format!("built-in panel '{name}'");
            }
            let new = builtin(name, Some(&ledger)).unwrap();
            assert_eq!(format!("{old:?}"), format!("{:?}", &*new));
            let actual = new
                .iter()
                .map(|primer| primer.configuration_bytes().unwrap())
                .sum::<u64>()
                + (new.capacity() * size_of::<PCRParams>()) as u64;
            assert_eq!(ledger.snapshot().reserved_memory_bytes, actual);
            drop(new);
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
            assert!(ledger.snapshot().valid);
        }
        let ledger = budget(64 << 20);
        let path = "tests/fixtures/test_panel.yaml";
        let mut old = load_panel_source(path).unwrap();
        for primer in &mut old {
            primer.source = format!("panel file '{path}'");
        }
        let new = source(path, Some(&ledger)).unwrap();
        assert_eq!(format!("{old:?}"), format!("{:?}", &*new));
        drop(new);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn local_panel_refusal_precedes_input_allocation_and_retries_on_original_ledger() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("panel.yaml");
        let text = format!(
            "{}\nnotes: \"{}\"\n",
            std::fs::read_to_string("tests/fixtures/test_panel.yaml").unwrap(),
            "A".repeat(65536)
        );
        std::fs::write(&path, text.as_bytes()).unwrap();
        let small = budget(4096);
        let error = source(path.to_str().unwrap(), Some(&small)).err().unwrap();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(std::fs::read(&path).unwrap(), text.as_bytes());
        assert_eq!(small.snapshot().reserved_memory_bytes, 0);
        let sufficient = budget(64 << 20);
        let params = source(path.to_str().unwrap(), Some(&sufficient)).unwrap();
        assert_eq!(params.len(), 1);
        // Top-level notes have dropped rather than remaining charged in PCR.
        assert!(sufficient.snapshot().reserved_memory_bytes < 8192);
        drop(params);
        assert_eq!(sufficient.snapshot().reserved_memory_bytes, 0);
    }
}

#[cfg(test)]
mod network_tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceLimits};
    use std::io::Write;
    use std::net::TcpListener;
    use std::time::Duration;

    fn budget() -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 64 << 20,
            process_overhead_bytes: 0,
            class_bytes: [64 << 20; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }
    fn serve(
        body: Vec<u8>,
        status: u16,
        gzip: bool,
        redirect: bool,
    ) -> (String, std::thread::JoinHandle<()>) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap();
        let handle = std::thread::spawn(move || {
            for index in 0..if redirect { 2 } else { 1 } {
                let (mut stream, _) = listener.accept().unwrap();
                stream
                    .set_read_timeout(Some(Duration::from_secs(5)))
                    .unwrap();
                let mut request = Vec::new();
                while !request.ends_with(b"\r\n\r\n") {
                    let mut byte = [0];
                    assert_eq!(stream.read(&mut byte).unwrap(), 1);
                    request.push(byte[0]);
                    assert!(request.len() <= 8192);
                }
                if redirect && index == 0 {
                    write!(stream, "HTTP/1.1 302 Found\r\nLocation: http://{address}/final\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
                } else {
                    write!(stream, "HTTP/1.1 {status} fixture\r\nContent-Length: {}\r\n{}Connection: close\r\n\r\n", body.len(), if gzip { "Content-Encoding: gzip\r\n" } else { "" }).unwrap();
                    stream.write_all(&body).unwrap();
                }
            }
        });
        (format!("http://{address}/panel"), handle)
    }

    #[test]
    fn owned_url_panel_preserves_gzip_redirect_fields_and_original_source_label() {
        let text = std::fs::read("tests/fixtures/test_panel.yaml").unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(&text).unwrap();
        let (url, thread) = serve(encoder.finish().unwrap(), 200, true, true);
        let ledger = budget();
        let params = source(&url, Some(&ledger)).unwrap();
        thread.join().unwrap();
        let mut old = load_panel_file("tests/fixtures/test_panel.yaml").unwrap();
        for primer in &mut old {
            primer.source = format!("panel URL '{url}'");
        }
        assert_eq!(format!("{old:?}"), format!("{:?}", &*params));
        let actual = params
            .iter()
            .map(|p| p.configuration_bytes().unwrap())
            .sum::<u64>()
            + (params.capacity() * size_of::<PCRParams>()) as u64;
        assert_eq!(ledger.snapshot().reserved_memory_bytes, actual);
        drop(params);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn url_panel_utf8_and_body_limit_failures_release_input_and_decoder() {
        for (body, message) in [
            (vec![0xff], "invalid utf-8"),
            (
                vec![b'A'; MAX_PANEL_YAML_BYTES as usize + 1],
                "exceeds maximum size",
            ),
        ] {
            let (url, thread) = serve(body, 200, false, false);
            let ledger = budget();
            let error = source(&url, Some(&ledger)).err().unwrap();
            thread.join().unwrap();
            assert!(format!("{error:#}").contains(message), "{error:#}");
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
            assert!(ledger.snapshot().valid);
        }
    }

    #[test]
    fn retained_url_status_error_propagates_without_charging_application_buffers() {
        let (url, thread) = serve(b"not found".to_vec(), 404, false, false);
        let ledger = budget();
        let error = source(&url, Some(&ledger)).err().unwrap();
        thread.join().unwrap();
        assert!(format!("{error:#}").contains("404"));
        assert!(matches!(
            error.downcast_ref::<ureq::Error>(),
            Some(ureq::Error::Status(404, _))
        ));
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        drop(error);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}

#[cfg(test)]
mod allocator_tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceLimits};

    #[test]
    fn actual_panel_models_cover_peak_and_every_retained_string_then_release() {
        const ENV: &str = "SHARKMER_PANEL_CONFIGURATION_ALLOCATOR_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact", "pcr::preconfigured::bounded::allocator_tests::actual_panel_models_cover_peak_and_every_retained_string_then_release", "--nocapture", "--test-threads=1"]).env(ENV,"1").status().unwrap();
            assert!(status.success());
            return;
        }
        let base = std::fs::read_to_string("tests/fixtures/test_panel.yaml").unwrap();
        drop(decode(base.as_bytes(), "fixture", "panel file", None).unwrap());
        for width in [1, 128, 4096] {
            let value = "Σßλ".repeat(width);
            let input = format!(
                "{}    region: V9\n    index: 2\n    compartment: \"{value}\"\n    gene_type: \"{value}\"\n    copy_number: \"{value}\"\n    deprecated_by: \"{value}\"\n    deprecated_reason: \"{value}\"\nnotes: \"{}\"\n",
                base.replace("Test citation", &value)
                    .replace("\"test\"", &format!("\"{value}\"")),
                "A".repeat(65536)
            );
            let ledger = ResourceBudget::new(ResourceLimits {
                memory_bytes: 64 << 20,
                process_overhead_bytes: 0,
                class_bytes: [64 << 20; MEMORY_CLASSES],
                temporary_disk_bytes: 0,
                temporary_open_files: 0,
            })
            .unwrap();
            let baseline = crate::PEAK_ALLOC.current_usage();
            crate::PEAK_ALLOC.reset_peak_usage();
            let params = decode(input.as_bytes(), "fixture", "panel file", Some(&ledger)).unwrap();
            let observed_peak = crate::PEAK_ALLOC.peak_usage().saturating_sub(baseline);
            let observed_retained = crate::PEAK_ALLOC.current_usage().saturating_sub(baseline);
            let snapshot = ledger.snapshot();
            assert!(observed_peak as u64 <= snapshot.peak_reserved_memory_bytes);
            assert!(
                observed_retained as u64 <= snapshot.reserved_memory_bytes,
                "width={width} retained={observed_retained} reserved={}",
                snapshot.reserved_memory_bytes
            );
            println!(
                "panel_configuration_allocator width={width} peak={observed_peak} reserved_peak={} retained={observed_retained} reserved_retained={}",
                snapshot.peak_reserved_memory_bytes, snapshot.reserved_memory_bytes
            );
            assert_eq!(params.len(), 1);
            drop(params);
            assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
            assert!(ledger.snapshot().valid);
        }
    }
}
