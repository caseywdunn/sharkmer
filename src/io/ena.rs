//! Owned ENA response and first-row selection, on the original run ledger.

use super::*;
use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::text::BudgetString;
use std::io::Read;

// Preserve ureq 2.12.1 into_string's raw-body limit and lossy UTF-8 semantics.
const BODY_LIMIT: usize = 10 * 1024 * 1024;
const CLASS: MemoryClass = MemoryClass::RetainedResults;

/// A preview starts no count plan and must ignore unused counting settings.
/// Bound its metadata within existing process headroom and half the supplied
/// envelope; this allowance never augments a real ingestion run's ledger.
pub(crate) fn preview_budget(process_ceiling: u64) -> Result<Arc<ResourceBudget>> {
    let memory = crate::process_memory::NORMAL_PROCESS_HEADROOM_BYTES.min(process_ceiling / 2);
    let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
    classes[CLASS as usize] = memory;
    ResourceBudget::new(ResourceLimits {
        memory_bytes: memory,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: 0,
        temporary_open_files: 0,
    })
}

/// Metadata owners survive transfers into the original-pass read-source plan.
pub(crate) struct EnaResult {
    pub(crate) urls: BudgetVec<BudgetString>,
    pub(crate) scientific_name: Option<BudgetString>,
}

/// Resolve the same first data row and URL order as the released parser.
pub(crate) fn get_ena_fastq_urls(
    accession: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<EnaResult> {
    let url = BudgetString::format(
        format_args!(
            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp,scientific_name"
        ),
        budget,
        CLASS,
    )?;
    info!("Querying ENA for accession {}...", accession);
    fetch(&url, accession, budget)
}

fn fetch(url: &str, accession: &str, budget: Option<&Arc<ResourceBudget>>) -> Result<EnaResult> {
    let response = crate::http::get(url)
        .with_context(|| format!("Failed to query ENA API for accession {accession}"))?;
    let result = read_response(response, accession, budget)?;
    // Avoid materializing a second complete URL list solely for a log message.
    info!(
        "Found {} FASTQ file(s) for {}",
        result.urls.len(),
        accession
    );
    if let Some(name) = &result.scientific_name {
        info!("Scientific name: {}", name);
    }
    Ok(result)
}

fn read_body<R: Read>(
    mut reader: R,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<u8>> {
    let mut bytes = BudgetVec::in_class(budget, CLASS);
    let mut scratch = [0; 8192];
    let mut reader = (&mut reader).take(BODY_LIMIT as u64 + 1);
    loop {
        let n = match reader.read(&mut scratch) {
            Ok(0) => break,
            Ok(n) => n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error).context("Unable to read ENA response body"),
        };
        ensure!(
            bytes.len().checked_add(n).is_some_and(|n| n <= BODY_LIMIT),
            "response too big for into_string"
        );
        bytes.extend_from_slice(&scratch[..n])?;
    }
    Ok(bytes)
}

fn read_response<R: Read>(
    reader: R,
    accession: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<EnaResult> {
    // Read the entire bounded body before selecting a row: later I/O errors or
    // oversize data must not silently become a successful first-row response.
    let bytes = read_body(reader, budget).context("Failed to read ENA API response")?;
    if let Ok(text) = std::str::from_utf8(&bytes) {
        parse(text, accession, budget)
    } else {
        let text = decode_lossy(&bytes, budget)?;
        drop(bytes);
        parse(&text, accession, budget)
    }
}

fn lossy_length(mut bytes: &[u8]) -> Result<usize> {
    let mut length = 0_usize;
    while !bytes.is_empty() {
        match std::str::from_utf8(bytes) {
            Ok(text) => {
                return length
                    .checked_add(text.len())
                    .context("ENA decoded length overflow");
            }
            Err(error) => {
                length = length
                    .checked_add(error.valid_up_to())
                    .and_then(|n| n.checked_add(3))
                    .context("ENA decoded length overflow")?;
                let end = error.valid_up_to()
                    + error
                        .error_len()
                        .unwrap_or(bytes.len() - error.valid_up_to());
                bytes = &bytes[end..];
            }
        }
    }
    Ok(length)
}

fn decode_lossy(mut bytes: &[u8], budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetString> {
    let mut text = BudgetString::with_capacity(lossy_length(bytes)?, budget, CLASS)?;
    while !bytes.is_empty() {
        match std::str::from_utf8(bytes) {
            Ok(valid) => {
                text.push_str(valid)?;
                break;
            }
            Err(error) => {
                text.push_str(std::str::from_utf8(&bytes[..error.valid_up_to()])?)?;
                text.push_char(char::REPLACEMENT_CHARACTER)?;
                let end = error.valid_up_to()
                    + error
                        .error_len()
                        .unwrap_or(bytes.len() - error.valid_up_to());
                bytes = &bytes[end..];
            }
        }
    }
    Ok(text)
}

pub(super) fn parse(
    body: &str,
    accession: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<EnaResult> {
    let mut lines = body.lines();
    let header = lines.next().unwrap_or("");
    let data = lines.next().with_context(|| format!("ENA returned no results for accession '{accession}'. Check that the accession is valid."))?;
    let mut ftp_index = None;
    let mut scientific_index = None;
    for (index, field) in header.split('\t').enumerate() {
        if field == "fastq_ftp" && ftp_index.is_none() {
            ftp_index = Some(index);
        }
        if field == "scientific_name" && scientific_index.is_none() {
            scientific_index = Some(index);
        }
    }
    let ftp_index = ftp_index.context("ENA response missing fastq_ftp column")?;
    let ftp = data.split('\t').nth(ftp_index).filter(|s| !s.is_empty()).with_context(|| format!("ENA returned no FASTQ URLs for accession '{accession}'. The run may not have public FASTQ files."))?;
    // One exact vector admission before allocating any nested URL strings.
    let mut urls = BudgetVec::with_capacity_in_class(ftp.split(';').count(), budget, CLASS)?;
    for source in ftp.split(';') {
        let prefix = if source.starts_with("ftp://") || source.starts_with("http") {
            ""
        } else {
            "http://"
        };
        let length = prefix
            .len()
            .checked_add(source.len())
            .context("ENA URL length overflow")?;
        let mut url = BudgetString::with_capacity(length, budget, CLASS)?;
        url.push_str(prefix)?;
        url.push_str(source)?;
        urls.push(url)?;
    }
    let scientific_name = scientific_index
        .and_then(|index| data.split('\t').nth(index))
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|name| BudgetString::from_str(name, budget, CLASS))
        .transpose()?;
    Ok(EnaResult {
        urls,
        scientific_name,
    })
}

/// Construct the same sample spelling without an intermediate replaced name.
pub(crate) fn sample_name(
    args: &Args,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<(BudgetString, Option<EnaResult>)> {
    let (sample, cached) = if let Some(sample) = &args.sample {
        (BudgetString::from_str(sample, budget, CLASS)?, None)
    } else if let Some(accession) = &args.ena {
        let result = get_ena_fastq_urls(accession, budget)?;
        let sample = derived_sample(accession, result.scientific_name.as_deref(), budget)?;
        warn!(
            "No --sample provided, using '{}' derived from ENA metadata",
            sample
        );
        (sample, Some(result))
    } else {
        bail!(
            "--sample is required. Provide a sample name as output file prefix.\nWhen using --ena, the sample name can be derived automatically from ENA metadata."
        );
    };
    ensure!(
        sample
            .chars()
            .all(|c| c.is_alphanumeric() || c == '_' || c == '-' || c == '.'),
        "Sample name '{}' contains characters that are unsafe for filenames. Use only alphanumeric characters, hyphens, underscores, and periods.",
        sample
    );
    Ok((sample, cached))
}

fn derived_sample(
    accession: &str,
    name: Option<&str>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetString> {
    let Some(name) = name else {
        return BudgetString::from_str(accession, budget, CLASS);
    };
    let length = name
        .len()
        .checked_add(accession.len())
        .and_then(|n| n.checked_add(1))
        .context("Derived sample length overflow")?;
    let mut sample = BudgetString::with_capacity(length, budget, CLASS)?;
    for c in name.chars() {
        sample.push_char(if c == ' ' { '_' } else { c })?;
    }
    sample.push_char('_')?;
    sample.push_str(accession)?;
    Ok(sample)
}

#[cfg(test)]
mod tests;
