//! Optional Linux address-space safeguard, installed before ordinary CLI parsing.
//! This bounds virtual mappings, not kernel caches or other processes. Owned
//! stage reservations remain responsible for controlled admission and spilling.
#[cfg(any(target_os = "linux", test))]
use anyhow::ensure;
use anyhow::{Context, Result};
use serde::Serialize;
use std::ffi::OsStr;

const ENV: &str = "SHARKMER_PROCESS_ADDRESS_SPACE_LIMIT_MIB";

/// Normal Linux envelope measured for the 16 GB hardware baseline.
pub(crate) const NORMAL_ADDRESS_SPACE_BYTES: u64 = 12 << 30;
/// Keep planning space beyond owned buffers for stacks, libraries and allocator.
pub(crate) const NORMAL_PROCESS_HEADROOM_BYTES: u64 = 2 << 30;

/// Actual installed process cap; omitted when no safeguard was installed.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct ProcessAddressSpace {
    requested_limit_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    inherited_soft_limit_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    inherited_hard_limit_bytes: Option<u64>,
    effective_limit_bytes: u64,
    observed_address_space_bytes: u64,
}

fn parse_request(value: &OsStr) -> Result<Option<u64>> {
    let mib: u64 = value
        .to_str()
        .with_context(|| format!("{ENV} must contain an unsigned integer"))?
        .parse()
        .with_context(|| format!("{ENV} must contain an unsigned integer"))?;
    if mib == 0 {
        return Ok(None);
    }
    let bytes = mib
        .checked_mul(1 << 20)
        .context("Process address-space limit overflow")?;
    Ok(Some(bytes))
}

#[cfg(any(target_os = "linux", test))]
impl ProcessAddressSpace {
    fn plan(
        requested: u64,
        soft: Option<u64>,
        hard: Option<u64>,
        observed: u64,
        page: u64,
    ) -> Result<Self> {
        ensure!(page > 0, "Invalid operating-system page size");
        let ceiling = requested
            .min(soft.unwrap_or(u64::MAX))
            .min(hard.unwrap_or(u64::MAX));
        let effective = ceiling / page * page;
        ensure!(effective > 0, "Process address-space cap is below one page");
        ensure!(
            observed <= effective,
            "Process image already uses {observed} address-space bytes, exceeding the requested or inherited cap of {effective} bytes"
        );
        Ok(Self {
            requested_limit_bytes: requested,
            inherited_soft_limit_bytes: soft,
            inherited_hard_limit_bytes: hard,
            effective_limit_bytes: effective,
            observed_address_space_bytes: observed,
        })
    }
}

impl ProcessAddressSpace {
    /// Actual kernel ceiling used by production resource preflight.
    pub(crate) fn effective_limit_bytes(&self) -> u64 {
        self.effective_limit_bytes
    }
}

/// Bound metadata previews by an installed or inherited Linux process limit.
/// This observes limits without installing a cap or evaluating count controls.
pub(crate) fn metadata_ceiling(existing: Option<&ProcessAddressSpace>) -> u64 {
    let ceiling = existing.map_or(
        NORMAL_ADDRESS_SPACE_BYTES,
        ProcessAddressSpace::effective_limit_bytes,
    );
    inherited_ceiling(ceiling)
}

/// Use a lower inherited Linux limit when planning an explicit memory budget.
pub(crate) fn inherited_ceiling(ceiling: u64) -> u64 {
    #[cfg(target_os = "linux")]
    {
        let inherited = rustix::process::getrlimit(rustix::process::Resource::As);
        ceiling
            .min(inherited.current.unwrap_or(u64::MAX))
            .min(inherited.maximum.unwrap_or(u64::MAX))
    }
    #[cfg(not(target_os = "linux"))]
    ceiling
}

/// Read the optional development setting before CLI, panel and ENA decoding.
/// Zero or absence preserves inherited limits; no production default is chosen.
pub(crate) fn from_environment() -> Result<Option<ProcessAddressSpace>> {
    let Some(value) = std::env::var_os(ENV) else {
        return Ok(None);
    };
    let Some(bytes) = parse_request(&value)? else {
        return Ok(None);
    };
    apply(bytes).map(Some)
}

/// Apply the normal Linux safeguard before panel/ENA decoding. An explicitly
/// disabled safeguard preserves inherited limits. An explicit CLI budget can
/// lower a cap already installed from the environment.
pub(crate) fn for_normal_run(
    existing: Option<ProcessAddressSpace>,
    normal_run: bool,
    requested: Option<u64>,
) -> Result<Option<ProcessAddressSpace>> {
    if !normal_run {
        return Ok(existing);
    }
    if let Some(cap) = existing {
        if requested.is_none_or(|bytes| bytes >= cap.effective_limit_bytes()) {
            return Ok(existing);
        }
    } else if std::env::var_os(ENV).is_some() {
        return Ok(None);
    }
    #[cfg(target_os = "linux")]
    {
        apply(requested.unwrap_or(NORMAL_ADDRESS_SPACE_BYTES)).map(Some)
    }
    #[cfg(not(target_os = "linux"))]
    {
        Ok(None)
    }
}

#[cfg(not(target_os = "linux"))]
fn apply(_requested: u64) -> Result<ProcessAddressSpace> {
    anyhow::bail!("The development process address-space cap is supported on Linux only")
}

#[cfg(target_os = "linux")]
fn apply(requested: u64) -> Result<ProcessAddressSpace> {
    use rustix::process::{Resource, Rlimit, getrlimit, setrlimit};
    let inherited = getrlimit(Resource::As);
    let page = u64::try_from(rustix::param::page_size())?;
    let observed = observed_bytes(page)?;
    let plan = ProcessAddressSpace::plan(
        requested,
        inherited.current,
        inherited.maximum,
        observed,
        page,
    )?;
    // Lower both limits. Never raise a stricter inherited soft or hard ceiling.
    setrlimit(
        Resource::As,
        Rlimit {
            current: Some(plan.effective_limit_bytes),
            maximum: Some(plan.effective_limit_bytes),
        },
    )
    .context("Unable to install the process address-space cap")?;
    let actual = getrlimit(Resource::As);
    ensure!(
        actual.current == Some(plan.effective_limit_bytes) && actual.maximum == actual.current,
        "Operating system did not retain the requested process address-space cap"
    );
    Ok(plan)
}

#[cfg(target_os = "linux")]
fn observed_bytes(page: u64) -> Result<u64> {
    use std::io::Read;
    let mut file = std::fs::File::open("/proc/self/statm")
        .context("Unable to inspect current process address space")?;
    let mut bytes = [0u8; 64];
    let mut used = 0;
    loop {
        ensure!(
            used < bytes.len(),
            "Invalid process address-space observation"
        );
        let read = match file.read(&mut bytes[used..]) {
            Ok(0) => anyhow::bail!("Incomplete process address-space observation"),
            Ok(read) => read,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error).context("Unable to read process address space"),
        };
        used += read;
        if let Some(end) = bytes[..used].iter().position(u8::is_ascii_whitespace) {
            let pages: u64 = std::str::from_utf8(&bytes[..end])?.parse()?;
            return pages
                .checked_mul(page)
                .context("Process address-space observation overflow");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn process_cap_values_preserve_omission_and_check_units_and_overflow() {
        assert_eq!(parse_request(OsStr::new("0")).unwrap(), None);
        assert_eq!(parse_request(OsStr::new("256")).unwrap(), Some(256 << 20));
        for value in ["", "-1", "1.5", "18446744073709551615"] {
            assert!(parse_request(OsStr::new(value)).is_err(), "{value}");
        }
    }

    #[test]
    fn process_cap_plans_never_raise_inherited_limits_or_admit_an_oversized_image() {
        for (soft, hard, expected) in [
            (None, None, 16384),
            (Some(8193), None, 8192),
            (None, Some(12289), 12288),
            (Some(4096), Some(8192), 4096),
        ] {
            let plan = ProcessAddressSpace::plan(16385, soft, hard, 4096, 4096).unwrap();
            assert_eq!(plan.effective_limit_bytes, expected);
            assert_eq!(plan.inherited_soft_limit_bytes, soft);
            assert_eq!(plan.inherited_hard_limit_bytes, hard);
        }
        assert!(ProcessAddressSpace::plan(4096, None, None, 4097, 4096).is_err());
        assert!(ProcessAddressSpace::plan(4095, None, None, 1, 4096).is_err());
        assert!(ProcessAddressSpace::plan(4096, None, None, 1, 0).is_err());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn process_cap_is_enforced_by_the_kernel_in_an_isolated_child() {
        const CHILD: &str = "SHARKMER_PROCESS_CAP_AUDIT_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let status=std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact","process_memory::tests::process_cap_is_enforced_by_the_kernel_in_an_isolated_child","--nocapture","--test-threads=1"])
                .env(CHILD,"1").env(ENV,"256").status().unwrap();
            assert!(status.success());
            return;
        }
        use rustix::process::{Resource, getrlimit};
        let before = getrlimit(Resource::As);
        assert!(apply(1).is_err());
        assert_eq!(getrlimit(Resource::As), before);
        let report = from_environment().unwrap().unwrap();
        assert!(report.effective_limit_bytes <= 256 << 20);
        let actual = getrlimit(Resource::As);
        assert_eq!(actual.current, Some(report.effective_limit_bytes));
        assert_eq!(actual.maximum, actual.current);
        let mut oversized = Vec::<u8>::new();
        assert!(oversized.try_reserve_exact(512 << 20).is_err());
        assert_eq!(oversized.capacity(), 0);
        let again = apply(512 << 20).unwrap();
        assert_eq!(again.effective_limit_bytes, report.effective_limit_bytes);
        println!(
            "process_address_space requested_bytes={} effective_bytes={} observed_bytes={} refused_allocation_bytes={}",
            report.requested_limit_bytes,
            report.effective_limit_bytes,
            report.observed_address_space_bytes,
            512u64 << 20
        );
    }
}
