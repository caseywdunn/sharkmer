//! Process-aware allowance for the shared count/replay descriptor ledger.
//!
//! This is a startup preflight, not interception of arbitrary library opens or
//! protection against another thread changing the process limit afterward.

#[cfg(target_os = "linux")]
use anyhow::Context;
use anyhow::Result;
use serde::Serialize;

/// Startup observations and the resulting count/replay descriptor ceiling.
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct ProcessFileLimits {
    /// Maximum count/replay ownership requested by the counting frontend.
    pub(crate) requested_temporary_file_limit: u64,
    /// Actual finite OS soft limit, absent when unlimited or unavailable.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) process_soft_limit: Option<u64>,
    /// Open process descriptors at startup, excluding the census directory.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) process_open_at_start: Option<u64>,
    /// Additional allowance for input, output and library opens after startup.
    /// These descriptors are not individually controlled by this ledger yet.
    pub(crate) other_file_allowance: u64,
    /// Output lease known to open after this census; absent when already open.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_lease_allowance: Option<u64>,
    /// Effective shared count/replay limit after the startup preflight.
    pub(crate) temporary_file_limit: u64,
}

impl ProcessFileLimits {
    /// Leave conservative additional headroom for other stages and libraries.
    pub(crate) fn detect(requested: u64, output_lease_allowance: u64) -> Result<Self> {
        #[cfg(unix)]
        let soft_limit = rustix::process::getrlimit(rustix::process::Resource::Nofile).current;
        #[cfg(not(unix))]
        let soft_limit = None;

        #[cfg(target_os = "linux")]
        let open_at_start = Some(open_file_count()?);
        #[cfg(not(target_os = "linux"))]
        let open_at_start = None;

        Ok(Self::from_observations(
            requested,
            soft_limit,
            open_at_start,
            16,
            output_lease_allowance,
        ))
    }

    fn from_observations(
        requested: u64,
        soft_limit: Option<u64>,
        open_at_start: Option<u64>,
        other_file_allowance: u64,
        output_lease_allowance: u64,
    ) -> Self {
        let temporary_file_limit = soft_limit.map_or(requested, |limit| {
            limit
                .saturating_sub(open_at_start.unwrap_or(0))
                .saturating_sub(other_file_allowance)
                .saturating_sub(output_lease_allowance)
                .min(requested)
        });
        Self {
            requested_temporary_file_limit: requested,
            process_soft_limit: soft_limit,
            process_open_at_start: open_at_start,
            other_file_allowance,
            output_lease_allowance: (output_lease_allowance > 0).then_some(output_lease_allowance),
            temporary_file_limit,
        }
    }
}

#[cfg(target_os = "linux")]
fn open_file_count() -> Result<u64> {
    let directory =
        std::fs::read_dir("/proc/self/fd").context("Unable to inspect open process descriptors")?;
    let mut count = 0u64;
    for entry in directory {
        let entry = entry.context("Unable to inspect open process descriptor")?;
        if entry
            .file_name()
            .to_str()
            .is_some_and(|name| name.parse::<u64>().is_ok())
        {
            count = count
                .checked_add(1)
                .context("Process descriptor count overflow")?;
        }
    }
    // read_dir opened its own descriptor, which was included in the census and
    // has now closed. Do not charge it as a descriptor surviving startup.
    count
        .checked_sub(1)
        .context("Process descriptor census omitted its own directory")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_limits_keep_existing_and_future_owners_out_of_the_temporary_pool() {
        for (soft, open, headroom, expected) in [
            (1024, 4, 16, 128),
            (64, 7, 16, 41),
            (32, 7, 16, 9),
            (16, 7, 16, 0),
            (16, 17, 0, 0),
            (u64::MAX, u64::MAX, u64::MAX, 0),
        ] {
            let limits =
                ProcessFileLimits::from_observations(128, Some(soft), Some(open), headroom, 0);
            assert_eq!(limits.temporary_file_limit, expected);
            assert_eq!(limits.process_soft_limit, Some(soft));
            assert_eq!(limits.process_open_at_start, Some(open));
        }
        let unlimited = ProcessFileLimits::from_observations(128, None, Some(7), 16, 0);
        assert_eq!(unlimited.temporary_file_limit, 128);
        let unavailable = ProcessFileLimits::from_observations(128, Some(32), None, 16, 0);
        assert_eq!(unavailable.temporary_file_limit, 16);
        assert!(unavailable.process_open_at_start.is_none());
    }

    #[test]
    fn planned_output_lease_preserves_the_pool_without_inventing_a_census() {
        for soft in [0, 3, 16, 32, 64, 1024, u64::MAX] {
            for open in [0, 3, 30, 48, u64::MAX - 1] {
                let earlier =
                    ProcessFileLimits::from_observations(128, Some(soft), Some(open), 16, 1);
                let later =
                    ProcessFileLimits::from_observations(128, Some(soft), Some(open + 1), 16, 0);
                assert_eq!(earlier.temporary_file_limit, later.temporary_file_limit);
                assert_eq!(earlier.process_open_at_start, Some(open));
                assert_eq!(earlier.other_file_allowance, 16);
                assert_eq!(earlier.output_lease_allowance, Some(1));
            }
        }
        let unavailable = ProcessFileLimits::from_observations(128, Some(32), None, 16, 1);
        assert_eq!(unavailable.temporary_file_limit, 15);
        assert!(unavailable.process_open_at_start.is_none());
        assert!(
            ProcessFileLimits::from_observations(128, None, None, 16, 0)
                .output_lease_allowance
                .is_none()
        );
    }
}
