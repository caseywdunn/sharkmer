//! Delegate remote requests and content decoding to the HTTP client.
//!
//! Application body/record buffers have their own limits. Client internals use
//! process headroom and platform limits; they are outside the application ledger.

use anyhow::Result;
use std::io::Read;
use std::time::Duration;

/// Preserve the default client's redirects and automatic content decoding.
pub(crate) fn get(url: &str) -> Result<Box<dyn Read + Send + Sync>> {
    Ok(ureq::get(url).call()?.into_reader())
}

/// Keep the existing primer-panel connect and read timeouts.
pub(crate) fn panel(url: &str) -> Result<Box<dyn Read + Send + Sync>> {
    let agent = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(10))
        .timeout_read(Duration::from_secs(30))
        .build();
    Ok(agent.get(url).call()?.into_reader())
}

#[cfg(test)]
mod tests;
