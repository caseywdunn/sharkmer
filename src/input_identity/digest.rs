//! Fallible SHA-256 state for the unchanged selected-input byte stream.

use anyhow::{Context, Result, ensure};
use ring::digest::{Context as DigestContext, SHA256};

// SHA-256 encodes its input length in 64 bits. Check before calling the
// provider, whose finalization rejects an unrepresentable bit length by panic.
const MAX_BYTES: u64 = u64::MAX / 8;

#[derive(Clone)]
pub(crate) struct SelectedDigest {
    context: DigestContext,
    bytes: u64,
    failed: bool,
}

impl SelectedDigest {
    pub(crate) fn new() -> Self {
        Self {
            context: DigestContext::new(&SHA256),
            bytes: 0,
            failed: false,
        }
    }

    pub(crate) fn update(&mut self, bytes: impl AsRef<[u8]>) -> Result<()> {
        ensure!(!self.failed, "Selected-input SHA-256 previously failed");
        let bytes = bytes.as_ref();
        let next = u64::try_from(bytes.len())
            .ok()
            .and_then(|length| self.bytes.checked_add(length))
            .filter(|total| *total <= MAX_BYTES);
        let Some(next) = next else {
            self.failed = true;
            anyhow::bail!("Selected-input SHA-256 byte length exceeds the supported limit");
        };
        self.context.update(bytes);
        self.bytes = next;
        Ok(())
    }

    pub(crate) fn finalize(self) -> Result<[u8; 32]> {
        ensure!(!self.failed, "Selected-input SHA-256 previously failed");
        self.context
            .finish()
            .as_ref()
            .try_into()
            .context("Selected-input SHA-256 digest has an unexpected size")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

    #[test]
    fn empty_chunked_and_cloned_streams_match_an_independent_provider() {
        for length in [0, 1, 55, 56, 63, 64, 65, 65537] {
            let bytes = vec![b'x'; length];
            let mut state = SelectedDigest::new();
            for chunk in bytes.chunks(17) {
                state.update(chunk).unwrap();
            }
            let expected = <[u8; 32]>::from(Sha256::digest(&bytes));
            assert_eq!(state.clone().finalize().unwrap(), expected);
            state.update(b"continued").unwrap();
            let mut reference = Sha256::new();
            reference.update(&bytes);
            reference.update(b"continued");
            assert_eq!(
                state.finalize().unwrap(),
                <[u8; 32]>::from(reference.finalize())
            );
        }
    }

    #[test]
    fn length_rejection_precedes_provider_update_and_poison_is_permanent() {
        let mut state = SelectedDigest::new();
        // Exercise the guard at the limit without allocating an enormous input.
        state.bytes = MAX_BYTES - 1;
        state.update(b"x").unwrap();
        state.update(b"").unwrap();
        assert!(state.update(b"y").is_err());
        assert_eq!(state.bytes, MAX_BYTES);
        assert!(state.update(b"").is_err());
        assert!(state.clone().finalize().is_err());
        state.bytes = u64::MAX;
        assert!(state.finalize().is_err());
    }
}
