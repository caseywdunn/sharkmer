//! Typed identities and positioned windows shared by v4 count-store boundaries.
//!
//! Store implementations may pack raw bits internally, but identities crossing
//! a store boundary retain their k so mismatched queries can be rejected.

use anyhow::{Result, bail, ensure};
use std::iter::FusedIterator;

/// Validated length for the current two-bit, u64 k-mer representation.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct KmerSize(u8);

impl KmerSize {
    /// Validate a supported k without truncating the caller's value.
    pub fn new(k: usize) -> Result<Self> {
        ensure!(
            (1..=31).contains(&k),
            "k must be between 1 and 31, got {}",
            k
        );
        Ok(Self(k as u8))
    }

    /// Return the number of encoded bases.
    pub fn get(self) -> usize {
        usize::from(self.0)
    }

    fn mask(self) -> u64 {
        (1u64 << (2 * self.get())) - 1
    }
}

/// An oriented graph/primer identity; reverse complements remain distinct.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct OrientedKey {
    bits: u64,
    size: KmerSize,
}

// These constructors form the persisted/query boundary used by later backends.
#[allow(dead_code)]
impl OrientedKey {
    /// Reject bits outside k instead of silently aliasing another sequence.
    pub fn new(bits: u64, size: KmerSize) -> Result<Self> {
        ensure!(
            bits & !size.mask() == 0,
            "Encoded k-mer exceeds k={}",
            size.get()
        );
        Ok(Self { bits, size })
    }

    /// Return the directional encoding.
    pub fn bits(self) -> u64 {
        self.bits
    }

    /// Return the identity's validated length.
    pub fn size(self) -> KmerSize {
        self.size
    }

    /// Reverse-complement this identity while preserving its k.
    pub fn reverse_complement(self) -> Self {
        let mut remaining = self.bits;
        let mut bits = 0;
        for _ in 0..self.size.get() {
            bits = (bits << 2) | (3 - (remaining & 3));
            remaining >>= 2;
        }
        Self {
            bits,
            size: self.size,
        }
    }

    /// Normalize orientation for exact count lookup.
    pub fn canonical(self) -> CanonicalKey {
        CanonicalKey {
            bits: self.bits.min(self.reverse_complement().bits),
            size: self.size,
        }
    }
}

/// Exact canonical count identity, including k; distinct k values cannot alias.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct CanonicalKey {
    bits: u64,
    size: KmerSize,
}

impl CanonicalKey {
    /// Decode and normalize an externally supplied oriented encoding.
    #[allow(dead_code)]
    pub fn new(bits: u64, size: KmerSize) -> Result<Self> {
        Ok(OrientedKey::new(bits, size)?.canonical())
    }

    /// Return the canonical encoding for internal packed storage.
    pub fn bits(self) -> u64 {
        self.bits
    }

    /// Return k for validating a query against a store's metadata.
    #[allow(dead_code)]
    pub fn size(self) -> KmerSize {
        self.size
    }
}

/// A window's original zero-based byte position and both identity roles.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KmerWindow {
    pub start: usize,
    pub oriented: OrientedKey,
    pub canonical: CanonicalKey,
}

/// Allocation-free rolling extraction under the released uppercase ACGTN policy.
///
/// N breaks windows without erasing original positions. Other bytes produce
/// one error and terminate the iterator. Consumers streaming into mutable stores
/// must abort publication on error; collecting consumers discard partial results.
pub struct KmerWindows<'a> {
    sequence: &'a [u8],
    size: KmerSize,
    mask: u64,
    next: usize,
    forward: u64,
    reverse: u64,
    valid: usize,
    valid_bases_seen: usize,
}

impl<'a> KmerWindows<'a> {
    /// Validate k before yielding any windows.
    pub fn new(sequence: &'a [u8], k: usize) -> Result<Self> {
        let size = KmerSize::new(k)?;
        Ok(Self {
            sequence,
            size,
            mask: size.mask(),
            next: 0,
            forward: 0,
            reverse: 0,
            valid: 0,
            valid_bases_seen: 0,
        })
    }

    /// Number of ACGT bytes processed so far, including runs shorter than k.
    pub fn valid_bases_seen(&self) -> usize {
        self.valid_bases_seen
    }
}

impl Iterator for KmerWindows<'_> {
    type Item = Result<KmerWindow>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(&byte) = self.sequence.get(self.next) {
            let position = self.next;
            self.next += 1;
            let base = match byte {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'N' => {
                    self.forward = 0;
                    self.reverse = 0;
                    self.valid = 0;
                    continue;
                }
                _ => {
                    self.next = self.sequence.len();
                    return Some(invalid_base(byte));
                }
            };
            self.forward = ((self.forward << 2) | base) & self.mask;
            self.reverse = (self.reverse >> 2) | ((3 - base) << (2 * (self.size.get() - 1)));
            self.valid += 1;
            self.valid_bases_seen += 1;
            if self.valid >= self.size.get() {
                return Some(Ok(KmerWindow {
                    start: position + 1 - self.size.get(),
                    oriented: OrientedKey {
                        bits: self.forward,
                        size: self.size,
                    },
                    canonical: CanonicalKey {
                        bits: self.forward.min(self.reverse),
                        size: self.size,
                    },
                }));
            }
        }
        None
    }
}

impl FusedIterator for KmerWindows<'_> {}

fn invalid_base(byte: u8) -> Result<KmerWindow> {
    bail!(
        "Invalid character '{}' in sequence. Only ACGTN allowed.",
        byte as char
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // Deliberately scalar: re-encode each literal window independently rather
    // than sharing rolling state or the production reverse-complement routine.
    fn oracle(sequence: &[u8], k: usize) -> Vec<(usize, u64, u64)> {
        sequence
            .windows(k)
            .enumerate()
            .filter_map(|(start, window)| {
                let bases = window
                    .iter()
                    .map(|base| match base {
                        b'A' => Some(0u64),
                        b'C' => Some(1),
                        b'G' => Some(2),
                        b'T' => Some(3),
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>()?;
                let forward = bases.iter().fold(0, |bits, base| (bits << 2) | base);
                let reverse = bases
                    .iter()
                    .rev()
                    .fold(0, |bits, base| (bits << 2) | (3 - base));
                Some((start, forward, forward.min(reverse)))
            })
            .collect()
    }

    #[test]
    fn rolling_windows_match_independent_oracle_with_gaps_and_boundary_lengths() {
        let mut state = 7u64;
        for length in [0, 1, 2, 3, 30, 31, 32, 62, 150] {
            for _ in 0..40 {
                let sequence = (0..length)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                        b"ACGTN"[(state >> 32) as usize % 5]
                    })
                    .collect::<Vec<_>>();
                for k in [1, 2, 3, 19, 30, 31] {
                    let observed = KmerWindows::new(&sequence, k)
                        .unwrap()
                        .map(|window| {
                            let window = window.unwrap();
                            assert_eq!(window.oriented.size().get(), k);
                            assert_eq!(window.canonical.size().get(), k);
                            assert_eq!(window.oriented.canonical(), window.canonical);
                            (
                                window.start,
                                window.oriented.bits(),
                                window.canonical.bits(),
                            )
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(
                        observed,
                        oracle(&sequence, k),
                        "k={k}, sequence={sequence:?}"
                    );
                }
            }
        }
        for sequence in [
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_slice(),
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            b"ACGTACGTACGTACGTACGTACGTACGTACGT",
        ] {
            for k in [1, 2, 30, 31] {
                let observed = KmerWindows::new(sequence, k)
                    .unwrap()
                    .map(|w| {
                        let w = w.unwrap();
                        (w.start, w.oriented.bits(), w.canonical.bits())
                    })
                    .collect::<Vec<_>>();
                assert_eq!(observed, oracle(sequence, k));
            }
        }
    }

    #[test]
    fn identity_rejects_invalid_lengths_bits_and_cross_length_aliasing() {
        for k in [0, 32, 256, usize::MAX] {
            assert!(KmerSize::new(k).is_err());
        }
        let size = KmerSize::new(3).unwrap();
        assert!(OrientedKey::new(64, size).is_err());
        assert!(CanonicalKey::new(u64::MAX, size).is_err());
        let key = OrientedKey::new(7, size).unwrap();
        assert_ne!(key, key.reverse_complement());
        assert_eq!(key.canonical(), key.reverse_complement().canonical());
        assert_ne!(
            CanonicalKey::new(0, size).unwrap(),
            CanonicalKey::new(0, KmerSize::new(2).unwrap()).unwrap()
        );
        assert_eq!(key.reverse_complement().reverse_complement(), key);
    }

    #[test]
    fn invalid_bytes_fail_once_and_never_bridge_gaps() {
        for invalid in [b'a', b'R', 0, 255] {
            let sequence = [b'A', b'C', invalid, b'T', b'G'];
            let mut windows = KmerWindows::new(&sequence, 2).unwrap();
            assert!(windows.next().unwrap().is_ok());
            assert!(windows.next().unwrap().is_err());
            assert!(windows.next().is_none());
            assert!(windows.next().is_none());
        }
        let starts = KmerWindows::new(b"ACNNTG", 2)
            .unwrap()
            .map(|w| w.unwrap().start)
            .collect::<Vec<_>>();
        assert_eq!(starts, vec![0, 4]);
    }
}
