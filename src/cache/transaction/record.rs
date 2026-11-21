//! Fixed-size, checksummed cache publication receipts; no decoded heap model.

use super::*;

pub(super) const CONTEXT_BYTES: usize = 169;
pub(super) const JOURNAL_BYTES: usize = 233;
const CONTEXT_MAGIC: &[u8; 8] = b"SHKCCX01";
const JOURNAL_MAGIC: &[u8; 8] = b"SHKCCJ01";

#[derive(Clone, Copy)]
pub(super) struct Record {
    pub(super) key: [u8; 64],
    pub(super) previous: Option<[[u8; 32]; 2]>,
    pub(super) next: Option<[[u8; 32]; 2]>,
}

impl Record {
    pub(super) fn initial(key: &str, public: &[BudgetPath; 2], previous: bool) -> Result<Self> {
        ensure!(valid_key(key.as_bytes()), "Invalid cache transaction key");
        let mut bytes = [0; 64];
        bytes.copy_from_slice(key.as_bytes());
        Ok(Self {
            key: bytes,
            previous: if previous {
                Some([hash(&public[0])?, hash(&public[1])?])
            } else {
                None
            },
            next: None,
        })
    }

    pub(super) fn encode<const N: usize>(&self) -> Result<[u8; N]> {
        ensure!(
            N == CONTEXT_BYTES || N == JOURNAL_BYTES,
            "Invalid cache receipt size"
        );
        ensure!(
            (N == JOURNAL_BYTES) == self.next.is_some(),
            "Incomplete cache receipt"
        );
        let mut bytes = [0; N];
        bytes[..8].copy_from_slice(if N == CONTEXT_BYTES {
            CONTEXT_MAGIC
        } else {
            JOURNAL_MAGIC
        });
        bytes[8..72].copy_from_slice(&self.key);
        if let Some(previous) = &self.previous {
            bytes[72] = 1;
            bytes[73..105].copy_from_slice(&previous[0]);
            bytes[105..137].copy_from_slice(&previous[1]);
        }
        if let Some(next) = &self.next {
            bytes[137..169].copy_from_slice(&next[0]);
            bytes[169..201].copy_from_slice(&next[1]);
        }
        let digest = Sha256::digest(&bytes[..N - 32]);
        bytes[N - 32..].copy_from_slice(&digest);
        Ok(bytes)
    }

    pub(super) fn read<const N: usize>(path: &Path) -> Result<Option<Self>> {
        let Some(metadata) = regular(path)? else {
            return Ok(None);
        };
        ensure!(
            metadata.len() == N as u64 && (N == CONTEXT_BYTES || N == JOURNAL_BYTES),
            "Invalid cache receipt extent; cache preserved"
        );
        let mut bytes = [0; N];
        let mut file = File::open(path)?;
        file.read_exact(&mut bytes)?;
        ensure!(
            file.read(&mut [0])? == 0,
            "Cache receipt changed while reading"
        );
        drop(file);
        ensure!(
            &bytes[..8]
                == if N == CONTEXT_BYTES {
                    CONTEXT_MAGIC
                } else {
                    JOURNAL_MAGIC
                }
                && valid_key(&bytes[8..72])
                && bytes[72] <= 1
                && Sha256::digest(&bytes[..N - 32])[..] == bytes[N - 32..],
            "Damaged cache receipt; cache preserved"
        );
        let mut key = [0; 64];
        key.copy_from_slice(&bytes[8..72]);
        let mut previous = [[0; 32]; 2];
        previous[0].copy_from_slice(&bytes[73..105]);
        previous[1].copy_from_slice(&bytes[105..137]);
        ensure!(
            bytes[72] == 1 || previous == [[0; 32]; 2],
            "Invalid absent cache generation"
        );
        let mut next = [[0; 32]; 2];
        if N == JOURNAL_BYTES {
            next[0].copy_from_slice(&bytes[137..169]);
            next[1].copy_from_slice(&bytes[169..201]);
        }
        Ok(Some(Self {
            key,
            previous: (bytes[72] == 1).then_some(previous),
            next: (N == JOURNAL_BYTES).then_some(next),
        }))
    }
}

pub(super) fn valid_key(key: &[u8]) -> bool {
    key.len() == 64
        && key
            .iter()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(b))
}

pub(super) fn hash(path: &Path) -> Result<[u8; 32]> {
    regular(path)?.context("Missing cache recovery file")?;
    let mut file = File::open(path)?;
    let mut digest = Sha256::new();
    let mut bytes = [0; 8192];
    loop {
        let n = file.read(&mut bytes)?;
        if n == 0 {
            break;
        }
        digest.update(&bytes[..n]);
    }
    Ok(digest.finalize().into())
}
