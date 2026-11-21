//! Workspace parent identity from kernel objects, independent of NFS endpoints.
//!
//! NFS object handles survive alternate server addresses and client mount/device
//! IDs. Local storage also binds the boot/device/inode and optional birth time,
//! so another host's temporary directory cannot stand in for the original.

use super::*;

#[cfg(target_os = "linux")]
mod handle;

/// Preadmitted pathname conversion workspace survives live and cold cleanup.
pub(super) struct Identity {
    digest: [u8; 32],
    _path_workspace: Option<crate::kmer::budget::MemoryReservation>,
}

impl Identity {
    pub(super) fn capture(parent: &BudgetPath) -> Result<Self> {
        let path_bytes = parent
            .as_os_str()
            .len()
            .checked_add(1)
            .context("Workspace path workspace overflow")?;
        let path_workspace = parent
            .budget()
            .map(|budget| budget.reserve(MemoryClass::RetainedResults, u64::try_from(path_bytes)?))
            .transpose()?;
        let digest = observe(parent)?;
        Ok(Self {
            digest,
            _path_workspace: path_workspace,
        })
    }

    pub(super) fn recorded(parent: &BudgetPath, digest: [u8; 32]) -> Result<Self> {
        let identity = Self::capture(parent)?;
        ensure!(
            identity.digest == digest,
            "Workspace storage identity changed; registration and private data preserved. Retry with the original temporary storage available."
        );
        Ok(identity)
    }

    pub(super) fn digest(&self) -> &[u8; 32] {
        &self.digest
    }

    pub(super) fn verify(&self, parent: &BudgetPath) -> Result<()> {
        ensure!(
            observe(parent)? == self.digest,
            "Workspace storage identity changed; registration and private data preserved. Retry with the original temporary storage available."
        );
        Ok(())
    }
}

#[cfg(target_os = "linux")]
fn observe(parent: &Path) -> Result<[u8; 32]> {
    use rustix::fs::{AtFlags, CWD, StatxFlags, statx};
    use std::os::unix::fs::MetadataExt;
    let metadata = fs::symlink_metadata(parent)
        .context("Workspace parent storage unavailable; registration and private data preserved")?;
    ensure!(
        metadata.file_type().is_dir(),
        "Workspace parent is not a directory; registration and private data preserved"
    );
    let kind = rustix::fs::statfs(parent)?.f_type;
    let digest = if kind == rustix::fs::NFS_SUPER_MAGIC {
        handle::fingerprint(parent).context(
            "Unable to identify NFS workspace object; registration and private data preserved",
        )?
    } else {
        let extended = match statx(
            CWD,
            parent,
            AtFlags::SYMLINK_NOFOLLOW,
            StatxFlags::INO | StatxFlags::BTIME,
        ) {
            Ok(value) => Some(value),
            Err(
                rustix::io::Errno::NOSYS | rustix::io::Errno::OPNOTSUPP | rustix::io::Errno::INVAL,
            ) => None,
            Err(error) => return Err(error).context("Unable to inspect workspace object identity"),
        };
        if let Some(value) = &extended {
            ensure!(
                value.stx_ino == metadata.ino(),
                "Workspace parent changed while observing identity"
            );
        }
        let mut boot = [0; 37];
        File::open("/proc/sys/kernel/random/boot_id")?.read_exact(&mut boot)?;
        ensure!(
            boot[36] == b'\n'
                && boot[..36]
                    .iter()
                    .enumerate()
                    .all(|(i, b)| if matches!(i, 8 | 13 | 18 | 23) {
                        *b == b'-'
                    } else {
                        b.is_ascii_hexdigit()
                    }),
            "Invalid kernel boot identity"
        );
        let mut hash = Sha256::new();
        hash.update(b"sharkmer-workspace-local-linux-v3\0");
        hash.update(boot);
        hash.update(metadata.dev().to_le_bytes());
        hash.update(metadata.ino().to_le_bytes());
        if let Some(value) = &extended
            && value.stx_mask & StatxFlags::BTIME.bits() != 0
        {
            hash.update(value.stx_btime.tv_sec.to_le_bytes());
            hash.update(value.stx_btime.tv_nsec.to_le_bytes());
        }
        hash.finalize().into()
    };
    let after = fs::symlink_metadata(parent)?;
    ensure!(
        after.file_type().is_dir()
            && after.ino() == metadata.ino()
            && after.dev() == metadata.dev(),
        "Workspace parent changed while observing identity"
    );
    Ok(digest)
}

#[cfg(not(target_os = "linux"))]
fn observe(parent: &Path) -> Result<[u8; 32]> {
    let metadata = fs::symlink_metadata(parent)
        .context("Workspace parent storage unavailable; registration and private data preserved")?;
    ensure!(
        metadata.file_type().is_dir(),
        "Workspace parent is not a directory; registration and private data preserved"
    );
    let mut digest = Sha256::new();
    digest.update(b"sharkmer-workspace-parent-other-v2\0");
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        digest.update(rustix::system::uname().nodename().to_bytes());
        digest.update(metadata.dev().to_le_bytes());
        digest.update(metadata.ino().to_le_bytes());
    }
    #[cfg(not(unix))]
    {
        let host = std::env::var("COMPUTERNAME").context("Unable to identify workspace host")?;
        digest.update(host.as_bytes());
    }
    // Creation time is stable across child creation/removal, unlike ctime.
    let created = metadata
        .created()
        .context("Unable to identify workspace directory creation time")?
        .duration_since(std::time::UNIX_EPOCH)?;
    digest.update(created.as_secs().to_le_bytes());
    digest.update(created.subsec_nanos().to_le_bytes());
    Ok(digest.finalize().into())
}

#[cfg(test)]
mod tests;
