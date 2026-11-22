//! One bounded binding for the Linux filesystem-object identity syscall.
//!
//! rustix 1.1.4 has no safe name_to_handle_at adapter. Endpoint names cannot
//! identify a multipath NFS object; statfs IDs may be zero. Keep this necessary
//! FFI boundary isolated rather than weakening recovery to pathname identity.

use super::*;
use std::ffi::CString;
use std::os::unix::ffi::OsStrExt;

const CAPACITY: usize = libc::MAX_HANDLE_SZ as usize;

#[repr(C)]
struct Buffer {
    header: libc::file_handle,
    data: [u8; CAPACITY],
}

// The C flexible-array payload must immediately follow the public libc header.
const _: () = assert!(std::mem::offset_of!(Buffer, data) == size_of::<libc::file_handle>());

fn frame(kind: i32, length: usize, bytes: &[u8]) -> Result<[u8; 32]> {
    ensure!(
        length > 0 && length <= CAPACITY,
        "Unsupported NFS object-handle size; registration preserved"
    );
    let bytes = bytes.get(..length).context("Truncated NFS object handle")?;
    let mut hash = Sha256::new();
    hash.update(b"sharkmer-workspace-nfs-object-v3\0");
    hash.update(kind.to_le_bytes());
    hash.update(u32::try_from(length)?.to_le_bytes());
    hash.update(bytes);
    Ok(hash.finalize().into())
}

/// Query only identity; this never opens by handle, changes files or raises privilege.
pub(super) fn fingerprint(parent: &Path) -> Result<[u8; 32]> {
    // The caller retains a reservation for this exact NUL-terminated copy.
    let path = CString::new(parent.as_os_str().as_bytes())?;
    let mut buffer = Buffer {
        header: libc::file_handle {
            handle_bytes: libc::MAX_HANDLE_SZ as libc::c_uint,
            handle_type: 0,
            f_handle: [],
        },
        data: [0; CAPACITY],
    };
    let mut mount_id: libc::c_int = 0;
    loop {
        // SAFETY: repr(C) plus the assertion above gives libc's header followed
        // by CAPACITY initialized writable bytes. The pointer derives from the
        // WHOLE buffer, not a header-only borrow; the kernel is told precisely
        // that capacity. path remains a valid C string throughout the call;
        // mount_id is an exclusive writable c_int. No pointer escapes. Oversize
        // handles fail without retrying/growing, and success is bounds-checked
        // before indexing. mount_id is deliberately excluded from identity.
        let result = unsafe {
            libc::name_to_handle_at(
                libc::AT_FDCWD,
                path.as_ptr(),
                std::ptr::from_mut(&mut buffer).cast(),
                &mut mount_id,
                0,
            )
        };
        if result == 0 {
            break;
        }
        let error = std::io::Error::last_os_error();
        if error.kind() == std::io::ErrorKind::Interrupted {
            continue;
        }
        return Err(error).context("Unable to query persistent NFS object handle");
    }
    frame(
        buffer.header.handle_type,
        usize::try_from(buffer.header.handle_bytes)?,
        &buffer.data,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn buffer_layout_and_frame_reject_bad_lengths_and_distinguish_objects() {
        assert_eq!(
            std::mem::offset_of!(Buffer, data),
            size_of::<libc::file_handle>()
        );
        assert!(std::mem::align_of::<Buffer>() >= std::mem::align_of::<libc::file_handle>());
        assert_eq!(
            size_of::<Buffer>(),
            size_of::<libc::file_handle>() + CAPACITY
        );
        for length in [0, CAPACITY + 1, usize::MAX] {
            assert!(frame(8, length, &[0; CAPACITY]).is_err());
        }
        assert!(frame(8, 4, &[0; 3]).is_err());
        let first = frame(8, 4, b"abcd").unwrap();
        assert_ne!(first, frame(9, 4, b"abcd").unwrap());
        assert_ne!(first, frame(8, 3, b"abcd").unwrap());
        assert_ne!(first, frame(8, 4, b"abce").unwrap());
        assert_eq!(first, frame(8, 4, b"abcd ignored capacity").unwrap());
    }
}
