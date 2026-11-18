//! Compare the streamed candidate with the pinned, unmodified decoder.
#[rustfmt::skip]
// The client supplies buffered transport; fixtures use in-memory readers.
#[allow(clippy::unbuffered_bytes)]
#[path = "../vendor/ureq-2.12.1/src/chunked/decoder.rs"]
mod candidate;
#[rustfmt::skip]
#[allow(clippy::all)]
#[path = "fixtures/ureq-2.12.1/decoder.rs"]
mod legacy;

use std::io::{self, Cursor, Read};

#[global_allocator]
static HEAP: peak_alloc::PeakAlloc = peak_alloc::PeakAlloc;

#[derive(Debug, PartialEq)]
struct Outcome {
    bytes: Vec<u8>,
    error: Option<(io::ErrorKind, String)>,
    consumed: u64,
}

fn compare(token: &[u8]) {
    let mut wire = token.to_vec();
    wire.extend_from_slice(b"\r\nABCD\r\n0\r\n\r\n");
    fn read(mut decoder: impl Read) -> (Vec<u8>, Option<(io::ErrorKind, String)>) {
        let mut buffer = [0; 8];
        match decoder.read(&mut buffer) {
            Ok(n) => (buffer[..n].to_vec(), None),
            Err(error) => (vec![], Some((error.kind(), error.to_string()))),
        }
    }
    let mut original = legacy::Decoder::new(Cursor::new(wire.as_slice()));
    let (bytes, error) = read(&mut original);
    let expected = Outcome {
        bytes,
        error,
        consumed: original.into_inner().position(),
    };
    let mut streamed = candidate::Decoder::new(Cursor::new(wire.as_slice()));
    let (bytes, error) = read(&mut streamed);
    let actual = Outcome {
        bytes,
        error,
        consumed: streamed.into_inner().position(),
    };
    assert_eq!(actual, expected, "token={token:?}");
}

#[test]
fn framing_and_errors_match_for_exhaustive_short_tokens() {
    let alphabet = b"01aFg+ -\t\n;\xc2\xa0\xff";
    compare(b"");
    for length in 1..=4u32 {
        for mut index in 0..alphabet.len().pow(length) {
            let mut token = [0; 4];
            for byte in &mut token[..length as usize] {
                *byte = alphabet[index % alphabet.len()];
                index /= alphabet.len();
            }
            compare(&token[..length as usize]);
        }
    }
    for byte in 0..=255 {
        compare(&[byte]);
    }
}

#[test]
fn unicode_whitespace_zeroes_signs_and_integer_boundaries_match() {
    for scalar in 0..=0x10ffff {
        if let Some(character) = char::from_u32(scalar)
            && character.is_whitespace()
            && character != '\r'
        {
            for middle in ["1", "+1", "0", "", "+", "1 1", " 1 "] {
                compare(format!("{character}{middle}{character}").as_bytes());
            }
        }
    }
    for length in [8, 16, 32, 63, 64, 65, 4096, 65536] {
        for value in ["1", "f", "+1", "", " ", "g"] {
            compare(format!("{}{value}", "0".repeat(length)).as_bytes());
        }
    }
    for token in [
        format!("{:x}", usize::MAX),
        format!("{:x}0", usize::MAX),
        format!("+{:x}", usize::MAX),
        "+".repeat(100),
    ] {
        compare(token.as_bytes());
    }
    compare(b"1;\xff\xfe");
    compare(b"\xc2;ignored");
}

#[test]
fn long_size_and_extension_stream_without_retaining_input() {
    const ENV: &str = "SHARKMER_TEST_CHUNK_FRAMING_ALLOCATIONS";
    if std::env::var_os(ENV).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "long_size_and_extension_stream_without_retaining_input",
                "--nocapture",
                "--test-threads=1",
            ])
            .env(ENV, "1")
            .status()
            .unwrap();
        assert!(status.success());
        return;
    }
    let baseline = HEAP.current_usage();
    HEAP.reset_peak_usage();
    let zeros = io::repeat(b'0').take(384 << 20);
    let input = zeros.chain(Cursor::new(b"1\r\nx\r\n0\r\n\r\n"));
    let mut decoder = candidate::Decoder::new(input);
    let mut byte = [0];
    assert_eq!(decoder.read(&mut byte).unwrap(), 1);
    assert_eq!(byte, *b"x");
    assert_eq!(decoder.read(&mut byte).unwrap(), 0);
    let input = Cursor::new(b"1;")
        .chain(io::repeat(b'x').take(384 << 20))
        .chain(Cursor::new(b"\r\nx\r\n0\r\n\r\n"));
    let mut decoder = candidate::Decoder::new(input);
    assert_eq!(decoder.read(&mut byte).unwrap(), 1);
    assert_eq!(byte, *b"x");
    let observed = HEAP.peak_usage() - baseline;
    assert_eq!(observed, 0, "valid chunk framing allocated heap storage");
    assert_eq!(HEAP.current_usage(), baseline);
    println!(
        "chunk_framing_allocator observed={observed} streamed_bytes={}",
        2 * (384u64 << 20)
    );
}

#[derive(Clone)]
struct Fallible<'a> {
    bytes: &'a [u8],
    position: usize,
    fail_at: usize,
    failed: bool,
    kind: io::ErrorKind,
}
impl Read for Fallible<'_> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if self.position == self.fail_at && !self.failed {
            self.failed = true;
            return Err(self.kind.into());
        }
        let before_error = if self.failed {
            usize::MAX
        } else {
            self.fail_at - self.position
        };
        let count = output
            .len()
            .min(self.bytes.len() - self.position)
            .min(before_error);
        output[..count].copy_from_slice(&self.bytes[self.position..self.position + count]);
        self.position += count;
        Ok(count)
    }
}

#[test]
fn malformed_tokens_preserve_consumption_and_source_error_precedence() {
    for token in [
        "1",
        "g",
        "+",
        "1 1",
        "\u{2003}+01\u{2003}",
        "\u{fffd}",
        "ffffffffffffffffffffffff",
    ] {
        let wire = format!("{token}\r\nABCD\r\n0\r\n\r\n");
        for fail_at in 0..=wire.len() {
            for kind in [
                io::ErrorKind::PermissionDenied,
                io::ErrorKind::ConnectionReset,
                io::ErrorKind::UnexpectedEof,
            ] {
                let input = Fallible {
                    bytes: wire.as_bytes(),
                    position: 0,
                    fail_at,
                    failed: false,
                    kind,
                };
                let mut original = legacy::Decoder::new(input.clone());
                let mut streamed = candidate::Decoder::new(input);
                let mut expected = [0; 8];
                let mut actual = [0; 8];
                let old = original
                    .read(&mut expected)
                    .map_err(|e| (e.kind(), e.to_string()));
                let new = streamed
                    .read(&mut actual)
                    .map_err(|e| (e.kind(), e.to_string()));
                assert_eq!(new, old, "token={token:?} fail_at={fail_at}");
                assert_eq!(actual, expected);
                assert_eq!(
                    streamed.into_inner().position,
                    original.into_inner().position
                );
            }
        }
    }
}
