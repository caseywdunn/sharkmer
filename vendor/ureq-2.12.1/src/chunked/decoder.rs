// Copyright 2015 The tiny-http Contributors
// Copyright 2015 The rust-chunked-transfer Contributors
// Forked into ureq, 2024, from https://github.com/frewsxcv/rust-chunked-transfer
// Forked under dual MIT and Apache 2.0 license (see adjacent LICENSE-MIT and LICENSE-APACHE file)

use std::error::Error;
use std::fmt;
use std::io::Error as IoError;
use std::io::ErrorKind;
use std::io::Read;
use std::io::Result as IoResult;

/// Reads HTTP chunks and sends back real data.
///
/// # Example
///
/// ```no_compile
/// use chunked_transfer::Decoder;
/// use std::io::Read;
///
/// let encoded = b"3\r\nhel\r\nb\r\nlo world!!!\r\n0\r\n\r\n";
/// let mut decoded = String::new();
///
/// let mut decoder = Decoder::new(encoded as &[u8]);
/// decoder.read_to_string(&mut decoded);
///
/// assert_eq!(decoded, "hello world!!!");
/// ```
pub struct Decoder<R> {
    // where the chunks come from
    source: R,

    // remaining size of the chunk being read
    // none if we are not in a chunk
    remaining_chunks_size: Option<usize>,
}

// Stream the size token instead of retaining a potentially unbounded line.
// Preserve the historical UTF-8/trim/from_str_radix grammar, including a plus
// sign, arbitrary leading zeroes and Unicode whitespace. Invalid tokens still
// drain through CRLF so framing and source-I/O error precedence stay unchanged.
#[derive(Clone, Copy)]
enum SizePhase {
    Leading,
    Plus,
    Digits,
    Trailing,
    Invalid,
}

struct ChunkSize {
    value: usize,
    phase: SizePhase,
    utf8: [u8; 4],
    utf8_len: usize,
}

impl ChunkSize {
    fn new() -> Self {
        Self {
            value: 0,
            phase: SizePhase::Leading,
            utf8: [0; 4],
            utf8_len: 0,
        }
    }

    fn push(&mut self, byte: u8) {
        if matches!(self.phase, SizePhase::Invalid) {
            return;
        }
        if self.utf8_len == 0 && byte.is_ascii() {
            self.character(char::from(byte));
            return;
        }
        self.utf8[self.utf8_len] = byte;
        self.utf8_len += 1;
        match std::str::from_utf8(&self.utf8[..self.utf8_len]) {
            Ok(text) => {
                if let Some(character) = text.chars().next() {
                    self.utf8_len = 0;
                    self.character(character);
                }
            }
            Err(error) if error.error_len().is_none() && self.utf8_len < 4 => {}
            Err(_) => self.phase = SizePhase::Invalid,
        }
    }

    fn character(&mut self, character: char) {
        if character.is_whitespace() {
            self.phase = match self.phase {
                SizePhase::Leading => SizePhase::Leading,
                SizePhase::Digits | SizePhase::Trailing => SizePhase::Trailing,
                _ => SizePhase::Invalid,
            };
            return;
        }
        if character == '+' && matches!(self.phase, SizePhase::Leading) {
            self.phase = SizePhase::Plus;
            return;
        }
        let digit = if character.is_ascii() {
            character.to_digit(16)
        } else {
            None
        };
        let next = digit
            .filter(|_| {
                matches!(
                    self.phase,
                    SizePhase::Leading | SizePhase::Plus | SizePhase::Digits
                )
            })
            .and_then(|digit| self.value.checked_mul(16)?.checked_add(digit as usize));
        if let Some(value) = next {
            self.value = value;
            self.phase = SizePhase::Digits;
            return;
        }
        self.phase = SizePhase::Invalid;
    }

    fn finish(self) -> Option<usize> {
        (self.utf8_len == 0 && matches!(self.phase, SizePhase::Digits | SizePhase::Trailing))
            .then_some(self.value)
    }
}

impl<R> Decoder<R>
where
    R: Read,
{
    pub fn new(source: R) -> Decoder<R> {
        Decoder {
            source,
            remaining_chunks_size: None,
        }
    }

    /// Unwraps the Decoder into its inner `Read` source.
    pub fn into_inner(self) -> R {
        self.source
    }

    fn read_chunk_size(&mut self) -> IoResult<usize> {
        let mut chunk_size = ChunkSize::new();
        let mut has_ext = false;

        loop {
            let byte = match self.source.by_ref().bytes().next() {
                Some(b) => b?,
                None => return Err(IoError::new(ErrorKind::InvalidInput, DecoderError)),
            };

            if byte == b'\r' {
                break;
            }

            if byte == b';' {
                has_ext = true;
                break;
            }

            chunk_size.push(byte);
        }

        // Ignore extensions for now
        if has_ext {
            loop {
                let byte = match self.source.by_ref().bytes().next() {
                    Some(b) => b?,
                    None => return Err(IoError::new(ErrorKind::InvalidInput, DecoderError)),
                };
                if byte == b'\r' {
                    break;
                }
            }
        }

        self.read_line_feed()?;

        let chunk_size = chunk_size
            .finish()
            .ok_or_else(|| IoError::new(ErrorKind::InvalidInput, DecoderError))?;

        Ok(chunk_size)
    }

    fn read_carriage_return(&mut self) -> IoResult<()> {
        match self.source.by_ref().bytes().next() {
            Some(Ok(b'\r')) => Ok(()),
            _ => Err(IoError::new(ErrorKind::InvalidInput, DecoderError)),
        }
    }

    fn read_line_feed(&mut self) -> IoResult<()> {
        match self.source.by_ref().bytes().next() {
            Some(Ok(b'\n')) => Ok(()),
            _ => Err(IoError::new(ErrorKind::InvalidInput, DecoderError)),
        }
    }

    // Sometimes the last \r\n is missing.
    fn read_end(&mut self) -> IoResult<()> {
        fn expect_or_end(
            bytes: &mut impl Iterator<Item = IoResult<u8>>,
            expected: u8,
        ) -> IoResult<()> {
            match bytes.next() {
                Some(Ok(c)) => {
                    if c == expected {
                        Ok(())
                    } else {
                        Err(IoError::new(ErrorKind::InvalidInput, DecoderError))
                    }
                }
                Some(Err(e)) => {
                    match e.kind() {
                        // Closed connections are ok.
                        ErrorKind::ConnectionReset | ErrorKind::ConnectionAborted => Ok(()),
                        _ => Err(IoError::new(ErrorKind::InvalidInput, DecoderError)),
                    }
                }
                None => Ok(()), // End of iterator is ok
            }
        }

        let mut bytes = self.source.by_ref().bytes();

        expect_or_end(&mut bytes, b'\r')?;
        expect_or_end(&mut bytes, b'\n')?;

        Ok(())
    }
}

impl<R> Read for Decoder<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> IoResult<usize> {
        let remaining_chunks_size = match self.remaining_chunks_size {
            Some(c) => c,
            None => {
                // first possibility: we are not in a chunk, so we'll attempt to determine
                // the chunks size
                let chunk_size = self.read_chunk_size()?;

                // if the chunk size is 0, we are at EOF
                if chunk_size == 0 {
                    self.read_end()?;
                    return Ok(0);
                }

                chunk_size
            }
        };

        // second possibility: we continue reading from a chunk
        if buf.len() < remaining_chunks_size {
            let read = self.source.read(buf)?;
            self.remaining_chunks_size = Some(remaining_chunks_size - read);
            return Ok(read);
        }

        // third possibility: the read request goes further than the current chunk
        // we simply read until the end of the chunk and return
        assert!(buf.len() >= remaining_chunks_size);

        let buf = &mut buf[..remaining_chunks_size];
        let read = self.source.read(buf)?;

        self.remaining_chunks_size = if read == remaining_chunks_size {
            self.read_carriage_return()?;
            self.read_line_feed()?;
            None
        } else {
            Some(remaining_chunks_size - read)
        };

        Ok(read)
    }
}

#[derive(Debug, Copy, Clone)]
struct DecoderError;

impl fmt::Display for DecoderError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(fmt, "Error while decoding chunks")
    }
}

impl Error for DecoderError {
    fn description(&self) -> &str {
        "Error while decoding chunks"
    }
}

#[cfg(test)]
mod test {
    use super::Decoder;
    use std::io;
    use std::io::Read;

    /// This unit test is taken from from Hyper
    /// https://github.com/hyperium/hyper
    /// Copyright (c) 2014 Sean McArthur
    #[test]
    fn test_read_chunk_size() {
        fn read(s: &str, expected: usize) {
            let mut decoded = Decoder::new(s.as_bytes());
            let actual = decoded.read_chunk_size().unwrap();
            assert_eq!(expected, actual);
        }

        fn read_err(s: &str) {
            let mut decoded = Decoder::new(s.as_bytes());
            let err_kind = decoded.read_chunk_size().unwrap_err().kind();
            assert_eq!(err_kind, io::ErrorKind::InvalidInput);
        }

        read("1\r\n", 1);
        read("01\r\n", 1);
        read("0\r\n", 0);
        read("00\r\n", 0);
        read("A\r\n", 10);
        read("a\r\n", 10);
        read("Ff\r\n", 255);
        read("Ff   \r\n", 255);
        // Missing LF or CRLF
        read_err("F\rF");
        read_err("F");
        // Invalid hex digit
        read_err("X\r\n");
        read_err("1X\r\n");
        read_err("-\r\n");
        read_err("-1\r\n");
        // Acceptable (if not fully valid) extensions do not influence the size
        read("1;extension\r\n", 1);
        read("a;ext name=value\r\n", 10);
        read("1;extension;extension2\r\n", 1);
        read("1;;;  ;\r\n", 1);
        read("2; extension...\r\n", 2);
        read("3   ; extension=123\r\n", 3);
        read("3   ;\r\n", 3);
        read("3   ;   \r\n", 3);
        // Invalid extensions cause an error
        read_err("1 invalid extension\r\n");
        read_err("1 A\r\n");
        read_err("1;no CRLF");
    }

    #[test]
    fn test_valid_chunk_decode() {
        let source = io::Cursor::new(
            "3\r\nhel\r\nb\r\nlo world!!!\r\n0\r\n\r\n"
                .to_string()
                .into_bytes(),
        );
        let mut decoded = Decoder::new(source);

        let mut string = String::new();
        decoded.read_to_string(&mut string).unwrap();

        assert_eq!(string, "hello world!!!");
    }

    #[test]
    fn test_decode_zero_length() {
        let mut decoder = Decoder::new(b"0\r\n\r\n" as &[u8]);

        let mut decoded = String::new();
        decoder.read_to_string(&mut decoded).unwrap();

        assert_eq!(decoded, "");
    }

    #[test]
    fn test_decode_invalid_chunk_length() {
        let mut decoder = Decoder::new(b"m\r\n\r\n" as &[u8]);

        let mut decoded = String::new();
        assert!(decoder.read_to_string(&mut decoded).is_err());
    }

    #[test]
    fn invalid_input1() {
        let source = io::Cursor::new(
            "2\r\nhel\r\nb\r\nlo world!!!\r\n0\r\n"
                .to_string()
                .into_bytes(),
        );
        let mut decoded = Decoder::new(source);

        let mut string = String::new();
        assert!(decoded.read_to_string(&mut string).is_err());
    }

    #[test]
    fn invalid_input2() {
        let source = io::Cursor::new(
            "3\rhel\r\nb\r\nlo world!!!\r\n0\r\n"
                .to_string()
                .into_bytes(),
        );
        let mut decoded = Decoder::new(source);

        let mut string = String::new();
        assert!(decoded.read_to_string(&mut string).is_err());
    }

    #[test]
    fn test_decode_end_missing_last_crlf() {
        // This has been observed in the wild.
        // See https://github.com/algesten/ureq/issues/325

        // Missing last \r\n
        let source = io::Cursor::new(
            "3\r\nhel\r\nb\r\nlo world!!!\r\n0\r\n"
                .to_string()
                .into_bytes(),
        );
        let mut decoded = Decoder::new(source);

        let mut string = String::new();
        decoded.read_to_string(&mut string).unwrap();

        assert_eq!(string, "hello world!!!");
    }
}
