use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceBudget, ResourceLimits};
use std::io::{self, BufRead, Write};
use std::net::TcpListener;
use std::sync::Arc;

pub(super) fn budget(bytes: u64) -> Arc<ResourceBudget> {
    ResourceBudget::new(ResourceLimits {
        memory_bytes: bytes.max(1),
        process_overhead_bytes: 0,
        class_bytes: [bytes; MEMORY_CLASSES],
        temporary_disk_bytes: 0,
        temporary_open_files: 0,
    })
    .unwrap()
}

fn gzip(data: &[u8], fields: usize, header_crc: bool) -> Vec<u8> {
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    encoder.write_all(data).unwrap();
    let compressed = encoder.finish().unwrap();
    if fields == 0 && !header_crc {
        return compressed;
    }
    let mut bytes = compressed[..10].to_vec();
    bytes[3] = 4 | 8 | 16 | if header_crc { 2 } else { 0 };
    bytes.extend((fields.min(u16::MAX as usize) as u16).to_le_bytes());
    bytes.extend(std::iter::repeat_n(0xfd, fields.min(u16::MAX as usize)));
    for _ in 0..2 {
        bytes.extend(std::iter::repeat_n(0xfe, fields));
        bytes.push(0);
    }
    if header_crc {
        let mut crc = flate2::Crc::new();
        crc.update(&bytes);
        bytes.extend((crc.sum() as u16).to_le_bytes());
    }
    bytes.extend_from_slice(&compressed[10..]);
    bytes
}

fn serve(responses: Vec<(u16, &'static str, Vec<u8>)>) -> (String, std::thread::JoinHandle<()>) {
    serve_owned(
        responses
            .into_iter()
            .map(|(s, h, b)| (s, h.to_owned(), b))
            .collect(),
    )
}

pub(super) fn serve_owned(
    responses: Vec<(u16, String, Vec<u8>)>,
) -> (String, std::thread::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}/reads", listener.local_addr().unwrap());
    let server = std::thread::spawn(move || {
        for (status, headers, body) in responses {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(std::time::Duration::from_secs(5)))
                .unwrap();
            stream
                .set_write_timeout(Some(std::time::Duration::from_secs(5)))
                .unwrap();
            let mut request = io::BufReader::new(stream.try_clone().unwrap());
            loop {
                let mut line = String::new();
                assert!(request.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
            }
            write!(stream, "HTTP/1.1 {status} Test\r\nContent-Length: {}\r\n{headers}Connection: close\r\n\r\n", body.len()).unwrap();
            // Refusal or a malformed eager header can close the connection
            // before all bytes are sent; these fixtures intentionally allow it.
            let _ = stream.write_all(&body);
        }
    });
    (url, server)
}

#[test]
fn automatic_gzip_preserves_members_maximum_headers_and_redirects() {
    let mut body = gzip(b"first\r\n", 65535, true);
    body.extend(gzip(b"second\n", 65535, false));
    let (url, server) = serve(vec![
        (
            302,
            "Location: /final\r\nContent-Encoding: gzip\r\n",
            gzip(b"ignored", 65535, true),
        ),
        (200, "Content-Encoding: x-gzip\r\n", body),
    ]);
    let mut reader = get(&url).unwrap();
    let mut actual = String::new();
    reader.read_to_string(&mut actual).unwrap();
    assert_eq!(actual, "first\r\nsecond\n");
    drop(reader);
    server.join().unwrap();
}

#[test]
fn http_status_and_late_gzip_errors_are_propagated() {
    let (url, server) = serve(vec![(
        404,
        "Content-Encoding: gzip\r\n",
        gzip(b"missing", 65535, true),
    )]);
    let error = get(&url).err().unwrap();
    server.join().unwrap();
    assert!(matches!(
        error.downcast_ref::<ureq::Error>(),
        Some(ureq::Error::Status(404, _))
    ));
    for too_long in [false, true] {
        let mut body = gzip(b"first\n", if too_long { 65536 } else { 65535 }, true);
        if !too_long {
            *body.last_mut().unwrap() ^= 1;
        }
        let (url, server) = serve(vec![(200, "Content-Encoding: gzip\r\n", body)]);
        let mut reader = get(&url).unwrap();
        assert!(io::copy(&mut reader, &mut io::sink()).is_err());
        drop(reader);
        server.join().unwrap();
    }
}

#[test]
fn http_content_decoding_preserves_application_gzip_and_its_buffer_budget() {
    let records = b"@r1\nACGT\n+\n!!!!\n@r2\nTGCA\n+\n####\n";
    let body = gzip(&gzip(records, 65535, true), 65535, true);
    let (url, server) = serve(vec![(200, "Content-Encoding: gzip\r\n", body)]);
    let ledger = budget(2 << 20);
    let reader = get(&url).unwrap();
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    let mut decoder = crate::codec::Decoder::new(reader, Some(&ledger)).unwrap();
    assert!(ledger.snapshot().reserved_memory_bytes > 0);
    let mut actual = Vec::new();
    decoder.read_to_end(&mut actual).unwrap();
    assert_eq!(actual, records);
    drop(decoder);
    server.join().unwrap();
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn default_stream_closes_at_eof_without_extending_the_agents_pool_lifetime() {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let input = format!("http://{}/stream", listener.local_addr().unwrap());
    let (sent, received) = std::sync::mpsc::channel();
    let server = std::thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(std::time::Duration::from_secs(10)))
            .unwrap();
        stream
            .set_write_timeout(Some(std::time::Duration::from_secs(10)))
            .unwrap();
        let mut request = io::BufReader::new(stream.try_clone().unwrap());
        loop {
            let mut line = String::new();
            assert!(request.read_line(&mut line).unwrap() > 0);
            if line == "\r\n" {
                break;
            }
        }
        // Length exceeds the client's 8-KiB eager buffer. HTTP/1.1 without
        // Connection: close leaves a pool-return Weak in the live body stream.
        stream
            .write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 32768\r\n\r\n")
            .unwrap();
        stream.write_all(&[b'x'; 32768]).unwrap();
        assert_eq!(stream.read(&mut [0u8; 1]).unwrap(), 0);
        sent.send(()).unwrap();
    });
    let mut reader = get(&input).unwrap();
    let mut body = Vec::new();
    reader.read_to_end(&mut body).unwrap();
    assert_eq!(body, [b'x'; 32768]);
    // A live body reader must not extend the shortcut Agent pool lifetime.
    received
        .recv_timeout(std::time::Duration::from_secs(10))
        .unwrap();
    drop(reader);
    server.join().unwrap();
}
