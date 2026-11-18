use super::*;
use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};
use clap::Parser;
use std::io::{BufRead, Cursor, Write};
use std::net::TcpListener;

fn budget(memory: u64) -> Arc<ResourceBudget> {
    ResourceBudget::new(ResourceLimits {
        memory_bytes: memory.max(1),
        process_overhead_bytes: 0,
        class_bytes: [memory; MEMORY_CLASSES],
        temporary_disk_bytes: 0,
        temporary_open_files: 0,
    })
    .unwrap()
}

fn values(result: &EnaResult) -> (Vec<String>, Option<String>) {
    (
        result.urls.iter().map(|s| s.as_str().to_owned()).collect(),
        result.scientific_name.as_deref().map(str::to_owned),
    )
}

// Frozen v3/v4-pre-migration selection: whole-response lines/column vectors,
// first data row only, first duplicate header, existing URL prefix semantics.
fn historical(body: &[u8]) -> Result<(Vec<String>, Option<String>)> {
    let body = String::from_utf8_lossy(body);
    let lines: Vec<_> = body.lines().collect();
    ensure!(lines.len() >= 2, "no results");
    let headers: Vec<_> = lines[0].split('\t').collect();
    let ftp = headers
        .iter()
        .position(|h| *h == "fastq_ftp")
        .context("missing ftp")?;
    let scientific = headers.iter().position(|h| *h == "scientific_name");
    let fields: Vec<_> = lines[1].split('\t').collect();
    ensure!(fields.len() > ftp && !fields[ftp].is_empty(), "no URLs");
    let urls = fields[ftp]
        .split(';')
        .map(|u| {
            if u.starts_with("ftp://") || u.starts_with("http") {
                u.to_owned()
            } else {
                format!("http://{u}")
            }
        })
        .collect();
    let name = scientific
        .and_then(|index| fields.get(index))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(str::to_owned);
    Ok((urls, name))
}

#[test]
fn iterator_selection_matches_historical_rows_columns_url_order_and_lossy_text() {
    let cases: &[&[u8]] = &[
        b"", b"fastq_ftp\n", b"unrecognized\nsource\n", b"fastq_ftp\tname\n\tvalue\n",
        b"run_accession\tfastq_ftp\tscientific_name\nSRR1\tftp.sra/a_1.gz;ftp.sra/a_2.gz\tHomo sapiens\nSRR2\tignored.gz\tDifferent species\n",
        b"scientific_name\tfastq_ftp\r\n  Species name  \tftp://host/a;https://host/b;httpx://odd;;host/c\r\n",
        b"fastq_ftp\tfastq_ftp\tscientific_name\tscientific_name\nhost/first\thost/ignored\tfirst\tignored\n",
        b"fastq_ftp\tscientific_name\nhost/a;\t \n", b"fastq_ftp\n;;;\n",
        b"fastq_ftp\tscientific_name\nhost/\xff\tSpecies\xf0\x9f\x92\nignored\xff\n",
        "scientific_name\tfastq_ftp\n\u{2003}Aurelia aurita\u{2003}\thost/a\n".as_bytes(),
        b"fastq_ftp\tscientific_name\n\nignored.gz\tignored\n",
    ];
    for bytes in cases {
        let ledger = budget(2 << 20);
        let actual = read_response(Cursor::new(bytes), "SRR1", Some(&ledger));
        let expected = historical(bytes);
        assert_eq!(actual.is_ok(), expected.is_ok(), "body={bytes:?}");
        if let Ok(result) = &actual {
            assert_eq!(values(result), expected.unwrap());
        }
        drop(actual);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}

#[test]
fn lossy_decoding_matches_standard_replacement_grouping_and_exact_lengths() {
    for first in 0..=255 {
        for suffix in [
            &b""[..],
            &b"\x80"[..],
            &b"\xbf\xff"[..],
            &b"\x82\xac"[..],
            &b"\x9f\x92\xa9"[..],
        ] {
            let bytes: Vec<u8> = std::iter::once(first)
                .chain(suffix.iter().copied())
                .collect();
            let expected = String::from_utf8_lossy(&bytes);
            let ledger = budget(1024);
            let actual = decode_lossy(&bytes, Some(&ledger)).unwrap();
            assert_eq!(actual.as_str(), expected);
            assert_eq!(lossy_length(&bytes).unwrap(), expected.len());
            assert_eq!(
                ledger.snapshot().reserved_memory_bytes,
                expected.len() as u64
            );
            drop(actual);
            assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        }
    }
}

#[test]
fn quota_refusal_releases_partial_models_and_transferred_urls_keep_their_owners() {
    let bytes = b"fastq_ftp\tscientific_name\nhost/first;host/second\tSpecies name\n";
    let mut refused = 0;
    let mut succeeded = 0;
    for memory in (0..=1024).step_by(8) {
        let ledger = budget(memory);
        match read_response(Cursor::new(bytes), "SRR1", Some(&ledger)) {
            Err(error) => {
                assert!(
                    error.downcast_ref::<ResourceExhausted>().is_some(),
                    "{error:#}"
                );
                refused += 1;
            }
            Ok(result) => {
                let expected = (result.urls.capacity() * size_of::<BudgetString>()
                    + result.urls.iter().map(|s| s.len()).sum::<usize>()
                    + result.scientific_name.as_ref().map_or(0, |s| s.len()))
                    as u64;
                assert_eq!(ledger.snapshot().reserved_memory_bytes, expected);
                let EnaResult {
                    urls,
                    scientific_name,
                } = result;
                drop(scientific_name);
                let plan = ReadSourcePlan::UncachedRemote(urls);
                assert!(ledger.snapshot().reserved_memory_bytes > 0);
                drop(plan);
                succeeded += 1;
            }
        }
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert_eq!(Arc::strong_count(&ledger), 1);
        assert!(ledger.snapshot().valid);
    }
    assert!(refused > 0 && succeeded > 0);
}

#[test]
fn whole_body_limit_and_tail_io_errors_are_enforced_before_first_row_success() {
    let mut bytes = b"fastq_ftp\nhost/a\n".to_vec();
    bytes.resize(BODY_LIMIT, b'\n');
    let ledger = budget(64 << 20);
    let result = read_response(Cursor::new(&bytes), "SRR1", Some(&ledger)).unwrap();
    assert_eq!(result.urls[0].as_str(), "http://host/a");
    drop(result);
    bytes.push(b'x');
    let error = read_response(Cursor::new(&bytes), "SRR1", Some(&ledger))
        .err()
        .unwrap();
    assert!(format!("{error:#}").contains("response too big for into_string"));
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    struct TailError {
        data: Cursor<&'static [u8]>,
        interrupted: bool,
    }
    impl Read for TailError {
        fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
            if !self.interrupted {
                self.interrupted = true;
                return Err(std::io::ErrorKind::Interrupted.into());
            }
            if self.data.position() == self.data.get_ref().len() as u64 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::PermissionDenied,
                    "late response error",
                ));
            }
            self.data.read(output)
        }
    }
    let error = read_response(
        TailError {
            data: Cursor::new(b"fastq_ftp\nhost/a\n"),
            interrupted: false,
        },
        "SRR1",
        Some(&ledger),
    )
    .err()
    .unwrap();
    assert_eq!(
        error.downcast_ref::<std::io::Error>().unwrap().kind(),
        std::io::ErrorKind::PermissionDenied
    );
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn sample_spelling_and_explicit_sample_query_bypass_are_preserved() {
    let ledger = budget(2 << 20);
    for name in [
        None,
        Some("Homo sapiens"),
        Some("Zoë  alpha"),
        Some("A\u{2003}B"),
    ] {
        let sample = derived_sample("SRR1", name, Some(&ledger)).unwrap();
        let expected = name.map_or_else(
            || "SRR1".to_owned(),
            |name| format!("{}_SRR1", name.replace(' ', "_")),
        );
        assert_eq!(sample.as_str(), expected);
    }
    let args = Args::parse_from([
        "sharkmer",
        "--ena",
        "invalid-accession-must-not-be-queried",
        "--sample",
        "explicit_name",
    ]);
    let (sample, cached) = sample_name(&args, Some(&ledger)).unwrap();
    assert_eq!(sample.as_str(), "explicit_name");
    assert!(cached.is_none());
    drop(sample);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
    let prepared = prepare_ingestion(&args).unwrap();
    assert!(Arc::ptr_eq(
        prepared.input_budget().unwrap(),
        prepared.result_budget().unwrap()
    ));
    let args = Args::parse_from([
        "sharkmer",
        "--dry-run",
        "--sample",
        "preview",
        "--count-table-budget-mib",
        "18446744073709551615",
    ]);
    let preview = preview_budget(512 << 20).unwrap();
    let (name, cached) = sample_name(&args, Some(&preview)).unwrap();
    assert_eq!(name.as_str(), "preview");
    assert!(cached.is_none());
    assert_eq!(preview.class_limit(CLASS), 256 << 20);
    drop(name);
    assert_eq!(preview.snapshot().reserved_memory_bytes, 0);
}

fn serve_once(body: Vec<u8>) -> (String, std::thread::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}/filereport", listener.local_addr().unwrap());
    let server = std::thread::spawn(move || {
        let (mut stream, _) = listener.accept().unwrap();
        stream
            .set_read_timeout(Some(std::time::Duration::from_secs(5)))
            .unwrap();
        let mut request = std::io::BufReader::new(stream.try_clone().unwrap());
        loop {
            let mut line = String::new();
            request.read_line(&mut line).unwrap();
            if line == "\r\n" {
                break;
            }
        }
        write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nContent-Type: text/plain; charset=iso-8859-1\r\nConnection: close\r\n\r\n", body.len()).unwrap();
        stream.write_all(&body).unwrap();
    });
    (url, server)
}

#[test]
fn real_http_response_uses_owned_body_and_preserves_first_row_and_lossy_decoding() {
    let body =
        b"fastq_ftp\tscientific_name\nhost/a_1.gz;host/a_2.gz\tSpecies\xff\nignored\tignored\n"
            .to_vec();
    let (historical_url, historical_server) = serve_once(body.clone());
    let historical_text = ureq::get(&historical_url)
        .call()
        .unwrap()
        .into_string()
        .unwrap();
    historical_server.join().unwrap();
    let expected = historical(historical_text.as_bytes()).unwrap();
    let (url, server) = serve_once(body);
    let ledger = budget(2 << 20);
    let result = fetch(&url, "SRR1", Some(&ledger)).unwrap();
    server.join().unwrap();
    assert_eq!(values(&result), expected);
    assert!(ledger.snapshot().reserved_memory_bytes > 0);
    drop(result);
    assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
}

#[test]
fn response_allocation_observations_release_body_and_surviving_model_storage() {
    const ENV: &str = "SHARKMER_TEST_ENA_ALLOCATIONS";
    if std::env::var_os(ENV).is_none() {
        let status = std::process::Command::new(std::env::current_exe().unwrap()).args(["--exact", "io::ena::tests::response_allocation_observations_release_body_and_surviving_model_storage", "--nocapture", "--test-threads=1"]).env(ENV, "1").status().unwrap();
        assert!(status.success());
        return;
    }
    let mut ignored = b"fastq_ftp\tscientific_name\nhost/a\tSpecies name\n".to_vec();
    ignored.resize(1 << 20, b'\n');
    let cases = [
        b"fastq_ftp\tscientific_name\nhost/a;host/b\tSpecies name\n".to_vec(),
        ignored,
        b"fastq_ftp\tscientific_name\nhost/\xff;host/b\tSpecies\xf0\x9f\x92\nignored\xff\n"
            .to_vec(),
    ];
    for (index, bytes) in cases.iter().enumerate() {
        let ledger = budget(64 << 20);
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let result = read_response(Cursor::new(bytes), "SRR1", Some(&ledger)).unwrap();
        let observed = crate::PEAK_ALLOC.peak_usage() - baseline;
        let reserved = ledger.snapshot().peak_reserved_memory_bytes;
        assert!(
            observed as u64 <= reserved,
            "case={index} observed={observed} reserved={reserved}"
        );
        drop(result);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        println!("ena_allocator case={index} observed={observed} reserved={reserved}");
    }
}
