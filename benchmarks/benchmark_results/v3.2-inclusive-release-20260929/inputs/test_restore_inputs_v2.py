"""Author synthetic controls; no network, biological prefix or assembler use."""

import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("restore_inputs", HERE / "restore_inputs_v2.py")
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def record(name, sequence="ACGT"):
    return f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()


def expected(raw):
    return {"records": len(raw.splitlines()) // 4,
            "bases": sum(len(line) for line in raw.splitlines()[1::4]),
            "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


class Response(io.BytesIO):
    status = 200
    headers = {}
    url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/test.fastq.gz"


def fixture(parts, eof=False):
    sources, compressed = [], []
    for index, raw in enumerate(parts):
        packed = gzip.compress(raw, mtime=0)
        compressed.append(packed)
        sources.append({"url": f"https://ftp.sra.ebi.ac.uk/vol1/fastq/test{index}.fastq.gz",
                        "records": expected(raw)["records"], "eof_required": eof,
                        "archived_compressed_sha256": hashlib.sha256(packed).hexdigest(),
                        "archived_full_compressed_md5": hashlib.md5(packed).hexdigest(),
                        "archived_full_compressed_size_bytes": len(packed)})
    sample = {"accession": "SRR1", "panel": "test", "expected_prefix": expected(b"".join(parts)),
              "current_file": None, "restoration_sources": sources}
    return sample, compressed


class RestoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)

    def restore(self, sample, compressed):
        responses = iter(compressed)
        return R.restore(sample, self.path / "prefix.partial", [],
                         opener=lambda *_args, **_kwargs: Response(next(responses)))

    def test_prefix_stops_before_additional_records(self):
        raw = record("first")
        sample, _ = fixture([raw])
        self.assertEqual(self.restore(sample, [gzip.compress(raw + record("extra"))]), expected(raw))

    def test_three_complete_sources_in_exact_order(self):
        parts = [record("a"), record("b", "AAC"), record("c", "AT")]
        sample, compressed = fixture(parts, eof=True)
        self.assertEqual(self.restore(sample, compressed), expected(b"".join(parts)))
        self.assertEqual((self.path / "prefix.partial").read_bytes(), b"".join(parts))

    def test_complete_first_source_then_partial_final_source(self):
        parts = [record("r1a") + record("r1b"), record("r2a")]
        sample, compressed = fixture(parts, eof=True)
        sample["restoration_sources"][-1]["eof_required"] = False
        compressed[-1] = gzip.compress(parts[-1] + record("r2extra"))
        R.validate_inventory({"schema_version": 1, "samples": [sample], "panels": [{"name": "test"}]})
        self.assertEqual(self.restore(sample, compressed), expected(b"".join(parts)))

    def test_complete_source_rejects_extra_record(self):
        sample, _ = fixture([record("first")], eof=True)
        with self.assertRaisesRegex(ValueError, "additional"):
            self.restore(sample, [gzip.compress(record("first") + record("extra"))])

    def test_complete_source_checks_compressed_identity(self):
        sample, compressed = fixture([record("first")], eof=True)
        sample["restoration_sources"][0]["archived_full_compressed_md5"] = "0" * 32
        with self.assertRaisesRegex(ValueError, "compressed source differs"):
            self.restore(sample, compressed)

    def test_distinct_local_cache_gzip_sha_is_not_upstream_pin(self):
        raw = record("first")
        sample, compressed = fixture([raw], eof=True)
        local_cache = gzip.compress(raw, mtime=12345)
        sample["restoration_sources"][0]["archived_compressed_sha256"] = hashlib.sha256(local_cache).hexdigest()
        self.assertNotEqual(sample["restoration_sources"][0]["archived_compressed_sha256"], hashlib.sha256(compressed[0]).hexdigest())
        transports = []
        observed = R.restore(sample, self.path / "prefix.partial", transports,
                             opener=lambda *a, **k: Response(compressed[0]))
        self.assertEqual(observed, expected(raw))
        self.assertTrue(transports[0]["upstream_md5_and_size_verified"])
        self.assertFalse(transports[0]["historical_cache_sha256_is_upstream_pin"])

    def test_actual_remote_size_mismatch_rejected(self):
        sample, compressed = fixture([record("first")], eof=True)
        sample["restoration_sources"][0]["archived_full_compressed_size_bytes"] += 1
        with self.assertRaisesRegex(ValueError, "ENA MD5/size"):
            self.restore(sample, compressed)

    def test_wrong_raw_digest_rejected(self):
        sample, compressed = fixture([record("first")])
        sample["expected_prefix"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "differs"):
            self.restore(sample, compressed)

    def test_truncated_and_invalid_quality_rejected(self):
        for raw in (b"@a\nAC\n+\nI", b"@a\nAC\n+\nI \n"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                R.consume(io.BytesIO(raw), 1, R.PrefixCounter(expected(record("a", "AC"))))

    def test_byte_count_bound_rejected(self):
        sample, compressed = fixture([record("first")])
        sample["expected_prefix"]["size_bytes"] -= 1
        with self.assertRaisesRegex(ValueError, "exceeds"):
            self.restore(sample, compressed)

    def test_compressed_read_limit_rejected(self):
        with mock.patch.object(R, "MAX_COMPRESSED_BYTES", 2), self.assertRaisesRegex(ValueError, "limit"):
            R.BoundedResponse(io.BytesIO(b"abcd")).read(10)

    def test_existing_prefix_exact_and_extra_bytes(self):
        raw = record("first")
        path = self.path / "existing.fastq"
        path.write_bytes(raw)
        self.assertEqual(R.verify_existing(path, expected(raw))[0], expected(raw))
        path.write_bytes(raw + b"\n")
        with self.assertRaisesRegex(ValueError, "byte count"):
            R.verify_existing(path, expected(raw))

    def test_symlink_rejected(self):
        target = self.path / "target"
        target.write_bytes(record("a"))
        link = self.path / "link"
        link.symlink_to(target)
        with self.assertRaisesRegex(ValueError, "non-symlink"):
            R.verify_existing(link, expected(target.read_bytes()))

    def test_failed_attempt_retains_partial_and_receipt(self):
        sample, compressed = fixture([record("first")])
        sample["expected_prefix"]["sha256"] = "0" * 64
        result = R.restore_one(sample, self.path, {}, opener=lambda *a, **k: Response(compressed[0]))
        self.assertEqual(result["status"], "failed")
        self.assertTrue((self.path / "SRR1.fastq.partial").is_file())
        self.assertTrue((self.path / "SRR1.restore.json").is_file())
        self.assertFalse((self.path / "SRR1.fastq").exists())
        with self.assertRaisesRegex(ValueError, "exists"):
            R.restore_one(sample, self.path, {})

    def test_successful_existing_attempt_source_binding(self):
        sample, _ = fixture([record("first")])
        path = self.path / "existing.fastq"
        path.write_bytes(record("first"))
        sample["current_file"] = {"path": str(path)}
        frozen = {"helper": R.fingerprint(R.__file__)}
        result = R.restore_one(sample, self.path, frozen)
        self.assertEqual(result["status"], "verified")
        self.assertEqual(result["sources_before"], result["sources_after"])
        self.assertEqual(result["available_path"], str(path))

    def test_changed_source_blocks_publication(self):
        sample, compressed = fixture([record("first")])
        wrong = R.fingerprint(R.__file__)
        wrong["sha256"] = "0" * 64
        result = R.restore_one(sample, self.path, {"helper": wrong},
                               opener=lambda *a, **k: Response(compressed[0]))
        self.assertEqual(result["status"], "failed")
        self.assertFalse((self.path / "SRR1.fastq").exists())

    def test_inventory_rejects_unbounded_or_ambiguous_sources(self):
        sample, _ = fixture([record("a"), record("b")], eof=False)
        document = {"schema_version": 1, "samples": [sample], "panels": [{"name": "test"}]}
        with self.assertRaisesRegex(ValueError, "Concatenated"):
            R.validate_inventory(document)
        sample["expected_prefix"]["records"] = True
        with self.assertRaisesRegex(ValueError, "count"):
            R.validate_inventory(document)

    def test_wrong_external_pin_writes_failure_receipt_without_network(self):
        manifest = self.path / "inventory.json"
        manifest.write_text("{}")
        output = self.path / "attempt"
        result = R.run(manifest, output, "0" * 64)
        self.assertEqual(result["status"], "failed")
        self.assertIn("external SHA256", result["error"])
        self.assertTrue((output / "availability.json").is_file())
        with self.assertRaises(FileExistsError):
            R.run(manifest, output, "0" * 64)


if __name__ == "__main__":
    unittest.main()
