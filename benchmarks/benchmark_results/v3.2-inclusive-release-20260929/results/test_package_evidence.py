"""Tiny packaging controls; no actual benchmark collection or archive publication."""

import io
import json
from pathlib import Path
import tarfile
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import package_evidence as P


class PackageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.scratch = self.root / "bundle", self.root / "scratch"
        self.bundle.mkdir()
        self.scratch.mkdir()
        helper, source = self.root / "prior.py", self.root / "packager.py"
        helper.write_text("# synthetic prior helper\n")
        source.write_text("# synthetic packager source\n")
        patch = mock.patch.multiple(P, ROOT=self.root, HERE=self.bundle, SCRATCH=self.scratch,
                                    PRIOR_HELPER=helper, __file__=str(source))
        patch.start()
        self.addCleanup(patch.stop)
        patch = mock.patch.object(P.BASE, "ROOT", self.root)
        patch.start()
        self.addCleanup(patch.stop)

    def write(self, relative, content):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(content) if isinstance(content, (dict, list)) else content)
        return path

    def fixture(self):
        versions = {}
        required_review = []
        for role in ("baseline", "candidate"):
            root = self.scratch / role / "source"
            cargo = self.write(f"scratch/{role}/source/Cargo.toml", "[package]\n")
            fastq = self.write(f"scratch/{role}/source/tests/tiny.fastq", "@x\nAC\n+\nII\n")
            files = {str(path.relative_to(root)): P.fingerprint(path) for path in (cargo, fastq)}
            snapshot = self.write(f"scratch/{role}/snapshot.json", {"source": str(root), "files": files})
            build = self.write(f"scratch/{role}/build.json", {"status": "complete", "returncode": 0,
                                "source_before": files, "source_after": files})
            versions[role] = {"snapshot": P.fingerprint(snapshot), "build": P.fingerprint(build)}
            required_review += [snapshot, build]
        builds = self.write("scratch/builds.json", {"status": "pass", "versions": versions})
        marker = P.fingerprint(self.write("scratch/restoration_source.py", "# synthetic\n"))
        availability = self.write("scratch/availability.json", {"status": "pass", "complete_inventory_selected": True,
            "sources_before": {"helper": marker}, "sources_after": {"helper": marker},
            "samples": [{"accession": f"SRR{i}"} for i in range(10)], "source_attempts": []})
        jobs = []
        for index in range(60):
            job = {"status": "pass", "invocation": {"invocation_id": str(index)},
                   "outputs": {"files": {}}, "execution": {"files": {}}}
            path = self.write(f"scratch/runs/{index}.json", job)
            job["receipt"] = P.fingerprint(path)
            jobs.append(job)
        summary = self.write("scratch/summary.json", {"status": "pass"})
        execution = self.write("scratch/execution.json", {"status": "pass", "jobs": jobs,
            "summary": P.fingerprint(summary), "bindings": {str(path): P.fingerprint(path) for path in (builds, availability)}})
        binding = {"execution": P.fingerprint(execution), "summary": P.fingerprint(summary)}
        report = self.write("bundle/report.json", {"status": "pass", "sources_before": binding, "sources_after": binding})
        diagnostic = self.write("scratch/diagnostic/source/lib.rs", "fn synthetic() {}\n")
        files = {"lib.rs": P.fingerprint(diagnostic)}
        frozen = {"source": str(diagnostic.parent), "files": {"lib.rs": {key: files["lib.rs"][key] for key in ("sha256", "size_bytes")}}}
        assessment = self.write("scratch/diagnostic/assessment.json", {"status": "pass"})
        dros = self.write("scratch/diagnostic/replay.receipt.json", {"status": "pass", "before": {}, "after": {}, "after_audit": {},
            "source_freeze": frozen, "assessment": P.fingerprint(assessment)})
        controls = self.write("scratch/diagnostic/controls.json", {"status": "pass", "source_before": files,
            "source_after": files, "steps": [{"status": "pass"}]})
        readme = self.write("bundle/README.md", "Synthetic documentation.\n")
        doc = self.write("results.md", "Synthetic result interpretation.\n")
        required_review += [builds, availability, execution, summary, report, dros, controls, assessment, readme, doc]
        review = self.write("bundle/outcome_review.json", {"status": "pass", "reviewed_artifacts": {
            "arbitrary-" + str(index): P.fingerprint(path) for index, path in enumerate(required_review)}})
        return SimpleNamespace(execution=execution, report=report, builds=builds, availability=availability,
                               dros_receipt=dros, dros_controls=controls, outcome_review=review, doc=[doc], evidence=[])

    def test_complete_collection_includes_exact_source_fixtures_and_failed_logs(self):
        args = self.fixture()
        self.write("scratch/failed-attempt/error.log", "synthetic failure\n")
        self.write("scratch/download.fastq.partial", "RAW EXCLUDED")
        self.write("scratch/source.tar", "ARCHIVE EXCLUDED")
        self.write("scratch/target/large.o", "TARGET EXCLUDED")
        self.write("scratch/sharkmer", "BINARY EXCLUDED")
        selection = P.collect(args, set())
        names = set(selection.files)
        self.assertIn("payload/scratch/baseline/source/tests/tiny.fastq", names)
        self.assertIn("payload/scratch/candidate/source/tests/tiny.fastq", names)
        self.assertIn("payload/scratch/failed-attempt/error.log", names)
        self.assertFalse(any("download.fastq" in name or name.endswith("source.tar") or "/target/" in name for name in names))

    def test_missing_reviewed_document_rejected(self):
        args = self.fixture()
        review = json.loads(args.outcome_review.read_text())
        review["reviewed_artifacts"] = {key: row for key, row in review["reviewed_artifacts"].items()
                                        if row["path"] != str(args.doc[0])}
        args.outcome_review.write_text(json.dumps(review))
        with self.assertRaisesRegex(ValueError, "does not bind required"):
            P.collect(args, set())

    def test_changed_bound_payload_rejected(self):
        args = self.fixture()
        self.write("scratch/baseline/source/Cargo.toml", "changed")
        with self.assertRaisesRegex(ValueError, "Bound identity differs"):
            P.collect(args, set())

    def test_incomplete_main_execution_rejected(self):
        args = self.fixture()
        execution = json.loads(args.execution.read_text())
        execution["jobs"].pop()
        args.execution.write_text(json.dumps(execution))
        with self.assertRaisesRegex(ValueError, "sixty"):
            P.collect(args, set())

    def test_tiny_archive_reopen_listing_and_no_clobber(self):
        path = self.write("payload.txt", "exact payload\n")
        selection = P.Selection()
        selection.add(path, "synthetic fixture")
        archive, receipt = self.root / "tiny.tar.gz", self.root / "receipt.json"
        result = P.publish(selection, archive, receipt, {})
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["member_count"], 3)
        with tarfile.open(archive) as stream:
            manifest = json.load(stream.extractfile("MANIFEST.json"))
            self.assertIn(P.fingerprint(path)["sha256"].encode(), stream.extractfile("SHA256SUMS").read())
        self.assertEqual(P.verify_archive(archive, manifest), result["archive"])
        before = archive.read_bytes()
        with self.assertRaisesRegex(ValueError, "fresh distinct"):
            P.publish(selection, archive, receipt, {})
        self.assertEqual(archive.read_bytes(), before)

    def test_reopen_rejects_corrupt_or_duplicate_member(self):
        path = self.write("payload.txt", "correct")
        selection = P.Selection()
        selection.add(path, "synthetic fixture")
        manifest = {"members": selection.files}
        name = next(iter(selection.files))
        for duplicate in (False, True):
            archive = self.root / f"bad-{duplicate}.tar.gz"
            with tarfile.open(archive, "w:gz") as stream:
                rows = [("MANIFEST.json", json.dumps(manifest).encode()), ("SHA256SUMS", P.sha_listing(manifest)),
                        (name, b"correct" if duplicate else b"corrupt")]
                if duplicate:
                    rows.append((name, b"correct"))
                for member_name, raw in rows:
                    member = tarfile.TarInfo(member_name)
                    member.size = len(raw)
                    stream.addfile(member, io.BytesIO(raw))
            with self.assertRaises(ValueError):
                P.verify_archive(archive, manifest)

    def test_symlink_rejected(self):
        path = self.write("real.txt", "payload")
        link = self.root / "link.txt"
        link.symlink_to(path)
        with self.assertRaisesRegex(ValueError, "symlink"):
            P.exclusion(link)


if __name__ == "__main__":
    unittest.main()
