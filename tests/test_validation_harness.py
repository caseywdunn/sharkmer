import sys
import subprocess
import json
import yaml
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from sharkmer_validate import blast_references, report, results, runner
import bootstrap_from_runs


class BlastClassificationTests(unittest.TestCase):
    def fixture(self, name):
        return (REPO_ROOT / "tests" / "fixtures" / "validation" / name).read_text()

    def test_wrong_gene_is_never_same_gene_success(self):
        match = blast_references._parse_blast_xml(
            self.fixture("wrong_gene.xml"), "target", "Taxon A"
        )
        self.assertEqual(match.status, "wrong_gene")
        self.assertFalse(match.on_target)
        score = report._score_gene(
            True,
            "target",
            "Taxon A",
            {**match.__dict__, "all_products_confirmed": False},
            {"target": {"Taxon A"}},
        )
        self.assertNotIn(score, {"+**", "+++"})

    def test_short_fragment_cannot_confirm_product(self):
        match = blast_references._parse_blast_xml(
            self.fixture("short_fragment.xml"), "target", "Taxon A"
        )
        self.assertEqual(match.status, "insufficient_alignment")
        self.assertEqual(match.query_coverage_pct, 15.0)
        self.assertFalse(match.on_target)

    def test_complementary_same_gene_references_are_split_evidence(self):
        match = blast_references._parse_blast_xml(
            self.fixture("same_gene_split.xml"), "target", "Taxon A"
        )
        self.assertEqual(match.status, "split_or_chimeric_alignment")
        self.assertTrue(match.split_alignment)
        self.assertFalse(match.on_target)
        self.assertEqual(len(match.alignment_evidence), 2)

    def test_every_product_is_evaluated(self):
        runs = [
            {
                "success": True,
                "genes": [
                    {
                        "gene": "target",
                        "products": [
                            {"product_index": 0, "sequence": "A" * 100},
                            {"product_index": 1, "sequence": "C" * 100},
                        ],
                    }
                ],
            }
        ]
        matches = [
            blast_references.RefBlastResult("confirmed_product", "target", "Taxon A"),
            blast_references.RefBlastResult("wrong_gene", "target", "Taxon A"),
        ]
        with mock.patch.object(
            blast_references, "blast_against_references", side_effect=matches
        ) as blast:
            blast_references.blast_all_products(
                runs, Path("db"), "Taxon A", reference_genes={"target"}
            )
        self.assertEqual(blast.call_count, 2)
        self.assertEqual(
            runs[0]["genes"][0]["products"][1]["reference_match"]["status"],
            "wrong_gene",
        )
        self.assertFalse(runs[0]["genes"][0]["reference_match"]["on_target"])

    def test_stronger_wrong_gene_beats_qualifying_expected_homolog(self):
        xml = """<BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-len>100</Iteration_query-len><Iteration_hits>
        <Hit><Hit_def>target|Taxon_A|T</Hit_def><Hit_hsps><Hsp><Hsp_identity>98</Hsp_identity><Hsp_align-len>100</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>100</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>100</Hsp_hit-to></Hsp></Hit_hsps></Hit>
        <Hit><Hit_def>other|Taxon_A|O</Hit_def><Hit_hsps><Hsp><Hsp_identity>99</Hsp_identity><Hsp_align-len>99</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>99</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>99</Hsp_hit-to></Hsp></Hit_hsps></Hit>
        </Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"""
        match = blast_references._parse_blast_xml(xml, "target", "Taxon A")
        self.assertEqual(match.status, "wrong_gene")

    def test_tied_taxa_are_not_species_confirmation(self):
        xml = """<BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-len>100</Iteration_query-len><Iteration_hits>
        <Hit><Hit_def>target|Taxon_A|A</Hit_def><Hit_hsps><Hsp><Hsp_identity>100</Hsp_identity><Hsp_align-len>100</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>100</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>100</Hsp_hit-to></Hsp></Hit_hsps></Hit>
        <Hit><Hit_def>target|Taxon_B|B</Hit_def><Hit_hsps><Hsp><Hsp_identity>100</Hsp_identity><Hsp_align-len>100</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>100</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>100</Hsp_hit-to></Hsp></Hit_hsps></Hit>
        </Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"""
        match = blast_references._parse_blast_xml(xml, "target", "Taxon A")
        self.assertEqual(match.status, "ambiguous_taxon")
        self.assertFalse(match.on_target)

    def test_overlapping_hsps_do_not_sum_to_full_coverage(self):
        xml = """<BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-len>400</Iteration_query-len><Iteration_hits><Hit><Hit_def>target|Taxon_A|A</Hit_def><Hit_hsps>
        <Hsp><Hsp_identity>200</Hsp_identity><Hsp_align-len>200</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>200</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>200</Hsp_hit-to></Hsp>
        <Hsp><Hsp_identity>200</Hsp_identity><Hsp_align-len>200</Hsp_align-len><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>200</Hsp_query-to><Hsp_hit-from>201</Hsp_hit-from><Hsp_hit-to>400</Hsp_hit-to></Hsp>
        </Hit_hsps></Hit></Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"""
        match = blast_references._parse_blast_xml(xml, "target", "Taxon A")
        self.assertEqual(match.status, "insufficient_alignment")
        self.assertFalse(match.on_target)

    def test_blast_requests_every_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "db"
            Path(f"{database}.reference_count").write_text("22")
            completed = subprocess.CompletedProcess(
                [], 0, self.fixture("wrong_gene.xml"), ""
            )
            with mock.patch.object(
                blast_references.subprocess, "run", return_value=completed
            ) as command:
                blast_references.blast_against_references(
                    "A" * 100, database, "target", "Taxon A"
                )
        arguments = command.call_args.args[0]
        self.assertEqual(arguments[arguments.index("-num_alignments") + 1], "22")
        self.assertNotIn("-max_target_seqs", arguments)


class CurrentRunManifestTests(unittest.TestCase):
    def test_execution_failures_are_detectable_after_artifact_generation(self):
        self.assertTrue(runner.has_execution_failures([({}, [{"success": False}])]))
        self.assertFalse(runner.has_execution_failures([({}, [{"success": True}])]))

    def test_manifest_file_list_excludes_stale_fasta_and_parses_each_header(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory)
            current = run_dir / "sample_panel_gene.fasta"
            stale = run_dir / "sample_panel_stale.fasta"
            current.write_text(
                ">p0 sample=s gene=panel_gene product=0 length=4 kmer_count_median=7\nAAAA\n"
                ">p1 sample=s gene=panel_gene product=1 length=5 kmer_count_median=11\nCCCCC\n"
            )
            stale.write_text(">stale\nGGGG\n")
            parsed = runner.parse_fasta_products(
                "sample", run_dir, ["sample_panel_gene.fasta"]
            )
        self.assertEqual(len(parsed), 1)
        self.assertEqual([7, 11], [item["kmer_count_median"] for item in parsed[0]["products"]])
        self.assertEqual([0, 1], [item["product_index"] for item in parsed[0]["products"]])

    def test_parameter_mismatch_is_rejected(self):
        stats = {"sample": "run", "kmer_length": 31, "pcr_results": []}
        with self.assertRaisesRegex(ValueError, "kmer_length"):
            runner._validate_stats_manifest(stats, "run", 19, set(), "counting-only")

    def test_counting_manifest_may_omit_pcr_results(self):
        stats = {"sample": "run", "kmer_length": 19}
        runner._validate_stats_manifest(stats, "run", 19, set(), "counting-only")

    def test_extra_args_cannot_add_unfingerprinted_input(self):
        with tempfile.NamedTemporaryFile(suffix=".fastq") as input_file:
            with self.assertRaisesRegex(ValueError, "untracked input"):
                runner._validate_extra_args([input_file.name])

    def test_failed_execution_preserves_stdout_and_stderr(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fastq"
            input_path.write_text("@r\nACGT\n+\n!!!!\n")
            (root / "panel.yaml").write_text("name: panel\nprimers: []\n")
            completed = subprocess.CompletedProcess([], 2, "partial stdout", "useful stderr")
            with mock.patch.object(runner, "_run_with_rss", return_value=(completed, None, "unavailable")):
                run = runner.run_sharkmer(
                    root / "panel.yaml", "panel", "sample", 1, root / "runs",
                    executable=root / "sharkmer", input_path=input_path,
                )
            self.assertEqual(Path(run["logs"]["stdout"]).read_text(), "partial stdout")
            self.assertEqual(Path(run["logs"]["stderr"]).read_text(), "useful stderr")
            self.assertEqual(run["failure"]["kind"], "failed_run")

    def test_success_without_current_stats_rejects_stale_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "reads.fastq"
            input_path.write_text("@r\nACGT\n+\n!!!!\n")
            (root / "panel.yaml").write_text("name: panel\nprimers: []\n")
            run_root = root / "runs"
            run_root.mkdir()
            (run_root / "panel_sample_0k_panel_gene.fasta").write_text(">stale\nAAAA\n")
            completed = subprocess.CompletedProcess([], 0, "", "")
            with mock.patch.object(runner, "_run_with_rss", return_value=(completed, None, "unavailable")):
                run = runner.run_sharkmer(
                    root / "panel.yaml", "panel", "sample", 1, run_root,
                    executable=root / "sharkmer", input_path=input_path,
                )
            self.assertFalse(run["success"])
            self.assertEqual(run["failure"]["kind"], "invalid_manifest")

    def test_remote_run_requires_checksummed_consumed_source_plan(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel_path = root / "panel.yaml"
            panel_path.write_text("name: panel\nprimers: []\n")

            def successful_without_source_plan(command, resource_path):
                output_dir = Path(command[command.index("-o") + 1])
                sample = command[command.index("-s") + 1]
                stats = {
                    "sample": sample,
                    "kmer_length": 19,
                    "command": " ".join(command),
                    "n_reads_read": 10,
                    "n_bases_read": 100,
                }
                (output_dir / f"{sample}.stats.yaml").write_text(yaml.safe_dump(stats))
                return subprocess.CompletedProcess([], 0, "", ""), None, "unavailable"

            with mock.patch.object(runner, "CACHE_DIR", root / "cache"), mock.patch.object(
                runner, "_run_with_rss", side_effect=successful_without_source_plan
            ):
                run = runner.run_sharkmer(
                    panel_path,
                    "panel",
                    "remote-sample",
                    10,
                    root / "runs",
                    executable=root / "sharkmer",
                    benchmark_scope="counting-only",
                )
            self.assertFalse(run["success"])
            self.assertIn("input-source provenance", run["failure"]["message"])


class OutputTransactionTests(unittest.TestCase):
    def write_transaction(self, directory, sample="sample", successful=True):
        stats = {
            "sharkmer_version": "3.2.0-dev",
            "sample": sample,
            "kmer_length": 19,
            "run_id": "unique-current-run",
            "run_status": "complete",
            "output_manifest": f"{sample}.manifest.yaml",
            "pcr_results": [],
        }
        files = []
        if successful:
            fasta_name = f"{sample}_panel_gene.fasta"
            (directory / fasta_name).write_text(">record gene=panel_gene product=0\nACGT\n")
            files.append(fasta_name)
            stats["pcr_results"].append({
                "gene_name": "panel_gene", "status": "success", "n_products": 1,
                "product_lengths": [4], "output_file": fasta_name,
            })
        stats_name = f"{sample}.stats.yaml"
        (directory / stats_name).write_text(yaml.safe_dump(stats))
        files.append(stats_name)
        manifest = {
            "schema_version": 1, "producer": "sharkmer", "sample": sample,
            "run_id": stats["run_id"], "status": "complete",
            "files": [
                {"path": file_name, "sha256": runner._sha256_file(directory / file_name)}
                for file_name in files
            ],
        }
        (directory / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
        return stats, manifest

    def test_complete_transaction_excludes_unrelated_stale_products(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_transaction(root)
            stale = root / "sample_panel_stale.fasta"
            stale.write_text(">stale\nAAAA\n")
            products = runner.parse_fasta_products("sample", root)
            self.assertEqual([product["gene"] for product in products], ["panel_gene"])
            self.assertTrue(stale.exists())

    def test_interrupted_and_failed_transactions_are_not_success(self):
        for status in ("in_progress", "publishing", "failed"):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stats, manifest = self.write_transaction(root)
                manifest["status"] = status
                (root / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
                with self.assertRaisesRegex(ValueError, "not complete"):
                    runner.parse_fasta_products("sample", root)
                with self.assertRaisesRegex(ValueError, "not complete"):
                    runner.parse_fasta_products("sample", root, ["sample_panel_gene.fasta"])

    def test_matching_completed_run_identity_is_required(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats, manifest = self.write_transaction(root)
            manifest["run_id"] = "previous-run"
            (root / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
            with self.assertRaisesRegex(ValueError, "different or incomplete"):
                runner._validate_output_manifest(stats, "sample", root)

    def test_current_stats_and_products_are_checksummed(self):
        for file_name in ("sample.stats.yaml", "sample_panel_gene.fasta"):
            with self.subTest(file_name=file_name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stats, _ = self.write_transaction(root)
                with (root / file_name).open("a") as artifact:
                    artifact.write("\n")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    runner._validate_output_manifest(stats, "sample", root)

    def test_receipt_set_must_match_current_stats_exactly(self):
        for mutation in ("missing", "extra", "duplicate", "traversal"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stats, manifest = self.write_transaction(root)
                if mutation == "missing":
                    manifest["files"].pop()
                elif mutation == "extra":
                    manifest["files"].append({"path": "unrelated.txt", "sha256": "0" * 64})
                elif mutation == "duplicate":
                    manifest["files"].append(manifest["files"][0])
                else:
                    manifest["files"][0]["path"] = "../outside.fasta"
                (root / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
                with self.assertRaises(ValueError):
                    runner._validate_output_manifest(stats, "sample", root)

    def test_successful_output_path_must_belong_to_its_sample_and_gene(self):
        for file_name in ("foreign.txt", "other_panel_gene.fasta", "sample_other.fasta"):
            with self.subTest(file_name=file_name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stats, manifest = self.write_transaction(root)
                (root / "sample_panel_gene.fasta").rename(root / file_name)
                stats["pcr_results"][0]["output_file"] = file_name
                (root / "sample.stats.yaml").write_text(yaml.safe_dump(stats))
                manifest["files"][0]["path"] = file_name
                manifest["files"][1]["sha256"] = runner._sha256_file(root / "sample.stats.yaml")
                (root / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
                with self.assertRaisesRegex(ValueError, "sample and gene"):
                    runner._validate_output_manifest(stats, "sample", root)

    def test_symlinks_are_not_owned_output_files(self):
        for file_name in ("sample.manifest.yaml", "sample.stats.yaml", "sample_panel_gene.fasta"):
            with self.subTest(file_name=file_name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                stats, _ = self.write_transaction(root)
                original = root / file_name
                target = root / "unrelated-preserved"
                original.rename(target)
                original.symlink_to(target)
                with self.assertRaises(ValueError):
                    runner._validate_output_manifest(stats, "sample", root)
                self.assertTrue(target.exists())

    def test_new_versions_require_a_transaction_and_old_snapshots_remain_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats = {"sharkmer_version": "3.2.0-dev", "sample": "sample"}
            with self.assertRaisesRegex(ValueError, "manifest is missing"):
                runner._validate_output_manifest(stats, "sample", root)
            stats["sharkmer_version"] = "3.1.0"
            self.assertIsNone(runner._validate_output_manifest(stats, "sample", root))

    def test_counting_only_transaction_has_only_stats_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stats, manifest = self.write_transaction(root, successful=False)
            self.assertEqual(runner._validate_output_manifest(stats, "sample", root), manifest)
            self.assertEqual(runner.parse_fasta_products("sample", root), [])

    def test_bootstrap_does_not_prefer_old_root_fasta_over_interrupted_invocation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = "panel_accession_1k"
            (root / f"{prefix}_panel_gene.fasta").write_text(">stale\nAAAA\n")
            invocation = root / f"{prefix}_new"
            invocation.mkdir()
            stats, manifest = self.write_transaction(invocation, sample=prefix)
            manifest["status"] = "publishing"
            (invocation / stats["output_manifest"]).write_text(yaml.safe_dump(manifest))
            panel = {
                "name": "panel", "primers": [],
                "validation": {"samples": [{"accession": "accession", "max_reads": [1000]}]},
            }
            with self.assertRaisesRegex(ValueError, "not complete"):
                bootstrap_from_runs.collect_amplicons_from_runs(panel, root, "panel")

    def test_non_mapping_stats_are_rejected_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.stats.yaml"
            path.write_text("- not-a-stats-mapping\n")
            with self.assertRaisesRegex(ValueError, "must be a mapping"):
                runner._parse_stats_yaml(path)

    def test_transaction_change_during_product_reading_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, manifest = self.write_transaction(root)
            changed_manifest = {**manifest, "run_id": "a-new-run"}
            with mock.patch.object(
                runner, "_validate_output_manifest", side_effect=[manifest, changed_manifest]
            ):
                with self.assertRaisesRegex(ValueError, "changed while products"):
                    runner.parse_fasta_products("sample", root)


class ExecutableProvenanceTests(unittest.TestCase):
    def test_same_second_runs_get_distinct_report_names(self):
        fixed_time = mock.Mock()
        fixed_time.strftime.return_value = "20260908_162116_000000"
        panel = {"name": "panel", "panel_version": "1"}
        with mock.patch.object(runner, "datetime") as datetime_class:
            datetime_class.now.return_value = fixed_time
            first_id = runner.unique_run_id()
            second_id = runner.unique_run_id()
        first_name = results.result_filename(panel, "3.2.0", first_id)
        second_name = results.result_filename(panel, "3.2.0", second_id)
        self.assertNotEqual(first_id, second_id)
        self.assertNotEqual(first_name, second_name)
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            results.write_result({"cache_mode": "warm"}, output_dir / first_name)
            results.write_result({"cache_mode": "cold"}, output_dir / second_name)
            written = sorted(output_dir.glob("*.yaml"))
            self.assertEqual(len(written), 2)
            self.assertEqual(
                {results.load_result(path)["cache_mode"] for path in written},
                {"warm", "cold"},
            )

    def test_explicit_executable_is_fingerprinted_without_build(self):
        with tempfile.NamedTemporaryFile() as executable:
            executable.write(b"binary")
            executable.flush()
            with mock.patch.object(runner, "get_git_commit", return_value="revision"), mock.patch.object(
                runner, "_git_dirty", return_value=True
            ), mock.patch.object(runner, "_source_tree_sha256", return_value="tree"):
                provenance = runner.build_sharkmer(Path(executable.name))
        self.assertTrue(provenance["selected_explicitly"])
        self.assertFalse(provenance["built_for_run"])
        self.assertEqual(provenance["workspace_source_observation"]["relationship_to_explicit_binary"], "unknown")

    def test_default_executable_is_built_even_when_it_exists(self):
        with tempfile.NamedTemporaryFile() as executable:
            executable.write(b"binary")
            executable.flush()
            artifact = {
                "reason": "compiler-artifact",
                "target": {"name": "sharkmer"},
                "executable": executable.name,
                "features": ["ahashmap", "default"],
            }
            completed = subprocess.CompletedProcess([], 0, json.dumps(artifact), "")
            with mock.patch.object(runner, "SHARKMER_BIN", Path(executable.name)), mock.patch.object(
                runner.subprocess, "run", return_value=completed
            ) as command, mock.patch.object(runner, "get_git_commit", return_value="revision"), mock.patch.object(
                runner, "_git_dirty", return_value=False
            ), mock.patch.object(runner, "_source_tree_sha256", return_value="tree"):
                provenance = runner.build_sharkmer()
        self.assertEqual(command.call_args.args[0][:3], ["cargo", "build", "--release"])
        self.assertTrue(provenance["built_for_run"])
        self.assertEqual(provenance["hash_backend"], "ahashmap")


class ResultStatusTests(unittest.TestCase):
    def test_failed_and_filtered_genes_are_distinct(self):
        panel = {"name": "panel", "primers": [{"gene": "a"}, {"gene": "b"}]}
        sample = {"accession": "X", "taxon": "Taxon A"}
        failed = {"max_reads": 10, "success": False, "genes": []}
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml") as panel_file:
            panel_file.write("name: panel\n")
            panel_file.flush()
            built = results.build_result(
                Path(panel_file.name), panel, [(sample, [failed])], "test", evaluated_genes=["a"]
            )
        statuses = {gene["gene"]: gene["evaluation_status"] for gene in built["samples"][0]["depths"][0]["genes"]}
        self.assertEqual(statuses, {"a": "failed_run", "b": "not_evaluated"})


if __name__ == "__main__":
    unittest.main()
