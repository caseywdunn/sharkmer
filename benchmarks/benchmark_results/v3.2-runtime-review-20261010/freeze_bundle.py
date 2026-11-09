#!/usr/bin/env python3
"""Seal the completed runtime evidence bundle without changing earlier artifacts."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fingerprint(path):
    data = path.read_bytes()
    return {"sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}


def main():
    required = ["saved_analysis.json", "saved_identity_deltas.json", "independent_saved_review.json",
                "independent_protocol_review.json", "results-v1/result.json", "analysis-v1/analysis.json",
                "independent_human_results_review.json"]
    for relative in required:
        expected_status = "approved_for_bounded_diagnostic" if relative == "independent_protocol_review.json" else "pass"
        if json.loads((HERE / relative).read_text())["status"] != expected_status:
            raise ValueError("Missing passing result/review: " + relative)
    protocol = json.loads((HERE / "protocol-v1.json").read_text())
    for binding in protocol["sources"]:
        expected = {key: binding[key] for key in ("sha256", "size_bytes")}
        if fingerprint(Path(binding["path"])) != expected:
            raise ValueError("Frozen execution source changed: " + binding["path"])
    files = {str(path.relative_to(HERE)): fingerprint(path) for path in sorted(HERE.rglob("*"))
             if path.is_file() and "__pycache__" not in path.parts and path.name != "bundle_manifest.json"}
    result = {"status": "pass", "scope": "Completed saved-run analysis, bounded Human diagnostic, full outputs/traces, independent reviews and findings",
              "required_passing_results": required, "execution_sources_reverified": len(protocol["sources"]),
              "files": files, "file_count": len(files), "bytes": sum(row["size_bytes"] for row in files.values()),
              "exclusions": ["This self-referential manifest", "Python __pycache__ files"],
              "no_release_acceptance_claim": True}
    with (HERE / "bundle_manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "pass", "files": len(files), "bytes": result["bytes"],
                      "manifest": fingerprint(HERE / "bundle_manifest.json")}))


if __name__ == "__main__":
    main()
