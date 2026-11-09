import argparse
import collections
import json
from pathlib import Path

import yaml


parser = argparse.ArgumentParser()
parser.add_argument("--execution", type=Path, required=True)
parser.add_argument("--analysis", type=Path, required=True)
parser.add_argument("--pre-fix-results", type=Path)
parser.add_argument("--output", type=Path)
arguments = parser.parse_args()
analysis = json.loads((arguments.analysis / "analysis.json").read_text())
protocol = json.loads((arguments.execution / "provenance.json").read_text())["protocol"]
scopes = {
    panel: {entry["gene"]: entry["scope"] for entry in entries}
    for panel, entries in protocol["target_metadata"].items()
}
selected = collections.defaultdict(dict)
for path in sorted((arguments.execution / "results").glob("*.json")):
    result = json.loads(path.read_text())
    invocation = result["signature"]["invocation"]
    assert result["status"] == "complete"
    if invocation["pair_index"] == 1:
        selected[invocation["cell"]][invocation["version"]] = result
totals = {role: collections.Counter() for role in ("baseline", "candidate")}
audits = []
pre_fix_checks = []
for cell in analysis["cells"]:
    versions = selected[cell["cell"]]
    panel = versions["candidate"]["signature"]["invocation"]["panel"]
    for role, result in versions.items():
        for gene in result["genes"]:
            totals[role][scopes[panel][gene["gene"]]] += len(gene["products"])
    lost = cell["pairs"][0]["lost_products"]
    high_copy_lost = [product for product in lost if product["scope"] == "high_copy_candidate"]
    stats = yaml.safe_load(Path(versions["candidate"]["stats_path"]).read_text())
    stats_by_gene = {gene["gene_name"]: gene for gene in stats["pcr_results"]}
    for gene_name in sorted({product["gene"] for product in high_copy_lost}):
        gene_stats = stats_by_gene[f"{panel}_{gene_name}"]
        audits.append({
            "cell": cell["cell"], "gene": gene_name,
            "lost_products": [product for product in high_copy_lost if product["gene"] == gene_name],
            "candidate_gene_status": gene_stats["status"],
            "candidate_product_count": gene_stats["n_products"],
            "candidate_failure_reason": gene_stats.get("failure_reason"),
            "threshold_diagnostics": gene_stats["threshold_diagnostics"],
        })
if arguments.pre_fix_results:
    for path in sorted(arguments.pre_fix_results.glob("*.json")):
        previous = json.loads(path.read_text())
        invocation = previous["signature"]["invocation"]
        if invocation["version"] != "candidate" or invocation["pair_index"] != 1 or invocation["cell"] not in selected:
            continue
        current = selected[invocation["cell"]]["candidate"]
        assert previous["status"] == "complete"
        for field in ("panel_sha256", "input_sha256", "input_subset"):
            assert previous["signature"][field] == current["signature"][field]
        memberships = [
            {(gene["gene"], product["sha256"]) for gene in result["genes"] for product in gene["products"]}
            for result in (previous, current)
        ]
        assert memberships[0] == memberships[1]
        pre_fix_checks.append({"cell": invocation["cell"], "products": len(memberships[1]), "unchanged": True})
    assert len(pre_fix_checks) == len(selected)
output_path = arguments.output or (arguments.analysis / "diagnostic-audit.json")
with output_path.open("x") as output:
    json.dump({
        "representation": "First pair per cell; biological loci/products are not multiplied by repeated timings.",
        "product_totals_by_scope": totals,
        "high_copy_lost_gene_audits": audits,
        "pre_fix_sequence_comparison": pre_fix_checks,
        "interpretation": "Markers and exhausted search do not establish sequence absence or exact repeat-copy truth; alternative product order is not an abundance estimate.",
    }, output, indent=2)
    output.write("\n")
print(json.dumps(totals))
