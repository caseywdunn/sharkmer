import copy
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


workspace = Path(__file__).resolve().parent
previous = Path('/tmp/sharkmer-high-copy-20260911')
original_path = previous / 'final-protocol.json'
protocol = json.loads(original_path.read_text())
protocol['schema_version'] = 1
protocol['registered_at'] = datetime.now(timezone.utc).isoformat()
protocol['purpose'] = 'Preregistered longer-k calibration under #153; no production/default changes or release'
protocol['allowed_k'] = [19, 23, 27, 31]
del protocol['settings']['k']
protocol['samples'] = [sample for sample in protocol['samples'] if sample['panel'] == 'insecta']
protocol['target_metadata'] = {'insecta': protocol['target_metadata']['insecta']}
protocol['ordering'] = {
    'outer_order': ['pair_index', 'protocol_sample_order'],
    'k_rotation': 'left_by_(pair_index_minus_1_plus_zero_based_sample_index)_modulo_k_count',
    'version_order': 'baseline_first_when_(pair_index_plus_zero_based_sample_index_plus_allowed_k_index)_modulo_2_equals_1',
    'same_k_version_pair_contiguous': True,
}
protocol['fixed_reference'] = {
    'execution': str(previous / 'final-execution'),
    'analysis': str(previous / 'final-analysis'),
    'protocol_sha256': hashlib.sha256(original_path.read_bytes()).hexdigest(),
    'diagnostic_audit_sha256': hashlib.sha256((previous / 'final-analysis/diagnostic-audit.json').read_bytes()).hexdigest(),
    'baseline_k': 19,
    'expected_lost_high_copy_products': 7,
    'expected_candidate_high_copy_products_all_ten_samples': 64,
    'exact_its2_sha256': '23163ca76f97f5d1685f0fcf4bda6ffc798bdf58c63953416833e78a480af68c',
}
protocol['confirmation_selection'] = {
    'status': 'preregistered_before_discovery',
    'eligible': 'k > 19; all discovery invocations and classifications valid and same-version sequence/classification stable; restores at least one of the seven exact missing high-copy products in every replicate; preserves every k19 candidate confirmed_product sequence in the affected samples in every replicate',
    'rank_descending': ['exact_reference_ITS2_restored', 'number_of_seven_exact_missing_products_restored', 'number_of_k19_released_high_copy_sequences_retained'],
    'tie_break': 'lowest k',
    'select_at_most': 1,
    'confirmation': 'Three alternating pairs at selected k for the seven remaining non-insect samples from the frozen ten-sample protocol; same inputs, panels, and settings',
    'if_none_eligible': 'Do not expand timings; report the result and uncertainties for user review',
    'release_gate': 'Selection is exploratory calibration, not a change to default k or release approval; report every other high-copy loss and resource increase',
}
protocol['comparison_contract'] = {
    'primary': 'Exact gene-plus-sequence membership, not product index or raw gene/product totals',
    'same_k': 'Both versions at each k, including aggregate read/base/kmer occurrence parity',
    'cross_k': 'Both versions against fixed k19 released and candidate memberships; read/base/prefix equality required, kmer occurrence totals may differ',
    'changed_endpoints': 'Keep exact sequence changes visible; only separately report boundary equivalence if demonstrated from frozen FASTA/panel evidence; never count as exact rescue',
    'mechanism': 'Retain seed/threshold/SCC/collision/node/DFS/path diagnostics; k changes seed context and effective support as well as topology, so cannot attribute every change solely to repeat resolution',
    'resources': 'Whole unchanged-panel wall/RSS; same-k version comparisons and within-version cross-k comparisons; three same-host repeats are descriptive',
}
protocol['limitations'].extend([
    'All inputs are previously inspected calibration data, not independent held-out truth.',
    'Larger k changes coverage, observed seed context and accepted kmer workload; no coverage, primer or search-budget retuning.',
    'Exact recovery is not independent read-supported repeat-copy truth; no release/default change is authorized.',
])
assert len(protocol['samples']) == 3
assert sum(sample['pairs'] for sample in protocol['samples']) * len(protocol['allowed_k']) * 2 == 72
with (workspace / 'discovery-protocol.json').open('x') as output:
    json.dump(protocol, output, indent=2)
    output.write('\n')
confirmation = copy.deepcopy(json.loads(original_path.read_text()))
confirmation['samples'] = [sample for sample in confirmation['samples'] if sample['panel'] != 'insecta']
confirmation['target_metadata'].pop('insecta')
assert len(confirmation['samples']) == 7
with (workspace / 'confirmation-template.json').open('x') as output:
    json.dump(confirmation, output, indent=2)
    output.write('\n')
print(workspace / 'discovery-protocol.json')
