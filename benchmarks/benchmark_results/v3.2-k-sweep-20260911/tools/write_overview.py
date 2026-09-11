import hashlib
import json
from pathlib import Path


workspace = Path(__file__).resolve().parent
analysis_path = workspace / 'discovery-analysis/analysis.json'
summary_path = workspace / 'discovery-analysis/summary.json'
analysis = json.loads(analysis_path.read_text())
summary = json.loads(summary_path.read_text())
options = analysis['confirmation_options_ranked']
assert len(options) == 3 and all(not option['eligible'] for option in options)
overview = {
    'status': 'Discovery complete; no eligible longer k; conditional confirmation not run',
    'release': 'Held under #153; no default or production changes',
    'source_files': {
        str(path.relative_to(workspace)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [analysis_path, summary_path, workspace / 'discovery-protocol.json']
    },
    'counts': {'samples': 3, 'k_values': [19, 23, 27, 31], 'versions': 2, 'pairs_per_cell': 3, 'timed_invocations': 72},
    'per_k': {
        key: {
            field: value[field]
            for field in [
                'high_copy_exact_products_total',
                'baseline_wall_sum_cell_medians_s',
                'candidate_wall_sum_cell_medians_s',
                'candidate_rss_median_cell_medians_bytes',
                'same_k_count_parity_every_pair',
                'same_k_sequence_parity',
                'repeat_sequence_and_classification_stable',
            ]
        }
        for key, value in summary['per_k'].items()
    },
    'options': [
        {
            'k': option['kmer_length'],
            'eligible': option['eligible'],
            'exact_missing_products_restored': len(option['exact_restored_known_losses_every_replicate']),
            'released_k19_high_copy_retained': len(option['k19_released_high_copy_exact_retained_every_replicate']),
            'confirmed_k19_candidate_sequences_lost': len(option['confirmed_k19_candidate_sequences_lost']),
        }
        for option in options
    ],
    'scope_note': 'Only the three insect samples; products are not multiplied by timing repetitions. Exact restoration differs from new gene recovery or reference classification.',
}
with (workspace / 'overview.json').open('x') as output:
    json.dump(overview, output, indent=2)
    output.write('\n')
print(workspace / 'overview.json')
