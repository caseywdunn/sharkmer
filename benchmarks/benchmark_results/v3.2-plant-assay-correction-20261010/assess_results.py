#!/usr/bin/env python3
"""Verify all eight plant results and join independent exact/partial evidence.

No new assembly, read scan or BLAST query. Previously assessed queries join by
sample, target and full sequence hash. Full literal reference occurrences are
checked against the pinned catalog. Outputs require a fresh directory.
"""
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/results-v1/plant'
PRIOR = ROOT / 'benchmarks/benchmark_results/v3.2-target-coverage-20261009'
TARGET = 'angiospermae_trnV-atpE'
PLANTS = {'SRR25378184', 'ERR14009273'}
spec = importlib.util.spec_from_file_location('frozen_plant_runner', HERE / 'run_comparison.py')
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)
from sharkmer_validate import reference_provenance


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seqsha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def dump(path, data):
    with path.open('x') as stream:
        stream.write(json.dumps(data, indent=2, sort_keys=True) + '\n')


def identity(row, accession=None):
    return (accession or row['accession'], row['gene'], row['sha256'])


def inventory(records):
    # Keep actual sequence and multiplicity: matching only a count/hash set is
    # insufficient for the requested historical reproduction comparison.
    return Counter((r['gene'], r['sequence'], r['sha256'], r['assembly_status']) for r in records)


def placements(sequence, sources, assay_intervals):
    result = []
    for accession, record in sources.items():
        source = record['sequence']
        n = len(source)
        if len(sequence) > n:
            continue
        for strand, oriented in [('+', source), ('-', reference_provenance.reverse_complement(source))]:
            search = oriented + oriented[:len(sequence) - 1] if record['topology'] == 'circular' else oriented
            offset = 0
            while True:
                index = search.find(sequence, offset)
                if index < 0 or index >= n:
                    break
                stop = index + len(sequence)
                if strand == '+':
                    start, end = index, stop if stop <= n else stop - n
                else:
                    start, end = (n - stop) % n, n - index
                annotation = assay_intervals.get(accession)
                assigned = bool(annotation and strand == '+' and [start, end] == [annotation['start'], annotation['end']])
                result.append({'accession': accession, 'source_taxon': record['organism'], 'start': start, 'end': end,
                               'strand': strand, 'wraps_origin': stop > n,
                               'annotation': 'exact_reviewed_corrected_assay_interval' if assigned else 'not_assigned_here',
                               'source_sequence_sha256': record['sequence_sha256']})
                offset = index + 1
    return result


def tsv(path, fields, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    output = HERE / 'assessment-v1'
    if output.exists():
        raise SystemExit('assessment-v1 must not already exist')
    summary = json.loads((RUNS / 'summary.json').read_text())
    plan = json.loads((RUNS / 'plan.json').read_text())
    mapping = json.loads((HERE / 'proposal-v2/mapping.json').read_text())
    previous = json.loads((PRIOR / 'reference/assessment-v2/report.json').read_text())
    saved_inventory = json.loads((PRIOR / 'matrix/final-v1/sequence_inventory.json').read_text())
    historical = json.loads((PRIOR / 'matrix/historical_inventory.json').read_text())
    matrix = json.loads((PRIOR / 'matrix/final-v1/target_matrix.json').read_text())
    assert summary['status'] == 'pass' and len(summary['jobs']) == 8
    assert {(j['scope'], j['sample_id'], j['variant']) for j in summary['jobs']} == {
        (scope, sample, variant) for scope, samples in [('real', PLANTS), ('synthetic', {'NC_008326.1', 'NC_056221.1'})]
        for sample in samples for variant in ('before', 'after')}
    fingerprints = {row['path']: row for row in [plan['binary'], *plan['panels'].values(), *plan['bindings'].values()]}
    input_hashes = {}
    run_rows = []
    for job in summary['jobs']:
        assert job['status'] == 'pass' and job['execution']['status'] == 'pass'
        receipt = Path(job['directory']) / 'receipt.json'
        assert json.loads(receipt.read_text()) == job
        input_hashes[str(receipt.relative_to(ROOT))] = sha(receipt)
        for fingerprint in [job['input_file'], *job['execution']['files'].values(), *job['outputs']['files'].values()]:
            fingerprints[fingerprint['path']] = fingerprint
        counts = Counter(r['assembly_status'] for r in job['outputs']['records'])
        target = [r for r in job['outputs']['records'] if r['gene'] == TARGET]
        target_counts = Counter(r['assembly_status'] for r in target)
        run_rows.append({'scope': job['scope'], 'sample_id': job['sample_id'], 'variant': job['variant'],
                         'total': sum(counts.values()), 'standard': counts['standard'], 'uncertain': counts['uncertain'],
                         'target_total': len(target), 'target_standard': target_counts['standard'], 'target_uncertain': target_counts['uncertain'],
                         'wall_seconds': job['execution']['metrics']['wall_seconds'],
                         'peak_rss_bytes': job['execution']['metrics']['peak_rss_bytes']})
        genes = {r['gene_name']: 'assay_correction_evaluation' for r in job['outputs']['pcr_results']}
        assert len(genes) == 10
        invocation = {'sample_prefix': job['sample_prefix'], 'version': 'candidate'}
        recaptured = producer.capture.capture(Path(job['outputs_directory']), invocation, job, job['command'], genes)
        assert recaptured == job['outputs']
    for fingerprint in fingerprints.values():
        producer.capture.verify(fingerprint)
    assert producer.assess(summary['jobs'], genes, mapping['references']) == summary['comparisons']
    assert all(c['aggregate_counts_unchanged'] and c['other_nine_targets_unchanged'] for c in summary['comparisons'])
    assert all(c['exact_corrected_source_emitted_after'] and not c['exact_corrected_source_emitted_before']
               for c in summary['comparisons'] if c['scope'] == 'synthetic')

    before_matches = []
    for sample in sorted(PLANTS):
        before = next(j for j in summary['jobs'] if j['scope'] == 'real' and j['sample_id'] == sample and j['variant'] == 'before')
        saved = [r for r in saved_inventory if r['accession'] == sample]
        assert inventory(before['outputs']['records']) == inventory(saved)
        before_matches.append({'sample_id': sample, 'records': len(saved),
                               'full_sequence_hash_class_multiplicity_matches_saved_matrix': True})
    catalog, catalog_receipt = reference_provenance.load_catalog()
    assert catalog_receipt['sha256'] == mapping['catalog']['sha256']
    expected_intervals = {r['accession']: r['corrected_reference_interval'] for r in mapping['records']}
    previous_index = {identity(r): r for r in previous['candidates']}
    rows, cache = [], {}
    for job in summary['jobs']:
        for record in job['outputs']['records']:
            if record['gene'] != TARGET:
                continue
            sequence = record['sequence']
            assert seqsha(sequence) == record['sha256'] and len(sequence) == record['length']
            if sequence not in cache:
                cache[sequence] = placements(sequence, catalog['records'], expected_intervals)
            exact = cache[sequence]
            row = {'scope': job['scope'], 'sample_id': job['sample_id'], 'variant': job['variant'],
                   'gene': TARGET, 'sha256': record['sha256'], 'length': record['length'],
                   'assembly_status': record['assembly_status'], 'exact_full_reference_placements': exact,
                   'sample_full_product_phase': 'not_established'}
            if exact:
                assert any(p['annotation'] == 'exact_reviewed_corrected_assay_interval' for p in exact)
                row.update(reference_status='full_exact_annotated_corrected_assay',
                           partial_reference_assessment='superseded_by_exact_full_sequence_evidence',
                           reference_taxon_relation=('same_species' if job['sample_id'] == 'SRR25378184' else
                                                      'same_genus_related_species_proxy' if job['sample_id'] == 'ERR14009273' else 'synthetic_source'))
            else:
                assert job['scope'] == 'real' and job['variant'] == 'before'
                old = previous_index[identity(record, job['sample_id'])]
                assert old['length'] == len(sequence)
                assert old['strong_public_single_hsp_count'] == old['expected_target_hsp_count'] == old['same_taxon_hsp_count'] == 0
                assert not old['trnV_mapped_assay_overlap_hsps']
                row.update(reference_status='non_target_or_ambiguous_partial_hits_no_full_support',
                           partial_reference_assessment='reused_hash_bound_saved_query', previous_query=old,
                           prior_evidence_path=str((PRIOR / 'reference/assessment-v2/report.json').relative_to(ROOT)))
            rows.append(row)
    assert len(rows) == 45 and sum(r['partial_reference_assessment'] == 'reused_hash_bound_saved_query' for r in rows) == 41

    # Explicitly derived combination: eight saved samples and two fresh plants.
    # This is not a fresh ten-sample run using the final binary.
    combined = []
    for record in saved_inventory:
        assert seqsha(record['sequence']) == record['sha256']
        if record['accession'] not in PLANTS:
            combined.append({'accession': record['accession'], 'gene': record['gene'], 'sha256': record['sha256'],
                             'sequence': record['sequence'], 'length': record['length'], 'assembly_status': record['assembly_status'],
                             'evidence_age': 'saved_previous_eight_sample_inventory', 'source': record['source']})
    for job in summary['jobs']:
        if job['scope'] == 'real' and job['variant'] == 'after':
            for record in job['outputs']['records']:
                combined.append({'accession': job['sample_id'], 'gene': record['gene'], 'sha256': record['sha256'],
                                 'sequence': record['sequence'], 'length': record['length'], 'assembly_status': record['assembly_status'],
                                 'evidence_age': 'fresh_final_binary_corrected_plant_panel', 'source': str(Path(job['directory']) / 'receipt.json')})
    new_keys = {identity(r) for r in combined}
    old_keys = {identity(r) for r in historical}
    for r in historical:
        assert seqsha(r['sequence']) == r['sha256']
    assert len(new_keys) == len(combined) == 287 and len(old_keys) == 86
    represented = {(r['accession'], r['gene']) for r in combined}
    all_targets = {(r['accession'], r['gene']) for r in matrix}
    assert represented <= all_targets
    class_counts = Counter(r['assembly_status'] for r in combined)
    derived = {'status': 'derived_mixed_age_inventory_not_fresh_ten_sample_benchmark',
               'total': len(combined), 'standard': class_counts['standard'], 'uncertain': class_counts['uncertain'],
               'release_identities_retained': len(old_keys & new_keys), 'release_identities_lost': len(old_keys - new_keys),
               'new_vs_release': len(new_keys - old_keys), 'represented_targets': len(represented), 'target_rows': len(all_targets),
               'source_composition': 'Eight samples from the saved final-v1 matrix; two newly run plants use final patched-dependency binary and plant panel2.0.0.',
               'assembly_policy_or_performance_comparison': 'not_established_by_this_derived_view',
               'target_identity_compatibility': 'trnV-atpE name is retained but its assay boundaries changed from panel1.x to2.0.0.'}
    assert [derived[k] for k in ('standard', 'uncertain', 'release_identities_retained', 'release_identities_lost', 'new_vs_release', 'represented_targets', 'target_rows')] == [66, 221, 86, 0, 201, 75, 118]

    paths = [Path(__file__), HERE / 'run_comparison.py', HERE / 'proposal-v2/mapping.json', RUNS / 'summary.json', RUNS / 'plan.json',
             PRIOR / 'reference/assessment-v2/report.json', PRIOR / 'matrix/final-v1/sequence_inventory.json',
             PRIOR / 'matrix/historical_inventory.json', PRIOR / 'matrix/final-v1/target_matrix.json', ROOT / 'panels/reference_sources.json.gz']
    input_hashes.update({str(p.relative_to(ROOT)): sha(p) for p in paths})
    report = {'status': 'pass', 'binary': plan['binary'], 'source_hashes': input_hashes,
              'verified_run_file_fingerprints': list(fingerprints.values()), 'run_rows': run_rows,
              'comparisons': summary['comparisons'], 'before_runs_reproduce_previous_matrix': before_matches,
              'reference_assessment_complete': True, 'reference_candidates': rows,
              'reference_method': 'All literal full-sequence occurrences in both orientations of all121 pinned source records, including circular origins. For41 old nonexact hypotheses, reuse prior full-catalog BLAST assessment by sample+gene+fullhash; no new alignment or read scan.',
              'catalog_receipt': catalog_receipt, 'derived_combined_inventory': derived,
              'interpretation': [
                  'Both synthetic sources reconstruct exactly as standard products only after the versioned primer correction.',
                  'Liriodendron replaces41 unsupported old uncertain strings with one909bp standard string exact in two independently annotated same-species plastomes.',
                  'Acer gains one929bp standard string exact in the related-species Acer opalus annotated plastome.',
                  'All nine other targets preserve sequence, assembly class, multiplicity, outcomes and search diagnostics; all five aggregate read/k-mer counts remain equal.',
                  'The change corrects an assay definition; it is not a retention-policy filter applied to the old41 hypotheses.',
                  'Exact reference agreement and standard assembly status do not independently establish the complete haplotype in either real sample.'
              ],
              'limitations': ['The catalog is bounded, not an exhaustive public search.', 'Acer opalus is a related-species proxy.',
                             'Historical Liriodendron input is656431 R1 plus343569 R2 records in a concatenated file, not paired-fragment analysis.',
                             'Per-run times are functionality-check observations, not collaborator-approved performance acceptance.',
                             'No new release is authorized or performed.']}
    # Replace pending runner annotations only in this new downstream assessment.
    for comparison in report['comparisons']:
        comparison['partial_reference_assessment'] = 'completed_in_this_downstream_assessment'
    output.mkdir()
    dump(output / 'report.json', report)
    dump(output / 'derived_combined_inventory.json', combined)
    tsv(output / 'run_summary.tsv', ['scope', 'sample_id', 'variant', 'total', 'standard', 'uncertain', 'target_total', 'target_standard', 'target_uncertain', 'wall_seconds', 'peak_rss_bytes'], run_rows)
    controls = [{'scope': c['scope'], 'sample_id': c['sample_id'], 'gene': t['gene'],
                 'assay_changed': t['assay_changed'], 'before_records': sum(r['count'] for r in t['before']),
                 'after_records': sum(r['count'] for r in t['after']),
                 'literal_outputs_unchanged': t['literal_identity_status_multiplicity_unchanged'],
                 'diagnostics_unchanged': t['outcome_and_search_diagnostics_unchanged']}
                for c in report['comparisons'] for t in c['targets']]
    tsv(output / 'target_comparison.tsv', list(controls[0]), controls)
    compact_refs = [{**{k: r[k] for k in ('scope', 'sample_id', 'variant', 'sha256', 'length', 'assembly_status', 'reference_status', 'partial_reference_assessment')},
                     'exact_sources': ';'.join(f"{p['accession']}:{p['start']}-{p['end']}:{p['strand']}" for p in r['exact_full_reference_placements'])}
                    for r in rows]
    tsv(output / 'target_reference_assessment.tsv', list(compact_refs[0]), compact_refs)
    dump(output / 'checksums.json', {p.name: sha(p) for p in sorted(output.iterdir())})
    print(json.dumps({'status': 'pass', 'output': str(output), 'reference_records': len(rows), 'derived': derived}))


if __name__ == '__main__':
    main()
