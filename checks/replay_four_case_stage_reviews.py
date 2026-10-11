"""Replay the retained independent reviews without changing executed sources.

The original audit scripts expect a two-branch staging layout. This wrapper
supplies their paths explicitly and copies only the small classical protocol
needed for the cross-method comparison. New audit outputs remain separate
from the frozen evidence. No production helper is imported by either review.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True,
                        help='New or empty directory for replayed audit outputs')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise SystemExit('Output directory must be new or empty')
    audits = [('classical', 'audit_final_classical_arithmetic.py',
               'final_classical_arithmetic_audit.json'),
              ('quantum', 'audit_final_quantum_algebra.py',
               'final_quantum_algebra_audit.json')]
    for _, source, result in audits:
        path = ROOT / 'reviews/four_case_stage_audits' / source
        retained = json.loads((ROOT / 'results' / result).read_text())
        if sha(path) != retained['auditor_sha256']:
            raise SystemExit('Original review source failed verification: ' + source)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'licence_review').mkdir()
    protocol_dir = output / 'classical/results'
    protocol_dir.mkdir(parents=True)
    protocol = ROOT / 'results/four_case_classical_protocol.json'
    (protocol_dir / protocol.name).write_bytes(protocol.read_bytes())
    manifest = {'schema': 'qbp-portable-stage-review-replay-v1',
                'wrapper_sha256': sha(__file__),
                'classical_protocol_sha256': sha(protocol), 'reviews': [],
                'scope': 'Original audit source bytes are unmodified. Only path globals are supplied for the installed layout. Fresh timing metadata differs from retained results.'}
    for kind, source, result in audits:
        path = ROOT / 'reviews/four_case_stage_audits' / source
        spec = importlib.util.spec_from_file_location('qbp_review_' + kind, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.STAGE = output
        module.ROOT = ROOT
        module.RESULTS = ROOT / 'results'
        module.REPO = ROOT
        module.main()
        replay = output / 'licence_review' / result
        manifest['reviews'].append({'source': source, 'source_sha256': sha(path),
                                    'replayed_result': str(replay.relative_to(output)),
                                    'replayed_result_sha256': sha(replay),
                                    'status': json.loads(replay.read_text())['status']})
    (output / 'replay_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
