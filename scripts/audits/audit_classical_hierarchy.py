"""Retain the independent monotone/Taylor audit and branch-regression results."""
from pathlib import Path
import argparse
import hashlib
import io
import json
import os
import sys
import unittest


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
    from audit_qpe_benchmark import taylor_audit
    mathematical=taylor_audit(args.root)
    os.environ['QBP_AUDIT_ROOT']=str(args.root)
    suite=unittest.defaultTestLoader.discover(str(HERE),pattern='test_classical_hierarchy.py')
    text=io.StringIO();tests=unittest.TextTestRunner(stream=text,verbosity=2).run(suite)
    assert tests.wasSuccessful()
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    report={'schema':'qbp-monotone-taylor-hierarchy-independent-audit-v1','status':'passed',
            'mathematical_and_result_audit':mathematical,'test_groups':tests.testsRun,
            'retained_test_output':text.getvalue(),
            'source_sha256':{n:sha(args.root/n)for n in
                ('taylor_order_benchmark.py','taylor_pricing.py','continuous_quadrature.py')},
            'input_sha256':{n:sha(args.root/n)for n in
                ('results/taylor_order_benchmark.json','results/taylor_pricing.json',
                 'results/continuous_quadrature.json','results/generalisation_protocol.json')},
            'audit_source_sha256':{p.name:sha(p)for p in
                (Path(__file__),HERE/'audit_qpe_benchmark.py',HERE/'test_classical_hierarchy.py')}}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','test_groups':tests.testsRun,
                      'family_policy_runs':mathematical['family_policy_runs'],
                      'cold_standalone_runs':mathematical['cold_standalone_runs'],
                      'counters':mathematical['counters']},indent=2))
