"""Retain the actual coherent and finite-arithmetic QPE audit outcome."""
from pathlib import Path
import hashlib
import io
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from qpe_estimator import plan,compiled_plan
from fractions import Fraction as F


def main():
    suite=unittest.defaultTestLoader.discover(str(ROOT/'checks'),pattern='test_qpe_estimator.py')
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    record={'schema':'qbp-canonical-qpe-estimator-audit-v1',
            'status':'passed' if result.wasSuccessful() else 'failed',
            'test_groups':result.testsRun,'failures':len(result.failures),
            'errors':len(result.errors),'retained_test_output':stream.getvalue(),
            'coherent_scope':{'ideal_inverse_DFT_matrix_widths':[1,2,3,4,5,6],
                'ideal_inverse_DFT_entries_compared':5460,
                'ideal_complete_QPE_distributions':15,
                'complete_QPE_probability_cases':['0','1/4','1/2','3/4','1'],
                'complete_QPE_phase_widths':[2,3,4],
                'literal_Clifford_T_inverse_DFT_widths':[1,2,3],
                'literal_Clifford_T_inverse_DFT_columns':14,
                'literal_Clifford_T_complete_QPE_distributions':13,
                'literal_Clifford_T_complete_QPE_phase_widths':[2,3],
                'nonpalindromic_bit_order_counterexample':'p=1/4, m=3: omitting the MSB-to-LSB Fourier relabelling changes the phase law by total variation greater than 0.1',
                'wrong_controlled_global_sign_counterexample':
                'For p=1 and m=3, omitting Z(control) moves the phase atom from 4 to 0',
                'published_good_probability_numerical_sweep_cases':505,
                'outward_phase_output_cases':1022,
                'finite_binary_plan_cases':16,
                'retained_emitted_fourier_profiles':32},
            'source_sha256':{},
            'scope':'coarse numerical coherent replay and exact finite-planner checks; no physical quantum execution; theorem and interval certificates supply guarantees',
            'examples':[],'actual_compiled_examples':[]}
    files=('qpe_estimator.py','checks/test_qpe_estimator.py',
           'checks/audit_qpe_estimator.py','qpe_estimator_proof.tex',
           'qpe_estimator_sources.md','certified_scalar.py','quantum_gates.py',
           'qft_synthesis.py','results/qft_rotation_tapes.json',
           'results/qft_rotation_certificates.json')
    for name in files:
        record['source_sha256'][name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    for epsilon in (F(1,200),F(1,2000),F(1,5000)):
        for prior in (F(1),F(1,10),F(1,100)):
            record['examples'].append(plan(epsilon,prior_upper=prior).record())
            setup,profile=compiled_plan(epsilon,prior_upper=prior)
            record['actual_compiled_examples'].append(setup.record())
    out=ROOT/'results'/'qpe_estimator_audit.json'
    out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'status':record['status'],'test_groups':result.testsRun,
                      'result_file':str(out)}))
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__=='__main__':
    main()
