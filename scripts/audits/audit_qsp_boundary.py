"""Exact endpoint/Boolean audits and diagnostic small projector-phase matrices.

The matrix diagnostics use floating arithmetic; the mathematical identity and
probability counts are proved separately. No QSP performance result is implied.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import importlib
import json
import math
import os
import sys

REPO = Path(os.environ.get('QBP_PROJECT_SOURCE', str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(REPO))
vo = importlib.import_module('value_oracle')


def exact_comparator_audit():
    rows = []
    total = 0
    for c in range(7):
        circuit = vo.ValueCircuit(mode='exact7t')
        a = circuit.alloc(c+1)
        j = circuit.alloc(c)
        flag = circuit.alloc(1)[0]
        vo.range_predicate_into(circuit, a, j, c, flag)
        basis_cases = 0
        successes = []
        for value in range((1 << c)+1):
            good = 0
            for J, dirty in product(range(1 << c), (0, 1)):
                state = circuit.run({a:value, j:J, (flag,):dirty})
                wanted = int(J < value)
                assert state[flag] == dirty ^ wanted
                assert circuit.value(state,a) == value
                assert circuit.value(state,j) == J
                assert all(bit==0 for w,bit in enumerate(state)
                           if w not in a+j+(flag,))
                basis_cases += 1
                if dirty == 0:
                    good += state[flag]
            assert good == value
            successes.append(str(F(good,1 << c)))
        literal_source_zero = [sum(not(value<J) for J in range(1 << c))
                               for value in range(1 << c)]
        assert literal_source_zero == list(range(1,(1 << c)+1))
        resources = circuit.resources()
        rows.append({'comparator_bits':c,'valid_dirty_basis_cases':basis_cases,
                     'all_valid_value_success_probabilities':successes,
                     'literal_source_zero_counts':literal_source_zero,
                     'logical_T_count':resources['logical_T_count'],
                     'literal_expanded_gate_counts':resources['clifford_T_expanded_gate_counts'],
                     'scratch_clean':True,'endpoint_one_exact':True})
        total += basis_cases
    return {'basis_cases':total,'profiles':rows}


def exact_loading_obstruction():
    data = json.loads((REPO/'results/generalisation_protocol.json').read_text())
    rows = []
    for model in data['models']:
        if model['assets'] not in (4,8):
            continue
        w = list(map(F,model['weights']))
        R = [list(map(F,r)) for r in model['loadings']]
        rg = [sum(wi*r[j] for wi,r in zip(w,R)) for j in range(len(R[0]))]
        v = [rg[1],-rg[0]]+[F(0)]*(len(rg)-2)
        assert sum(a*b for a,b in zip(rg,v)) == 0
        exposures = [sum(a*b for a,b in zip(r,v)) for r in R]
        assert any(z!=0 for z in exposures)
        rows.append({'model_id':model['model_id'],'assets':model['assets'],
                     'v':list(map(str,v)), 'weighted_loading_dot_v':'0',
                     'asset_loading_dot_v':list(map(str,exposures)),
                     'strict_convexity_criterion':True})
    assert len(rows) == 6
    return rows


def signed_phase_diagnostic():
    maximum = 0.0
    for phi in (0.,.37,-1.1,math.pi/2,math.pi):
        for b,s in product((0,1),repeat=2):
            # CX, Rz(2phi), CX: selector eigenvalue changes under XOR.
            physical = complex(math.cos(-phi*(1-2*(s^b))),
                               math.sin(-phi*(1-2*(s^b))))
            angle = phi*(1-2*s)*(2*b-1)
            desired = complex(math.cos(angle),math.sin(angle))
            maximum = max(maximum,abs(physical-desired))
    assert maximum < 1e-14
    return {'basis_cases':20,'maximum_absolute_matrix_entry_discrepancy':maximum,
            'scope':'floating diagnostic; identity proved algebraically'}


if __name__ == '__main__':
    paths = ['value_oracle.py','architecture.py','reversible.py',
             'results/generalisation_protocol.json']
    result = {'status':'passed','comparator':exact_comparator_audit(),
              'higher_dimensional_obstruction':exact_loading_obstruction(),
              'signed_projector_phase':signed_phase_diagnostic(),
              'source_hashes':{p:hashlib.sha256((REPO/p).read_bytes()).hexdigest() for p in paths},
              'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'no_nontrivial_QSP_baseline_implemented':True}
    output = Path(__file__).resolve().parents[1]/'results/qsp_boundary_audit.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'basis_cases':result['comparator']['basis_cases'],
                      'obstruction_models':len(result['higher_dimensional_obstruction']),
                      'output':str(output)}))
