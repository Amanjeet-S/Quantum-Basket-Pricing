"""Literal compiled-circuit diagnostics, including charged amplification.

The Acb rotation certificate, not this floating replay, supplies gate-error
rigour. The complete remainder observable and estimator are not implemented.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from tilt_gate_replay import replay
import quantum_gates as qg
data=json.loads((ROOT/'results/coherent_tilt_toy_compiled.json').read_text())
ideal=json.loads((ROOT/'results/coherent_tilt_toy.json').read_text())
cert=json.loads((ROOT/'results/coherent_tilt_toy_rotation_certificates.json').read_text())
A=[tuple(g) for g in data['chronological_gates']]
state=np.zeros(256,complex);state[0]=1
out=replay(A,state);indices=np.array([j<<2 for j in range(16)])
wanted=np.array(list(map(float,ideal['unnormalised_success_amplitudes'])),complex)
discrepancy=float(np.linalg.norm(out[indices]-wanted))
zeta=F(cert['LCU_preparation_operator_norm_upper'])
assert discrepancy<=float(zeta)+1e-12
assert sum(abs(out[i])**2 for i in range(256) if i>>6)<1e-20

r=ideal['standard_AA_iterations']
good=qg.expand_toffoli(qg.zero_phase((0,1),()))
zero=qg.expand_toffoli(qg.zero_phase(tuple(range(6)),tuple(range(6,10))))
Q=good+qg.adjoint(A)+zero+A
state=np.zeros(1024,complex);state[0]=1
state=replay(A,state)
for _ in range(r):state=replay(Q,state)
# Q omits the irrelevant uncontrolled global minus. Keep it in the comparison.
state*=(-1)**r
kappa=sum(abs(state[j<<2])**2 for j in range(16))
ideal_kappa=float(ideal['standard_AA_success_probability'])
aa_error=(2*r+1)*zeta
assert abs(kappa-ideal_kappa)<=2*float(aa_error)+1e-10
scratch_probability=sum(abs(state[i])**2 for i in range(1024) if i>>6)
assert scratch_probability<1e-19
counts=Counter(g[0] for g in A)
reflection_counts=Counter(g[0] for g in good+zero)
total=Counter({g:(2*r+1)*n for g,n in counts.items()});total.update({g:r*n for g,n in reflection_counts.items()})
T=total['T']+total['Tdg']
assert T==(2*r+1)*data['actual_T_count']+r*56

# Equally informed base-state opponent: the same product Gaussian restriction
# is prepared directly, with no selector and no cancellation postselection.
base=[]
angles=[math.pi/2,2*math.atan(math.exp(.5)),2*math.atan(math.exp(-.5))]
def physical_ry(angle,wire):
    key=min(data['rotation_tapes'],key=lambda key:abs(float(key)-abs(angle)))
    assert abs(float(key)-abs(angle))<1e-12
    tape=[(g,wire) for g in data['rotation_tapes'][key]['chronological_gates']]
    return tape if angle>=0 else qg.adjoint(tape)
for lo,hi in ((0,1),(2,3)):
    # Ry(pi/2)=H Z, chronological Z then H, exactly.
    base+=[('Z',hi),('H',hi)]
    for branch,theta in [(0,angles[1]),(1,angles[2])]:
        if branch==0:base.append(('X',hi))
        base+=physical_ry(theta/2,lo)+[('CX',hi,lo)]+physical_ry(-theta/2,lo)+[('CX',hi,lo)]
        if branch==0:base.append(('X',hi))
base_state=np.zeros(16,complex);base_state[0]=1;base_out=replay(base,base_state)
weights=list(map(float,ideal['coordinate_probabilities']))
expected=np.array([math.sqrt(weights[j1]*weights[j2]) for j2 in range(4) for j1 in range(4)],complex)
base_discrepancy=float(np.linalg.norm(base_out-expected))
assert base_discrepancy<1e-7
base_counts=Counter(g[0] for g in base)
base_T=base_counts['T']+base_counts['Tdg']
assert base_T==736
toy_Umax=math.cosh(.75)-1
range_to_importance=float(toy_Umax**2/float(ideal['M2']))
out={'status':'passed','LCU_success_vector_discrepancy_diagnostic':discrepancy,
     'LCU_operator_norm_upper':str(zeta),'amplification_iterations':r,
     'amplification_A_or_Aadjoint_calls':2*r+1,'amplification_logical_qubits':10,
     'amplification_gate_counts':dict(total),'amplification_actual_T_count':T,
     'amplification_operator_norm_upper':str(aa_error),
     'amplified_success_probability_diagnostic':float(kappa),
     'amplified_ideal_success_probability_diagnostic':ideal_kappa,
     'amplified_dirty_scratch_probability_diagnostic':float(scratch_probability),
     'direct_base_preparation':{'logical_qubits':4,'gate_counts':dict(base_counts),
                                'actual_T_count':base_T,'state_discrepancy_diagnostic':base_discrepancy},
     'matched_range_weight_ratio_Umax_squared_to_M2_diagnostic':range_to_importance,
     'ideal_A_query_ratio_importance_to_direct_range_diagnostic':(2*r+1)/range_to_importance,
     'signal_T_ratio_importance_to_direct_range_diagnostic':T/(base_T*range_to_importance),
     'compiled_sha256':hashlib.sha256((ROOT/'results/coherent_tilt_toy_compiled.json').read_bytes()).hexdigest(),
     'certificates_sha256':hashlib.sha256((ROOT/'results/coherent_tilt_toy_rotation_certificates.json').read_bytes()).hexdigest(),
     'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'dependency_source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['tilt_gate_replay.py','quantum_gates.py']},
     'scope':'Signal preparation and its amplification. Monetary observable, confidence postprocessing, sampling/reset and end-to-end pricing remain unimplemented.'}
(ROOT/'results/coherent_tilt_toy_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'passed','LCU_discrepancy':discrepancy,'AA_T':T,'AA_qubits':10,'AA_success':float(kappa),'AA_error_upper':float(aa_error)}))
