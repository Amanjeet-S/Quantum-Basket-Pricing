"""Compile every toy Ry into actual retained Clifford+T gates.

An existing isolated official pygridsynth1.2.0 runtime is used. No install,
licence purchase, or user environment change is performed.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent
VENDOR=Path('/private/tmp/qbp-qft-synthesis/synthesis_vendor')
sys.path.insert(0,str(VENDOR))
sys.path.insert(0,'/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
import mpmath as mp
from pygridsynth.gridsynth import gridsynth_gates
import quantum_gates as qg
mp.mp.dps=100
data=json.loads((ROOT/'results/coherent_tilt_toy_tape.json').read_text())
rotations={};compiled=[];global_phase=F(0)
for gate in data['gates']:
    if gate[0]!='RY':compiled.append(tuple(gate));continue
    theta=mp.mpf(gate[2]);sign=1 if theta>=0 else -1
    key=mp.nstr(abs(theta),90)
    if key not in rotations:
        raw=gridsynth_gates(theta=abs(theta),epsilon=mp.mpf(1)/(2**30),seed=1729,dps=100,up_to_phase=False)
        assert set(raw)<=set('HTSXW')
        rz=list(reversed(raw.replace('W','')))
        rz=[h for g in rz for h in (['H','S','S','H'] if g=='X' else [g])]
        # C=S H maps Z to Y. Chronological C†, Rz, C.
        ry=['Sdg','H']+rz+['H','S']
        phase=F(raw.count('W'),4)
        rotations[key]={'theta_decimal':key,'phase_pi':str(phase),'chronological_gates':ry,
                        'gate_counts':dict(Counter(ry)),'T_count':sum(g in ('T','Tdg') for g in ry)}
        print('new rotation',len(rotations),'T',rotations[key]['T_count'],flush=True)
    rotation=rotations[key];tape=[(g,gate[1]) for g in rotation['chronological_gates']]
    if sign<0:tape=qg.adjoint(tape)
    compiled+=tape;global_phase+=sign*F(rotation['phase_pi'])
assert global_phase==0,global_phase
compiled=qg.expand_toffoli(compiled)
counts=Counter(g[0] for g in compiled)
out={'schema':'qbp-signed-lcu-toy-compiled-v1','logical_qubits':8,'chronological_gates':compiled,
     'rotation_tapes':rotations,'aggregate_phase_pi':str(global_phase),
     'gate_counts':dict(counts),'actual_T_count':counts['T']+counts['Tdg'],
     'metrics':qg.resources(compiled,8),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'input_tape_sha256':hashlib.sha256((ROOT/'results/coherent_tilt_toy_tape.json').read_bytes()).hexdigest(),
     'synthesis':{'pygridsynth':'1.2.0','mpmath':'1.4.1','epsilon':'1/1073741824','up_to_phase':False,
                  'precision_decimal_digits':100,'seed':1729},
     'scope':'LCU signal preparation only. Financial observable oracle, final amplification phases, sampling/reset and estimator remain unimplemented.'}
(ROOT/'results/coherent_tilt_toy_compiled.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'rotations':len(rotations),'actual_T_count':out['actual_T_count'],'primitive_gates':sum(counts.values()),'aggregate_phase_pi':str(global_phase)}))
