"""Independent Acb384 certificates of actual toy rotation tapes.

Targets are rebuilt from the exact finite Gaussian-density model, not copied
from the numerical synthesis input. State replay is a separate diagnostic.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib,json,sys
from flint import arb,acb,acb_mat,fmpq,ctx
ROOT=Path(__file__).resolve().parent;ctx.prec=384
data=json.loads((ROOT/'results/coherent_tilt_toy_compiled.json').read_text())


def rat(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))


def phase(x):return acb(0,arb.pi()*rat(x)).exp()


def mat(g):
    if g=='H':return acb_mat([[1,1],[1,-1]])/arb(2).sqrt()
    angle={'S':F(1,2),'Sdg':F(-1,2),'T':F(1,4),'Tdg':F(-1,4)}[g]
    return acb_mat([[1,0],[0,phase(angle)]])


def endpoint(x,upper):
    a,b=(x.upper() if upper else x.lower()).man_exp();a,b=int(a),int(b)
    return F(a*(1<<b)) if b>=0 else F(a,1<<(-b))


points=list(map(rat,['-3/2','-1/2','1/2','3/2']))
weights=[(-y*y/2).exp() for y in points];total=sum(weights,arb(0));weights=[p/total for p in weights]
delta=rat('1/4');z=sum((p*(2*delta*y).exp() for p,y in zip(weights,points)),arb(0));B=z/2
targets=[{'role':'selector_uncontrolled_top','angle':2*(B/(1+B)).sqrt().atan()},
         {'role':'selector_controlled_half','angle':B.sqrt().atan()}]
for v in [arb(0),delta,-delta]:
    numerators=[(-y*y/2+2*v*y).exp() for y in points]
    top=((numerators[2]+numerators[3])/(numerators[0]+numerators[1])).sqrt().atan()
    low0=(v*(points[1]-points[0])-(points[1]*points[1]-points[0]*points[0])/4).exp().atan()
    low1=(v*(points[3]-points[2])-(points[3]*points[3]-points[2]*points[2])/4).exp().atan()
    for role,angle in [('coordinate_top_half',top),('coordinate_lower0_half',low0),('coordinate_lower1_half',low1)]:
        targets.append({'role':role,'tilt':str(v),'angle':angle})
certificates=[];error_by_key={}
for key,record in data['rotation_tapes'].items():
    near=rat(key);candidate=min(targets,key=lambda r:float(abs(r['angle']-near).upper()))
    assert abs(candidate['angle']-near).upper()<rat('1e-85')
    theta=candidate['angle'];c=(theta/2).cos();s=(theta/2).sin()
    target=acb_mat([[c,-s],[s,c]])
    U=acb_mat([[1,0],[0,1]])
    for g in record['chronological_gates']:U=mat(g)*U
    U=phase(record['phase_pi'])*U
    diff=U-target
    norm=sum((abs(a).upper()*abs(a).upper() for a in diff.entries()),arb(0)).nonnegative_part().sqrt().upper()
    upper=endpoint(norm,True)
    assert upper<F(1,1<<27)
    assert Counter(record['chronological_gates'])==Counter(record['gate_counts'])
    assert record['T_count']==sum(g in ('T','Tdg') for g in record['chronological_gates'])
    error_by_key[key]=upper
    certificates.append({'rotation_key':key,'target_role':candidate['role'],
                         'operator_norm_upper':str(upper),'target_angle_lower':str(endpoint(theta,False)),
                         'target_angle_upper':str(endpoint(theta,True)),'precision_bits':384,'passed':True})

# Sum every actual use, including both signs and repetitions.
logical=json.loads((ROOT/'results/coherent_tilt_toy_tape.json').read_text())
import mpmath as mp
mp.mp.dps=100
use_count=Counter(mp.nstr(abs(mp.mpf(g[2])),90) for g in logical['gates'] if g[0]=='RY')
zeta=sum((n*error_by_key[key] for key,n in use_count.items()),F(0))
assert data['aggregate_phase_pi']=='0'
counts=Counter(g[0] for g in data['chronological_gates'])
assert counts==Counter(data['gate_counts']) and data['actual_T_count']==counts['T']+counts['Tdg']
out={'status':'passed','rotation_certificates':certificates,'actual_rotation_uses':sum(use_count.values()),
     'LCU_preparation_operator_norm_upper':str(zeta),'actual_T_count':data['actual_T_count'],
     'global_phase_correction_scope':'Positive and negative rotation uses pair as exact adjoints; selector PREP and PREP† phases also cancel. Aggregate retained phase is exactly zero.',
     'compiled_sha256':hashlib.sha256((ROOT/'results/coherent_tilt_toy_compiled.json').read_bytes()).hexdigest(),
     'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'scope':'Certifies signal preparation only. No remainder evaluator or end-to-end quantum pricing estimator is supplied.'}
(ROOT/'results/coherent_tilt_toy_rotation_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'passed','rotations':len(certificates),'uses':sum(use_count.values()),'LCU_norm_upper':float(zeta),'T':data['actual_T_count']}))
