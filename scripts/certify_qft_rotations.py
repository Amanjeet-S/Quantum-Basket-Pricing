"""Independent Acb product certificates for retained chronological gate tapes."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import hashlib,json
from flint import arb,acb,acb_mat,fmpq,ctx
ctx.prec=384
ROOT=Path(__file__).resolve().parent
P=ROOT/'results'/'qft_rotation_tapes.json'
data=json.loads(P.read_text())
def rat(s):
    f=Fraction(s);return arb(fmpq(f.numerator,f.denominator))
def phase(s):return acb(0,arb.pi()*rat(s)).exp()
def mat(g):
    if g=='H':return acb_mat([[1,1],[1,-1]])/arb(2).sqrt()
    if g in ('S','Sdg','T','Tdg'):
        p={'S':'1/2','Sdg':'-1/2','T':'1/4','Tdg':'-1/4'}[g]
        return acb_mat([[1,0],[0,phase(p)]])
    raise ValueError(g)
def identity(n):return acb_mat([[int(i==j)for j in range(n)]for i in range(n)])
def frobenius(B):return sum((abs(z).upper()*abs(z).upper() for z in B.entries()),arb(0)).nonnegative_part().sqrt().upper().fmpq()
records=[]
for key,r in data['rotations'].items():
    U=identity(2)
    for g in r['chronological_gates']: U=mat(g)*U
    U=phase(r['phase_pi'])*U
    theta=arb.pi()*rat(r['theta_pi'])
    target=acb_mat([[acb(0,-theta/2).exp(),0],[0,acb(0,theta/2).exp()]])
    upper=frobenius(U-target)
    assert upper<fmpq(1,131072),(key,str(upper))
    assert Counter(r['chronological_gates'])==Counter(r['gate_counts'])
    assert sum(g in ('T','Tdg')for g in r['chronological_gates'])==r['T_count']
    records.append({'key':key,'operator_norm_upper_rational':str(upper),'certificate':'operator norm <= Frobenius norm of the full Acb-enclosed matrix difference','precision_bits':ctx.prec,'passed':True})
out={'schema':'qbp-qft-rotation-certificate-v1','tapes_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'arithmetic':'python-flint 0.9.0, Acb at 384 bits','independence':'Only the retained literal tapes are read. No synthesis library or project circuit/resource API is imported. The target pi and each H/S/T matrix are enclosed afresh.','certificates':records}
(ROOT/'results'/'qft_rotation_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
print(len(records),'passed')
