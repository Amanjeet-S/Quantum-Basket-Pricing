"""Fresh small-register coherent matrices and literal tape/count audits."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import hashlib,json,math
from flint import arb,acb,acb_mat,fmpq,ctx
import qft_synthesis as api
ctx.prec=384
ROOT=Path(__file__).resolve().parent

def rat(f):
    f=Fraction(str(f));return arb(fmpq(f.numerator,f.denominator))
def exp_pi(f):return acb(0,arb.pi()*rat(f)).exp()
def I(n):return acb_mat([[int(i==j)for j in range(n)]for i in range(n)])
def gate(g):
    if g=='H':return acb_mat([[1,1],[1,-1]])/arb(2).sqrt()
    phases={'S':Fraction(1,2),'Sdg':Fraction(-1,2),'T':Fraction(1,4),'Tdg':Fraction(-1,4)}
    return acb_mat([[1,0],[0,exp_pi(phases[g])]])
def norm(B):
    squared=sum((abs(z).upper()*abs(z).upper() for z in B.entries()),arb(0))
    return squared.upper().nonnegative_part().sqrt().upper().fmpq()
def apply_single(U,g,w,n):
    stride=1<<(n-1-w)
    for x in range(U.nrows()):
        if x&stride:continue
        y=x|stride
        for c in range(U.ncols()):
            a,b=U[x,c],U[y,c];U[x,c]=g[0,0]*a+g[0,1]*b;U[y,c]=g[1,0]*a+g[1,1]*b

def apply_cx(U,control,target,n):
    cb=1<<(n-1-control);tb=1<<(n-1-target)
    for x in range(U.nrows()):
        if x&cb and not x&tb:
            y=x|tb
            for c in range(U.ncols()):U[x,c],U[y,c]=U[y,c],U[x,c]

def apply_two(U,V,control,target,n):
    cb=1<<(n-1-control);tb=1<<(n-1-target)
    for x in range(U.nrows()):
        if x&(cb|tb):continue
        ix=[x,x|tb,x|cb,x|cb|tb]
        for c in range(U.ncols()):
            col=[U[i,c]for i in ix]
            for i in range(4):U[ix[i],c]=sum((V[i,j]*col[j]for j in range(4)),acb(0))

def product(tape,n):
    U=I(1<<n)
    for op in tape:
        if op[0]=='CX':apply_cx(U,op[1],op[2],n)
        else:apply_single(U,gate(op[0]),op[1],n)
    return U

cp=[]; cp_cache={}
for d in range(1,32):
    for sign in (-1,1):
        record=api.controlled_phase(d,0,1,sign)
        U=product(record['gates'],2)*exp_pi(record['global_phase_correction_pi'])
        target=I(4);target[3,3]=exp_pi(Fraction(sign,2**d))
        error=norm(U-target)
        # Frobenius upper for two qubits may be sqrt(2) times the embedding operator upper.
        assert rat(error)<rat(record['norm_error_upper'])*arb(2).sqrt()+arb('1e-65'),d
        cp.append({'distance':d,'sign':sign,'Frobenius_error_upper_rational':str(error),'operator_error_upper_from_telescope':str(record['norm_error_upper']),'literal_gate_count':len(record['gates']),'passed':True})
        if sign==-1:cp_cache[d]=U

full=[]
for n in range(1,6):
    tape=api.inverse_qft_tape(n)
    # Use already independently multiplied literal two-qubit gate tapes as exact blocks.
    U=I(1<<n)
    for j in range(n//2):
        k=n-1-j;apply_cx(U,j,k,n);apply_cx(U,k,j,n);apply_cx(U,j,k,n)
    for j in range(n-1,-1,-1):
        for k in range(n-1,j,-1):apply_two(U,cp_cache[k-j],k,j,n)
        apply_single(U,gate('H'),j,n)
    N=1<<n
    target=acb_mat([[exp_pi(Fraction(-2*x*y,N))/arb(N).sqrt()for y in range(N)]for x in range(N)])
    error=norm(U-target)
    assert rat(error)<rat(tape['norm_error_upper'])*arb(N).sqrt()+arb('1e-65'),(n,str(error),str(tape['norm_error_upper']))
    # Independent literal full-tape replay for m<=3, including the recorded single scalar.
    if n<=3:
        literal=product(tape['gates'],n)*exp_pi(tape['global_phase_correction_pi'])
        assert rat(norm(literal-U))<arb('1e-65')
    full.append({'phase_bits':n,'basis_columns':N,'matrix_entries':N*N,'operator_target':'inverse DFT, U[x,y]=exp(-2*pi*i*x*y/2^m)/sqrt(2^m), wire 0 most significant','Frobenius_error_upper_rational':str(error),'literal_full_gate_replay':n<=3,'passed':True})
    print('QFT',n,'passed',flush=True)

counts=[]
for m in range(1,33):
    tape=api.inverse_qft_tape(m);p=api.inverse_qft_profile(m)
    c=Counter(op[0]for op in tape['gates'])
    assert dict(c)==p['gate_counts']
    assert c['T']+c['Tdg']==p['T_count']==p['serial_T_depth_upper']
    assert sum(c.values())==p['metrics']['expanded_primitive_operations']
    assert p['norm_error_upper']<=Fraction(1,64)
    counts.append({'phase_bits':m,'literal_operations':sum(c.values()),'T':c['T']+c['Tdg'],'passed':True})

out={'schema':'qbp-qft-coherent-independent-audit-v1','arithmetic':'Acb 384 bits, fresh exact gate matrices','source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [ROOT/'qft_synthesis.py',ROOT/'results'/'qft_rotation_tapes.json',ROOT/'results'/'qft_rotation_certificates.json']},'scope':'All 62 signed controlled phase matrices from literal tapes. Full inverse DFT all computational-basis columns m=1..5, thus every coherent input by linearity. Literal complete primitive replay m=1..3. All 32 production-width resource counts. No controlled Grover operation is audited here.','controlled_phase_records':cp,'full_qft_records':full,'count_records':counts}
(ROOT/'results'/'qft_synthesis_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('62 CP, 5 coherent QFT, 32 count audits passed')
