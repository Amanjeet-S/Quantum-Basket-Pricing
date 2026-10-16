"""Explicit signed LCU preparation on a finite product Gaussian toy model.

Classical terminal distribution: two independent four-point Gaussian-density
restrictions, equal-weight two-asset basket, common Gaussian smoothing.
Actual signal state uses two selector, four coordinate and two clean work
qubits. Arbitrary rotations are retained explicitly for separate synthesis.
"""
from pathlib import Path
from collections import Counter
from itertools import product
import hashlib,json,math,sys
ROOT=Path(__file__).resolve().parent
VENDOR=Path('/private/tmp/qbp-qft-synthesis/synthesis_vendor')
sys.path.insert(0,str(VENDOR))
import mpmath as mp
mp.mp.dps=100
Y=list(map(mp.mpf,['-1.5','-.5','.5','1.5']))
weights=[mp.exp(-y*y/2) for y in Y];normal=sum(weights)
weights=[p/normal for p in weights]


def z(v):return sum(p*mp.exp(2*v*y) for p,y in zip(weights,Y))


def coord_amplitudes(v):
    Z=z(v)
    return [mp.sqrt(p)*mp.exp(v*y)/mp.sqrt(Z) for p,y in zip(weights,Y)]


def coord_angles(v):
    a=coord_amplitudes(v)
    return [2*mp.atan(mp.sqrt((a[2]**2+a[3]**2)/(a[0]**2+a[1]**2))),
            2*mp.atan(a[1]/a[0]),2*mp.atan(a[3]/a[2])]


def adjoint(gates):
    output=[]
    for g in reversed(gates):
        if g[0]=='RY':output.append(('RY',g[1],mp.nstr(-mp.mpf(g[2]),105)))
        else:output.append(g)
    return output


def controlled_ry(controls,target,theta,work=(6,7)):
    controls=list(controls);g=[]
    if len(controls)==1:bit=controls[0];body=[]
    else:
        body=[('CCX',controls[0],controls[1],work[0])]
        for k in range(2,len(controls)):body.append(('CCX',work[k-2],controls[k],work[k-1]))
        bit=work[len(controls)-2]
    half=theta/2
    return body+[('RY',target,mp.nstr(half,105)),('CX',bit,target),
                 ('RY',target,mp.nstr(-half,105)),('CX',bit,target)]+body[::-1]


def build(delta):
    delta=mp.mpf(delta);B=z(delta)/2;L=1+2*B
    prep=[('RY',1,mp.nstr(2*mp.atan(mp.sqrt(B/(1+B))),105)),('X',1)]
    prep+=controlled_ry([1],0,2*mp.atan(mp.sqrt(B)))
    prep+=[('X',1)]
    select=[]
    for index,v in enumerate([mp.mpf(0),delta,-delta]):
        negate=[wire for wire in (0,1) if not (index>>wire)&1]
        select+=[('X',w) for w in negate]
        if index==0:select.append(('CZ',0,1))
        angles=coord_angles(v)
        for lo,hi in ((2,3),(4,5)):
            select+=controlled_ry([0,1],hi,angles[0])
            select+=[('X',hi)]+controlled_ry([0,1,hi],lo,angles[1])+[('X',hi)]
            select+=controlled_ry([0,1,hi],lo,angles[2])
        select+=[('X',w) for w in reversed(negate)]
    tape=prep+select+adjoint(prep)
    states=[];M=mp.mpf(0)
    for j2,j1 in product(range(4),repeat=2):
        p=weights[j1]*weights[j2];u=mp.cosh(delta*(Y[j1]+Y[j2]))-1
        states.append(mp.sqrt(p)*u/L);M+=p*u*u
    eta=M/(L*L)
    theta=mp.asin(mp.sqrt(eta));aa=max(0,int(mp.floor(mp.pi/(4*theta)-mp.mpf('.5'))))
    return tape,{'delta':mp.nstr(delta,105),'coordinate_points':list(map(str,Y)),
                 'coordinate_probabilities':list(map(lambda x:mp.nstr(x,105),weights)),
                 'component_coefficients':list(map(lambda x:mp.nstr(x,105),[-1,B,B])),
                 'L':mp.nstr(L,105),'M2':mp.nstr(M,105),'postselection_probability':mp.nstr(eta,105),
                 'unnormalised_success_amplitudes':list(map(lambda x:mp.nstr(x,105),states)),
                 'success_state_amplitudes':list(map(lambda x:mp.nstr(x/mp.sqrt(eta),105),states)),
                 'standard_AA_iterations':aa,'standard_AA_success_probability':mp.nstr(mp.sin((2*aa+1)*theta)**2,105),
                 'standard_AA_A_or_Aadjoint_calls':2*aa+1}


def replay(tape,state):
    import numpy as np
    state=np.asarray(state,dtype=complex).copy();N=len(state);indices=np.arange(N)
    for g in tape:
        kind=g[0]
        if kind=='RY':
            wire=g[1];theta=float(g[2]);a=indices[(indices&(1<<wire))==0];b=a|(1<<wire)
            x,y=state[a].copy(),state[b].copy();c,s=math.cos(theta/2),math.sin(theta/2)
            state[a],state[b]=c*x-s*y,s*x+c*y
        elif kind in ('X','CX','CCX'):
            *controls,target=g[1:];a=indices[(indices&(1<<target))==0]
            for wire in controls:a=a[(a&(1<<wire))!=0]
            b=a|(1<<target);state[a],state[b]=state[b].copy(),state[a].copy()
        elif kind in ('Z','CZ'):
            a=indices
            for wire in g[1:]:a=a[(a&(1<<wire))!=0]
            state[a]*=-1
        elif kind=='H':
            wire=g[1];a=indices[(indices&(1<<wire))==0];b=a|(1<<wire)
            x,y=state[a].copy(),state[b].copy();state[a],state[b]=(x+y)/math.sqrt(2),(x-y)/math.sqrt(2)
        elif kind in ('S','Sdg','T','Tdg'):
            phase={'S':1j,'Sdg':-1j,'T':complex(math.cos(math.pi/4),math.sin(math.pi/4)),'Tdg':complex(math.cos(-math.pi/4),math.sin(-math.pi/4))}[kind]
            state[(indices&(1<<g[1]))!=0]*=phase
        else:raise ValueError(kind)
    return state


def financial_identity(delta):
    delta=mp.mpf(delta);K=mp.mpf(100);H=mp.mpf(10);lam=mp.mpf('.2');forward=K*mp.exp(mp.mpf('.01'))
    # Exact declared finite model, with asset means matched to the forward.
    z1=sum(p*mp.exp(delta*y) for p,y in zip(weights,Y))
    mu=mp.log(forward/K)-lam*lam/2-2*mp.log(z1)
    Phi=lambda a:(1+mp.erf(a/mp.sqrt(2)))/2
    def call(m,strike):
        a=(m-mp.log(strike/K))/lam
        return K*mp.exp(m+lam*lam/2)*Phi(a+lam)-strike*Phi(a)
    def g(t):
        m=mu+mp.log(t)
        return call(m,K)-call(m,K+H)
    g0=g(1);slope=mp.diff(g,mp.mpf(1));EU=mp.mpf(0);price=mp.mpf(0);M2=mp.mpf(0);rem=mp.mpf(0)
    for j1,j2 in product(range(4),repeat=2):
        p=weights[j1]*weights[j2];u=mp.cosh(delta*(Y[j1]+Y[j2]))-1
        R=g(1+u)-g0-slope*u
        EU+=p*u;price+=p*g(1+u);rem+=p*R;M2+=p*u*u
    controlled=g0+slope*EU
    assert abs(price-controlled-rem)<mp.mpf('1e-90')
    return {'mu':mp.nstr(mu,105),'forward':mp.nstr(forward,105),
            'model_asset_forward_check':mp.nstr(K*mp.exp(mu+lam*lam/2)*z1*z1,105),
            'undiscounted_price':mp.nstr(price,105),'classical_order1_price':mp.nstr(controlled,105),
            'correction_mean':mp.nstr(rem,105),'M2':mp.nstr(M2,105),
            'payoff_bounds':{'lower':'0','upper':'10'},
            'no_arbitrage_scope':'Exact positive terminal measure and forward normalisation; no claim that independently rounded cross-contract estimates satisfy every joint constraint.'}


def main():
    import numpy as np
    tape,record=build('1/4');state=np.zeros(256,complex);state[0]=1
    out=replay(tape,state);indices=[j<<2 for j in range(16)]
    wanted=np.array(list(map(float,record['unnormalised_success_amplitudes'])),complex)
    discrepancy=float(np.linalg.norm(out[indices]-wanted));eta=sum(abs(out[i])**2 for i in indices)
    dirty=sum(abs(out[i])**2 for i in range(256) if i>>6)
    assert discrepancy<1e-13 and dirty<1e-26 and abs(eta-float(record['postselection_probability']))<1e-13
    counts=Counter(g[0] for g in tape)
    assert counts['CCX']==60 and counts['RY']==42
    record.update({'logical_qubits':8,'high_level_gate_counts':dict(counts),
                   'exact_7T_Toffoli_T_count':7*counts['CCX'],'logical_state_norm_diagnostic':float(np.linalg.norm(out)),
                   'success_amplitude_discrepancy_diagnostic':discrepancy,'dirty_scratch_probability_diagnostic':float(dirty),
                   'financial_identity':financial_identity('1/4'),
                   'weak_dispersion_sweep':[build(v)[1] for v in ['1/4','1/8','1/16','1/32','1/64']]})
    (ROOT/'results/coherent_tilt_toy_tape.json').write_text(json.dumps({'schema':'qbp-signed-lcu-toy-v1','gates':tape,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    (ROOT/'results/coherent_tilt_toy.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'logical_qubits':8,'counts':dict(counts),'eta':float(eta),'amplitude_discrepancy':discrepancy,'dirty_probability':float(dirty)}))


if __name__=='__main__':main()
