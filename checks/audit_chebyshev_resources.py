"""Independent complete-resource replay for bounded Chebyshev histories.

The retained exact QPE and comparator formulas below reproduce the existing
audit, specialised without importing either producer or charging API. New
parity samplers, actual histories, full reset, tape counts and random words
are checked independently. Full financial table transcendental certification
remains the separately retained four-survivor audit, not repeated here.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from math import comb
from functools import lru_cache
import gzip,hashlib,json
HERE=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def lo(iv):return F(iv['lower_exact_fraction'])
def hi(iv):return F(iv['upper_exact_fraction'])
def clean(c):return {k:v for k,v in c.items() if v}
def failure(n,s):return sum((F(comb(n,j))*(1-s)**j*s**(n-j) for j in range((n+1)//2,n+1)),F(0))
PI=F(355,113); EPS=F(1,100); ALPHA=F(1,3600)

def fourier(bits,tapes,bounds):
    seq=[];error=F(0);phase=F(0);retained=omitted=0
    for j in range(bits//2):
        k=bits-1-j;seq.extend([('CX',j,k),('CX',k,j),('CX',j,k)])
    for j in range(bits-1,-1,-1):
        for k in range(bits-1,j,-1):
            d=k-j
            if d>31:error+=PI/F(1<<d);omitted+=1;continue
            rm=tapes['rotations'][f'-1:{d+1}'];rp=tapes['rotations'][f'1:{d+1}']
            minus=rm['chronological_gates'];plus=rp['chronological_gates']
            if minus or plus:
                seq.extend((g,k) for g in minus);seq.extend((g,j) for g in minus)
                seq.append(('CX',k,j));seq.extend((g,j) for g in plus);seq.append(('CX',k,j))
            error+=2*bounds[f'-1:{d+1}']+bounds[f'1:{d+1}'];retained+=1
            phase+=-F(1,1<<(d+2))+2*F(rm['phase_pi'])+F(rp['phase_pi'])
        seq.append(('H',j))
    return seq,error,phase,retained,omitted

def reflection(core):
    # Reconstruct the conjunction ladder and the literal exact CCX expansion.
    ctrl=core;work=list(range(core+1,2*core));body=[]
    body.append(('CCX',0,1,work[0]))
    for j in range(2,core):body.append(('CCX',work[j-2],j,work[j-1]))
    logical=[('CZ',ctrl,core-1)]+[('X',j) for j in range(core)]+body
    logical+=[('CZ',work[core-2],ctrl)]+body[::-1]+[('X',j) for j in range(core)]+[('Z',ctrl)]
    seq=[]
    for g in logical:
        if g[0]!='CCX':seq.append(g);continue
        _,a,b,t=g
        seq.extend([('H',t),('CX',b,t),('Tdg',t),('CX',a,t),('T',t),
            ('CX',b,t),('Tdg',t),('CX',a,t),('T',b),('T',t),
            ('H',t),('CX',a,b),('T',a),('Tdg',b),('CX',a,b)])
    return seq,2*(core-1)

def verify_profile(profile,known_tables):
    forward=Counter()
    for b in profile['blocks']:
        c=Counter(b['primitive_counts']);assert b['logical_T_count']==c['T']+c['Tdg']
        assert not set(c)-{'X','H','CX','CZ','Z','S','Sdg','T','Tdg'}
        if b['table_sha256'] is not None:assert b['table_sha256'] in known_tables
        forward.update(c)
    assert clean(forward)==clean(Counter(profile['primitive_counts_forward_body']))
    complete=Counter({g:2*n for g,n in forward.items()})
    complete.update(profile['comparator_primitive_counts']);complete['H']+=profile['uniform_H']
    assert clean(complete)==clean(Counter(profile['primitive_counts']))
    assert complete['T']+complete['Tdg']==profile['logical_T_count']
    assert profile['explicit_T_depth_serial_upper']==2*sum(b['explicit_T_depth_serial_upper'] for b in profile['blocks'])+profile_comparator_depth(profile,'explicit')
    assert profile['conservative_5layer_T_depth_serial_upper']==2*sum(b['conservative_5layer_T_depth_serial_upper'] for b in profile['blocks'])+profile_comparator_depth(profile,'five')
    return complete

@lru_cache(None)
def comparator(width):
    # Construct the subtraction ripple from Boolean operations independently.
    a=list(range(width));b=list(range(width,2*width));flag=2*width
    next_wire=flag+1;carry=next_wire;next_wire+=1;body=[]
    for ai,bi in zip(a,b):
        si=next_wire;next_wire+=1
        body.extend(('CX',x,si) for x in (ai,bi,carry))
        na=next_wire;next_wire+=1;body.extend([('X',na),('CX',ai,na)])
        nc=next_wire;next_wire+=1
        body.extend([('CCX',na,bi,nc),('CCX',na,carry,nc),('CCX',bi,carry,nc)])
        carry=nc
    logical=body+[('CX',carry,flag)]+body[::-1];sequence=[]
    five=[0]*next_wire
    for g in logical:
        level=max(five[w] for w in g[1:])+(5 if g[0]=='CCX' else 0)
        for w in g[1:]:five[w]=level
        if g[0]!='CCX':sequence.append(g);continue
        _,aa,bb,t=g
        sequence.extend([('H',t),('CX',bb,t),('Tdg',t),('CX',aa,t),('T',t),
            ('CX',bb,t),('Tdg',t),('CX',aa,t),('T',bb),('T',t),
            ('H',t),('CX',aa,bb),('T',aa),('Tdg',bb),('CX',aa,bb)])
    explicit=[0]*next_wire
    for g in sequence:
        level=max(explicit[w] for w in g[1:])+(g[0] in ('T','Tdg'))
        for w in g[1:]:explicit[w]=level
    return Counter(g[0] for g in sequence),max(explicit),max(five)

def profile_comparator_depth(profile,mode):
    counts,explicit,five=comparator(profile['comparison_bits']+1)
    assert clean(counts)==clean(Counter(profile['comparator_primitive_counts']))
    return explicit if mode=='explicit' else five

def verify_plan(flag,profile,tapes,bounds):
    plan=flag['plan'];R=plan['repetitions'];M=plan['M'];m=plan['phase_register_bits']
    assert M==1<<m and R>0 and R%2==1
    epsilon=F(flag['probability_tolerance']);prior=F(flag['prior_upper']);alpha=F(plan['alpha'])
    assert F(plan['epsilon'])==epsilon and F(plan['prior_upper'])==prior
    variance=prior*(1-prior) if prior<=F(1,2) else F(1,4)
    assert F(plan['variance_upper'])==variance and prior>2*epsilon
    seq,error,phase,retained,omitted=fourier(m,tapes,bounds);counts=Counter(g[0] for g in seq)
    qft=flag['fourier_profile'];assert digest(profile)==flag['financial_profile_sha256']
    assert qft['chronological_tape_sha256']==hashlib.sha256(json.dumps(seq,separators=(',',':')).encode()).hexdigest()
    assert error==F(qft['norm_error_upper']) and phase==F(qft['global_phase_correction_pi'])
    assert clean(counts)==clean(Counter(qft['gate_counts']))
    assert qft['T_count']==counts['T']+counts['Tdg']
    assert (retained,omitted)==(qft['controlled_phases_retained'],qft['controlled_phases_omitted'])
    scale=1<<80;scaled=error*scale;rounded=F(-(-scaled.numerator//scaled.denominator),scale)
    assert F(plan['qft_norm_error_allowance'])==rounded and rounded<=F(1,64)
    success=8/PI**2-rounded;assert success==F(plan['good_probability_lower'])
    fail=failure(R,success);assert fail==F(plan['binomial_failure_upper']) and fail<=alpha
    if R>1:assert failure(R-2,success)>alpha
    ae=F(plan['amplitude_error_upper']);root=(ae-PI**2/M**2)*M/(2*PI)
    assert root>=0 and root*root>=variance
    assert ae+F(plan['postprocess_allowance'])<=epsilon
    if M>2:assert 2*PI*root/(M//2)+PI**2/(M//2)**2+F(plan['postprocess_allowance'])>epsilon
    ap=Counter(profile['primitive_counts']);adj=Counter()
    swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    for op,n in ap.items():adj[swap.get(op,op)]+=n
    refl,ccx=reflection(profile['reflection_core_qubits']);rc=Counter(g[0] for g in refl)
    qp=ap+adj+rc;per=ap+Counter({op:(M-1)*n for op,n in qp.items()})+counts;per['H']+=m
    charge=flag['complete_charge'];ref=charge['controlled_reflection']
    assert clean(rc)==clean(Counter(ref['gate_counts']))
    assert ref['chronological_tape_sha256']==hashlib.sha256(json.dumps(refl,separators=(',',':')).encode()).hexdigest()
    assert ref['T']==14*(profile['reflection_core_qubits']-1)
    peak=max(profile['logical_qubits_upper'],2*profile['reflection_core_qubits'])+m
    assert charge['peak_logical_qubits']==peak
    expected={'T':per['T']+per['Tdg'],'measurements':peak,'feedback_CZ':0,'reset_X':peak,
        'expanded_primitive_operations':sum(per.values())+2*peak,
        'uniform_H':profile['uniform_H']+(M-1)*2*profile['uniform_H']+m,
        'explicit_T_depth_serial_upper':profile['explicit_T_depth_serial_upper']+(M-1)*(2*profile['explicit_T_depth_serial_upper']+4*ccx)+qft['T_count'],
        'conservative_5layer_T_depth_serial_upper':profile['conservative_5layer_T_depth_serial_upper']+(M-1)*(2*profile['conservative_5layer_T_depth_serial_upper']+5*ccx)+qft['T_count']}
    expected.update({f'gate_{op}':n for op,n in per.items()})
    for op in ('H','S','Sdg','T','Tdg','CX'):expected.setdefault(f'gate_{op}',0)
    assert expected==charge['per_repetition_metrics']
    assert {k:R*v for k,v in expected.items()}==charge['metrics']
    assert charge['controlled_Q_calls']==R*(M-1)
    assert charge['complete_reinitialisations']==R and charge['phase_register_measurements']==R*m
    assert charge['complete_oracle_T']==ap['T']+ap['Tdg'] and charge['controlled_Q_T']==qp['T']+qp['Tdg']
    return charge['metrics']['T']


import argparse, cmath, math, random, statistics


def ccx_tape(a,b,t):
    return [('H',t),('CX',b,t),('Tdg',t),('CX',a,t),('T',t),
            ('CX',b,t),('Tdg',t),('CX',a,t),('T',b),('T',t),
            ('H',t),('CX',a,b),('T',a),('Tdg',b),('CX',a,b)]


def uncontrolled_reflection(core):
    body=[]
    if core>2:
        body=[('CCX',0,1,core)]
        for j in range(2,core-1):
            body.append(('CCX',core+j-2,j,core+j-1))
    tape=[('X',j) for j in range(core)]+body
    tape += [('CZ',0,1)] if core==2 else [('CZ',2*core-3,core-1)]
    tape+=body[::-1]+[('X',j) for j in range(core)]
    expanded=[]
    for g in tape:
        expanded+=ccx_tape(*g[1:]) if g[0]=='CCX' else [g]
    return tape,expanded


def multiply(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def vector(a,b):
    return [sum(a[i][j]*b[j] for j in range(2)) for i in range(2)]


def power(a,k):
    result=[[1.,0.],[0.,1.]]
    while k:
        if k&1:result=multiply(result,a)
        a=multiply(a,a);k//=2
    return result


def two_dimensional_checks():
    checked=0
    for theta in [0.,.01,.2,math.pi/4,1.5,math.pi/2]:
        c,s=math.cos(theta),math.sin(theta)
        A=[[c,-s],[s,c]];adj=[[c,s],[-s,c]]
        sg=[[1.,0.],[0.,-1.]];sz=[[-1.,0.],[0.,1.]]
        Q=multiply(multiply(multiply(A,sz),adj),sg)
        Q=[[-v for v in row] for row in Q]
        a=[c,s]
        for d in list(range(1,33))+[2047]:
            k=d//2
            if d%2:
                v=vector(power(Q,k),a);actual=abs(v[1])**2
                target=math.sin(d*theta)**2
            else:
                v=vector(adj,vector(sg,vector(power(Q,k-1),a)))
                actual=abs(v[0])**2;target=math.cos(d*theta)**2
            assert abs(actual-target)<1e-9
            assert abs(actual-math.cos(d*math.acos(s))**2)<1e-9
            checked+=1
    return checked


def ccx_expansion_checks():
    # Independent numerical replay of all basis inputs for one actual exact
    # seven-T expansion, including phase rather than Boolean counts alone.
    def replay(seq,state):
        state=state[:]
        for g in seq:
            op,*w=g
            if op=='H':
                bit=1<<w[0]
                for i in range(8):
                    if i&bit:continue
                    a,b=state[i],state[i|bit]
                    state[i]=(a+b)/math.sqrt(2);state[i|bit]=(a-b)/math.sqrt(2)
            elif op=='CX':
                target=1<<w[1]
                for i in range(8):
                    if i&target or not ((i>>w[0])&1):continue
                    state[i],state[i|target]=state[i|target],state[i]
            else:
                assert op in ('T','Tdg')
                phase=cmath.exp((1 if op=='T' else -1)*1j*math.pi/4)
                for i in range(8):
                    if (i>>w[0])&1:state[i]*=phase
        return state
    for b in range(8):
        state=[0j]*8;state[b]=1
        actual=replay(ccx_tape(0,1,2),state)
        target=b^(4 if b&1 and b&2 else 0)
        assert all(abs(v-(1 if i==target else 0))<1e-14 for i,v in enumerate(actual))
    checks=0
    for n in range(2,9):
        ideal,_=uncontrolled_reflection(n)
        for b in range(1<<n):
            bits=[(b>>j)&1 for j in range(n)]+[0]*(n-2);phase=1
            for op,*w in ideal:
                if op=='X':bits[w[0]]^=1
                elif op=='CCX':bits[w[2]]^=bits[w[0]]&bits[w[1]]
                else:
                    assert op=='CZ'
                    if bits[w[0]]&bits[w[1]]:phase=-phase
            assert bits[:n]==[(b>>j)&1 for j in range(n)] and not any(bits[n:])
            assert phase==(-1 if b==0 else 1)
            checks+=1
    return checks


def verify_sampler(p,profile,ref):
    d=p['degree'];k=d//2
    forward=Counter(profile['primitive_counts']);adj=Counter()
    for op,n in forward.items():adj[{'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}.get(op,op)]+=n
    f,a,z,g=(k+1,k,k,k) if d%2 else (k,k,k-1,k)
    assert (p['A_forward_calls'],p['A_adjoint_calls'],p['zero_reflection_calls'],p['good_Z_calls'])==(f,a,z,g)
    counts=Counter({op:f*n for op,n in forward.items()})
    counts.update({op:a*n for op,n in adj.items()})
    counts.update({op:z*n for op,n in ref['primitive_counts'].items()})
    counts['Z']+=g
    assert clean(counts)==p['primitive_counts']
    wires=max(profile['logical_qubits_upper'],2*profile['reflection_core_qubits']-2)
    assert p['peak_logical_qubits']==wires
    metrics={'T':counts['T']+counts['Tdg'],'measurements':wires,'reset_X':wires,'feedback_CZ':0,
        'expanded_primitive_operations':sum(counts.values())+2*wires,
        'uniform_H':d*profile['uniform_H'],
        'explicit_T_depth_serial_upper':d*profile['explicit_T_depth_serial_upper']+z*4*2*(profile['reflection_core_qubits']-2),
        'conservative_5layer_T_depth_serial_upper':d*profile['conservative_5layer_T_depth_serial_upper']+z*5*2*(profile['reflection_core_qubits']-2)}
    metrics.update({f'gate_{op}':n for op,n in counts.items() if n})
    assert metrics==p['per_shot_metrics']
    assert p['financial_preparation_sha256']==digest(profile)
    return metrics['T']


def main(repo,result_dir):
    from flint import arb,ctx
    ctx.prec=512
    protocol=json.loads((result_dir/'chebyshev_protocol.json').read_text())
    result=json.loads((result_dir/'chebyshev_resources.json').read_text())
    summary=json.loads((result_dir/'chebyshev_resource_summary.json').read_text())
    assert result['protocol_sha256']==sha(result_dir/'chebyshev_protocol.json')
    assert summary['result_sha256']==sha(result_dir/'chebyshev_resources.json')
    for name,h in protocol['source_sha256'].items():
        path=HERE/name if name.startswith('chebyshev_') else repo/name
        assert sha(path)==h
    for name,h in protocol['input_sha256'].items():assert sha(repo/name)==h
    retained=json.loads((repo/'results/four_survivor_quantum_resources.json').read_text())
    tapes=json.loads((repo/'results/qft_rotation_tapes.json').read_text())
    qft_cert=json.loads((repo/'results/qft_rotation_certificates.json').read_text())
    bounds={r['key']:F(r['operator_norm_upper_rational']) for r in qft_cert['certificates']}
    assert qft_cert['tapes_sha256']==sha(repo/'results/qft_rotation_tapes.json')
    history_count=batch_count=word_count=profile_count=signal_count=0
    for ci,case in enumerate(result['cases']):
        source=retained['selected_prices'][ci]['methods']['positive_r1']
        row=retained['records'][source['candidate_index']]
        original=source['selected']['result'];profile=row['profiles']['unitary8t']['centred']
        info=case['retained_financial_information'];prior=original['correction_prior'];original_flag=original['flags'][0]
        assert info==protocol['cases'][ci]
        assert info['financial_preparation_profile']==profile
        assert info['probability_lower']==prior['probability_lower'] and info['probability_upper']==prior['probability_upper']
        assert info['probability_tolerance']==original_flag['probability_tolerance']
        verify_profile(profile,{b['table_sha256'] for b in profile['blocks'] if b['table_sha256'] is not None})
        # Reuse checked synthetic financial preparations, without substituting
        # a sampled price or reconstructing new financial output tables.
        assert F(info['deterministic_price_error_upper'])+F(info['monetary_signal_factor_upper'])*F(info['probability_tolerance'])<=F(1,100)
        ideal,expanded=uncontrolled_reflection(profile['reflection_core_qubits'])
        ref=case['uncontrolled_zero_reflection']
        assert ref['ideal_gate_list']==[list(g) for g in ideal]
        assert ref['chronological_primitive_gate_list']==[list(g) for g in expanded]
        assert ref['chronological_tape_sha256']==digest(expanded)
        assert ref['primitive_counts']==clean(Counter(g[0] for g in expanded))
        assert ref['T_count']==14*(profile['reflection_core_qubits']-2)
        fallback=case['half_alpha_QPE_fallback'];plan=fallback['plan']
        assert F(plan['alpha'])==F(1,7200)
        assert F(plan['epsilon'])==F(info['probability_tolerance'])
        assert F(plan['prior_upper'])==F(info['probability_upper'])
        flag={'plan':plan,'probability_tolerance':info['probability_tolerance'],
            'prior_upper':info['probability_upper'],'financial_profile_sha256':digest(profile),
            'fourier_profile':fallback['fourier_profile'],'complete_charge':fallback['complete_charge']}
        assert verify_plan(flag,profile,tapes,bounds)==fallback['T_count']
        assert case['pilot_T_cap']==fallback['T_count']//10
        assert case['whole_history_T_cap']==case['pilot_T_cap']+fallback['T_count']
        costs={}
        at=profile['logical_T_count'];rt=ref['T_count']
        for d in range(1,2048):
            z=(d//2 if d%2 else d//2-1)
            costs[d]=100*(d*at+z*rt)
        assert digest(costs)==case['batch_cost_by_degree_sha256']
        assert len(case['diagnostics'])==100
        for key,p in case['sampler_profiles'].items():
            assert int(key)==p['degree'];assert verify_sampler(p,profile,ref)*100==costs[int(key)]
            profile_count+=1
        for hi,diag in enumerate(case['diagnostics']):
            history_count+=1
            mi,si=divmod(hi,20)
            assert diag['mesh_index']==mi and diag['seed_index']==si
            assert diag['seed']==1729+1000000*ci+1000*mi+si
            truth=F(info['probability_lower'])+F(mi,4)*(F(info['probability_upper'])-F(info['probability_lower']))
            assert F(diag['truth_probability'])==truth
            truth_ball=arb(truth.numerator)/truth.denominator;theta=truth_ball.sqrt().asin()
            rng=random.Random(diag['seed']);controller=diag['controller_result'];batches=diag['raw_batches'];charged=0
            assert len(controller['history'])==len(batches)<=12
            assert controller['pilot_shots']==100*len(batches)
            for b,h in zip(batches,controller['history']):
                batch_count+=1;d=b['degree'];words=b['seeded_u64_words'];signal=b['signal']
                assert len(words)==b['shots']==h['shots']==100
                assert words==[rng.getrandbits(64) for _ in range(100)]
                word_count+=100
                value=(theta*d).sin()**2 if d%2 else (theta*d).cos()**2
                lo=(value*(1<<64)).lower().floor().unique_fmpz();up=(value*(1<<64)).upper().floor().unique_fmpz()
                assert lo==up==signal['probability_q64_numerator'];signal_count+=1
                assert signal['probability_q64_denominator']==1<<64
                assert sum(w<int(lo) for w in words)==b['successes']==h['successes']
                assert b['batch_T_count']==h['complete_batch_T']==costs[d]
                charged+=costs[d];assert charged==h['cumulative_pilot_T']
                assert F(h['conditional_failure_upper'])==F(1,86400)
            assert charged==controller['pilot_T_count']<=case['pilot_T_cap']
            assert F(controller['pilot_failure_upper'])==len(batches)*F(1,86400)<=F(1,7200)
            needs=controller['status']!='certified'
            assert needs==diag['fallback_invoked_in_resource_history']
            assert diag['actual_history_T_count']==charged+(fallback['T_count'] if needs else 0)<=case['whole_history_T_cap']
            if not needs:
                assert F(controller['probability_radius'])<=F(info['probability_tolerance'])
            interval=controller['probability_interval']
            covered=F(interval[0])<=truth<=F(interval[1])
            assert covered==diag['pilot_interval_covers_truth']
        totals=[d['actual_history_T_count'] for d in case['diagnostics']]
        stops=sum(d['controller_result']['status']=='certified' for d in case['diagnostics'])
        s=case['diagnostic_summary']
        assert s['histories']==100 and s['pilot_certified']==stops and s['fallback_charged']==100-stops
        assert (s['observed_history_T_min'],s['observed_history_T_median'],s['observed_history_T_max'])==(min(totals),statistics.median(totals),max(totals))
    output={'schema':'qbp-chebyshev-resource-independent-audit-v1','status':'passed',
        'source_sha256':sha(__file__),'protocol_sha256':sha(result_dir/'chebyshev_protocol.json'),
        'result_sha256':sha(result_dir/'chebyshev_resources.json'),
        'retained_four_survivor_result_sha256':sha(repo/'results/four_survivor_quantum_resources.json'),
        'registered_input_hashes':protocol['input_sha256'],
        'two_dimensional_parity_checks':two_dimensional_checks(),
        'uncontrolled_reflection_basis_checks':ccx_expansion_checks(),
        'seven_T_CCX_phase_basis_checks':8,
        'financial_preparation_block_and_comparator_profiles':4,'fresh_half_alpha_QPE_plans':4,
        'sampler_profiles_reconstructed':profile_count,'diagnostic_histories':history_count,
        'complete_fresh_batches':batch_count,'pseudorandom_u64_words_replayed':word_count,
        'signal_q64_floors_reconstructed_at_512_bits':signal_count,
        'financial_mean_or_price_outputs_used':False,'full_wire_observation_and_reset_charged':True,
        'scope':'Independent parity, reflection, complete opcode counts, full QPE fallback, monetary decoding and every seeded history. Fixed-binomial inversion and controller scheduling are independently checked in the controller audit. Financial ordinate/financial-law proofs rely on their separately retained audits; no physical quantum execution or advantage is certified.'}
    (result_dir/'chebyshev_resource_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=HERE)
    parser.add_argument('--results',type=Path,default=HERE/'results')
    args=parser.parse_args();main(args.repo.resolve(),args.results.resolve())
