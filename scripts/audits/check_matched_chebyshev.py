"""Independent literal-resource, financial-interface and paired-history audit.

No matched producer, controller or charging API is imported. Existing
independent reflection, sampler and QPE formula checks are reused unchanged.
Finite confidence root brackets are checked by exact integer comparisons.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
from math import comb
import argparse, gzip, hashlib, importlib.util, json, random, statistics, sys
sys.dont_write_bytecode = True

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v): return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p): return json.loads(gzip.decompress(Path(p).read_bytes()) if str(p).endswith('.gz') else Path(p).read_bytes())
def load(p):
    spec=importlib.util.spec_from_file_location('retained_independent_sampler',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

@lru_cache(None)
def tail_sign(p, n, start, end, target):
    p,target=F(p),F(target);a,D=p.numerator,p.denominator;b=D-a
    num=sum(comb(n,j)*a**j*b**(n-j) for j in range(start,end+1))
    diff=num*target.denominator-D**n*target.numerator
    return (diff>0)-(diff<0)

def cp_check(iv, successes):
    assert iv['method']=='fixed_sample_Clopper_Pearson_exact_rational'
    assert F(iv['tail_target'])==F(1,172800) and iv['inverse_bits']==64
    lo,hi=map(F,iv['lower_root_bracket']);u,v=map(F,iv['upper_root_bracket'])
    assert F(iv['lower'])==lo and F(iv['upper'])==v
    assert 0<=lo<=hi<=u<=v<=1 and hi-lo<=F(1,1<<64) and v-u<=F(1,1<<64)
    if successes:
        assert tail_sign(str(lo),100,successes,100,str(F(1,172800)))<=0
        assert tail_sign(str(hi),100,successes,100,str(F(1,172800)))>=0
    else: assert lo==hi==0
    if successes<100:
        assert tail_sign(str(u),100,0,successes,str(F(1,172800)))>=0
        assert tail_sign(str(v),100,0,successes,str(F(1,172800)))<=0
    else: assert u==v==1

def main(args):
    repo,folder=Path(args.repo),Path(args.results)
    independent=load(repo/'checks/audit_chebyshev_resources.py')
    protocol=read(folder/'matched_chebyshev_protocol.json')
    result=read(folder/'matched_chebyshev_resources.json.gz')
    assert result['status']=='completed' and result['protocol_sha256']==sha(folder/'matched_chebyshev_protocol.json')
    assert F(protocol['alpha_per_alternative_price'])==F(1,3600)
    assert F(protocol['all_twelve_alternative_union_failure_upper'])==F(1,300)
    for rel,h in protocol['evidence_sha256'].items():
        assert sha(Path(rel) if Path(rel).is_absolute() else repo/rel)==h
    tapes=read(repo/'results/qft_rotation_tapes.json')
    qft=read(repo/'results/qft_rotation_certificates.json')
    bounds={x['key']:F(x['operator_norm_upper_rational']) for x in qft['certificates']}
    from flint import arb,ctx
    ctx.prec=512
    histories=batches_count=profiles_count=words_count=0
    dominance=[];coefficients={};signal_cache={};summary=[]
    for ci,(case,pcase) in enumerate(zip(result['cases'],protocol['cases'])):
        assert case['model_id']==pcase['model_id']
        lower,upper=map(F,pcase['common_signal_interval'])
        assert [str(lower),str(upper)]==case['common_signal_interval']
        for name,m in case['methods'].items():
            info=pcase['methods'][name];p=info['financial_preparation_profile']
            assert digest(p)==info['financial_preparation_profile_sha256']==m['financial_preparation_profile_sha256']
            assert F(info['probability_lower'])<=lower<=upper<=F(info['probability_upper'])
            assert F(info['finite_ideal_mean_lower'])<=lower<=upper<=F(info['finite_ideal_mean_upper'])
            monetary=F(info['deterministic_price_error_upper'])+F(info['monetary_signal_factor_upper'])*F(info['probability_tolerance'])
            assert monetary<=F(1,100) and monetary==F(m['financial_monetary_error_check'])
            body=Counter()
            for block in p['blocks']:
                count=Counter(block['primitive_counts']);assert block['logical_T_count']==count['T']+count['Tdg'];body.update(count)
            assert body==Counter(p['primitive_counts_forward_body'])
            complete=Counter({k:2*v for k,v in body.items()})
            complete.update(p.get('comparator_primitive_counts',p.get('predicate_primitive_counts',{})))
            complete['H']+=p['uniform_H'];assert complete==Counter(p['primitive_counts'])
            at=complete['T']+complete['Tdg'];assert at==p['logical_T_count']
            ideal,tape=independent.uncontrolled_reflection(p['reflection_core_qubits'])
            ref=m['uncontrolled_zero_reflection']
            assert ref['ideal_gate_list']==[list(x) for x in ideal]
            assert ref['chronological_primitive_gate_list']==[list(x) for x in tape]
            assert ref['primitive_counts']==independent.clean(Counter(x[0] for x in tape))
            assert ref['chronological_tape_sha256']==digest(tape)
            rt=14*(p['reflection_core_qubits']-2);assert rt==ref['T_count']
            coefficients[ci,name]=(at,rt)
            fb=m['half_alpha_QPE_fallback'];plan=fb['plan']
            assert F(plan['alpha'])==F(1,7200)
            flag={'plan':plan,'probability_tolerance':info['probability_tolerance'],
                  'prior_upper':info['probability_upper'],'financial_profile_sha256':digest(p),
                  'fourier_profile':fb['fourier_profile'],'complete_charge':fb['complete_charge']}
            assert independent.verify_plan(flag,p,tapes,bounds)==fb['T_count']
            assert m['pilot_T_cap']==fb['T_count']//10
            assert m['whole_history_T_cap']==m['pilot_T_cap']+fb['T_count']
            costs={d:100*(d*at+(d//2 if d%2 else d//2-1)*rt) for d in range(1,2048)}
            assert digest(costs)==m['batch_cost_by_degree_sha256']
            for degree,profile in m['sampler_profiles'].items():
                assert independent.verify_sampler(profile,p,ref)*100==costs[int(degree)];profiles_count+=1
            assert len(m['diagnostics'])==100
            for idx,h in enumerate(m['diagnostics']):
                histories+=1;mi,si=divmod(idx,20)
                truth=lower+F(mi,4)*(upper-lower)
                assert h['mesh_index']==mi and h['seed_index']==si and F(h['truth_probability'])==truth
                assert h['seed']==1729+1000000*ci+1000*mi+si
                rng=random.Random(h['seed']);words=[[rng.getrandbits(64) for _ in range(100)] for _ in range(12)]
                ctrl=h['controller_result'];trace=ctrl['history'];raw=h['raw_batches'];paid=0
                assert len(trace)==len(raw)<=12 and ctrl['pilot_shots']==100*len(raw)
                for look,(b,t) in enumerate(zip(raw,trace)):
                    batches_count+=1;words_count+=100;d=b['degree']
                    assert b['look']==look+1 and b['shots']==t['shots']==100
                    assert b['paired_seeded_words_sha256']==digest(words[look])
                    key=ci,mi,d
                    if key not in signal_cache:
                        theta=(arb(truth.numerator)/truth.denominator).sqrt().asin()
                        signal=(theta*d).sin()**2 if d%2 else (theta*d).cos()**2
                        l=(signal*(1<<64)).lower().floor().unique_fmpz();u=(signal*(1<<64)).upper().floor().unique_fmpz()
                        assert l==u;signal_cache[key]=int(l)
                    threshold=signal_cache[key];assert threshold==b['signal']['probability_q64_numerator']
                    successes=sum(w<threshold for w in words[look])
                    assert successes==b['successes']==t['successes'];cp_check(t['signal_confidence'],successes)
                    assert costs[d]==b['batch_T_count']==t['complete_batch_T']
                    paid+=costs[d];assert paid==t['cumulative_pilot_T']
                    assert F(t['conditional_failure_upper'])==F(1,86400)
                assert paid==ctrl['pilot_T_count']<=m['pilot_T_cap']
                assert F(ctrl['pilot_failure_upper'])==F(len(raw),86400)<=F(1,7200)
                needs=ctrl['status']!='certified';assert needs==h['fallback_invoked_in_resource_history']
                assert h['actual_history_T_count']==paid+(fb['T_count'] if needs else 0)<=m['whole_history_T_cap']
                interval=list(map(F,ctrl['probability_interval']))
                assert h['pilot_interval_covers_stipulated_truth']==(interval[0]<=truth<=interval[1])
                if not needs: assert F(ctrl['probability_radius'])<=F(info['probability_tolerance'])
            totals=[h['actual_history_T_count'] for h in m['diagnostics']]
            stops=sum(h['controller_result']['status']=='certified' for h in m['diagnostics'])
            ss=m['diagnostic_summary'];assert ss['histories']==100 and ss['pilot_certified']==stops
            ordered=sorted(totals)
            assert F(ss['observed_history_T_median'])==F(ordered[49]+ordered[50],2)
            assert F(ss['observed_history_T_mean'])==F(sum(totals),100)
            summary.append({'model_id':case['model_id'],'method':name,'canonical_T':info['full_alpha_canonical_T'],
                            'half_alpha_fallback_T':fb['T_count'],'whole_history_cap':m['whole_history_T_cap'],
                            'pilot_certified':stops,'diagnostic_median_T':ss['observed_history_T_median']})
        ac,sc=coefficients[ci,'shared_threshold'];qc=2*ac+sc
        for rival in ['raw_payoff','numerical_conditional']:
            ab,sb=coefficients[ci,rival];qb=2*ab+sb
            assert ac<ab and qc<qb
            floor=min(F(ab-ac,ab),F(qb-qc,qb))
            for d in [1,2,3,4,2047,2048,10**6,10**6+1,2**40,2**40+1]:
                z=d//2 if d%2 else d//2-1;tc=d*ac+z*sc;tb=d*ab+z*sb
                assert tc<tb and F(tb-tc,tb)>=floor
            dominance.append({'model_id':case['model_id'],'rival':rival,'candidate_A_T':ac,'rival_A_T':ab,
                              'candidate_S0_T':sc,'rival_S0_T':sb,'candidate_pair_T':qc,'rival_pair_T':qb,
                              'uniform_matched_degree_saving_lower':str(floor),
                              'uniform_matched_degree_saving_per_cent':float(100*floor)})
    assert histories==1200
    output={'schema':'qbp-matched-chebyshev-independent-audit-v1','status':'passed','auditor_sha256':sha(__file__),
            'retained_independent_formula_source_sha256':sha(repo/'checks/audit_chebyshev_resources.py'),
            'protocol_sha256':sha(folder/'matched_chebyshev_protocol.json'),'result_sha256':sha(folder/'matched_chebyshev_resources.json.gz'),
            'histories':histories,'paired_mesh_seed_histories':400,'fresh_batches':batches_count,'replayed_u64_words':words_count,
            'sampler_profiles':profiles_count,'512_bit_signal_floors':len(signal_cache),'fresh_half_alpha_plans':12,
            'all_degree_dominance':dominance,'summary':summary,'imports_matched_producer_or_controller':False,
            'scope':'Exact preparation/profile bindings, independently reconstructed parity/reflection/opcode costs, complete QPE fallback, exact binomial root brackets, monetary interfaces and every paired synthetic history. All-degree claim follows the separate affine proof, not finite tests. Controller scheduling/inverse guarantees reuse unchanged previously audited controller mathematics. Diagnostics are not financial prices or certified expected costs.'}
    path=Path(args.output);assert not path.exists();path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,indent=2)+'\n')
    print('Independent matched audit passed:',histories,'histories;',batches_count,'batches;',len(dominance),'all-degree comparisons.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--results',required=True);p.add_argument('--output',required=True)
    main(p.parse_args())
