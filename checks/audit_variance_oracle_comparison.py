"""Independent retained-result accounting and rational classical audit.

No production comparator or quantum confidence-planner module is imported.
The audited Boolean gate emitters are reused only for fresh predicates; body
multisets, composition, reflection depth and all aggregate accounting are
reconstructed separately.
This checks recorded allocations, support/moment sharing, linear charges,
financial allowances, candidate minima and strengthened classical caps.
Quantum DP tightness remains covered by its separate scalar planner audit.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,gzip,hashlib,json
from collections import Counter
from functools import lru_cache
import value_oracle as vo
import quantum_gates as qg


def ceil_f(x):return -((-x.numerator)//x.denominator)
def decimal(x):x=F(x);return Decimal(x.numerator)/Decimal(x.denominator)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ranked(x,q):return x['T'],x['shots'],x.get('c',q)
def minimum(seq,q):return min((s for s in seq if 'T' in s),key=lambda s:ranked(s,q))



def expanded_multiset(counts):
    out=Counter()
    for kind,n in counts.items():
        if kind=='CCX':out.update({'H':2*n,'CX':6*n,'T':4*n,'Tdg':3*n})
        else:out[{'TD':'Tdg','SD':'Sdg'}.get(kind,kind)]+=n
    return out


@lru_cache(None)
def predicate(q,kind,coordinate):
    # Emit only the independently audited Boolean predicate, without a
    # production financial-body/profile factory or confidence controller.
    c=vo.ValueCircuit(store_gates=False,mode='measured4t');value=c.alloc(q+1)
    width=(coordinate+1 if coordinate<q else 0) if kind=='dyadic' else coordinate
    comparison=c.alloc(width);flag=c.alloc(1)[0]
    if kind=='dyadic':vo.bucket_predicate_into(c,value,comparison,coordinate,flag)
    elif kind=='range':vo.range_predicate_into(c,value,comparison,coordinate,flag)
    else:vo.clipped_range_predicate_into(c,value,comparison,coordinate,flag)
    return c.resources(),width,c.wire_count


@lru_cache(None)
def reflection(core):
    scratch=tuple(range(core,2*core-2));tape=qg.expand_toffoli(qg.zero_phase(tuple(range(core)),scratch))
    count=Counter(g[0] for g in tape);clocks=[0]*max(2*core-2,core)
    for gate in tape:
        kind,wires=gate[0],gate[1:]
        tick=max((clocks[w] for w in wires),default=0)+(kind in ('T','Tdg'))
        for wire in wires:clocks[wire]=tick
    return count,max(clocks,default=0),max(2*core-2,core),14*(core-2)


def independent_profile(old,q,kind,coordinate):
    body=Counter()
    for block in old['blocks']:body.update(block['primitive_gate_counts'])
    expanded=expanded_multiset(body);pm,width,wires=predicate(q,kind,coordinate)
    counts=Counter({k:2*n for k,n in expanded.items()});counts.update(pm['clifford_T_expanded_gate_counts'])
    coherent=old['input_qubits']+width;rc,rd,rq,rt=reflection(coherent+1)
    A_counts=counts+Counter({'H':coherent})
    Q_counts=Counter({k:2*n for k,n in A_counts.items()})+rc+Counter({'Z':1})
    T=A_counts['T']+A_counts['Tdg'];qt=Q_counts['T']+Q_counts['Tdg']
    assert qt==2*T+rt
    tphase=body['T']+body['TD'];ccx=body['CCX']
    four=2*(tphase+4*ccx)+pm['explicit_7T_sequence_T_depth']
    five=2*(tphase+5*ccx)+pm['conservative_5layer_CCX_T_depth_upper']
    A={'T':T,'measurements':counts['MZ'],'feedback_CZ':counts['CZ_IF'],'reset_X':counts['X_IF'],
       'explicit_T_depth_serial_upper':four,'conservative_5layer_T_depth_serial_upper':five,
       'expanded_primitive_operations':sum(A_counts.values()),'uniform_H':coherent}
    Q={'T':qt,'measurements':2*counts['MZ'],'feedback_CZ':2*counts['CZ_IF'],'reset_X':2*counts['X_IF'],
       'explicit_T_depth_serial_upper':2*four+rd,'conservative_5layer_T_depth_serial_upper':2*five+5*2*(coherent-1),
       'expanded_primitive_operations':sum(Q_counts.values()),'uniform_H':2*coherent}
    A.update({'gate_'+k:n for k,n in A_counts.items()});Q.update({'gate_'+k:n for k,n in Q_counts.items()})
    body_peak=old['logical_qubits_serial_upper_bound']-(q+2)
    private=wires-(q+1+width+1);peak=body_peak+width+1+private
    return A,Q,max(peak,rq)


def verify_charges(record,profile):
    A,Q,qubits=profile;n,k=record['shots'],record['Q_calls']
    expected={name:n*A.get(name,0)+k*Q.get(name,0) for name in set(A)|set(Q)}
    expected['T']=record['T'];expected['final_flag_measurements']=n;expected['complete_reinitialisations']=n
    if 'reinitialisation_additional_Z_measurements_upper' in record['charges']:
        z=n*A['uniform_H'];expected.update(reinitialisation_additional_Z_measurements_upper=z,
            reinitialisation_conditional_X_upper=z+n,all_measurements_including_final_and_reinitialisation_upper=expected['measurements']+n+z)
    assert record['charges']==expected
    if n:assert record['qubits']==qubits
    else:assert record['qubits']==0
    return A,Q

def audit(path,source_root=None):
    root=source_root or path.parents[1];container_bytes=path.read_bytes();raw_bytes=gzip.decompress(container_bytes) if path.suffix=='.gz' else container_bytes;data=json.loads(raw_bytes);counts={}
    if path.suffix=='.gz':
        manifest=json.loads((root/'results/variance_oracle_archive.json').read_text())
        assert manifest['compressed_sha256']==hashlib.sha256(container_bytes).hexdigest()
        assert manifest['uncompressed_sha256']==hashlib.sha256(raw_bytes).hexdigest()
        assert manifest['compressed_bytes']==len(container_bytes) and manifest['uncompressed_bytes']==len(raw_bytes)
        assert manifest['driver_sha256']==digest(root/'build_variance_archive.py')
    partial='configurations' not in data
    if partial:
        assert data['source_sha256']==digest(root/'variance_oracle_comparison.py')
        data={'configurations':data['rows'],'source_hashes':{'variance_oracle_comparison.py':data['source_sha256']},'input_hashes':{}}
    def hit(name,n=1):counts[name]=counts.get(name,0)+n
    with localcontext() as ctx:
        ctx.prec=300;ln_ref=F(Decimal(200).ln()+Decimal('1e-220'))
    for n,h in data['source_hashes'].items():assert digest(root/n)==h;hit('source_hashes')
    for n,h in data['input_hashes'].items():assert digest(root/n)==h;hit('input_hashes')
    original=json.loads((root/'results/architecture_resource_candidates.json').read_text())['configurations']
    original={(r['candidate_id'],r['architecture']):r['resources'] for r in original}
    for row in data['configurations']:
        old=original[(row['candidate_id'],row['architecture'])]
        q=row['q'];e=F(row['normalised_estimation_budget']);rho=F(row['range']);u=F(row['mean_upper']);v=F(row['second_moment_upper'])
        assert rho==F(row['max_integer'],1 << q)
        assert 0<=u<=rho<=1 and 0<=v<=rho*u
        assert F(row['deterministic_price_error_upper'])<=F(1,200)
        assert row['alpha']=='1/100';hit('shared_moment_support_rows')
        cb=row['classical_range_budget'];m=min(u,rho/2,v/rho) if rho else F(0)
        variance=m*(rho-m)
        assert F(cb['variance_upper'])==variance
        assert variance<=min(v,rho*rho/4)
        if rho:
            probability=m/rho
            assert 0<=probability<=1 and probability*rho<=u and probability*rho*rho<=v
            assert probability*rho*rho-(probability*rho)**2==variance
        b_ref=ceil_f((2*variance+2*rho*e/3)*ln_ref/(e*e))
        h_ref=ceil_f(rho*rho*ln_ref/(2*e*e))
        assert b_ref<=cb['Bernstein_sample_cap']<=b_ref+1
        assert h_ref<=cb['Hoeffding_sample_cap']<=h_ref+1
        assert cb['sample_cap']==min(cb['Bernstein_sample_cap'],cb['Hoeffding_sample_cap'])
        with localcontext() as ctx:
            ctx.prec=300;n=cb['sample_cap']
            den=2*variance+2*rho*e/3
            if n and den:
                if cb['Bernstein_sample_cap']<=cb['Hoeffding_sample_cap']:
                    failure=2*(-decimal(n*e*e/den)).exp()+Decimal('1e-220')
                else:failure=2*(-decimal(2*n*e*e/(rho*rho))).exp()+Decimal('1e-220')
                assert failure<=Decimal('.01')
        hit('classical_certificates')
        for controller,record in row['controllers'].items():
            scalar_candidates=[record['plain']]+record['range_candidates']+record['clipping_candidates']
            for c in scalar_candidates:
                if 'T' not in c:continue
                level=F(c['range']);tail=F(c['tail_upper']);epsilon=F(c['epsilon']);prior=F(c['prior'])
                assert level==F(1 << c['c'],1 << q)
                if c['kind']=='range':assert row['max_integer']<=1 << c['c'] and tail==0
                assert 0<=prior<=1 and tail>=0 and level*epsilon+tail/2<=e
                assert F(c['total_mean_error_upper'])==level*epsilon+tail/2
                assert c['moment_bounds_shared_with_all_opponents']
                assert c['charges']['T']==c['T']
                assert c['charges']['final_flag_measurements']==c['shots']
                assert c['charges']['complete_reinitialisations']==c['shots']
                if c['status']=='estimated':
                    assert c['scalar_plan_input']=={'epsilon':c['epsilon'],'alpha':'1/100','prior':c['prior']}
                    product=c['shots']*c['A_T']+c['Q_calls']*c['Q_T']
                    assert c['T']<=product
                    if controller=='three_five':assert c['T']==product
                else:assert c['T']==c['shots']==c['Q_calls']==0 and prior<=2*epsilon
                A,Q=verify_charges(c,independent_profile(old,q,c['kind'],c['c']))
                if c['status']=='estimated':assert c['A_T']==A['T'] and c['Q_T']==Q['T']
                hit('all_scalar_metric_charges');hit('scalar_candidates')
            for d in record['dyadic_candidates']:
                if 'T' not in d:continue
                assert F(d['total_mean_error_upper'])<=e
                assert F(d['second_moment_upper'])==v and d['max_integer']==row['max_integer']
                alpha_sum=F(0);error=F(0);T=shots=queries=peak=reset_z=0;allcharges={}
                for b in d['buckets']:
                    i=b['index'];weight=F(b['weight']);prior=F(b['prior_upper']);eps=F(b['probability_tolerance']);err=F(b['mean_error_allowance'])
                    assert weight==(F(1) if i==q else F(1 << (i+1),1 << q))
                    assert err==weight*eps and err<=weight*prior/2
                    if i==q and row['max_integer']<1 << q:assert prior==0
                    if i<q and row['max_integer']<1 << i:assert prior==0
                    if b['status']=='estimated':
                        assert eps>=F(1,10**6) and prior>2*eps
                        alpha_sum+=F(b['alpha'])
                    else:assert b['maximum_shots']==b['maximum_Q_calls']==0
                    A,Q,qubits=independent_profile(old,q,'dyadic',i)
                    n,k=b['maximum_shots'],b['maximum_Q_calls']
                    if n:peak=max(peak,qubits)
                    reset_z+=n*A['uniform_H']
                    for name in set(A)|set(Q):allcharges[name]=allcharges.get(name,0)+n*A.get(name,0)+k*Q.get(name,0)
                    error+=err;T+=b['metric_caps'].get('T',0)
                    shots+=b['maximum_shots'];queries+=b['maximum_Q_calls'];hit('dyadic_buckets')
                assert alpha_sum<=F(1,100)
                if d['omitted_mean_interval'] is not None:
                    a,b=map(F,d['omitted_mean_interval']);assert a==0 and b>=0
                    assert error+b/2==F(d['total_mean_error_upper'])
                assert T==d['T']==d['charges'].get('T',0)
                allcharges['T']=T;allcharges['final_flag_measurements']=shots;allcharges['complete_reinitialisations']=shots
                if 'reinitialisation_additional_Z_measurements_upper' in d['charges']:
                    allcharges.update(reinitialisation_additional_Z_measurements_upper=reset_z,
                        reinitialisation_conditional_X_upper=reset_z+shots,all_measurements_including_final_and_reinitialisation_upper=allcharges['measurements']+shots+reset_z)
                assert allcharges==d['charges']
                assert peak==d['qubits'];hit('all_dyadic_metric_charges')
                assert shots==d['shots'] and queries==d['Q_calls']
                assert d['charges'].get('final_flag_measurements',0)==shots
                assert d['charges'].get('complete_reinitialisations',0)==shots
                hit('dyadic_candidates')
            for name,peers in (
                ('best_range',record['range_candidates']),('best_clipped',record['clipping_candidates']),
                ('best_dyadic',record['dyadic_candidates']),
                ('best_single_flag',[record['plain'],record['best_range'],record['best_clipped']]),
                ('best_all',[record['best_single_flag'],record['best_dyadic']])):
                assert record[name]==minimum(peers,q);hit('profile_minima')
    if partial:
        return {'status':'passed partial checkpoint','input_path':str(path),'input_sha256':hashlib.sha256(container_bytes).hexdigest(),'uncompressed_sha256':hashlib.sha256(raw_bytes).hexdigest(),
                'audit_source_sha256':digest(Path(__file__)),'rows':len(data['configurations']),'checks':counts,
                'scope':'partial stored accounting/minima, support/moment/error transfer and independent 300-digit classical confidence; final all-case summary not yet checked'}
    selected=[]
    for case in sorted({r['case_id'] for r in data['configurations']}):
        group=[r for r in data['configurations'] if r['case_id']==case]
        for family in ('plain','best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
            pool=[(r,c,r['controllers'][c][family]) for r in group for c in ('three_five','miqae')]
            r,c,s=min(pool,key=lambda z:(z[2]['T'],z[2]['shots'],z[0]['candidate_id'],z[1]))
            selected.append({'case_id':case,'family':family,'candidate_id':r['candidate_id'],
                'method':r['method'],'architecture':r['architecture'],'controller':c,
                'T':s['T'],'shots':s['shots'],'Q_calls':s['Q_calls'],'qubits':s['qubits'],
                'c':s.get('c'),'allocation':s.get('allocation'),
                'deterministic_price_error_upper':r['deterministic_price_error_upper']})
            hit('financial_case_minima')
    assert selected==data['selected_cases']
    reinitialisations=[]
    row_index={(r['candidate_id'],r['architecture']):r for r in data['configurations']}
    for record in selected:
        row=row_index[(record['candidate_id'],record['architecture'])]
        plan=row['controllers'][record['controller']][record['family']]
        L=original[(record['candidate_id'],record['architecture'])]['input_qubits'];q=row['q']
        if plan['kind']=='dyadic':
            Z=sum(b['maximum_shots']*(L+(b['index']+1 if b['index']<q else 0)) for b in plan['buckets'])
            shots=sum(b['maximum_shots'] for b in plan['buckets'])
        else:Z=plan['shots']*(L+plan['c']);shots=plan['shots']
        reinitialisations.append({k:record[k] for k in ('case_id','family','candidate_id','architecture','controller')})
        reinitialisations[-1].update(nonflag_Z_measurements_upper=Z,conditional_Clifford_X_resets_upper=Z+shots,additional_T_count=0)
    assert len(data['configurations'])==data['summary']['profiles']==819
    assert len(selected)==data['summary']['selected_records']==162
    indexed={(r['case_id'],r['family']):r for r in selected};cases=sorted({r['case_id'] for r in selected})
    for family,report in data['summary']['family_ratios_to_plain'].items():
        ratios=[F(indexed[(c,family)]['T'],indexed[(c,'plain')]['T']) for c in cases]
        assert report=={'strict_improvements':sum(x<1 for x in ratios),'ties':sum(x==1 for x in ratios),
            'worse':sum(x>1 for x in ratios),'minimum_ratio':float(min(ratios)),'maximum_ratio':float(max(ratios))}
        hit('summary_families')
    ratios=[F(indexed[(c,'best_dyadic')]['T'],indexed[(c,'best_single_flag')]['T']) for c in cases]
    assert data['summary']['dyadic_vs_best_single_flag']=={'wins':sum(x<1 for x in ratios),
        'ties':sum(x==1 for x in ratios),'losses':sum(x>1 for x in ratios),
        'minimum_ratio':float(min(ratios)),'maximum_ratio':float(max(ratios))}
    hit('summary_families')
    return {'status':'passed','input_path':str(path),'input_sha256':hashlib.sha256(container_bytes).hexdigest(),'uncompressed_sha256':hashlib.sha256(raw_bytes).hexdigest(),
        'audit_source_sha256':digest(Path(__file__)),'checks':counts,
        'gate_emitter_source_hashes':{n:digest(root/n) for n in ('value_oracle.py','quantum_gates.py','architecture.py','reversible.py')},
        'postshot_reinitialisation_protocol_upper_caps':reinitialisations,
        'scope':'exact stored accounting/minima, support/moment/error transfer and independent 300-digit classical confidence; no production comparator/planner imports; separate scalar audit remains required for DP tightness'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path);parser.add_argument('--source-root',type=Path)
    args=parser.parse_args();out=audit(args.input,args.source_root)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','input_sha256','uncompressed_sha256','audit_source_sha256','checks')},indent=2))
