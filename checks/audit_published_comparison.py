"""Independent scheduling, charging and selection audit for all comparators.

The published coefficients are conditional source premises reported by the
source-reading review. This audit verifies their mathematical postprocessing
and complete charge. It does not imply that I consulted the journal text.
"""
from collections import Counter
from decimal import Decimal,localcontext,ROUND_CEILING
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
D=Decimal
R=F(69,1000)
ETA=F(1,1 << 70)
PAD=F(1,10**60)


def ceil(x):return -((-x.numerator)//x.denominator)
def floor(x):return x.numerator//x.denominator
def unpack(x):return F(int(x['numerator']),int(x['denominator']))
def dec(x):return D(x.numerator)/D(x.denominator)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    sys.modules[name]=result;spec.loader.exec_module(result);return result


def log_bound(x,sign=1):
    with localcontext() as ctx:
        ctx.prec=80
        return F(dec(x).ln()+sign*D('1e-60'))


def upper_sqrt(x):
    n=ceil(x*(1 << 200));z=math.isqrt(n)
    return F(z if z*z==n else z+1,1 << 100)


def positive_pi_reference():
    with localcontext() as ctx:
        ctx.prec=220
        term,total,n=D('.5'),D('.5'),0
        while term>D('1e-200'):
            term *= D((2*n+1)**2)/D(8*(n+1)*(2*n+3))
            total += term;n += 1
        return F(6*total-D('1e-180')),F(6*total+D('1e-180'))


def fenwick_cost(table,initial,last,a_cost,q_cost):
    """Independent prefix query tree, rather than the constructor's array."""
    tree=[None]*(last+2)
    values={};parents={}
    def query(index):
        best=None
        while index>0:
            if tree[index] is not None and (best is None or tree[index]>best):
                best=tree[index]
            index-=index&(-index)
        return best
    def update(index,value):
        while index<=last+1:
            if tree[index] is None or value>tree[index]:tree[index]=value
            index+=index&(-index)
    for m,n in table:
        local=n*(a_cost+q_cost*((m-1)//2))
        if m==initial:previous=None;value=local
        else:
            previous=query(m//3)
            if previous is None:continue
            value=previous[0]+local
        values[m]=value;parents[m]=None if previous is None else previous[1]
        update(m,(value,m))
    best=query(last)
    chain=[];cursor=None if best is None else best[1]
    while cursor is not None:chain.append(cursor);cursor=parents[cursor]
    chain.reverse()
    return (0 if best is None else best[0]),chain


@lru_cache(maxsize=2000)
def own_plan(e,alpha,u):
    guard=min(F(1),u+F(1,10**18))
    var=guard*(1-guard) if guard<=F(1,2) else F(1,4)
    gamma=2*upper_sqrt(var)
    if u<=2*e:
        return {'epsilon':e,'alpha':alpha,'prior_upper':u,
            'guarded_prior_upper':guard,'gamma_upper':gamma,
            'angle_tolerance':F(0),'normalised_angle_width_target':F(0),
            'multiplier_upper':F(0),'maximum_odd_multiplier':0,
            'initial_multiplier':1,'maximum_rounds':0,'maximum_shots':0,
            'maximum_Q_calls':0},()
    angle=(e-PAD)/gamma;target=4*angle/PI_UP;bound=1/target
    last=floor(bound);last-=1 if last%2==0 else 0
    upper=SCALAR.inverse_probability_interval(u)[1]
    assert SCALAR.probability_interval(upper)[1]<=guard
    initial=floor(1/upper);initial-=1 if initial%2==0 else 0
    initial=max(1,initial)
    assert initial*upper<=1 and (initial+2)*upper>1 and initial<=last
    rounds=1;v=initial
    while 3*v<=bound:v*=3;rounds+=1
    table=[]
    with localcontext() as ctx:
        ctx.prec=120
        for m in range(initial,last+1,2):
            if m!=initial and m<3*initial:continue
            failure=F(2,3)*alpha*m/bound
            n=ceil(log_bound(2/failure)/(2*R*R))
            # Independently check the actual integer satisfies Hoeffding.
            independent_failure=2*(-2*dec(R*R)*D(n)).exp()+D('1e-90')
            assert independent_failure<=dec(failure)
            table.append((m,n))
    shots,_=fenwick_cost(table,initial,last,1,0)
    calls,_=fenwick_cost(table,initial,last,0,1)
    values={'epsilon':e,'alpha':alpha,'prior_upper':u,
        'guarded_prior_upper':guard,'gamma_upper':gamma,
        'angle_tolerance':angle,'normalised_angle_width_target':target,
        'multiplier_upper':bound,'maximum_odd_multiplier':last,
        'initial_multiplier':initial,'maximum_rounds':rounds,
        'maximum_shots':shots,'maximum_Q_calls':calls}
    return values,tuple(table)


@lru_cache(maxsize=2000)
def source_caps(e,alpha,u):
    if u<=2*e:
        empty={'total_shots_upper':0,'total_grover_calls_upper':0,
               'grover_power_upper':0,'maximum_rounds':0}
        return {'iqae_published':empty,
            'miqae_published_continuous_shot_term':empty,
            'miqae_full_batch_integer_corollary':empty},F(0),F(0)
    variance=u*(1-u) if u<=F(1,2) else F(1,4)
    gamma=2*upper_sqrt(variance);angle=e/gamma
    assert angle<F(1,4)
    bound=PI_UP/(4*angle)
    log2up=log_bound(bound)/log_bound(F(2),-1)
    logarithm=log_bound(2*log2up/alpha)
    sine=upper_sqrt(SCALAR.probability_interval(F(1,7))[1])
    nmax=32*logarithm/(1-2*sine)**2
    rounds=1;m=F(1)
    while 2*m<=bound:m*=2;rounds+=1
    iqae={'total_shots_upper':rounds*(ceil(nmax)+100),
        'total_grover_calls_upper':ceil(50*logarithm/angle),
        'grover_power_upper':ceil(PI_UP/(8*angle)),
        'maximum_rounds':rounds,
        'published_maximum_rounds':max(1,ceil(log_bound(PI_UP/(8*angle))/log_bound(F(2),-1))),
        'base_shots':100,'maximum_shots_per_round_real_upper':nmax,
        'published_reported_Q_term_upper':ceil(50*logarithm/angle),
        'independent_structural_Q_cap_upper':ceil((nmax+100)*bound)}
    assert 32*log_bound((2/alpha)*(log_bound(SCALAR.pi_rational_interval()[0]/(4*angle),-1)/log_bound(F(2))),-1)>=100
    rounds=1;m=1
    while 3*m<=bound:m*=3;rounds+=1
    delta=SCALAR.probability_interval(F(3,7))[0]-F(1,4)
    initial=2/(delta*delta)*log_bound(3*bound/alpha)
    assert initial<F(1000000) and delta>F(13,100)
    # The declared numerical inverse guard is far below one extra shot's
    # conservative guaranteed span decrease in the retained domain.
    assert delta/(2*(F(1000000)+2))>F(6,10**8)
    assert 8*ETA<F(6,10**8)
    qterm=ceil(62*log_bound(6/alpha)/angle)
    real={'total_shots_upper':ceil(rounds*initial),
        'total_grover_calls_upper':qterm,
        'grover_power_upper':ceil((bound-1)/2),
        'maximum_rounds_including_k0':rounds,
        'initial_maximum_real_shots_upper':initial}
    integer={**real,'total_shots_upper':rounds*ceil(initial+3),
        'total_grover_calls_upper':qterm+ceil(9*bound/4),
        'integer_padding_extra_queries_upper':ceil(9*bound/4)}
    return {'iqae_published':iqae,
        'miqae_published_continuous_shot_term':real,
        'miqae_full_batch_integer_corollary':integer},gamma,angle


def charge_check(record,row,shots,queries,power):
    r=row['resources'];calls=shots+2*queries
    reflection=r['grover_zero_reflection'];ccx_a=r['expanded_gate_counts'].get('CCX',0)
    ccx_q=2*ccx_a+reflection['gate_counts'].get('CCX',0)
    ccx=calls*ccx_a+queries*reflection['gate_counts'].get('CCX',0)
    td_a=r['logical_T_depth_serial_upper'];td_q=r['grover_iterate_T_depth_serial_upper']
    a_depth=r['serial_all_gate_depth_upper']+1
    q_depth=2*a_depth+sum(reflection['gate_counts'].values())+1
    expected={'T_count_upper':shots*r['logical_T_count']+queries*r['grover_iterate_T_count'],
        'clean_A_and_adjoint_invocations_upper':calls,'grover_invocations_upper':queries,
        'flag_measurements_upper':shots,'lookup_measurements_upper':calls*r['measurement_count'],
        'all_measurements_upper':calls*r['measurement_count']+shots,
        'conditional_CZ_feedback_upper':calls*r['classical_feedback_CZ_count'],
        'conditional_lookup_resets_upper':calls*r['classical_reset_X_count'],
        'fresh_shot_register_resets_upper':shots,
        'legacy_assigned_three_layer_serial_T_depth':shots*td_a+queries*td_q,
        'expanded_five_layer_CCX_serial_T_depth_upper':shots*td_a+queries*td_q+2*ccx,
        'legacy_assigned_three_layer_maximum_shot_T_depth':td_a+power*td_q if shots else 0,
        'expanded_five_layer_CCX_maximum_shot_T_depth_upper':td_a+power*td_q+2*(ccx_a+power*ccx_q) if shots else 0,
        'serial_lookup_measurement_depth_upper':calls*r['measurement_depth_serial_upper'],
        'logical_macro_serial_dependency_depth_upper':shots*a_depth+queries*q_depth,
        'expanded_CCX_serial_all_gate_dependency_depth_upper':shots*a_depth+queries*q_depth+14*ccx,
        'logical_macro_maximum_shot_dependency_depth_upper':a_depth+power*q_depth if shots else 0,
        'peak_logical_qubits_upper':r['grover_iterate_logical_qubits_upper'] if shots else 0,
        'uniform_Hadamards_upper':calls*r['includes_finite_uniform_preparation_H_gates'],
        'rotation_synthesis_bias':0}
    for name,value in expected.items():assert record[name]==value,(row['candidate_id'],name,record[name],value)
    gates=Counter({k:calls*v for k,v in r['expanded_gate_counts'].items()})
    for k,v in reflection['gate_counts'].items():gates[{'CNOT':'CX'}.get(k,k)]+=queries*v
    gates['Z']+=queries;gates['H']+=calls*r['includes_finite_uniform_preparation_H_gates']
    n=gates.pop('CCX',0);assert n==ccx
    gates['CX']+=6*n;gates['H']+=2*n;gates['T']+=4*n;gates['TD']+=3*n
    assert dict(gates)==record['assigned_Clifford_T_and_instrument_operation_counts_upper']
    assert gates['T']+gates['TD']==record['T_count_upper']
    assert gates['MZ']==record['lookup_measurements_upper']


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo',type=Path,default=ROOT)
    p.add_argument('--comparator-root',type=Path,default=ROOT)
    p.add_argument('--config',type=Path)
    p.add_argument('--comparison',type=Path)
    p.add_argument('--output',type=Path,default=ROOT/'results/published_comparison_audit.json')
    a=p.parse_args();repo=a.repo.resolve();comp=a.comparator_root.resolve()
    config=a.config or repo/'results/architecture_resource_candidates.json'
    comparison=a.comparison or comp/'results/published_comparison.json'
    data=json.loads(config.read_text());retained=json.loads(comparison.read_text())
    for name,value in retained['source_hashes'].items():assert digest(repo/name)==value,(name,'financial source changed')
    for name,value in retained['own_source_hashes'].items():assert digest(comp/name)==value,(name,'comparator source changed')
    global SCALAR,PI_UP
    SCALAR=module(comp/'certified_scalar.py','independent_published_scalar')
    PI_UP=SCALAR.pi_rational_interval()[1]
    low,high=positive_pi_reference()
    assert SCALAR.pi_rational_interval()[0]<=low<=high<=PI_UP
    # Critical spans are reconstructed from the retained scalar enclosure
    # premise rather than invoking the comparator cover constructor.
    critical={s:[SCALAR.probability_interval(F(j,s)) for j in range(1,s)] for s in (3,5,7)}
    spans=[max(x[0],y[0],z[0])-min(x[1],y[1],z[1])
           for x in critical[3] for y in critical[5] for z in critical[7]]
    margin=min(spans)-2*R-8*ETA
    assert len(spans)==48 and margin>0
    sys.path.insert(0,str(repo))
    prior_module=module(repo/'estimation_prior.py','independent_published_financial_prior')
    models={m['model_id']:m for m in data['models']}
    key=lambda row:(row['architecture'],row['candidate_id'])
    sources={key(row):row for row in data['configurations']}
    assert len(sources)==len(data['configurations'])==819
    assert {key(row) for row in retained['configurations']}==set(sources)
    assert len(retained['configurations'])==819
    assert len({r['case_id'] for r in retained['configurations']})==27
    labels={'one_price_99_per_cent':F(1,100),'selected_27_simultaneous_99_per_cent':F(1,2700)}
    started=time.perf_counter();records=[];profile_count=0
    for index,row in enumerate(retained['configurations']):
        source=sources[key(row)]
        assert all(row[k]==source[k] for k in ('case_id','model_id','method','cap'))
        assert row['prefix']=={k:source[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits')}
        assert F(str(source['price_tolerance']))==F(1,100)
        assert F.from_float(source['reserved_statistical_price_allowance'])>=F(1,200)
        assert F(source['independent_deterministic_bias_upper_decimal'])<=F(1,200)
        prior=F(prior_module.prior_for_row(models[row['model_id']],source)['finite_mean_probability_upper_decimal'])
        assert unpack(row['proved_finite_probability_upper'])==prior
        assert set(row['confidence_profiles'])==set(labels)
        e=F(1,200*row['cap'])
        for scope,alpha in labels.items():
            profiles=row['confidence_profiles'][scope]
            assert set(profiles)=={'financial_prior','full_range'}
            for label,u in [('financial_prior',prior),('full_range',F(1))]:
                profile=profiles[label];own=profile['certified_full_batch_variant']
                values,table=own_plan(e,alpha,u)
                for name,value in values.items():
                    actual=unpack(own['plan'][name]) if isinstance(value,F) else own['plan'][name]
                    assert actual==value,(row['candidate_id'],scope,label,name)
                assert own['plan']['shot_table_entries']==len(table)
                assert own['plan']['shot_table_sha256']==hashlib.sha256(repr(table).encode()).hexdigest()
                power=max(0,(values['maximum_odd_multiplier']-1)//2)
                charge_check(own,source,values['maximum_shots'],values['maximum_Q_calls'],power)
                coupled=own['coupled_linear_T_cost_upper']
                cost,_=fenwick_cost(table,values['initial_multiplier'],values['maximum_odd_multiplier'],
                                  source['resources']['logical_T_count'],source['resources']['grover_iterate_T_count'])
                assert coupled['cost_upper']==cost<=own['T_count_upper']
                witness=coupled['witness_multipliers'];shot_table=dict(table)
                assert not witness or witness[0]==values['initial_multiplier']
                assert all(b>=3*x for x,b in zip(witness,witness[1:]))
                assert not witness or witness[-1]<=values['maximum_odd_multiplier']
                assert len(witness)<=values['maximum_rounds']
                assert coupled['witness_shots']==sum(shot_table[m] for m in witness)
                assert coupled['witness_Q_calls']==sum(shot_table[m]*((m-1)//2) for m in witness)
                assert cost==coupled['witness_shots']*source['resources']['logical_T_count']+coupled['witness_Q_calls']*source['resources']['grover_iterate_T_count']
                if witness:assert sum(F(2,3)*alpha*m/values['multiplier_upper'] for m in witness)<=alpha
                envelopes,gamma,angle=source_caps(e,alpha,u)
                post=profile['published_financial_postprocessing']
                assert post['zero_shots']==(u<=2*e)
                assert unpack(post['gamma'])==gamma and unpack(post['epsilon_angle'])==angle
                charged=profile['published_theorem_envelopes']
                assert set(charged)==set(envelopes)
                for name,caps in envelopes.items():
                    for field,value in caps.items():
                        actual=unpack(charged[name]['plan'][field]) if isinstance(value,F) else charged[name]['plan'][field]
                        assert actual==value,(row['candidate_id'],scope,label,name,field)
                    charge_check(charged[name],source,caps['total_shots_upper'],
                                 caps['total_grover_calls_upper'],caps['grover_power_upper'])
                    profile_count+=1
                profile_count+=1
        records.append({'architecture':row['architecture'],'candidate_id':row['candidate_id'],
            'case_id':row['case_id'],'method':row['method'],'matched_financial_target_prior_and_alpha':True,
            'independent_Fenwick_chain_caps':True,'all_preparation_query_and_instrument_charges':True})
        if (index+1)%200==0:print('Independent published candidate audit',index+1,flush=True)
    selected=retained['selected_cases']
    assert len(selected)==81 and len({(r['architecture'],r['case_id']) for r in selected})==81
    for selection in selected:
        pool=[r for r in retained['configurations'] if r['architecture']==selection['architecture'] and r['case_id']==selection['case_id']]
        expected_labels=set()
        for scope in labels:
            for comparator in ('certified_full_batch_variant','iqae_published','miqae_full_batch_integer_corollary'):
                profile_label=scope+'__'+comparator;expected_labels.add(profile_label)
                def value(r):
                    block=r['confidence_profiles'][scope]['financial_prior']
                    return block[comparator]['T_count_upper'] if comparator=='certified_full_batch_variant' else block['published_theorem_envelopes'][comparator]['T_count_upper']
                best=min(value(r) for r in pool);wins=[r for r in pool if value(r)==best]
                retained_selection=selection['profiles'][profile_label]
                assert retained_selection['minimum_T_count_upper']==best
                assert retained_selection['tied_methods']==sorted({r['method'] for r in wins})
                assert retained_selection['representative_prefix'] in [r['prefix'] for r in wins]
        assert set(selection['profiles'])==expected_labels
    out={'status':'independent819-candidate scheduling, full logical charging and81selection audit passed',
        'configuration_sha256':digest(config),'comparison_sha256':digest(comparison),
        'script_sha256':digest(Path(__file__)),'candidate_count':len(records),
        'complete_resource_profile_count':profile_count,'selected_case_architecture_count':len(selected),
        'cover_critical_triples':len(spans),'cover_margin_lower_numerator':str(margin.numerator),
        'cover_margin_lower_denominator':str(margin.denominator),
        'independent_DP':'Fenwick prefix maxima with independently checked witnesses',
        'confidence_scopes':['single price99 per cent','one preselected implementation per27cases jointly99 per cent'],
        'source_scope':'published IQAE50 and MIQAE62 coefficients are conditional source theorem premises, not independently rederived or transferred to the executable modified schedule; this auditor did not consult their source text',
        'scalar_scope':'shared outward scalar enclosure machinery is checked against an independent positive-series pi reference; its earlier scalar proof and audits remain prerequisites',
        'configurations':records,'elapsed_seconds':time.perf_counter()-started}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
    print('Published comparison audit passed:',len(records),'candidates;',profile_count,'resource profiles.')


if __name__=='__main__':main()
