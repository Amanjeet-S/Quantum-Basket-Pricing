"""Complete original-profile confidence budgets and labelled ideal diagnostics.

The main output retains all 108 method/contract profiles. Different confidence
scopes are explicit. No device, market data or physical gate runtime is used.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json
import importlib.util
import math
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from iterative_estimator import original_constructor_resources,decimal_record,estimate

def load_prior(prior):
    spec=importlib.util.spec_from_file_location('qbp_estimation_prior_for_study',prior/'estimation_prior.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(prior));spec.loader.exec_module(module)
    return module.prior_for_row


def budgets(repo=ROOT,prior=ROOT,output=ROOT/'results'):
    prior_for_row=load_prior(prior)
    output.mkdir(parents=True,exist_ok=True)
    source=repo/'results/generalisation_resources.json'
    raw=json.loads(source.read_text());models={m['model_id']:m for m in raw['models']}
    rows=[]
    for row in raw['configurations']:
        prior_record=prior_for_row(models[row['model_id']],row)
        upper=F(Decimal(prior_record['finite_mean_probability_upper_decimal']))
        profiles={}
        for name,alpha in [('one_price_99_per_cent',F(1,100)),
                           ('selected_27_simultaneous_99_per_cent',F(1,2700))]:
            improved=original_constructor_resources(row,alpha,upper)
            full=original_constructor_resources(row,alpha,F(1))
            with localcontext() as ctx:
                ctx.prec=220;ctx.rounding=ROUND_FLOOR
                alpha_lower=Decimal(alpha.numerator)/Decimal(alpha.denominator)
            classical=prior_for_row(models[row['model_id']],row,alpha_lower)
            profiles[name]={
                'prior_aware':decimal_record(improved),'full_probability_range':decimal_record(full),
                'prior_T_count_improvement_ratio':full['T_count_upper']/improved['T_count_upper'] if improved['T_count_upper'] else None,
                'classical_certified_samples':classical['classical_certified_samples'],
                'classical_selected_concentration_bound':classical['classical_selected_bound'],
                'comparison_units':'classical function evaluations versus logical T gates are different units; no runtime-speedup conclusion'}
            assert sum(improved['Clifford_T_gate_counts_upper'].get(k,0) for k in ('T','Tdg'))==improved['T_count_upper']
            assert improved['explicit_7T_sequence_all_gate_depth_serial_upper']==sum(improved['Clifford_T_gate_counts_upper'].values())
        rows.append({'case_id':row['case_id'],'model_id':row['model_id'],'method':row['method'],
            'cap':row['cap'],'q':row['q'],'finite_probability_upper':decimal_record(upper),
            'financial_prior':prior_record,'profiles':profiles})
    selected=[]
    for case in sorted(set(r['case_id'] for r in rows)):
        chosen=min((r for r in rows if r['case_id']==case),key=lambda r:r['profiles']['one_price_99_per_cent']['prior_aware']['T_count_upper'])
        selected.append({'case_id':case,'method':chosen['method'],
            'T_count_one_price_upper':chosen['profiles']['one_price_99_per_cent']['prior_aware']['T_count_upper'],
            'T_count_simultaneous_27_upper':chosen['profiles']['selected_27_simultaneous_99_per_cent']['prior_aware']['T_count_upper'],
            'selection_scope':'minimum among four original seven-T method profiles under complete prior-aware confidence budget'})
    summary={
        'original_profiles':len(rows),'selected_contracts':len(selected),
        'prior_aware_T_range_one_price':[min(r['profiles']['one_price_99_per_cent']['prior_aware']['T_count_upper'] for r in rows),max(r['profiles']['one_price_99_per_cent']['prior_aware']['T_count_upper'] for r in rows)],
        'finite_prior_range':[min(float(r['finite_probability_upper']['float']) for r in rows),max(float(r['finite_probability_upper']['float']) for r in rows)],
        'prior_query_improvement_range':[min(r['profiles']['one_price_99_per_cent']['prior_T_count_improvement_ratio'] for r in rows),max(r['profiles']['one_price_99_per_cent']['prior_T_count_improvement_ratio'] for r in rows)],
        'selected_method_counts':{method:sum(r['method']==method for r in selected) for method in sorted(set(r['method'] for r in rows))}}
    out={'status':'complete original exact-unitary seven-T profile budgets with rigorous prior-aware adaptive confidence controller',
        'configuration_sha256':sha256(source.read_bytes()).hexdigest(),
        'script_hashes':{p.name:sha256(p.read_bytes()).hexdigest() for p in [ROOT/'iterative_estimator.py',ROOT/'quantum_gates.py',ROOT/'estimator_proof.md',ROOT/'build_estimator_study.py',prior/'estimation_prior.py']},
        'confidence_scope':'each per-price profile is99 per cent; selected27 profile jointly99 per cent by union bound; all108 simultaneously would require alpha=.01/108',
        'logical_assumptions':'exact coherent finite A and inverse; independent fresh Born measurements; no physical noise or hardware runtime inference',
        'summary':summary,'configurations':rows,'selected_27':selected}
    (output/'estimator_resources.json').write_text(json.dumps(out,indent=2)+'\n')
    return summary


def simulation(output=ROOT/'results'):
    output.mkdir(parents=True,exist_ok=True)
    # Explicit synthetic p/valid-prior pairs; no financial price is inferred.
    import numpy as np
    rng=np.random.default_rng(20261005)
    points=[(0.,F(1)),(1.,F(1)),(.0001,F(1,1000)),(.01,F(1,50)),(.025,F(1,20)),
            (.25,F(1)),(.5,F(1)),(.75,F(1)),(.9999,F(1))]
    records=[]
    for cap in (1,10,25):
        for probability,prior in points:
            theta=math.asin(math.sqrt(probability))
            totals=[];covered=0;flagged=0
            for repeat in range(20):
                def backend(k,n):return int(rng.binomial(n,math.sin((2*k+1)*theta)**2))
                result=estimate(backend,F(1,200*cap),F(1,100),prior)
                lo,hi=result['amplitude_interval']
                covered+=float(lo)-1e-14<=probability<=float(hi)+1e-14
                flagged+=result['status'].startswith('confidence_failure')
                assert result['total_shots']<=result['plan'].total_shots_upper
                assert result['total_grover_calls']<=result['plan'].total_grover_calls_upper
                totals.append({'rounds':len(result['rounds']),'shots':result['total_shots'],'Q_calls':result['total_grover_calls']})
            records.append({'cap':cap,'synthetic_probability':probability,'proved_prior_upper':decimal_record(prior),
                'repetitions':20,'covered':int(covered),'flagged_failures':int(flagged),
                'shots_range':[min(t['shots'] for t in totals),max(t['shots'] for t in totals)],
                'Q_calls_range':[min(t['Q_calls'] for t in totals),max(t['Q_calls'] for t in totals)]})
    out={'scope':'ideal-binomial software diagnostic using declared synthetic amplitudes; no quantum device, market data or empirical coverage proof',
        'seed':20261005,'failure_probability_each_run':.01,'total_runs':540,'records':records}
    (output/'estimator_ideal_simulation.json').write_text(json.dumps(out,indent=2)+'\n')
    return {'diagnostic_runs':540,'covered':sum(r['covered'] for r in records),'flagged':sum(r['flagged_failures'] for r in records)}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--simulation',action='store_true')
    parser.add_argument('--repo',type=Path,default=ROOT,help='Financial source and results directory; default this project root')
    parser.add_argument('--prior',type=Path,default=ROOT,help='Directory containing estimation_prior.py; default this project root')
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    args=parser.parse_args()
    print(json.dumps(budgets(args.repo.resolve(),args.prior.resolve(),args.output.resolve()),indent=2))
    if args.simulation:print(json.dumps(simulation(args.output.resolve()),indent=2))
