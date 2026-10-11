"""Independent replay of every retained randomised digital-shift output.

This audit does not import the new engine. It uses the earlier independent
integer vector interpreter, a sequential author-recurrence implementation,
Python integer moments and 300-digit monetary endpoint reconstruction.
Seeded replay verifies computation, not the statistical randomness premise.
"""
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import generalisation as gen
from audit_generalisation import vector_integer


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def key(row):return tuple(row[f] for f in ('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits'))
def build(model,row):
    return gen.build_graph(model,row['method'],*[row[f] for f in
        ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],row['cap'])


def sequential_net(net_bits,dimension,address_bits):
    """Sequential trailing-one recurrence, independent of engine Gray vector."""
    lines=(ROOT/'public_data/new-joe-kuo-6.21201').read_text().splitlines()
    directions=[]
    for dimension_index in range(dimension):
        if dimension_index==0:
            vector=[1<<(31-j) for j in range(32)]
        else:
            dd,degree,polynomial,*initial=map(int,lines[dimension_index].split())
            assert dd==dimension_index+1
            vector=[initial[j]<<(31-j) for j in range(degree)]
            for j in range(degree,32):
                value=vector[j-degree]^(vector[j-degree]>>degree)
                for k in range(1,degree):
                    if polynomial&(1<<(degree-1-k)):value^=vector[j-k]
                vector.append(value)
        directions.append(vector)
    points=np.zeros((1<<net_bits,dimension),dtype=np.uint64)
    current=[0]*dimension
    for i in range(1,1<<net_bits):
        v=i-1;column=0
        while v&1:v>>=1;column+=1
        for j in range(dimension):
            current[j]^=directions[j][column]
            points[i,j]=current[j]>>(32-address_bits[j])
    return points


def audit_selection(protocol,resources,ranges):
    originals={key(r):r for r in resources['configurations'] if r['architecture']=='measured4t'}
    certificates={key(r):r for r in ranges['configurations']}
    profiles=protocol['all_candidate_profiles']
    assert len(originals)==273 and set(originals)==set(certificates)
    expected={k for k,r in originals.items() if r['method'] in ('full_global','full_cap')}
    assert Counter(tuple(p['key']) for p in profiles)==Counter({k:1 for k in expected})
    with localcontext() as ctx:
        ctx.prec=300
        logarithm=Decimal(4400).ln()
        for p in profiles:
            k=tuple(p['key']);row=originals[k];cert=certificates[k]
            assert p['row']=={x:y for x,y in row.items() if x!='resources'}
            assert p['range_certificate']==cert
            rho=Fraction(cert['finite_lattice_range_fraction'])
            u=Fraction(cert['fair_mean_upper_fraction']);v=Fraction(cert['fair_second_moment_upper_fraction'])
            m=min(u,rho/2,v/rho) if rho else Fraction(0)
            bound=m*(rho-m);epsilon=Fraction(1,200*row['cap'])
            def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
            reps=max(1,int(((2*dec(bound)+2*dec(rho)*dec(epsilon)/3)*logarithm/dec(epsilon)**2).to_integral_value(rounding=ROUND_CEILING)))
            assert Fraction(p['batch_variance_prior_exact_fraction'])==bound
            assert p['prospective_independent_shift_count']==reps
            assert p['selection_scalar_lane_work_upper']==reps*64*p['evaluator_profile']['scalar_lane_operation_upper_total']
    for selected in protocol['selected']:
        pool=[p for p in profiles if p['row']['case_id']==selected['row']['case_id']]
        minimum=min(p['selection_scalar_lane_work_upper'] for p in pool)
        ties=[p for p in pool if p['selection_scalar_lane_work_upper']==minimum]
        choice=min(ties,key=lambda p:(p['evaluator_profile']['prepared_table_bytes'],tuple(p['key'])))
        assert choice['key']==selected['key']
        assert {tuple(p['key']) for p in ties}=={tuple(k) for k in selected['tied_prefixes']}
    return {'conditioned_candidate_prefixes':len(profiles),'selected_prices':len(protocol['selected']),
            'independent_300_digit_prospective_count_and_finite_pool_selection':True}


def audit_moments(sums,record,q):
    r=len(sums);d=64*(1<<q);total=sum(sums);squares=sum(x*x for x in sums)
    pairs=r//2;difference=sum((sums[2*j]-sums[2*j+1])**2 for j in range(pairs))
    assert record['replicates']==r and int(record['batch_denominator'])==d
    assert int(record['integer_batch_sum'])==total and int(record['integer_batch_square_sum'])==squares
    assert record['independent_pairs']==pairs and int(record['integer_paired_difference_square_sum'])==difference
    assert Fraction(record['finite_mean_exact_fraction'])==Fraction(total,r*d)
    assert Fraction(record['paired_variance_estimate_exact_fraction'])==Fraction(difference,2*pairs*d*d)
    assert Fraction(record['ordinary_replicate_sample_variance_exact_fraction'])==Fraction(r*squares-total*total,r*(r-1)*d*d)


def audit_price(model,row,result):
    interval=result['confidence_interval'];confidence=result['checkpoints'][-1]['confidence']
    with localcontext() as ctx:
        ctx.prec=300
        correction=gen.geometric_correction(model,row['cap'])
        centre=Decimal(interval['geometric_correction']['undiscounted_value_decimal'])
        assert abs(centre-Decimal(correction['undiscounted_value_decimal']))<Decimal('1e-60')
        cerror=max(Decimal('1e-12'),Decimal.from_float(correction['absolute_allowance']))
        beta=Decimal.from_float(row['arithmetic_price_error_upper'])+Decimal.from_float(row['distribution_price_bound'])+cerror
        assert beta<=Decimal('.005')
        lower=Fraction(confidence['mean_lower_fraction']);upper=Fraction(confidence['mean_upper_fraction'])
        lo=Decimal(lower.numerator)/Decimal(lower.denominator)
        hi=Decimal(upper.numerator)/Decimal(upper.denominator)
        discount=(-Decimal(model['risk_free_rate'])).exp()
        reference_lo=max(Decimal(0),centre+row['cap']*lo-beta)*discount
        reference_hi=min(Decimal(row['cap']),centre+row['cap']*hi+beta)*discount
        assert Decimal(interval['price_interval_lower_decimal'])<=reference_lo
        assert Decimal(interval['price_interval_upper_decimal'])>=reference_hi
    lo=Fraction(Decimal(interval['price_interval_lower_decimal']));hi=Fraction(Decimal(interval['price_interval_upper_decimal']))
    midpoint=Fraction(Decimal(interval['display_midpoint_decimal']));radius=Fraction(Decimal(interval['maximum_radius_decimal']))
    assert lo<=midpoint<=hi and max(midpoint-lo,hi-midpoint)<=radius
    if result['status']=='certified per-price tolerance met':assert radius<=Fraction(1,100)
    return {'independent_300_digit_monetary_endpoints_contained':True,'display_midpoint_rounding_charged':True,
            'radius_exact_fraction':str(radius)}


def main():
    started=time.perf_counter()
    protocol_path=ROOT/'results/rqmc_pricing_protocol.json';price_path=ROOT/'results/rqmc_pricing_results.json'
    protocol=json.loads(protocol_path.read_text());prices=json.loads(price_path.read_text())
    resources=json.loads((ROOT/'results/architecture_resource_candidates.json').read_text())
    ranges=json.loads((ROOT/'results/finite_range_certificates.json').read_text())
    assert prices['protocol_sha256']==sha(protocol_path)
    assert protocol['configuration_sha256']==sha(ROOT/'results/architecture_resource_candidates.json')
    assert protocol['range_certificate_sha256']==sha(ROOT/'results/finite_range_certificates.json')
    for filename,h in protocol['source_hashes'].items():assert sha(ROOT/filename)==h,filename
    assert protocol['checkpoint_replicates']==[1<<k for k in range(6,17)]
    assert Fraction(protocol['alpha_per_event_exact_fraction'])*22==Fraction(1,100)
    selection=audit_selection(protocol,resources,ranges)
    models={m['model_id']:m for m in protocol['models']}; tables={};records=[]
    for case_index,(selected,result) in enumerate(zip(protocol['selected'],prices['configurations'])):
        row=selected['row'];graph=build(models[row['model_id']],row)
        assert result['key']==selected['key'] and result['seed_sequence']==[protocol['production_seed'],case_index]
        gen.make_tables(graph,tables)
        names=list(graph.inputs);bits=list(graph.inputs.values())
        points=sequential_net(protocol['net_bits'],len(names),bits)
        rng=np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(result['seed_sequence'])))
        expected_sums=list(map(int,result['batch_integer_sums']));seen=[]
        shash=hashlib.sha256();ohash=hashlib.sha256();replay_index={}
        for r in result['independent_scalar_replay_records']:
            replay_index.setdefault((r['replicate_index'],r['point_index']),[]).append(r)
        previous=0;checked=0
        for checkpoint in result['checkpoints']:
            target=checkpoint['replicates'];assert target in protocol['checkpoint_replicates'] and previous<target
            remaining=target-previous
            while remaining:
                reps=min(remaining,(1<<20)//64)
                shifts={name:rng.integers(0,1<<b,size=reps,dtype=np.uint64) for name,b in zip(names,bits)}
                for name,a in shifts.items():shash.update(name.encode());shash.update(a.astype('<u8',copy=False).tobytes())
                # Smaller audit chunks keep the earlier interpreter's full
                # node retention separate from the new prepared evaluator.
                for rep_start in range(0,reps,4096):
                    count=min(4096,reps-rep_start)
                    inputs={name:(shifts[name][rep_start:rep_start+count,None]^points[None,:,j]).reshape(-1)
                            for j,name in enumerate(names)}
                    outputs=vector_integer(graph,inputs,tables)
                    assert np.all(outputs<=selected['range_certificate']['max_integer'])
                    ohash.update(outputs.astype('<u8',copy=False).tobytes())
                    sums=[sum(map(int,batch)) for batch in outputs.reshape(count,64)]
                    seen.extend(sums);checked+=len(outputs)
                    for local_rep in range(count):
                        absolute_rep=previous+rep_start+local_rep
                        for point in (0,1,2,3,7,15,31,63):
                            if (absolute_rep,point) not in replay_index:continue
                            lane=local_rep*64+point
                            addresses={name:int(a[lane]) for name,a in inputs.items()}
                            scalar=graph.evaluate(addresses,tables)
                            assert scalar==int(outputs[lane])
                            for r in replay_index[(absolute_rep,point)]:
                                assert r['inputs']==addresses and r['integer_output']==scalar
                previous+=reps;remaining-=reps
            assert seen==expected_sums[:target]
            audit_moments(seen,checkpoint['exact_batch_statistics'],row['q'])
            if checkpoint['confidence']['status']=='nonempty':
                radius=Fraction(checkpoint['confidence']['halfwidth_upper_fraction'])
                assert checkpoint['tolerance_met']==(radius<=Fraction(1,200*row['cap']))
            else:
                assert checkpoint['confidence']['status']=='empty_intersection'
                assert checkpoint['confidence']['halfwidth_upper_fraction'] is None
                assert not checkpoint['tolerance_met']
                assert checkpoint==result['checkpoints'][-1]
                assert result['status']=='explicit empty confidence-intersection outcome'
            if checkpoint!=result['checkpoints'][-1]:assert not checkpoint['tolerance_met']
        assert len(seen)==result['replicates'] and seen==expected_sums
        assert shash.hexdigest()==result['shift_stream_sha256']
        assert ohash.hexdigest()==result['output_stream_sha256']
        assert checked==result['financial_evaluations']==result['replicates']*64
        monetary=(audit_price(models[row['model_id']],row,result)
                  if result['confidence_interval'] is not None else {'explicit_confidence_failure_no_price_interval':True})
        records.append({'case_id':result['case_id'],'financial_outputs_replayed':checked,
            'python_integer_batch_sums_and_all_checkpoint_moments':True,
            'scalar_financial_checks':len(result['independent_scalar_replay_records']),
            'monetary_interval':monetary})
        print('Independent digital-shift replay',result['case_id'],checked,'outputs',flush=True)
    output={'status':'complete independent digital-shift replay passed',
        'protocol_sha256':sha(protocol_path),'price_sha256':sha(price_path),'audit_source_sha256':sha(__file__),
        'independent_prior_interpreter_sha256':sha(ROOT/'checks/audit_generalisation.py'),
        'selection':selection,'configurations':records,'wall_seconds':time.perf_counter()-started,
        'summary':{'financial_outputs_replayed':sum(r['financial_outputs_replayed'] for r in records),
            'python_scalar_financial_checks':sum(r['scalar_financial_checks'] for r in records)},
        'qualification':'Every retained output, batch integer sum, checkpoint moment, seed stream and monetary interval replayed; confidence proof and ideal shift independence are separate assumptions, not established by numerical replay'}
    (ROOT/'results/rqmc_pricing_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['summary'],indent=2),flush=True)


if __name__=='__main__':main()
