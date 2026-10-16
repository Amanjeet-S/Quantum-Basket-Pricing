"""Independent replay and monetary-interval audit of complete finite prices.

The replay uses the earlier independent integer interpreter and Python object
moment reductions. It does not infer statistical coverage from a seeded run
or from comparison with an uncertified continuous numerical reference.
"""
from collections import Counter
from decimal import Decimal,localcontext,ROUND_CEILING,ROUND_HALF_EVEN
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import generalisation as gen
from classical_engine import FiniteEvaluator
from audit_generalisation import vector_integer

D = Decimal
METHODS = ('direct','interval','full_global','full_cap')
FIELDS = ('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits')


def key(row):return tuple(row[f] for f in FIELDS)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(model,row):
    return gen.build_graph(model,row['method'],*[row[f] for f in
        ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],row['cap'])


def independent_budget(certificate):
    with localcontext() as ctx:
        ctx.prec = 300
        v = D(certificate['finite_output_variance_upper_decimal'])
        t = D('.005')/D(certificate['cap']);alpha = D('.01');pad = D('1e-260')
        assert 0 <= v <= D('.25')
        log = (2/alpha).ln()+pad
        n_b = int(((2*v+2*t/3)*log/(t*t)+pad).to_integral_value(rounding=ROUND_CEILING))
        n_h = int((log/(2*t*t)+pad).to_integral_value(rounding=ROUND_CEILING))
        record = certificate['confidence_budgets']['single_price']
        assert D(record['failure_allocation_decimal']) == alpha
        assert record['Bernstein_samples'] == n_b and record['Hoeffding_samples'] == n_h
        assert record['certified_samples'] == min(n_b,n_h)
        expected = 'Bernstein' if n_b <= n_h else 'Hoeffding'
        assert record['selected_bound'] == expected
        n = record['certified_samples']
        failure = 2*(-D(n)*t*t/(2*v+2*t/3)).exp() if expected == 'Bernstein' else 2*(-2*D(n)*t*t).exp()
        emitted = D(record['certified_failure_upper_decimal'])
        assert failure+pad <= emitted <= alpha
        return n


def interval_audit(model,row,result):
    interval = result['confidence_interval'];moments = result['exact_moments']
    mean = Fraction(int(moments['integer_sum']),moments['count']*(1 << row['q']))
    assert Fraction(D(interval['normalised_sample_mean_lower_decimal'])) <= mean
    assert mean <= Fraction(D(interval['normalised_sample_mean_upper_decimal']))
    with localcontext() as ctx:
        ctx.prec = 220;ctx.rounding = ROUND_CEILING
        producer_correction = gen.geometric_correction(model,row['cap'])
    assert interval['geometric_correction'] == producer_correction
    with localcontext() as ctx:
        ctx.prec = 220;ctx.rounding = ROUND_HALF_EVEN
        default_correction = gen.geometric_correction(model,row['cap'])
    correction = interval['geometric_correction']
    assert D(interval['statistical_currency_allowance_exact']) == D('.005')
    c_error = max(D('1e-12'),D.from_float(correction['absolute_allowance']))
    with localcontext() as ctx:
        ctx.prec = 300
        a = D.from_float(row['arithmetic_price_error_upper'])
        b = D.from_float(row['distribution_price_bound'])
        beta = a+b+c_error
        assert beta <= D('.005')
        assert D(interval['deterministic_currency_allowance_decimal']) >= beta
        centre = D(correction['undiscounted_value_decimal'])
        centre_difference = abs(centre-D(default_correction['undiscounted_value_decimal']))
        assert centre_difference < D('1e-60')
        exact_mean = D(mean.numerator)/D(mean.denominator)
        low = max(D(0),centre+row['cap']*exact_mean-beta-D('.005'))
        high = min(D(row['cap']),centre+row['cap']*exact_mean+beta+D('.005'))
        discount = (-D(model['risk_free_rate'])).exp()
        reference_error = D('1e-250')
        p_low,p_high = D(interval['price_interval_lower_decimal']),D(interval['price_interval_upper_decimal'])
        assert p_low <= low*discount-reference_error
        assert p_high >= high*discount+reference_error
    # Rational comparisons avoid treating a displayed Decimal midpoint as
    # an exact unrounded average of the displayed endpoint decimals.
    lo = Fraction(D(interval['price_interval_lower_decimal']))
    hi = Fraction(D(interval['price_interval_upper_decimal']))
    midpoint = Fraction(D(interval['display_midpoint_decimal']))
    radius = Fraction(D(interval['maximum_radius_decimal']))
    assert lo <= midpoint <= hi
    assert max(midpoint-lo,hi-midpoint) <= radius <= Fraction(1,100)
    return {'outward_mean_enclosure':True,'producer_rounding_geometric_centre_reproduced':True,
        'default_rounding_geometric_centre_difference_decimal':str(centre_difference),
        'arithmetic_distribution_and_correction_budget_at_most_0_005':True,
        'independent_300_digit_endpoint_formula_contained':True,
        'display_midpoint_error_accounted_in_reported_radius':True,
        'radius_exact_fraction':str(radius)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--protocol',type=Path,default=ROOT/'results/classical_pricing_protocol.json')
    p.add_argument('--prices',type=Path,default=ROOT/'results/classical_pricing_results.json')
    p.add_argument('--resources',type=Path,default=ROOT/'results/architecture_resource_candidates.json')
    p.add_argument('--variance',type=Path,default=ROOT/'results/classical_variance_certificates.json')
    p.add_argument('--source-dir',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,default=ROOT/'results/classical_pricing_audit.json')
    args = p.parse_args();started = time.perf_counter()
    protocol = json.loads(args.protocol.read_text());prices = json.loads(args.prices.read_text())
    resources = json.loads(args.resources.read_text());variance = json.loads(args.variance.read_text())
    assert prices['protocol_sha256'] == sha(args.protocol)
    assert protocol['configuration_sha256'] == sha(args.resources)
    assert protocol['variance_certificate_sha256'] == sha(args.variance)
    assert protocol['source_hashes'] == prices['source_hashes']
    for f,h in protocol['source_hashes'].items():assert sha(args.source_dir/f) == h,f
    for f in ('generalisation.py','general_kernels.py','encoding_resources.py','classical_engine.py'):
        assert sha(ROOT/f) == protocol['source_hashes'][f],('independent staged dependency changed',f)
    assert protocol['numpy_version'] == np.__version__
    models = {m['model_id']:m for m in protocol['models']}
    assert models == {m['model_id']:m for m in resources['models']}
    originals = {key(r):r for r in resources['configurations']}
    certificates = {key(r):r for r in variance['configurations']}
    profiles = protocol['all_candidate_profiles'];selected = protocol['selected']
    assert len(profiles) == len(originals) == len(certificates) == 273
    assert Counter(tuple(x['key']) for x in profiles) == Counter({k:1 for k in originals})
    cases = sorted({r['case_id'] for r in originals.values()})
    assert len(cases) == 27 and len(selected) == len(prices['configurations']) == 108
    expected_cases = Counter((c,m) for c in cases for m in METHODS)
    assert Counter((x['row']['case_id'],x['row']['method']) for x in selected) == expected_cases
    assert Counter((x['case_id'],x['method']) for x in prices['configurations']) == expected_cases
    assert 'per price' in protocol['confidence_scope'] and 'not all108' in protocol['confidence_scope']
    bank = {};table_time = preparation_time = 0;profile_checks = []
    table_records = {(x['kernel'],x['address_bits'],x['q']):x for x in protocol['first_generation_table_records']}
    for profile in profiles:
        row = profile['row'];k = key(row)
        assert tuple(profile['key']) == k
        for field in ('model_id','cap','q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits',
                      'arithmetic_price_error_upper','distribution_price_bound'):
            assert row[field] == originals[k][field],(k,field)
        assert profile['certificate'] == certificates[k]
        n = independent_budget(profile['certificate'])
        assert profile['certified_samples'] == n
        graph = build(models[row['model_id']],row)
        start = time.perf_counter();gen.make_tables(graph,bank);table_time += time.perf_counter()-start
        start = time.perf_counter();evaluator = FiniteEvaluator(graph,bank);preparation_time += time.perf_counter()-start
        actual = evaluator.operation_profile();stored = profile['evaluator_profile']
        for field in ('q','width','inputs','nodes_by_kind','scalar_lane_operations_upper','scalar_lane_operation_upper_total',
                      'table_manifest','prepared_table_bytes','interpolation_product_branches','constant_product_branches'):
            assert stored[field] == actual[field],(k,field)
        assert profile['sample_evaluator_scalar_lane_work_upper'] == n*actual['scalar_lane_operation_upper_total']
        profile_checks.append({'key':list(k),'budget_and_operation_profile_equal':True})
        if len(profile_checks)%25 == 0:print('Independently checked classical candidate profiles',len(profile_checks),flush=True)
    assert set(bank) == set(table_records)
    for k,raw in bank.items():
        assert len(raw) == table_records[k]['entries']
        assert hashlib.sha256(json.dumps(raw,separators=(',',':')).encode()).hexdigest() == table_records[k]['sha256']
    for winner in selected:
        row = winner['row']
        pool = [x for x in profiles if (x['row']['case_id'],x['row']['method']) == (row['case_id'],row['method'])]
        objective = min(x['sample_evaluator_scalar_lane_work_upper'] for x in pool)
        ties = [x for x in pool if x['sample_evaluator_scalar_lane_work_upper'] == objective]
        expected = min(ties,key=lambda x:(x['evaluator_profile']['prepared_table_bytes'],tuple(x['key'])))
        assert winner['key'] == expected['key']
        assert winner['tied_prefixes'] == [x['key'] for x in ties]
    replay = [];replay_time = moment_time = generation_time = 0;scalar_states = 0
    for result in prices['configurations']:
        cid,method = result['case_id'],result['method']
        winner = next(x for x in selected if (x['row']['case_id'],x['row']['method']) == (cid,method))
        row = winner['row'];model = models[row['model_id']]
        assert result['key'] == winner['key']
        n = winner['certified_samples'];assert result['samples'] == n
        expected_seed = [protocol['seed'],cases.index(cid),METHODS.index(method)]
        assert result['seed_sequence'] == expected_seed
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(expected_seed)))
        graph = build(model,row);total = square_total = remaining = 0
        remaining = n;chunks = 0;input_hash = hashlib.sha256();output_hash = hashlib.sha256()
        while remaining:
            size = min(protocol['batch_size'],remaining)
            start = time.perf_counter()
            inputs = {name:rng.integers(0,1 << bits,size=size,dtype=np.uint64) for name,bits in graph.inputs.items()}
            generation_time += time.perf_counter()-start
            for name,array in inputs.items():
                input_hash.update(name.encode());input_hash.update(array.astype('<u8',copy=False).tobytes())
            start = time.perf_counter();out = vector_integer(graph,inputs,bank);replay_time += time.perf_counter()-start
            output_hash.update(out.astype('<u8',copy=False).tobytes())
            start = time.perf_counter();objects = out.astype(object)
            total += int(objects.sum());square_total += int((objects*objects).sum())
            moment_time += time.perf_counter()-start
            if chunks == 0:
                for index in sorted(set([0,min(1,size-1),size//2,size-1])):
                    assert int(out[index]) == graph.evaluate({name:int(array[index]) for name,array in inputs.items()},bank)
                    scalar_states += 1
            remaining -= size;chunks += 1
        assert result['chunks'] == chunks
        assert result['input_stream_sha256'] == input_hash.hexdigest()
        assert result['output_stream_sha256'] == output_hash.hexdigest()
        m = result['exact_moments'];unit = 1 << row['q']
        assert m['q'] == row['q'] and m['count'] == n
        assert int(m['integer_sum']) == total and int(m['integer_sum_of_squares']) == square_total
        numerator = n*square_total-total*total
        assert numerator >= 0 and int(m['sample_variance_integer_numerator']) == numerator
        assert Fraction(m['normalised_mean_exact_fraction']) == Fraction(total,n*unit)
        assert Fraction(m['normalised_sample_variance_exact_fraction']) == Fraction(numerator,n*(n-1)*unit*unit)
        assert result['confidence_failure_upper_decimal'] == winner['certificate']['confidence_budgets']['single_price']['certified_failure_upper_decimal']
        assert result['work']['financial_evaluations'] == n
        assert result['work']['evaluator_scalar_lane_operations_upper'] == winner['sample_evaluator_scalar_lane_work_upper']
        assert result['work']['random_finite_address_draws'] == n*len(graph.inputs)
        assert result['work']['random_address_bits_requested'] == n*sum(graph.inputs.values())
        assert result['work']['exact_moment_reduction_samples'] == n
        assert result['work']['prepared_table_bytes'] == winner['evaluator_profile']['prepared_table_bytes']
        keys = {(x['kernel'],x['address_bits'],x['q']) for x in winner['evaluator_profile']['table_manifest']}
        reconstructed = sum(table_records[k]['first_generation_wall_seconds'] for k in keys)
        assert abs(reconstructed-result['timing']['chosen_table_generation_wall_seconds_reconstructed']) <= 1e-10
        ci = interval_audit(model,row,result)
        replay.append({'case_id':cid,'method':method,'samples':n,'chunks':chunks,
            'full_stream_hashes_equal':True,'python_object_moments_equal':True,'monetary_interval':ci})
        if len(replay)%12 == 0:print('Independently replayed completed classical prices',len(replay),flush=True)
    samples = sum(x['samples'] for x in replay)
    assert prices['summary']['prices'] == len(replay) == 108
    assert prices['summary']['samples'] == samples
    assert prices['summary']['all_sample_price_interval_radii_at_most_0.01'] is True
    out = {'status':'passed','candidate_profiles':len(profile_checks),'deterministically_selected_prices':len(replay),
        'fully_replayed_samples':samples,'full_scalar_cross_checks':scalar_states,
        'table_records_independently_regenerated':len(bank),'candidate_budget_and_profile_checks':profile_checks,
        'selected_price_replays':replay,'protocol_sha256':sha(args.protocol),'prices_sha256':sha(args.prices),
        'reviewed_source_hashes':protocol['source_hashes'],
        'cold_shared_table_generation_seconds_diagnostic':table_time,
        'candidate_preparation_seconds_diagnostic':preparation_time,'independent_integer_replay_seconds_diagnostic':replay_time,
        'python_object_moment_reduction_seconds_diagnostic':moment_time,
        'random_stream_regeneration_seconds_diagnostic':generation_time,'elapsed_seconds_diagnostic':time.perf_counter()-started,
        'numpy_version':np.__version__,
        'claim_limit':'per-price analytical99 per cent coverage assumes idealIID finite uniforms; full pseudorandom replay checks execution only; geometric correction accuracy rests on its separate proved allowance, not uncertified quadrature'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('candidate_budget_and_profile_checks','selected_price_replays','reviewed_source_hashes')},indent=2))


if __name__ == '__main__':main()
