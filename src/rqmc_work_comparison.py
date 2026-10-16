"""Reconstruct ancillary work and a fixed IID plan from retained priors.

The IID counts are deterministic guaranteed plans, not new executed prices.
No realised batch mean or variance enters their prefix or sample selection.
All costs below are abstract operations with explicitly separated units.
"""
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def square_reduction_counts(words):
    # These counts cover the explicit limb arithmetic and its three reductions.
    # Array validation, memory allocation and arbitrary-precision limb costs
    # are kept outside the scalar-lane arithmetic unit.
    return Counter({'uint64_AND':words,'uint64_shift':words,'uint64_multiply':3*words,
        'uint64_reduction_add':3*max(0,words-1),
        'Python_integer_shift_calls':2,'Python_integer_add_calls':2})


def main():
    protocol_path=ROOT/'results/rqmc_pricing_protocol.json'
    price_path=ROOT/'results/rqmc_pricing_results.json'
    protocol=json.loads(protocol_path.read_text());prices=json.loads(price_path.read_text())
    assert prices['protocol_sha256']==sha(protocol_path)
    profiles=protocol['all_candidate_profiles'];rows=[]
    with localcontext() as ctx:
        ctx.prec=300
        logarithm=Decimal(200).ln()+Decimal('1e-250')
        def dec(text):
            x=Fraction(text);return Decimal(x.numerator)/Decimal(x.denominator)
        for profile in profiles:
            variance=dec(profile['batch_variance_prior_exact_fraction'])
            rho=dec(profile['finite_output_range_exact_fraction'])
            epsilon=dec(profile['estimation_tolerance_exact_fraction'])
            nb=max(1,int(((2*variance+2*rho*epsilon/3)*logarithm/epsilon**2+Decimal('1e-250')).to_integral_value(rounding=ROUND_CEILING)))
            nh=max(1,int((rho*rho*logarithm/(2*epsilon**2)+Decimal('1e-250')).to_integral_value(rounding=ROUND_CEILING)))
            samples=min(nb,nh);kind='Bernstein' if nb<=nh else 'Hoeffding'
            n=Decimal(samples)
            failure=(2*(-n*epsilon**2/(2*variance+2*rho*epsilon/3)).exp()
                     if kind=='Bernstein' else 2*(-2*n*epsilon**2/(rho*rho)).exp())
            assert failure+Decimal('1e-250')<Decimal('.01')
            rows.append({'key':profile['key'],'Bernstein_samples':nb,'Hoeffding_samples':nh,
                'selected_bound':kind,'certified_samples':samples,
                'failure_upper_decimal':str(failure+Decimal('1e-250')),
                'scalar_lane_operation_count':profile['evaluator_profile']['scalar_lane_operation_upper_total'],
                'sample_evaluator_scalar_lane_work_upper':samples*profile['evaluator_profile']['scalar_lane_operation_upper_total'],
                'separate_IID_generation_and_mean_reduction_work':{
                    'independent_random_address_word_draws':samples*len(profile['evaluator_profile']['inputs']),
                    'random_address_bits_requested':samples*sum(profile['evaluator_profile']['inputs'].values()),
                    'uint64_mean_reduction_additions_upper':samples,
                    'Python_integer_chunk_merge_additions_upper':(samples+(1<<20)-1)//(1<<20),
                    'maximum_uint64_mean_reduction_chunk':1<<20,
                    'normalised_mean_Fraction_constructor_calls':1,
                    'scope':'Known variance/range priors need only an integer output sum for this fixed-budget mean estimator; optional sample-square diagnostics are not required. At most2^20 outputs of at most2^30 are reduced per uint64 chunk.'},
                'finite_range_exact_fraction':profile['finite_output_range_exact_fraction'],
                'variance_prior_upper_exact_fraction':profile['batch_variance_prior_exact_fraction']})
    records=[]
    for selected,result in zip(protocol['selected'],prices['configurations']):
        pool=[p for p in rows if p['key'][0]==result['case_id']]
        minimum=min(p['sample_evaluator_scalar_lane_work_upper'] for p in pool)
        ties=[p for p in pool if p['sample_evaluator_scalar_lane_work_upper']==minimum]
        # Reproduce the graph selection's deterministic setup tie-breaker.
        original={tuple(p['key']):p for p in profiles}
        iid=min(ties,key=lambda p:(original[tuple(p['key'])]['evaluator_profile']['prepared_table_bytes'],tuple(p['key'])))
        work=Counter()
        for checkpoint in result['checkpoints']:
            r=checkpoint['replicates'];pairs=r//2
            work.update(square_reduction_counts(r))
            work.update(square_reduction_counts(pairs))
            work.update({'uint64_pair_maximum':pairs,'uint64_pair_minimum':pairs,
                'uint64_pair_subtraction':pairs,'uint64_reduction_add':r-1,
                'Python_fraction_constructor_calls':3,
                'Python_integer_multiply_calls_sample_variance':4,
                'Python_integer_subtraction_calls_sample_variance':1})
        ratio=Fraction(result['work']['evaluator_scalar_lane_operations_upper'],iid['sample_evaluator_scalar_lane_work_upper'])
        interval=result['confidence_interval']
        actual_radius=Fraction(Decimal(interval['maximum_radius_decimal'])) if interval else None
        checkpoints=result['checkpoints']
        records.append({'case_id':result['case_id'],'fixed_IID_plan':iid,
            'fixed_IID_tied_prefixes':[p['key'] for p in ties],
            'RQMC_to_fixed_IID_evaluator_work_ratio_exact_fraction':str(ratio),
            'RQMC_to_fixed_IID_evaluator_work_ratio_float_diagnostic':float(ratio),
            'RQMC_statistical_reservation_met':result['status']=='certified per-price tolerance met',
            'RQMC_actual_monetary_radius_at_most_0_01':actual_radius is not None and actual_radius<=Fraction(1,100),
            'explicit_checkpoint_moment_arithmetic_counts':dict(work),
            'checkpoint_log_upper_calls':3*len(checkpoints),
            'checkpoint_enclosed_square_root_calls':2*len(checkpoints),
            'retained_batch_sum_words':result['replicates'],
            'retained_batch_sum_storage_bits_bound':result['replicates']*(selected['row']['q']+7),
            'net_word_coordinates_stored':64*len(selected['evaluator_profile']['inputs']),
            'net_word_storage_bytes':64*len(selected['evaluator_profile']['inputs'])*8,
            'net_preparation_shared_scalar_lane_operations':{
                'Gray_code_shift':64,'Gray_code_XOR':64,'Gray_component_shift':6*64,
                'Gray_component_AND':6*64,'direction_component_multiply':6*64*len(selected['evaluator_profile']['inputs']),
                'point_component_XOR':6*64*len(selected['evaluator_profile']['inputs']),
                'address_truncation_shift':64*len(selected['evaluator_profile']['inputs'])}})
    ratios=[Fraction(r['RQMC_to_fixed_IID_evaluator_work_ratio_exact_fraction']) for r in records]
    output={'status':'deterministic ancillary-work and fixed-IID planning comparison complete',
        'source_sha256':sha(__file__),'protocol_sha256':sha(protocol_path),'prices_sha256':sha(price_path),
        'fixed_IID_failure_allocation_exact_fraction':'1/100',
        'fixed_IID_selection_scope':'Same137 fully conditioned finite prefixes, global deterministic range and moment priors only; computed after the experiment but without its observed outputs; no new IID execution is implied',
        'all_fixed_IID_candidate_plans':rows,'configurations':records,
        'summary':{'planned_IID_evaluator_work_lower_than_realised_RQMC_cases':sum(x>1 for x in ratios),
            'work_ratio_minimum_exact_fraction':str(min(ratios)),'work_ratio_maximum_exact_fraction':str(max(ratios)),
            'RQMC_statistical_reservation_successes':sum(r['RQMC_statistical_reservation_met'] for r in records),
            'RQMC_actual_monetary_tolerance_successes':sum(r['RQMC_actual_monetary_radius_at_most_0_01'] for r in records)},
        'qualification':'Financial scalar-lane arithmetic, uint64 reduction lanes, random word requests, table memory, special-function calls and Python arbitrary-precision calls have distinct units. These counts do not measure CPU instructions, wall-clock lower bounds or physical quantum speedup. Memory allocation, validation, direction-file parsing and variable-size integer limb costs are not folded into the scalar-lane metric.'}
    (ROOT/'results/rqmc_work_comparison.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['summary'],indent=2))


if __name__=='__main__':main()
