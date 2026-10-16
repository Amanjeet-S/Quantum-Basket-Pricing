"""Exhaustive deterministic prices for the retained finite uniform law.

Only the pre-existing inherited candidate pool is inspected. The financial
interpreter uses exact integer arithmetic; confidence has no random-bit premise.
Operation counts are scalar-lane algorithm counts, never processor instructions.
"""
from collections import Counter
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
import classical_engine as ce
import generalisation as gen


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def display_fraction(value, upwards=False):
    with localcontext() as ctx:
        ctx.prec = max(200, len(str(abs(value.numerator)))+len(str(value.denominator))+40)
        ctx.rounding = ROUND_CEILING if upwards else ROUND_FLOOR
        return str((Decimal(value.numerator)/Decimal(value.denominator)).quantize(Decimal('1e-30')))


def graph_for(model, row):
    return gen.build_graph(model, row['method'], row['q'], row['b_y'], row['b_z'], row['b_u'],
        row['exp_address_bits'], row['log_address_bits'], row['probability_address_bits'], row['cap'])


def finite_space(model, row):
    exponent = model['outer_dimension']*row['b_y']
    if row['method'] == 'direct': exponent += row['b_z']
    elif row['method'] == 'interval': exponent += row['b_u']
    return exponent, 1 << exponent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--protocol', type=Path, default=ROOT/'results/classical_pricing_protocol.json')
    parser.add_argument('--output', type=Path, default=ROOT/'results/deterministic_pricing.json')
    parser.add_argument('--maximum-states', type=int, default=1 << 17)
    args = parser.parse_args()
    if args.maximum_states != 1 << 17:
        raise ValueError('This retained study prespecifies maximum_states=131072')
    started = time.perf_counter()
    protocol = json.loads(args.protocol.read_text())
    models = {m['model_id']:m for m in protocol['models']}
    inventory=[]; records=[]; tables={}; generation_records=[]; corrections={}
    for item in protocol['all_candidate_profiles']:
        row=item['row']; model=models[row['model_id']]
        exponent, states = finite_space(model, row)
        retained_scope = model['assets']==2 and row['method'] in ('full_global','full_cap')
        inventory.append({'candidate_id':row['candidate_id'],'case_id':row['case_id'],
            'method':row['method'],'q':row['q'],'outer_dimension':model['outer_dimension'],
            'grid_address_bits':exponent,'tensor_states':str(states),
            'certified_MC_samples':item['certified_samples'],
            'tensor_states_to_certified_MC_samples_exact_fraction':str(Fraction(states,item['certified_samples'])),
            'enumerated':retained_scope and states<=args.maximum_states,
            'reason':('retained two-asset fully conditioned scope' if retained_scope and states<=args.maximum_states else
                'not enumerated under prespecified maximum-states cap' if states>args.maximum_states else
                'outside prespecified two-asset full-conditioning enumeration scope')})
        if not retained_scope or states>args.maximum_states: continue
        graph=graph_for(model,row)
        assert dict(graph.inputs)=={'iy0':row['b_y']}
        begin=time.perf_counter(); before=set(tables); gen.make_tables(graph,tables)
        generation_seconds=time.perf_counter()-begin
        new_tables=[]
        for key in sorted(set(tables)-before):
            kind,bits,q=key; n=len(tables[key])
            new_tables.append({'kernel':kind,'address_bits':bits,'q':q,'entries':n,
                'scalar_function_endpoint_evaluations':n+1 if kind!='quantile' else 0,
                'normal_inverse_seed_evaluations':n if kind=='quantile' else 0,
                'sha256':hashlib.sha256(json.dumps(tables[key],separators=(',',':')).encode()).hexdigest()})
        generation_records.extend(new_tables)
        engine=ce.FiniteEvaluator(graph,tables)
        begin=time.perf_counter(); addresses=np.arange(states,dtype=np.uint64)
        outputs=engine.evaluate({'iy0':addresses})
        # states <=2**17 and output <=2**30 prove this uint64 sum cannot overflow.
        assert states*(1 << row['q']) < 1 << 64
        integer_sum=int(np.sum(outputs,dtype=np.uint64))
        exact_mean=Fraction(integer_sum,states*(1 << row['q']))
        evaluation_seconds=time.perf_counter()-begin
        profile=engine.operation_profile()
        fn_nodes=profile['nodes_by_kind'].get('function',0)
        validation_predicates=len(graph.inputs)+fn_nodes+1
        components={'financial_evaluator_scalar_lane_operations_upper':states*profile['scalar_lane_operation_upper_total'],
            'input_address_assignments_upper':states,'input_address_increments_upper':states-1,
            'domain_comparisons_and_boolean_reductions_upper':validation_predicates*(2*states-1),
            'exact_mean_integer_additions':states-1,'final_rational_mean_construction':1}
        total_work=sum(components.values())
        ckey=(row['model_id'],row['cap'])
        if ckey not in corrections:corrections[ckey]=gen.geometric_correction(model,row['cap'])
        correction=corrections[ckey]
        point=Fraction(correction['undiscounted_value_decimal'])+row['cap']*exact_mean
        beta=Fraction(row['independent_deterministic_bias_upper_decimal'])
        assert beta<Fraction(1,200) and beta<Fraction(str(row['price_tolerance']))
        lower,upper=point-beta,point+beta
        record={'candidate_id':row['candidate_id'],'case_id':row['case_id'],'model_id':row['model_id'],
            'method':row['method'],'row':row,'enumerated_states':states,
            'integer_sum':str(integer_sum),'normalised_finite_mean_exact_fraction':str(exact_mean),
            'output_little_endian_u64_sha256':hashlib.sha256(outputs.astype('<u8',copy=False).tobytes()).hexdigest(),
            'undiscounted_price_point_exact_fraction':str(point),
            'undiscounted_price_interval_lower_exact_fraction':str(lower),
            'undiscounted_price_interval_upper_exact_fraction':str(upper),
            'undiscounted_price_interval_lower_decimal':display_fraction(lower),
            'undiscounted_price_interval_upper_decimal':display_fraction(upper,True),
            'continuous_target_error_bound_decimal':row['independent_deterministic_bias_upper_decimal'],
            'confidence':1,'statistical_error':0,'requested_price_tolerance':row['price_tolerance'],
            'retained_statistical_allowance_unused':row['reserved_statistical_price_allowance'],
            'geometric_correction':correction,'evaluator_profile':profile,
            'online_algorithm_operation_components_upper':components,
            'online_algorithm_operation_upper_total':total_work,
            'certified_MC_samples_same_prefix':item['certified_samples'],
            'finite_states_to_certified_MC_samples_exact_fraction':str(Fraction(states,item['certified_samples'])),
            'financial_evaluator_work_ratio_exact_fraction':str(Fraction(states,item['certified_samples'])),
            'setup':{'new_tables_generated_in_shared_run':new_tables,
                'table_generation_seconds_diagnostic':generation_seconds,
                'interpreter_preparation_seconds_diagnostic':engine.preparation_seconds,
                'graph_nodes':len(graph.nodes),
                'table_entries_prepared':sum(t['entries'] for t in profile['table_manifest']),
                'packed_table_word_reads':sum(t['entries'] for t in profile['table_manifest']),
                'unpacked_uint64_words_written':sum(t['entries']*(1 if t['kernel']=='quantile' else 2) for t in profile['table_manifest']),
                'qualification':'Complete retained table manifest and setup phase; special-function scalar arithmetic is not reduced to a CPU instruction count. Shared cache is disclosed, never charged as zero one-off preparation.'},
            'price_evaluation_and_sum_seconds_diagnostic':evaluation_seconds,
            'claim_limit':'Exact finite-grid average plus previously independently checked deterministic bias certificate. Neither empirical pricing validation nor a physical quantum/runtime comparison.'}
        records.append(record)
        if len(records)%12==0:print('Enumerated full prefixes',len(records),flush=True)
    selected=[]
    for case in sorted(set(r['case_id'] for r in records)):
        candidates=[r for r in records if r['case_id']==case]
        chosen=min(candidates,key=lambda r:(r['online_algorithm_operation_upper_total'],r['setup']['table_entries_prepared'],
            r['row']['q'],r['method'],r['candidate_id']))
        selected.append({'case_id':case,'candidate_id':chosen['candidate_id'],
            'method':chosen['method'],'q':chosen['row']['q'],'enumerated_states':chosen['enumerated_states'],
            'online_algorithm_operation_upper_total':chosen['online_algorithm_operation_upper_total'],
            'continuous_target_error_bound_decimal':chosen['continuous_target_error_bound_decimal'],
            'undiscounted_price_interval_lower_decimal':chosen['undiscounted_price_interval_lower_decimal'],
            'undiscounted_price_interval_upper_decimal':chosen['undiscounted_price_interval_upper_decimal'],
            'confidence':1,'selection':'Least retained online scalar-lane operation upper count; separately reported setup cells break ties; no global precision optimisation.'})
    out={'status':'completed exact finite-grid deterministic comparator','data':'Frozen synthetic models only',
        'scope':'All retained two-asset full_global and full_cap prefixes. No newly optimised financial prefix.',
        'protocol_sha256':sha(args.protocol),'source_hashes':{name:sha(ROOT/name) for name in
            ('deterministic_pricing.py','classical_engine.py','generalisation.py','general_kernels.py','encoding_resources.py','resource_study.py')},
        'maximum_enumerated_states':args.maximum_states,'tensor_inventory':inventory,
        'tensor_inventory_configurations':len(inventory),'enumerated_configurations':len(records),
        'total_enumerated_outputs':sum(r['enumerated_states'] for r in records),'configurations':records,'selected':selected,
        'all_first_generation_table_records':generation_records,
        'elapsed_seconds_diagnostic':time.perf_counter()-started,'numpy_version':np.__version__,
        'python_version':platform.python_version(),'platform':platform.platform(),
        'confidence_premise':'Deterministic integer enumeration and valid retained numerical-bias certificates; no pseudorandomness, sample variance or concentration inequality.',
        'comparison_limit':'Enumeration is an upper-cost constructive comparator, not a classical lower bound. Large tensor grids are not enumerated under the prespecified cap; no assertion that every alternative deterministic method is infeasible.',
        'novelty_boundary':'Exhaustive finite summation is established numerical practice. Its role here is a stronger matched comparator and an adverse check on quantum-advantage claims.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'configurations':len(records),'outputs':out['total_enumerated_outputs'],'selected':selected,
        'elapsed_seconds_diagnostic':out['elapsed_seconds_diagnostic']},indent=2))


if __name__=='__main__':main()
