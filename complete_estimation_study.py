"""Attach certified estimation budgets to each architecture's selected graph.

I retain every feasible prefix, method and confidence scope. Complete-cost
selection is restricted to the certified finite opponent-inherited pool.
"""
from collections import Counter
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_FLOOR
from pathlib import Path
import hashlib
import json
import argparse
from estimation_prior import prior_for_row
from iterative_estimator import plan, decimal_record

ROOT = Path(__file__).resolve().parent


def architecture_budget(row, alpha, prior_upper):
    setup = plan(Fraction(1, 200*int(row['cap'])), alpha, prior_upper)
    rr = row['resources']
    shots, queries = setup.total_shots_upper, setup.total_grover_calls_upper
    invocations = shots+2*queries
    reflection = rr['grover_zero_reflection']
    gates = Counter({k:invocations*v for k,v in rr['expanded_gate_counts'].items()})
    for key,value in reflection['gate_counts'].items():
        gates[{'CNOT':'CX'}.get(key,key)] += queries*value
    gates['Z'] += queries
    gates['H'] += invocations*rr['includes_finite_uniform_preparation_H_gates']
    ccx = gates.pop('CCX',0)
    gates['CX'] += 6*ccx
    gates['H'] += 2*ccx
    gates['T'] += 4*ccx
    gates['TD'] += 3*ccx
    total_t = shots*rr['logical_T_count']+queries*rr['grover_iterate_T_count']
    assert gates['T']+gates['TD'] == total_t
    oracle_measures = invocations*rr['measurement_count']
    assert gates['MZ'] == oracle_measures
    td_a = rr['logical_T_depth_serial_upper']
    td_q = rr['grover_iterate_T_depth_serial_upper']
    a_depth = rr['serial_all_gate_depth_upper']+1
    q_depth = 2*a_depth+sum(reflection['gate_counts'].values())+1
    ccx_a=rr['expanded_gate_counts'].get('CCX',0)
    ccx_q=2*ccx_a+reflection['gate_counts'].get('CCX',0)
    shot_ccx=ccx_a+setup.grover_power_upper*ccx_q
    return {'plan':decimal_record(setup),
        'T_count_upper':total_t,
        'clean_A_and_adjoint_invocations_upper':invocations,
        'grover_invocations_upper':queries,
        'flag_measurements_upper':shots,
        'lookup_measurements_upper':oracle_measures,
        'all_measurements_upper':oracle_measures+shots,
        'conditional_CZ_feedback_upper':invocations*rr['classical_feedback_CZ_count'],
        'conditional_lookup_resets_upper':invocations*rr['classical_reset_X_count'],
        'fresh_shot_register_resets_upper':shots,
        'assigned_Clifford_T_and_instrument_operation_counts_upper':dict(gates),
        'legacy_assigned_three_layer_serial_T_depth':shots*td_a+queries*td_q,
        'expanded_five_layer_CCX_serial_T_depth_upper':shots*td_a+queries*td_q+2*ccx,
        'legacy_assigned_three_layer_maximum_shot_T_depth':td_a+setup.grover_power_upper*td_q if shots else 0,
        'expanded_five_layer_CCX_maximum_shot_T_depth_upper':td_a+setup.grover_power_upper*td_q+2*shot_ccx if shots else 0,
        'serial_lookup_measurement_depth_upper':invocations*rr['measurement_depth_serial_upper'],
        'logical_macro_serial_dependency_depth_upper':shots*a_depth+queries*q_depth,
        'expanded_CCX_serial_all_gate_dependency_depth_upper':shots*a_depth+queries*q_depth+14*ccx,
        'logical_macro_maximum_shot_dependency_depth_upper':a_depth+setup.grover_power_upper*q_depth if shots else 0,
        'peak_logical_qubits_upper':rr['grover_iterate_logical_qubits_upper'] if shots else 0,
        'uniform_Hadamards_upper':invocations*rr['includes_finite_uniform_preparation_H_gates'],
        'rotation_synthesis_bias':0,
        'depth_convention':'legacy three-layer assigned CCX profile is not an audited bound for the emitted sequence; separately retained five-layer15-gate explicit CCX expansion supplies conservative serial bounds',
        'interpretation':'ideal exact logical channel and fresh independent shots; sequential register reuse; no physical timing, fault tolerance, noise or classical runtime comparison'}


def generate(data):
    models = {m['model_id']:m for m in data['models']}
    records = []
    for row in data['configurations']:
        financial = prior_for_row(models[row['model_id']], row)
        single_alpha = Fraction(1,100)
        joint_alpha = Fraction(1,2700)
        prior = Fraction(financial['finite_mean_probability_upper_decimal'])
        with localcontext() as ctx:
            ctx.prec=220;ctx.rounding=ROUND_FLOOR
            joint_classical_alpha = Decimal(1)/Decimal(2700)
        joint_classical = prior_for_row(models[row['model_id']], row, joint_classical_alpha)
        records.append({**{k:row[k] for k in ('model_id','case_id','method','cap','q','exp_address_bits','log_address_bits','probability_address_bits')},
            'architecture':row['resources']['architecture'],
            'financial_prior_and_classical_single_price':financial,
            'classical_simultaneous_27_prices':joint_classical,
            'quantum_single_price_prior':architecture_budget(row,single_alpha,prior),
            'quantum_single_price_full_range':architecture_budget(row,single_alpha,Fraction(1)),
            'quantum_simultaneous_27_prices_prior':architecture_budget(row,joint_alpha,prior),
            'quantum_simultaneous_27_prices_full_range':architecture_budget(row,joint_alpha,Fraction(1)),
            'confidence_scope':'single-price99%; or one preselected architecture/method/prefix for each of27cases jointly99%; not all candidate outputs jointly99%'} )
    selected = []
    for architecture in ('exact7t','unitary8t','measured4t'):
        for case in sorted({r['case_id'] for r in records}):
            pool = [r for r in records if r['architecture']==architecture and r['case_id']==case]
            best = min(r['quantum_single_price_prior']['T_count_upper'] for r in pool)
            ties = sorted({r['method'] for r in pool if r['quantum_single_price_prior']['T_count_upper']==best})
            best_joint=min(r['quantum_simultaneous_27_prices_prior']['T_count_upper'] for r in pool)
            method_minima=[]
            for method in ('direct','interval','full_global','full_cap'):
                arm=[r for r in pool if r['method']==method]
                def arm_min(field):
                    minimum=min(r[field]['T_count_upper'] for r in arm)
                    winning=[r for r in arm if r[field]['T_count_upper']==minimum]
                    representative=min(winning,key=lambda r:(r[field]['all_measurements_upper'],r[field]['peak_logical_qubits_upper'],r['q'],r['exp_address_bits'],r['log_address_bits'],r['probability_address_bits']))
                    return {'T_upper':minimum,'tied_prefixes':[{k:r[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits')} for r in winning],
                        'representative_prefix':{k:representative[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits')}}
                method_minima.append({'method':method,'single_price_prior':arm_min('quantum_single_price_prior'),'simultaneous_27_prior':arm_min('quantum_simultaneous_27_prices_prior')})
            selected.append({'case_id':case,'architecture':architecture,'minimum_complete_T_upper':best,
                'tied_methods':ties,
                'minimum_simultaneous_27_complete_T_upper':best_joint,
                'simultaneous_27_tied_methods':sorted({r['method'] for r in pool if r['quantum_simultaneous_27_prices_prior']['T_count_upper']==best_joint}),
                'method_minima':method_minima,
                'scope':'complete confidence budgets minimised over every retained certified finite-pool prefix; not global architecture/precision optimisation'})
    return records, selected


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--config',type=Path,default=ROOT/'results/architecture_resource_candidates.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/complete_estimation_resources.json')
    a=p.parse_args()
    data=json.loads(a.config.read_text());rows,selected=generate(data)
    scripts=('complete_estimation_study.py','estimation_prior.py','iterative_estimator.py','architecture.py')
    out={'status':'complete finite-logical confidence and monetary budgets',
        'configuration_sha256':hashlib.sha256(a.config.read_bytes()).hexdigest(),
        'source_hashes':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in scripts},
        'configurations':rows,'method_comparison':selected,
        'selection_warning':'complete-cost selection within finite opponent-inherited pool; no global circuit/precision optimum',
        'contribution_boundary':'locally derived financial prior and estimator/circuit guarantees; theoretical novelty and quantum advantage unestablished'}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'configurations':len(rows),'method_comparisons':len(selected),
        'prior_T_upper_range':[min(r['quantum_single_price_prior']['T_count_upper'] for r in rows),max(r['quantum_single_price_prior']['T_count_upper'] for r in rows)]},indent=2))


if __name__=='__main__':main()
