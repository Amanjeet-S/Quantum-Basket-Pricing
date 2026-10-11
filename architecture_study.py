"""Count the complete finite inherited pool under all three QROM modes.

All 819 rows are retained for later end-to-end estimator selection. The
324 oracle/Grover minima and 81 method comparisons are intermediate
diagnostics, with exact T ties preserved, not global precision optima.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import sys
import time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
import architecture as arch
import precision_candidates as pool
import generalisation as gen


def representative_key(row):
    r=row['resources']
    return (r['measurement_count'],r['logical_qubits_serial_upper_bound'],
            r['logical_T_depth_serial_upper'],row['q'],row['exp_address_bits'],
            row['log_address_bits'],row['probability_address_bits'])


def select_intermediate(rows):
    groups=defaultdict(list)
    for row in rows:groups[(row['case_id'],row['method'],row['architecture'])].append(row)
    selected=[]
    for key,options in sorted(groups.items()):
        least=min(r['resources']['grover_iterate_T_count'] for r in options)
        ties=[r for r in options if r['resources']['grover_iterate_T_count']==least]
        chosen=min(ties,key=representative_key).copy()
        chosen['oracle_selection_objective']='minimum emitted Grover-iterate T count in the independently feasible inherited pool'
        chosen['oracle_exact_T_ties']=[r['candidate_id'] for r in ties]
        chosen['feasible_candidates_compared']=len(options)
        selected.append(chosen)
    best=[];bycase=defaultdict(list)
    for row in selected:bycase[(row['case_id'],row['architecture'])].append(row)
    for key,options in sorted(bycase.items()):
        least=min(r['resources']['grover_iterate_T_count'] for r in options)
        ties=[r for r in options if r['resources']['grover_iterate_T_count']==least]
        chosen=min(ties,key=lambda r:(*representative_key(r),r['method'])).copy()
        chosen['winning_methods']=[r['method'] for r in ties]
        chosen['fair_tie_statement']='all methods attaining the exact minimum retained; representative chosen only for display'
        best.append(chosen)
    return selected,best


def main():
    start=time.perf_counter();candidate_path=ROOT/'results/inherited_precision_candidates.json'
    candidates=json.loads(candidate_path.read_text());historical=json.loads((ROOT/'results/generalisation_resources.json').read_text())
    files=['architecture.py','precision_candidates.py','architecture_study.py',*candidates['base_source_hashes']]
    before={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
    for f,h in candidates['base_source_hashes'].items():assert before[f]==h
    models={m['model_id']:m for m in candidates['models']};bank={};blocks={};rows=[];regressions=[]
    bycase=defaultdict(list)
    for row in candidates['feasible_candidates']:bycase[row['case_id']].append(row)
    modes=('measured4t','unitary8t','exact7t')
    for cid,case in bycase.items():
        for cfg in case:
            graph=pool.build_candidate_graph(models[cfg['model_id']],cfg);gen.make_tables(graph,bank)
            identifier=hashlib.sha256(json.dumps({k:cfg[k] for k in ('case_id','method','q','b_y','b_z','b_u',
                'exp_address_bits','log_address_bits','probability_address_bits')},sort_keys=True).encode()).hexdigest()[:20]
            for mode in modes:
                resources=arch.counted_composition(graph,bank,blocks,mode)
                row={**cfg,'candidate_id':identifier,'architecture':mode,'resources':resources};rows.append(row)
                if mode=='exact7t':
                    old=next(r for r in historical['configurations'] if r['case_id']==cid and r['method']==cfg['method'])
                    fields=('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')
                    if all(cfg[k]==old[k] for k in fields):
                        assert resources['logical_T_count']==old['resources']['logical_T_count_exact_7T_Toffoli']
                        assert resources['grover_iterate_T_count']==old['resources']['grover_iterate_T_count_exact_7T_Toffoli']
                        assert resources['expanded_gate_counts']['X']==old['resources']['expanded_gate_counts']['X']
                        assert resources['expanded_gate_counts']['CX']==old['resources']['expanded_gate_counts']['CNOT']
                        assert resources['expanded_gate_counts']['CCX']==old['resources']['expanded_gate_counts']['Toffoli']
                        regressions.append({'case_id':cid,'method':cfg['method'],'same_prefix_T_and_Clifford_counts_match':True})
        print('Counted inherited architecture pool',cid,len(case),'prefixes',3*len(case),'resource rows',flush=True)
    assert len(rows)==3*len(candidates['feasible_candidates'])==819
    assert len(regressions)==108
    selected,best=select_intermediate(rows);assert len(selected)==324 and len(best)==81
    after={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files};assert before==after,'Sources changed during counting'
    metadata={'source_hashes':before,'inherited_pool_sha256':hashlib.sha256(candidate_path.read_bytes()).hexdigest(),
        'historical_resource_sha256':hashlib.sha256((ROOT/'results/generalisation_resources.json').read_bytes()).hexdigest(),
        'models':candidates['models'],'scope':'finite inherited pool, no global optimum; historical results remain unchanged',
        'architecture_modes':list(modes),'quantum_model':'ideal logical gates; exact arithmetic Toffolis cost seven T; all QROM phase gates, measurements and worst-case conditional Cliffords emitted',
        'measurement_assumption':'corrected Kraus branch proof on clean selector work; inverse uses a fresh instrument, no reversal of measurement',
        'source_paths':'relative to the directory containing architecture_study.py',
        'elapsed_seconds':time.perf_counter()-start}
    allout={**metadata,'status':'all independently feasible inherited candidates counted',
        'configurations':rows,'same_prefix_historical_regressions':regressions,
        'summary':{'financial_cases':27,'feasible_prefix_rows':273,'architecture_resource_rows':len(rows),
            'rejected_financial_transfers':len(candidates['rejected_candidates']),'historical_same_prefix_regressions':len(regressions)}}
    output=ROOT/'results/architecture_resource_candidates.json';output.write_text(json.dumps(allout,indent=2)+'\n')
    chosen={**metadata,'status':'intermediate oracle/Grover T-count minima, before estimator confidence-budget selection',
        'configurations':selected,'best_by_case_and_architecture':best,
        'all_candidate_resources_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'fair_tie_policy':'retain every method attaining the exact minimum; no forced preference for cap-aware or generic labels',
        'summary':{'selected_case_method_architecture_rows':len(selected),'best_case_architecture_rows':len(best),
            'cases_with_interval_attaining_the_minimum':sum('interval' in r['winning_methods'] for r in best),
            'cap_aware_candidate_set_contains_generic_candidate_set':True}}
    target=ROOT/'results/architecture_selected_resources.json';target.write_text(json.dumps(chosen,indent=2)+'\n')
    print(output);print(target)


if __name__=='__main__':main()
