"""Complete independent Python-integer replay of every enumerated grid."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import generalisation as gen


def main():
    start=time.perf_counter()
    path=ROOT/'results/deterministic_pricing.json';out=json.loads(path.read_text())
    protocol_path=ROOT/'results/classical_pricing_protocol.json'
    assert out['protocol_sha256']==hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    protocol=json.loads(protocol_path.read_text());models={m['model_id']:m for m in protocol['models']}
    continuous=json.loads((ROOT/'results/continuous_quadrature.json').read_text())
    continuous_prices={r['case_id']:r['undiscounted_price_interval'] for r in continuous['configurations']}
    assert len(out['tensor_inventory'])==len(protocol['all_candidate_profiles'])==273
    inventory={x['candidate_id']:x for x in out['tensor_inventory']}
    for item in protocol['all_candidate_profiles']:
        row=item['row'];m=models[row['model_id']];exponent=m['outer_dimension']*row['b_y']
        if row['method']=='direct':exponent+=row['b_z']
        elif row['method']=='interval':exponent+=row['b_u']
        inv=inventory[row['candidate_id']]
        assert inv['grid_address_bits']==exponent and int(inv['tensor_states'])==1 << exponent
        assert Fraction(inv['tensor_states_to_certified_MC_samples_exact_fraction'])==Fraction(1 << exponent,item['certified_samples'])
    tables={};checks=[]
    for r in out['configurations']:
        row=r['row'];m=models[row['model_id']]
        assert m['assets']==2 and m['outer_dimension']==1 and row['method'] in ('full_global','full_cap')
        graph=gen.build_graph(m,row['method'],row['q'],row['b_y'],row['b_z'],row['b_u'],
            row['exp_address_bits'],row['log_address_bits'],row['probability_address_bits'],row['cap'])
        gen.make_tables(graph,tables)
        assert graph.inputs=={'iy0':row['b_y']}
        n=1 << row['b_y'];assert n==r['enumerated_states']<=out['maximum_enumerated_states']
        # Graph.evaluate uses arbitrary-precision Python multiplication and
        # signed clamps, independently of the vector limb-product interpreter.
        values=[graph.evaluate({'iy0':j},tables) for j in range(n)]
        assert min(values)>=0 and max(values)<=1 << row['q']
        assert sum(values)==int(r['integer_sum'])
        arr=np.asarray(values,dtype='<u8')
        assert hashlib.sha256(arr.tobytes()).hexdigest()==r['output_little_endian_u64_sha256']
        mean=Fraction(sum(values),n*(1 << row['q']))
        assert mean==Fraction(r['normalised_finite_mean_exact_fraction'])
        point=Fraction(r['geometric_correction']['undiscounted_value_decimal'])+row['cap']*mean
        beta=Fraction(row['independent_deterministic_bias_upper_decimal'])
        assert beta<Fraction(1,200)
        lo,hi=Fraction(r['undiscounted_price_interval_lower_exact_fraction']),Fraction(r['undiscounted_price_interval_upper_exact_fraction'])
        assert lo==point-beta and hi==point+beta
        assert Fraction(r['undiscounted_price_interval_lower_decimal'])<=lo
        assert Fraction(r['undiscounted_price_interval_upper_decimal'])>=hi
        assert r['confidence']==1 and r['statistical_error']==0
        c=continuous_prices[r['case_id']]
        assert lo<=Fraction(c['lower_exact_fraction'])<=Fraction(c['upper_exact_fraction'])<=hi
        work=r['online_algorithm_operation_components_upper']
        assert sum(work.values())==r['online_algorithm_operation_upper_total']
        nodes=r['evaluator_profile']['nodes_by_kind']; predicates=1+nodes['function']+1
        assert work['domain_comparisons_and_boolean_reductions_upper']==predicates*(2*n-1)
        assert work['financial_evaluator_scalar_lane_operations_upper']==n*r['evaluator_profile']['scalar_lane_operation_upper_total']
        checks.append({'candidate_id':r['candidate_id'],'case_id':r['case_id'],'states':n,
            'all_outputs_match_independent_Python_integer_interpreter':True,
            'exact_sum_matches':True,'continuous_Arb_enclosure_contained_in_finite_grid_certificate':True})
        if len(checks)%12==0:print('Independently replayed grids',len(checks),flush=True)
    for s in out['selected']:
        eligible=[r for r in out['configurations'] if r['case_id']==s['case_id']]
        winner=min(eligible,key=lambda r:(r['online_algorithm_operation_upper_total'],r['setup']['table_entries_prepared'],
            r['row']['q'],r['method'],r['candidate_id']))
        assert winner['candidate_id']==s['candidate_id']
    audit={'status':'passed','complete_grid_configurations':len(checks),
        'independently_replayed_outputs':sum(x['states'] for x in checks),
        'exact_tensor_inventory_configurations':len(inventory),'selected_price_certificates':len(out['selected']),
        'checks':checks,'result_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'continuous_result_sha256':hashlib.sha256((ROOT/'results/continuous_quadrature.json').read_bytes()).hexdigest(),
        'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'claim_limit':'Independent execution and exact summation check. Continuous accuracy remains conditional on the retained bias proof and validated Arb integration assumptions.'}
    (ROOT/'results/deterministic_pricing_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
