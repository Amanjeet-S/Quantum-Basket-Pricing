"""Independent Python-integer audit of the prepared finite classical engine."""
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import classical_engine as ce
import generalisation as gen


def product_audit(rng):
    checked = 0
    for q in range(18,31):
        w = q+7;mask = (1 << w)-1
        edges = [0,1,(1 << q)-1,1 << q,(1 << q)+1,
                 (1 << (w-1))-1,1 << (w-1),mask]
        pairs = list(product(edges,repeat=2))
        pairs += [(int(x),int(y)) for x,y in rng.integers(0,1 << w,(2048,2),dtype=np.uint64)]
        a,b = np.asarray(pairs,dtype=np.uint64).T
        observed = ce.fixed_multiply_u64(a,b,q,w)
        expected = np.asarray([((int(x)*int(y)) >> q)&mask for x,y in zip(a,b)],dtype=np.uint64)
        assert np.array_equal(observed,expected),(q,w)
        checked += len(pairs)
    full = rng.bit_generator.random_raw((2048,2))
    edges = [0,1,(1 << 32)-1,1 << 32,(1 << 32)+1,(1 << 63)-1,1 << 63,(1 << 64)-1]
    pairs = np.concatenate((np.asarray(list(product(edges,repeat=2)),dtype=np.uint64),full))
    a,b = pairs.T
    generic = 0
    for drop in (0,1,18,30,32,37,63,64,65,96,127):
        for width in (1,25,37,64):
            actual = ce.multiply_shift_u64(a,b,drop,width)
            expected = np.asarray([((int(x)*int(y)) >> drop)&((1 << width)-1) for x,y in zip(a,b)],dtype=np.uint64)
            assert np.array_equal(actual,expected),(drop,width)
            generic += len(a)
    return {'fixed_q_limb_products':checked,'generic_128_bit_window_products':generic,
        'q_range':[18,30],'window_drop_range':[0,127],'maximum_input_width':64,
        'unshifted_products_exceeding_uint64_included':True}


def moment_audit(rng):
    cases = 0
    for q in (18,26,30):
        unit = 1 << q
        a = np.concatenate((np.asarray([0,1,unit-1,unit],np.uint64),rng.integers(0,unit+1,size=513,dtype=np.uint64)))
        single = ce.ExactMoments(q).add(a)
        joined = ce.ExactMoments(q).add(a[:113]).merge(ce.ExactMoments(q).add(a[113:]))
        assert single.as_dict() == joined.as_dict()
        assert single.total == sum(map(int,a))
        assert single.total_squares == sum(int(x)**2 for x in a)
        cases += len(a)
    q = 30;unit = 1 << q
    # The sum of squares exceeds 2**64 and spans the reduction chunk boundary.
    a = np.full((1 << 20)+3,unit,dtype=np.uint64)
    m = ce.ExactMoments(q).add(a)
    assert m.total_squares == len(a)*unit*unit and m.total == len(a)*unit
    assert int(m.as_dict()['sample_variance_integer_numerator']) == 0
    return {'ordinary_outputs':cases,'all_one_chunk_boundary_outputs':len(a),
            'merge_identity_and_overflowing_squared_sum':True}


def systematic_inputs(graph,rng,random_count):
    names = list(graph.inputs);bits = list(graph.inputs.values())
    rows = [[p*((1 << b)-1) for p,b in zip(pattern,bits)] for pattern in product((0,1),repeat=len(names))]
    centre = [1 << (b-1) for b in bits]
    rows.append(centre)
    for j,b in enumerate(bits):
        for value in sorted(set([0,1,(1 << b)-2,(1 << b)-1,(1 << (b-1))-1,1 << (b-1),(1 << (b-1))+1])):
            row = centre.copy();row[j] = value;rows.append(row)
    rows.extend(np.column_stack([rng.integers(0,1 << b,size=random_count,dtype=np.uint64) for b in bits]).tolist())
    matrix = np.asarray(rows,dtype=np.uint64)
    return {name:matrix[:,j] for j,name in enumerate(names)}


def primitive_table_audit(tables,widths):
    states = 0;records = []
    for key,raw in tables.items():
        kind,bits,q = key
        if kind == 'quantile':continue
        w = widths[key];mask = (1 << w)-1
        from encoding_resources import function_span
        span = function_span(kind);f = q+span.bit_length()-1-bits
        spacing = 1 << f;limit = span*(1 << q)
        edges = np.arange((1 << bits)+1,dtype=np.uint64)*np.uint64(spacing)
        candidates = []
        for step in (-1,0,1):
            # Signed index creation avoids unsigned wrap at the first endpoint.
            x = np.clip(edges.astype(np.int64)+step,0,limit-1).astype(np.uint64)
            candidates.append(x)
        x = np.unique(np.concatenate(candidates))
        import encoding_resources as er
        g = er.Graph(q);g.width = w
        inp = g.inp('x',q+span.bit_length()-1)
        g.output = g.function(kind,inp,bits)
        engine = ce.FiniteEvaluator(g,{key:raw},require_normalised_output=False)
        actual = engine.evaluate({'x':x})
        expected = np.fromiter(((raw[int(v) >> f]&mask)+
            (((raw[int(v) >> f] >> w)*(int(v)&(spacing-1))) >> f) for v in x),np.uint64,count=len(x))&np.uint64(mask)
        assert np.array_equal(actual,expected),key
        assert np.all(actual[:-1] <= actual[1:]),key
        states += len(x)
        records.append({'kernel':kind,'bits':bits,'q':q,'states':len(x)})
    return {'all_table_cell_boundaries_and_adjacent_states':states,'tables':records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resources',type=Path,default=ROOT/'results/architecture_selected_resources.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/classical_engine_audit.json')
    parser.add_argument('--random-count',type=int,default=257)
    args = parser.parse_args()
    started = time.perf_counter();rng = np.random.default_rng(20261006)
    products = product_audit(rng);moments = moment_audit(rng)
    retained = json.loads(args.resources.read_text())
    models = {m['model_id']:m for m in retained['models']}
    unique = {}
    for row in retained['configurations']:
        identity = tuple(row[k] for k in ('model_id','method','cap','q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits'))
        unique[identity] = row
    tables = {};widths = {};records = []
    table_seconds = 0;eval_seconds = 0;scalar_seconds = 0
    for row in unique.values():
        g = gen.build_graph(models[row['model_id']],row['method'],row['q'],row['b_y'],row['b_z'],row['b_u'],
            row['exp_address_bits'],row['log_address_bits'],row['probability_address_bits'],row['cap'])
        start = time.perf_counter();gen.make_tables(g,tables);table_seconds += time.perf_counter()-start
        for node in g.nodes.values():
            if node.kind in ('quantile','function'):
                key = ('quantile' if node.kind == 'quantile' else node.meta,node.bits,g.q)
                assert key not in widths or widths[key] == g.width
                widths[key] = g.width
        inputs = systematic_inputs(g,rng,args.random_count)
        engine = ce.FiniteEvaluator(g,tables)
        start = time.perf_counter();actual = engine.evaluate(inputs);eval_seconds += time.perf_counter()-start
        start = time.perf_counter()
        expected = np.fromiter((g.evaluate({name:int(v[i]) for name,v in inputs.items()},tables)
            for i in range(len(actual))),np.uint64,count=len(actual))
        scalar_seconds += time.perf_counter()-start
        assert np.array_equal(actual,expected),(row['case_id'],row['method'],row['q'])
        mm = ce.ExactMoments(g.q).add(actual)
        assert mm.total == sum(map(int,expected)) and mm.total_squares == sum(int(x)**2 for x in expected)
        records.append({'case_id':row['case_id'],'method':row['method'],'q':row['q'],
            'configuration_identity':list(tuple(row[k] for k in ('b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits'))),
            'states':len(actual),'matches_scalar_python_integer_graph':True,
            'exact_moments':mm.as_dict()})
        if len(records)%12 == 0:print('Audited finite configurations',len(records),flush=True)
    primitive = primitive_table_audit(tables,widths)
    out = {'status':'passed','data':'frozen synthetic models only','finite_configurations':len(records),
        'production_and_systematic_states':sum(r['states'] for r in records),
        'products':products,'exact_moments':moments,'interpolation':primitive,
        'configuration_checks':records,'resources_sha256':hashlib.sha256(args.resources.read_bytes()).hexdigest(),
        'source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('classical_engine.py','generalisation.py','general_kernels.py','encoding_resources.py','reversible.py')},
        'table_generation_seconds_diagnostic':table_seconds,'vector_evaluation_seconds_diagnostic':eval_seconds,
        'scalar_reference_seconds_diagnostic':scalar_seconds,'elapsed_seconds_diagnostic':time.perf_counter()-started,
        'numpy_version':np.__version__,'python_version':platform.python_version(),
        'claim_limit':'exact finite interpreter audit; seeded pseudorandom samples do not prove independence or continuous-law accuracy'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('configuration_checks','interpolation','source_hashes')},indent=2))


if __name__ == '__main__':main()
