"""Finite-oracle prices from independently shifted deterministic digital nets.

Each shift, not each net point, is one independent observation in the proved
confidence model. Integer financial outputs and batch sums are exact. NumPy's
PCG64DXSM execution is reproducible; its seed does not prove independence.

The Joe--Kuo direction recurrence below is adapted from the retained author
source, copyright (c) 2008 Frances Y. Kuo and Stephen Joe. The complete BSD
notice and disclaimer are retained in public_data/licence and sobol.cc.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

import generalisation as gen
from classical_engine import FiniteEvaluator
from classical_pricing_study import build, prepare_tables, row_key

ROOT = Path(__file__).resolve().parent
NET_BITS = 6
POINTS_PER_SHIFT = 1 << NET_BITS
CHECKPOINTS = tuple(1 << k for k in range(6, 17))
PRODUCTION_SEED = 2026100707
ALPHA = Fraction(1, 100)
ALPHA_EVENT = ALPHA / (2 * len(CHECKPOINTS))
MAX_FINANCIAL_LANES_PER_CHUNK = 1 << 20
METHODS = ('full_global', 'full_cap')
SOURCE_FILES = ('rqmc_pricing.py', 'rqmc_confidence.py', 'classical_engine.py',
    'classical_pricing_study.py', 'generalisation.py', 'general_kernels.py',
    'encoding_resources.py', 'public_data/new-joe-kuo-6.21201',
    'public_data/sobol.cc', 'public_data/licence')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serialise_fraction(value):
    return str(Fraction(value))


def sobol_words(net_bits, dimension, direction_path=None):
    """Author-distributed 32-bit recurrence, evaluated in Gray-code order."""
    if not 1 <= net_bits <= 30 or not 1 <= dimension <= 8:
        raise ValueError('Outside the retained direction and study dimensions')
    path = direction_path or ROOT / 'public_data/new-joe-kuo-6.21201'
    lines = Path(path).read_text().splitlines()
    directions = np.zeros((dimension, 32), dtype=np.uint32)
    directions[0] = np.asarray([1 << (31-i) for i in range(32)], dtype=np.uint32)
    for dim in range(1, dimension):
        values = list(map(int, lines[dim].split()))
        dd, degree, polynomial = values[:3]
        initial = values[3:]
        if dd != dim + 1 or len(initial) != degree:
            raise ValueError('Direction file dimensions are inconsistent')
        for j in range(degree):
            directions[dim, j] = np.uint32(initial[j] << (31-j))
        for j in range(degree, 32):
            word = int(directions[dim, j-degree]) ^ (int(directions[dim, j-degree]) >> degree)
            for k in range(1, degree):
                if (polynomial >> (degree-1-k)) & 1:
                    word ^= int(directions[dim, j-k])
            directions[dim, j] = np.uint32(word)
    indices = np.arange(1 << net_bits, dtype=np.uint32)
    gray = indices ^ (indices >> np.uint32(1))
    points = np.zeros((1 << net_bits, dimension), dtype=np.uint32)
    for j in range(net_bits):
        points ^= ((gray >> np.uint32(j)) & np.uint32(1))[:, None] * directions[None, :, j]
    return points


def digital_points(inputs, net_bits=NET_BITS):
    names = tuple(inputs)
    words = sobol_words(net_bits, len(names))
    points = {}
    for j, (name, address_bits) in enumerate(inputs.items()):
        if not net_bits <= address_bits <= 32:
            raise ValueError('Digital net needs the retained finite address width')
        points[name] = (words[:, j] >> np.uint32(32-address_bits)).astype(np.uint64)
    return points


def draw_shifts(inputs, count, rng):
    return {name: rng.integers(0, 1 << bits, size=count, dtype=np.uint64)
            for name, bits in inputs.items()}


def shifted_inputs(points, shifts):
    """Point-major within each replicate, with no floating-point uniforms."""
    if set(points) != set(shifts):
        raise ValueError('Point and shift coordinates differ')
    return {name: np.bitwise_xor(shifts[name][:, None], points[name][None, :]).reshape(-1)
            for name in points}


def exact_u64_square_sum(values):
    """Exact integer square sum for <=2^20 words, each at most 2^36.

    With a=low+2^16 high, low^2 sums have at most 52 bits, high^2
    sums at most 61 bits and low*high sums at most 57 bits. Assembly is
    in Python arbitrary precision. The bounds include the endpoint 2^36.
    """
    a = np.asarray(values, dtype=np.uint64)
    if a.ndim != 1 or len(a) > 1 << 20 or np.any(a > 1 << 36):
        raise ValueError('Outside the proved exact batch-square range')
    low = a & np.uint64(65535)
    high = a >> np.uint64(16)
    return (int(np.sum(low * low, dtype=np.uint64))
            + (int(np.sum(low * high, dtype=np.uint64)) << 17)
            + (int(np.sum(high * high, dtype=np.uint64)) << 32))


def batch_statistics(sums, q, points_per_shift=POINTS_PER_SHIFT):
    a = np.asarray(sums, dtype=np.uint64)
    if not len(a) or len(a) > 1 << 20 or not 18 <= q <= 30:
        raise ValueError('Outside the retained production replicate family')
    denominator = points_per_shift * (1 << q)
    if denominator > 1 << 36 or np.any(a > denominator):
        raise ValueError('Invalid batch output')
    pairs = len(a) // 2
    left, right = a[:2*pairs:2], a[1:2*pairs:2]
    differences = np.maximum(left, right) - np.minimum(left, right)
    total = int(np.sum(a, dtype=np.uint64))
    squares = exact_u64_square_sum(a)
    paired_squares = exact_u64_square_sum(differences)
    return {'replicates': len(a), 'batch_denominator': str(denominator),
        'integer_batch_sum': str(total), 'integer_batch_square_sum': str(squares),
        'independent_pairs': pairs, 'integer_paired_difference_square_sum': str(paired_squares),
        'finite_mean_exact_fraction': str(Fraction(total, len(a) * denominator)),
        'paired_variance_estimate_exact_fraction': str(Fraction(paired_squares, 2 * pairs * denominator**2)),
        'ordinary_replicate_sample_variance_exact_fraction': str(Fraction(len(a)*squares-total**2,
            len(a)*(len(a)-1)*denominator**2)) if len(a) > 1 else None}


def graph_manifest(graph):
    return {'q': graph.q, 'width': graph.width, 'inputs': dict(graph.inputs), 'output': graph.output,
        'nodes': [{'name': name, 'kind': node.kind, 'args': list(node.args),
                   'meta': node.meta, 'bits': node.bits} for name, node in graph.nodes.items()]}


def pure_variance_bound(certificate):
    rho = Fraction(certificate['finite_lattice_range_fraction'])
    upper_mean = Fraction(certificate['fair_mean_upper_fraction'])
    second_moment = Fraction(certificate['fair_second_moment_upper_fraction'])
    if rho == 0:
        return rho, upper_mean, second_moment, Fraction(0)
    mean = min(upper_mean, rho / 2, second_moment / rho)
    return rho, upper_mean, second_moment, mean * (rho-mean)


def prospective_replicates(variance, rho, epsilon):
    # Jensen preserves the global second-moment bound for a batch, with no
    # assumed variance reduction. This prospective rule chooses only a graph.
    with localcontext() as ctx:
        ctx.prec = 100
        def dec(x): return Decimal(x.numerator) / Decimal(x.denominator)
        logarithm = (Decimal(2) / dec(ALPHA_EVENT)).ln()
        budget = (2*dec(variance) + 2*dec(rho)*dec(epsilon)/3) * logarithm / dec(epsilon)**2
        return max(1, int(budget.to_integral_value(rounding=ROUND_CEILING)))


def prepare_protocol(include_controls=True):
    start, cpu = time.perf_counter(), time.process_time()
    path = ROOT / 'results/architecture_resource_candidates.json'
    range_path = ROOT / 'results/finite_range_certificates.json'
    raw, ranges = json.loads(path.read_text()), json.loads(range_path.read_text())
    if ranges['configuration_sha256'] != digest(path):
        raise ValueError('Range certificate configuration hash differs')
    rows = {row_key(row): row for row in raw['configurations'] if row['architecture'] == 'measured4t'}
    if len(rows) != 273:
        raise ValueError('The inherited finite candidate pool changed')
    certs = {row_key(row): row for row in ranges['configurations']}
    if set(certs) != set(rows):
        raise ValueError('Range certificates and candidate pool differ')
    models = {model['model_id']: model for model in raw['models']}
    bank, table_records, profiles = {}, {}, []
    for key, row in sorted(rows.items()):
        if row['method'] not in METHODS:
            continue
        model = models[row['model_id']]
        if model['assets'] == 2 and not include_controls:
            continue
        graph = build(model, row)
        prepare_tables(graph, bank, table_records)
        evaluator = FiniteEvaluator(graph, bank)
        profile = evaluator.operation_profile()
        rho, mean, second, variance = pure_variance_bound(certs[key])
        epsilon = Fraction(1, 200 * row['cap'])
        plan = prospective_replicates(variance, rho, epsilon)
        stripped = {k:v for k,v in row.items() if k != 'resources'}
        profiles.append({'key': list(key), 'row': stripped, 'range_certificate': certs[key],
            'finite_output_range_exact_fraction': str(rho), 'batch_variance_prior_exact_fraction': str(variance),
            'batch_second_moment_prior_exact_fraction': str(second),
            'batch_mean_prior_exact_fraction': str(mean),
            'estimation_tolerance_exact_fraction': str(epsilon),
            'prospective_independent_shift_count': plan, 'evaluator_profile': profile,
            'selection_scalar_lane_work_upper': plan * POINTS_PER_SHIFT * profile['scalar_lane_operation_upper_total']})
    selected = []
    for case_id in sorted({p['row']['case_id'] for p in profiles}):
        pool = [p for p in profiles if p['row']['case_id'] == case_id]
        best = min(p['selection_scalar_lane_work_upper'] for p in pool)
        tied = [p for p in pool if p['selection_scalar_lane_work_upper'] == best]
        winner = min(tied, key=lambda p:(p['evaluator_profile']['prepared_table_bytes'], tuple(p['key']))).copy()
        winner['tied_prefixes'] = [p['key'] for p in tied]
        graph = build(models[winner['row']['model_id']], winner['row'])
        winner['graph_manifest'] = graph_manifest(graph)
        selected.append(winner)
    expected = 27 if include_controls else 18
    if len(selected) != expected:
        raise ValueError(('Unexpected contract count', len(selected)))
    sources = {name:digest(ROOT/name) for name in SOURCE_FILES if (ROOT/name).exists()}
    protocol = {'status':'frozen before production random shifts', 'configuration_sha256':digest(path),
        'range_certificate_sha256':digest(range_path), 'source_hashes':sources,
        'models': raw['models'], 'net_bits':NET_BITS, 'points_per_shift':POINTS_PER_SHIFT,
        'checkpoint_replicates':list(CHECKPOINTS), 'maximum_replicates':CHECKPOINTS[-1],
        'alpha_per_price_exact_fraction':str(ALPHA), 'alpha_per_event_exact_fraction':str(ALPHA_EVENT),
        'confidence_allocation':'Two events at each of eleven prespecified checkpoints; union bound across all checkpoints; adaptive stopping within this finite family only',
        'production_seed':PRODUCTION_SEED, 'random_generator':'NumPy PCG64DXSM SeedSequence([production_seed, sorted_case_index])',
        'randomness_premise':'Independent exactly uniform shift vectors in the finite product input group; reproducible PCG execution does not establish this premise',
        'selection_rule':'Minimise prospective Bernstein replicate count times 64 times exact-evaluator conservative scalar-lane count over both fully conditioned methods in all inherited prefixes. Global finite moments and range enter before data. No 1/N variance reduction is assumed.',
        'all_candidate_profiles':profiles, 'selected':selected,
        'table_first_generation_records':list(table_records.values()),
        'all_candidate_setup_wall_seconds':time.perf_counter()-start,
        'all_candidate_setup_process_cpu_seconds':time.process_time()-cpu,
        'python_version':sys.version, 'numpy_version':np.__version__, 'platform':platform.platform(),
        'direction_source':{'author_page':'https://web.maths.unsw.edu.au/~fkuo/sobol/',
            'consulted_version':'Previously retained author-distributed direction file and 21 October 2008 sobol.cc; associated journal papers are not implied read',
            'licence':'public_data/licence; copyright Frances Y. Kuo and Stephen Joe 2008',
            'randomisation':'Finite XOR digital shift; not Owen scrambling'},
        'scope':'18 higher-dimensional frozen synthetic prices plus nine optional two-asset validation controls; no market data, novelty or physical quantum comparison'}
    outfile = ROOT / 'results/rqmc_pricing_protocol.json'
    outfile.write_text(json.dumps(protocol, indent=2)+'\n')
    return protocol, bank, table_records


def price_interval(model, row, mean_low, mean_high):
    """Outward monetary interval using the actual finite-mean enclosure."""
    with localcontext() as ctx:
        ctx.prec = 220
        def low_decimal(x):
            ctx.rounding=ROUND_FLOOR; x=Fraction(x)
            return Decimal(x.numerator)/Decimal(x.denominator)
        def high_decimal(x):
            ctx.rounding=ROUND_CEILING; x=Fraction(x)
            return Decimal(x.numerator)/Decimal(x.denominator)
        pad=Decimal('1e-100')
        correction=gen.geometric_correction(model,row['cap'])
        centre=Decimal(correction['undiscounted_value_decimal'])
        error=max(Decimal('1e-12'),Decimal.from_float(correction['absolute_allowance']))
        arithmetic=Decimal.from_float(row['arithmetic_price_error_upper'])
        distribution=Decimal.from_float(row['distribution_price_bound'])
        deterministic=arithmetic+distribution+error+pad
        if deterministic > Decimal('.005'):
            raise ValueError('Deterministic bias reservation exceeded')
        ctx.rounding=ROUND_FLOOR
        undiscounted_low=max(Decimal(0),centre+Decimal(row['cap'])*low_decimal(mean_low)-deterministic)
        ctx.rounding=ROUND_CEILING
        undiscounted_high=min(Decimal(row['cap']),centre+Decimal(row['cap'])*high_decimal(mean_high)+deterministic)
        ctx.rounding=ROUND_FLOOR
        discount=(-Decimal(model['risk_free_rate'])).exp()
        low=undiscounted_low*(discount-pad)-pad
        ctx.rounding=ROUND_CEILING
        high=undiscounted_high*(discount+pad)+pad
        midpoint=(low+high)/2
        return {'price_interval_lower_decimal':str(low),'price_interval_upper_decimal':str(high),
            'display_midpoint_decimal':str(midpoint),'maximum_radius_decimal':str((high-low)/2+pad),
            'geometric_correction':correction,'deterministic_currency_allowance_decimal':str(deterministic),
            'numerical_premise':'Correctly rounded Decimal elementary functions at220 digits, with retained1e-12 correction allowance and1e-100 enclosing increment'}


def run_production(protocol, bank, table_records):
    from rqmc_confidence import confidence_from_batch_statistics
    protocol_path=ROOT/'results/rqmc_pricing_protocol.json'
    protocol_sha=digest(protocol_path)
    if 'rqmc_confidence.py' not in protocol['source_hashes']:
        raise ValueError('Confidence source must exist before freezing production protocol')
    models={m['model_id']:m for m in protocol['models']}
    results=[]
    for case_index, selected in enumerate(protocol['selected']):
        row=selected['row']; model=models[row['model_id']]
        graph=build(model,row); evaluator=FiniteEvaluator(graph,bank)
        points=digital_points(graph.inputs)
        rng=np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence([PRODUCTION_SEED,case_index])))
        rho=Fraction(selected['finite_output_range_exact_fraction'])
        epsilon=Fraction(selected['estimation_tolerance_exact_fraction'])
        retained_sums=[]; shift_hash=hashlib.sha256(); output_hash=hashlib.sha256()
        checkpoints=[]; replay=[]
        start,cpu=time.perf_counter(),time.process_time()
        generation_time=evaluation_time=reduction_time=confidence_time=0.0
        previous=0
        for target in CHECKPOINTS:
            remaining=target-previous
            while remaining:
                reps=min(remaining,MAX_FINANCIAL_LANES_PER_CHUNK//POINTS_PER_SHIFT)
                t=time.perf_counter(); shifts=draw_shifts(graph.inputs,reps,rng)
                inputs=shifted_inputs(points,shifts);generation_time+=time.perf_counter()-t
                for name, a in shifts.items():
                    shift_hash.update(name.encode());shift_hash.update(a.astype('<u8',copy=False).tobytes())
                t=time.perf_counter();outputs=evaluator.evaluate(inputs);evaluation_time+=time.perf_counter()-t
                if np.any(outputs > selected['range_certificate']['max_integer']):
                    raise ValueError('Financial output violates its certified finite range')
                output_hash.update(outputs.astype('<u8',copy=False).tobytes())
                t=time.perf_counter()
                batch_sums=np.sum(outputs.reshape(reps,POINTS_PER_SHIFT),axis=1,dtype=np.uint64)
                # 64 summands each <=2^30 give a sum <=2^36, exactly safe.
                retained_sums.extend(map(int,batch_sums))
                reduction_time+=time.perf_counter()-t
                for local_rep in range(min(8,reps)) if previous==0 else ():
                    for point_index in (0,1,2,3,7,15,31,63):
                        lane=local_rep*POINTS_PER_SHIFT+point_index
                        replay.append({'replicate_index':previous+local_rep,'point_index':point_index,
                            'inputs':{name:int(a[lane]) for name,a in inputs.items()},'integer_output':int(outputs[lane])})
                if remaining==reps:
                    local_rep=reps-1
                    for point_index in (0,1,2,3,7,15,31,63):
                        lane=local_rep*POINTS_PER_SHIFT+point_index
                        replay.append({'replicate_index':target-1,'point_index':point_index,
                            'inputs':{name:int(a[lane]) for name,a in inputs.items()},'integer_output':int(outputs[lane])})
                previous+=reps;remaining-=reps
            t=time.perf_counter()
            stats=batch_statistics(retained_sums,row['q'])
            confidence=confidence_from_batch_statistics(stats,rho=rho,
                alpha_variance=ALPHA_EVENT,alpha_mean=ALPHA_EVENT,
                variance_prior=Fraction(selected['batch_variance_prior_exact_fraction']),
                mean_prior=Fraction(selected['batch_mean_prior_exact_fraction']))
            confidence_time+=time.perf_counter()-t
            radius=(Fraction(confidence['halfwidth_upper_fraction'])
                    if confidence['status']=='nonempty' else None)
            checkpoint={'replicates':target,'exact_batch_statistics':stats,'confidence':confidence,
                        'tolerance_met':radius is not None and radius<=epsilon}
            checkpoints.append(checkpoint)
            print('RQMC checkpoint',row['case_id'],target,'radius',float(radius) if radius is not None else None,
                  'target',float(epsilon),'status',confidence['status'],flush=True)
            if radius is None or radius<=epsilon:
                break
        count=len(retained_sums); last=checkpoints[-1]
        if last['confidence']['status']=='nonempty':
            interval=price_interval(model,row,last['confidence']['mean_lower_fraction'],
                                    last['confidence']['mean_upper_fraction'])
            status='certified per-price tolerance met' if last['tolerance_met'] else 'prespecified maximum reached without requested certificate'
        else:
            interval=None
            status='explicit empty confidence-intersection outcome'
        keys={(r['kernel'],r['address_bits'],r['q']) for r in evaluator.table_manifest}
        result={'case_id':row['case_id'],'model_id':row['model_id'],'method':row['method'],'key':selected['key'],
            'assets':model['assets'],'outer_dimension':model['outer_dimension'],'status':status,
            'seed_sequence':[PRODUCTION_SEED,case_index],'replicates':count,
            'points_per_shift':POINTS_PER_SHIFT,'financial_evaluations':count*POINTS_PER_SHIFT,
            'batch_integer_sums':[str(x) for x in retained_sums],
            'shift_stream_sha256':shift_hash.hexdigest(),'output_stream_sha256':output_hash.hexdigest(),
            'checkpoints':checkpoints,'confidence_interval':interval,'independent_scalar_replay_records':replay,
            'timing':{'execution_wall_seconds':time.perf_counter()-start,
                'execution_process_cpu_seconds':time.process_time()-cpu,
                'shift_and_point_generation_wall_seconds':generation_time,
                'finite_graph_evaluation_wall_seconds':evaluation_time,
                'batch_reduction_wall_seconds':reduction_time,'confidence_construction_wall_seconds':confidence_time,
                'evaluator_preparation_wall_seconds':evaluator.preparation_seconds,
                'standalone_table_setup_wall_seconds_reconstructed':sum(table_records[k]['first_generation_wall_seconds'] for k in keys)},
            'work':{'financial_evaluations':count*POINTS_PER_SHIFT,
                'evaluator_scalar_lane_operations_upper':count*POINTS_PER_SHIFT*evaluator.operation_profile()['scalar_lane_operation_upper_total'],
                'random_shift_word_draws':count*len(graph.inputs),
                'random_shift_bits_requested':count*sum(graph.inputs.values()),
                'shifted_coordinate_XOR_operations':count*POINTS_PER_SHIFT*len(graph.inputs),
                'digital_point_coordinate_memory_reads':count*POINTS_PER_SHIFT*len(graph.inputs),
                'financial_batch_sum_additions':count*(POINTS_PER_SHIFT-1),
                'checkpoint_moment_reduction_input_words':sum(c['replicates'] for c in checkpoints),
                'checkpoint_pair_difference_input_pairs':sum(c['replicates']//2 for c in checkpoints),
                'prepared_table_bytes':evaluator.operation_profile()['prepared_table_bytes']},
            'qualification':'Independent uniform digital shifts are the mathematical premise; financial points within a shift are dependent; PCG seed and empirical replicate variation do not establish coverage or a 1/N rate'}
        results.append(result)
        (ROOT/'results/rqmc_pricing_checkpoint.json').write_text(json.dumps({'status':'in progress',
            'protocol_sha256':protocol_sha,'completed':len(results),'configurations':results},indent=2)+'\n')
    if digest(protocol_path)!=protocol_sha or any(digest(ROOT/name)!=sha for name,sha in protocol['source_hashes'].items()):
        raise ValueError('Frozen protocol or source changed during production')
    output={'status':'prespecified randomised digital-shift finite-oracle experiment complete',
        'protocol_sha256':protocol_sha,'source_hashes':protocol['source_hashes'],'configurations':results,
        'summary':{'prices':len(results),'higher_dimensional_prices':sum(r['assets']>2 for r in results),
            'two_asset_controls':sum(r['assets']==2 for r in results),
            'certified_tolerance_successes':sum(r['status']=='certified per-price tolerance met' for r in results),
            'financial_evaluations':sum(r['financial_evaluations'] for r in results),
            'random_shifts':sum(r['replicates'] for r in results),
            'execution_wall_seconds':sum(r['timing']['execution_wall_seconds'] for r in results)},
        'coverage_scope':'99 per cent per selected price under ideal independent exact-uniform shift vectors, simultaneous over that price\'s eleven checkpoints; not99 per cent for all27 prices jointly',
        'scope':'Synthetic finite financial functions, complete deterministic bias charged; no actual market sample, novelty or physical quantum-runtime assertion'}
    (ROOT/'results/rqmc_pricing_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['summary'],indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--omit-two-asset-controls',action='store_true')
    args=parser.parse_args()
    protocol,bank,records=prepare_protocol(not args.omit_two_asset_controls)
    print('Frozen digital-shift protocol',len(protocol['selected']),'selected prices',flush=True)
    if not args.prepare_only:
        run_production(protocol,bank,records)


if __name__=='__main__':
    main()
