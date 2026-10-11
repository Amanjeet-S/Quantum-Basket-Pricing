"""Complete classical pricing of the same certified finite financial targets.

Precision choice uses deterministic bounds and an operation-count objective.
Fresh samples never influence the chosen prefix or sample size. Timings have
their own units and do not imply physical quantum advantage.
"""
from collections import Counter
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import numpy as np
import generalisation as gen
import general_kernels as kernels
import encoding_resources as er
from classical_engine import FiniteEvaluator, ExactMoments

ROOT=Path(__file__).resolve().parent
D=Decimal
SEED=2026100503
BATCH_SIZE=65536
METHODS=('direct','interval','full_global','full_cap')
FIELDS=('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits')


def row_key(row):
    return tuple(row[k] for k in FIELDS)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(model,row):
    return gen.build_graph(model,row['method'],*[row[k] for k in
        ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],row['cap'])


def prepare_tables(graph,bank,records):
    for node in graph.nodes.values():
        if node.kind not in ('quantile','function'):continue
        key=('quantile' if node.kind=='quantile' else node.meta,node.bits,graph.q)
        if key in bank:continue
        start,cpustart=time.perf_counter(),time.process_time()
        bank[key]=(er.quantile_table(node.bits,graph.q)[0] if node.kind=='quantile'
            else kernels.function_table(node.meta,node.bits,graph.q,graph.width))
        records[key]={'kernel':key[0],'address_bits':key[1],'q':key[2],
            'entries':len(bank[key]),'sha256':hashlib.sha256(json.dumps(bank[key],separators=(',',':')).encode()).hexdigest(),
            'first_generation_wall_seconds':time.perf_counter()-start,
            'first_generation_process_cpu_seconds':time.process_time()-cpustart}


def price_interval(model,row,moments):
    """Outward price interval from exact finite sum and proved tolerance.

    Coverage concerns ideal independent uniform bits. Decimal elementary
    functions follow the explicitly stated correctly rounded premise.
    """
    with localcontext() as ctx:
        ctx.prec=220
        count=moments['count'];q=row['q'];total=int(moments['integer_sum'])
        denominator=count*(1 << q)
        ctx.rounding=ROUND_FLOOR;mean_lo=D(total)/D(denominator)
        ctx.rounding=ROUND_CEILING;mean_hi=D(total)/D(denominator)
        correction=gen.geometric_correction(model,row['cap'])
        c=D(correction['undiscounted_value_decimal']);c_error=max(D('1e-12'),D.from_float(correction['absolute_allowance']))
        arithmetic=D.from_float(row['arithmetic_price_error_upper'])
        distribution=D.from_float(row['distribution_price_bound'])
        statistical=D('.005');pad=D('1e-100')
        assert arithmetic+distribution+c_error<=D('.005')
        ctx.rounding=ROUND_FLOOR
        low=c-c_error+D(row['cap'])*mean_lo-arithmetic-distribution-statistical-pad
        ctx.rounding=ROUND_CEILING
        high=c+c_error+D(row['cap'])*mean_hi+arithmetic+distribution+statistical+pad
        # The payoff itself is in [0,H]. Intersect before applying discount.
        low=max(D(0),low);high=min(D(row['cap']),high)
        assert low<=high
        ctx.rounding=ROUND_FLOOR
        disc=(-D(model['risk_free_rate'])).exp()
        disc_lo=disc-pad;disc_hi=disc+pad
        plo=low*disc_lo-pad
        ctx.rounding=ROUND_CEILING;phi=high*disc_hi+pad
        # An exact midpoint of the enclosing endpoints delivers radius<.01.
        midpoint=(plo+phi)/2
        assert phi-plo<=D('.02')
        return {'normalised_sample_mean_lower_decimal':str(mean_lo),
            'normalised_sample_mean_upper_decimal':str(mean_hi),
            'geometric_correction':correction,
            'price_interval_lower_decimal':str(plo),'price_interval_upper_decimal':str(phi),
            'display_midpoint_decimal':str(midpoint),
            'maximum_radius_decimal':str((phi-plo)/2+pad),
            'statistical_currency_allowance_exact':'0.005',
            'deterministic_currency_allowance_decimal':str(arithmetic+distribution+c_error+pad),
            'numerical_allowance_decimal':str(pad),
            'numerical_premise':'correctly rounded Decimal division/exp at220 digits; geometric correction centre with retained1e-12 allowance'}


def main():
    start=time.perf_counter();config_path=ROOT/'results/architecture_resource_candidates.json'
    variance_path=ROOT/'results/classical_variance_certificates.json'
    data=json.loads(config_path.read_text());variance=json.loads(variance_path.read_text())
    assert variance['configuration_sha256']==digest(config_path)
    for f,h in variance['source_hashes'].items():assert digest(ROOT/f)==h,f
    models={m['model_id']:m for m in data['models']}
    original={row_key(r):r for r in data['configurations']}
    assert len(original)==273
    certs={row_key(r):r for r in variance['configurations']}
    assert set(original)==set(certs) and len(certs)==len(variance['configurations'])
    sources=('classical_pricing_study.py','classical_engine.py','variance_certificate.py','generalisation.py','general_kernels.py','encoding_resources.py')
    hashes={f:digest(ROOT/f) for f in sources}
    bank={};table_records={};profiles=[]
    for key,row in sorted(original.items()):
        graph=build(models[row['model_id']],row)
        prepare_tables(graph,bank,table_records)
        evaluator=FiniteEvaluator(graph,bank);profile=evaluator.operation_profile()
        certificate=certs[key]
        n=certificate['confidence_budgets']['single_price']['certified_samples']
        objective=n*profile['scalar_lane_operation_upper_total']
        # Keep setup separate: Decimal endpoint generation and NumPy lanes
        # have different abstract operation costs.
        profiles.append({'key':list(key),'row':{k:v for k,v in row.items() if k!='resources' and k!='architecture'},
            'certificate':certificate,'evaluator_profile':profile,'certified_samples':n,
            'sample_evaluator_scalar_lane_work_upper':objective,
            'selection_objective':'minimum certified sample count times conservative scalar-lane evaluator operations within finite inherited pool; not wall-clock or global optimisation'})
    chosen=[]
    for cid in sorted({r['case_id'] for r in original.values()}):
        for method in METHODS:
            pool=[p for p in profiles if p['row']['case_id']==cid and p['row']['method']==method]
            minimum=min(p['sample_evaluator_scalar_lane_work_upper'] for p in pool)
            ties=[p for p in pool if p['sample_evaluator_scalar_lane_work_upper']==minimum]
            winner=min(ties,key=lambda p:(p['evaluator_profile']['prepared_table_bytes'],tuple(p['key']))).copy()
            winner['tied_prefixes']=[p['key'] for p in ties]
            chosen.append(winner)
    assert len(chosen)==108
    protocol={'status':'frozen deterministic complete classical sampling protocol',
        'configuration_sha256':digest(config_path),'variance_certificate_sha256':digest(variance_path),
        'source_hashes':hashes,'seed':SEED,'batch_size':BATCH_SIZE,
        'random_generator':'NumPy PCG64DXSM with SeedSequence([seed,case_index,method_index])',
        'randomness_premise':'ideal independent uniform finite addresses for the mathematical guarantee; pseudorandom execution is reproducible and does not prove this premise',
        'confidence_scope':'99per cent per price; not all108 outputs jointly99per cent',
        'models':data['models'],'all_candidate_profiles':profiles,'selected':chosen,
        'first_generation_table_records':list(table_records.values()),
        'selection_and_all_candidate_setup_wall_seconds':time.perf_counter()-start,
        'setup_scope':'cache generation once per unique table; chosen standalone generation cost reconstructed from first-generation records; neither reconstructed setup nor scalar-lane work is a physical quantum runtime',
        'numpy_version':np.__version__,'python_version':sys.version,'platform':platform.platform()}
    protocol_path=ROOT/'results/classical_pricing_protocol.json'
    protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
    protocol_sha=digest(protocol_path);results=[]
    cases=sorted({p['row']['case_id'] for p in chosen})
    for ci,cid in enumerate(cases):
        for mi,method in enumerate(METHODS):
            selected=next(p for p in chosen if p['row']['case_id']==cid and p['row']['method']==method)
            row=selected['row'];model=models[row['model_id']]
            graph=build(model,row);evaluator=FiniteEvaluator(graph,bank)
            rng=np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence([SEED,ci,mi])))
            moments=ExactMoments(row['q']);input_hash=hashlib.sha256();output_hash=hashlib.sha256()
            generation_time=evaluation_time=reduction_time=0.0
            wall,cpu=time.perf_counter(),time.process_time()
            remaining=selected['certified_samples'];chunks=0
            while remaining:
                n=min(BATCH_SIZE,remaining);t=time.perf_counter()
                inputs={name:rng.integers(0,1 << bits,size=n,dtype=np.uint64) for name,bits in graph.inputs.items()}
                generation_time+=time.perf_counter()-t
                for name,array in inputs.items():
                    input_hash.update(name.encode());input_hash.update(array.astype('<u8',copy=False).tobytes())
                t=time.perf_counter();out=evaluator.evaluate(inputs);evaluation_time+=time.perf_counter()-t
                output_hash.update(out.astype('<u8',copy=False).tobytes())
                t=time.perf_counter();moments.add(out);reduction_time+=time.perf_counter()-t
                remaining-=n;chunks+=1
            executed=moments.as_dict();interval=price_interval(model,row,executed)
            keys={(x['kernel'],x['address_bits'],x['q']) for x in evaluator.table_manifest}
            standalone_setup=sum(table_records[k]['first_generation_wall_seconds'] for k in keys)
            result={'case_id':cid,'method':method,'key':selected['key'],'seed_sequence':[SEED,ci,mi],
                'samples':selected['certified_samples'],'chunks':chunks,'input_stream_sha256':input_hash.hexdigest(),
                'output_stream_sha256':output_hash.hexdigest(),'exact_moments':executed,'confidence_interval':interval,
                'confidence_failure_upper_decimal':selected['certificate']['confidence_budgets']['single_price']['certified_failure_upper_decimal'],
                'timing':{'sampling_execution_wall_seconds':time.perf_counter()-wall,
                    'sampling_execution_process_cpu_seconds':time.process_time()-cpu,
                    'random_address_generation_wall_seconds':generation_time,'integer_evaluation_wall_seconds':evaluation_time,
                    'exact_moment_reduction_wall_seconds':reduction_time,
                    'chosen_table_generation_wall_seconds_reconstructed':standalone_setup,
                    'chosen_evaluator_preparation_wall_seconds':evaluator.preparation_seconds},
                'work':{'financial_evaluations':executed['count'],
                    'evaluator_scalar_lane_operations_upper':executed['count']*evaluator.operation_profile()['scalar_lane_operation_upper_total'],
                    'random_finite_address_draws':executed['count']*len(graph.inputs),
                    'random_address_bits_requested':executed['count']*sum(graph.inputs.values()),
                    'exact_moment_reduction_samples':executed['count'],'prepared_table_bytes':evaluator.operation_profile()['prepared_table_bytes']},
                'qualification':'certified finite-output sampling under idealIID-bit premise; actual NumPy pseudorandom stream; no physical quantum or market execution'}
            results.append(result)
        print('Completed classical prices',cid,'all four methods',flush=True)
        checkpoint={'status':'in progress','completed':len(results),'protocol_sha256':protocol_sha,'configurations':results}
        (ROOT/'results/classical_pricing_checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
    assert {f:digest(ROOT/f) for f in sources}==hashes
    assert digest(protocol_path)==protocol_sha
    output={'status':'108 complete classical finite-output prices executed',
        'protocol_sha256':protocol_sha,'source_hashes':hashes,'configurations':results,
        'summary':{'prices':len(results),'samples':sum(r['samples'] for r in results),
            'sampling_execution_wall_seconds':sum(r['timing']['sampling_execution_wall_seconds'] for r in results),
            'all_sample_price_interval_radii_at_most_0.01':all(D(r['confidence_interval']['maximum_radius_decimal'])<=D('.01') for r in results)},
        'scope':'per-price mathematical99per cent underidealindependentfiniteuniformbits; exact integer financial arithmetic; runtime diagnostics only; no novelty or quantum speedup claim'}
    (ROOT/'results/classical_pricing_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['summary'],indent=2))


if __name__=='__main__':main()
