"""Prospective complete n2 financial profiles, retaining every variant.

These emitted circuits are input to the independent monetary certificate and
complete estimator planner. A word width alone never implies price accuracy.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,time
import integral_value_oracle as oracle
ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--table-bits',type=int,default=12)
    parser.add_argument('--all-models',action='store_true');parser.add_argument('--output',type=Path,default=ROOT/'results/integral_production_profiles.json')
    args=parser.parse_args();start=time.monotonic()
    source=ROOT/'integral_value_oracle.py';source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    protocol=json.loads((ROOT/'results/integral_oracle_protocol.json').read_text())
    family=json.loads((ROOT/'results/generalisation_protocol.json').read_text())
    domains={row['model_id']:row for row in json.loads((ROOT/'results/integral_finite_input_domains.json').read_text())['domains']}
    norm=json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())['records']
    D={(F(row['lambda']),row['derivative_order']):F(row['norm_upper']) for row in norm}
    records=[];cache={};table_cache={};graphs=0
    for model in family['models']:
        if not args.all_models and model['assets']!=2:continue
        domain=domains[model['model_id']];lam=F(model['lambda']);M=max(abs(F(x)) for x in domain['coordinate_bounds'])
        for cap in protocol['caps']:
            for kind in protocol['methods']:
                for encoding in (['centred'] if kind=='positive_r1' else ['centred','positive_part','negative_part']):
                    graph=oracle.build_graph(model,cap,q=28,b_y=10,beta_bits=8,strike_bits=8,
                         table_bits=args.table_bits,outer_bound=6,U_star=F(domain['U_star']),D2=D[lam,2],D1=D[lam,1],
                         kind=kind,retained_bound=M,output_encoding=encoding)
                    tables=oracle.make_tables(graph,table_cache);graphs+=1
                    for mode in protocol['coherent_lookup_architectures']:
                        profile=oracle.counted_composition(graph,tables,cache,mode)
                        profile.update(case_id=f"{model['model_id']}_cap{cap}",accuracy_status='Requires independent complete monetary certificate',
                             graph_nodes=len(graph.nodes),scalar_specs=[node.meta.as_dict() for node in graph.nodes.values() if node.kind=='function'])
                        records.append(profile)
                    print(model['model_id'],cap,kind,encoding,'graphs',graphs,'seconds',round(time.monotonic()-start,2),flush=True)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash,'Financial source changed during profile generation; rerun after source freeze'
    result={'schema':'qbp-integral-production-profiles-v1','status':'complete financial gate profiles;monetary accuracy separately certified',
      'primary_table_bits':args.table_bits,'profiles':records,'graph_variants':graphs,
      'architectures':protocol['coherent_lookup_architectures'],'source_sha256':source_hash,
      'driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'input_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
          [ROOT/'results/integral_oracle_protocol.json',ROOT/'results/generalisation_protocol.json',ROOT/'results/integral_finite_input_domains.json',ROOT/'results/gaussian_derivative_norms.json']},
      'dependency_source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['value_oracle.py','architecture.py','reversible.py','quantum_gates.py','integral_value_certificate.py']},
      'elapsed_seconds_diagnostic':time.monotonic()-start,'scope':'Matched finite support and architectures across integralR2/directR2,centred/split forms,positiveR1. No quantum execution or price-accuracy claim is inferred.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'profiles':len(records),'graph_variants':graphs,'seconds':result['elapsed_seconds_diagnostic'],'output':str(args.output)}))


if __name__=='__main__':main()
