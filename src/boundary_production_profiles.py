"""Retain the complete common-table seed family and complete circuit charges.

All n2 payoff cases are retained. Higher-dimensional graphs remain executable,
but complete financial-transfer certificates are separately required.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import gzip,hashlib,json,time
import boundary_value_oracle as o
import boundary_value_certificate as c
from generalisation_models import models

ROOT=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(output='results/boundary_production_profiles.json.gz'):
    start=time.perf_counter();table_cache={};block_cache={};profiles=[];certificates=[];scalar_archive={};graphs=[]
    dependencies=['boundary_value_oracle.py','boundary_value_certificate.py','boundary_derivative_certificate.py',
        'integral_value_oracle.py','integral_value_certificate.py','value_oracle.py','architecture.py',
        'quantum_gates.py','reversible.py','generalisation_models.py','results/boundary_derivative_norms.json',
        'results/generalisation_protocol.json','results/integral_finite_input_domains.json']
    provenance={name:sha(ROOT/name) for name in dependencies}
    for model in models():
        if model['assets']!=2:continue
        domain=c.finite_input_domain(model,6,10,28);M=max(abs(F(x)) for x in domain['coordinate_bounds'])
        for payoff,cap in [('capped',1),('capped',10),('capped',25),('vanilla',None)]:
            for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
                for encoding in (('centred','positive_part','negative_part') if payoff=='capped' and kind!='positive_r1' else ('centred',)):
                    graph=o.build_graph(model,cap,kind=kind,payoff=payoff,U_star=domain['U_star'],retained_bound=M,output_encoding=encoding)
                    tables,scalars=o.make_tables(graph,table_cache,True)
                    certified=c.certify_graph(graph,tables,scalars)
                    key=f"{model['model_id']}:{payoff}:{cap}:{kind}:{encoding}"
                    certified['record_id']=key;certificates.append(certified)
                    nodes=[]
                    for name,node in graph.nodes.items():
                        nodes.append({'name':name,'kind':node.kind,'args':list(node.args),
                            'meta':node.meta.as_dict() if isinstance(node.meta,o.ScalarSpec) else str(node.meta),
                            'bits':node.bits,'output_bits':node.output_bits,'bounds':list(map(str,graph.bounds[name]))})
                        if name in scalars:
                            record=scalars[name];digest=hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
                            scalar_archive.setdefault(digest,record)
                    graphs.append({'record_id':key,'configuration':graph.configuration,'width':graph.width,'inputs':graph.inputs,
                        'output':graph.output,'semantic_nodes':graph.semantic_nodes,'nodes':nodes})
                    for mode in ('exact7t','unitary8t'):
                        profile=o.counted_composition(graph,tables,block_cache,mode);profile['record_id']=key+':'+mode;profiles.append(profile)
                    print(key,'W',graph.width,'money',float(F(certified['monetary_value_error_upper'])),
                          'T',profiles[-1]['logical_T_count'],flush=True)
    now={name:sha(ROOT/name) for name in dependencies}
    if now!=provenance:raise RuntimeError('A dependency changed during production; rerun under one source snapshot')
    result={'schema':'qbp-boundary-production-family-v1','status':'passed' if all(x['status']=='passed' for x in certificates) else 'contains rejected cases',
        'created_utc':datetime.now(timezone.utc).isoformat(),'dependency_sha256':provenance,
        'source_sha256':sha(__file__),'profiles':profiles,'graph_certificates':certificates,'graphs':graphs,
        'scalar_certificate_archive':scalar_archive,'unique_scalar_certificates':len(scalar_archive),
        'unique_counted_blocks':len(block_cache),'profile_count':len(profiles),'financial_graph_count':len(graphs),
        'configuration':'All nine n2 capped cases and all three n2 vanilla cases, all four methods and all capped R2 flag encodings, both coherent lookup architectures, q28/b_y10/Beta8/table12.',
        'scope':'Complete emitted financial A preparations and arithmetic certificates. Continuous transfer, coefficients, auxiliary discrepancy and complete mean estimation are separate required certificates.',
        'runtime_seconds_diagnostic':time.perf_counter()-start}
    target=ROOT/output
    with gzip.open(target,'wt') as f:json.dump(result,f,separators=(',',':'));f.write('\n')
    summary={k:v for k,v in result.items() if k not in ('profiles','graph_certificates','graphs','scalar_certificate_archive')}
    summary['archive_sha256']=sha(target)
    (ROOT/'results/boundary_production_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Retained',len(profiles),'profiles;',len(graphs),'financial graphs;',len(block_cache),'unique blocks',flush=True)
    return result

if __name__=='__main__':run()
