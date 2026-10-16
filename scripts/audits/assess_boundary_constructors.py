"""All36 financial cases, with no continuous-transfer claim for high dimension."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from generalisation_models import models
import boundary_value_oracle as o
import boundary_value_certificate as c

ROOT=Path(__file__).resolve().parents[1]

def main():
    cache={};rows=[]
    for model in models():
        domain=c.finite_input_domain(model,6,10,28);M=max(abs(F(x)) for x in domain['coordinate_bounds'])
        for payoff,cap in [('capped',1),('capped',10),('capped',25),('vanilla',None)]:
            for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
                graph=o.build_graph(model,cap,kind=kind,payoff=payoff,U_star=domain['U_star'],retained_bound=M)
                tables,sc=o.make_tables(graph,cache,True);certificate=c.certify_graph(graph,tables,sc)
                rows.append({'model_id':model['model_id'],'payoff':payoff,'cap':cap,'kind':kind,
                    'outer_dimension':model['outer_dimension'],'configuration':graph.configuration,
                    'node_count':len(graph.nodes),'word_width':graph.width,'input_bits':sum(graph.inputs.values()),
                    'scalar_lookup_nodes':sum(n.kind=='function' for n in graph.nodes.values()),
                    'arithmetic_status':certificate['status'],'monetary_value_error_upper':certificate['monetary_value_error_upper'],
                    'scope':'Finite value arithmetic certificate only; no continuous-target transfer, complete estimator, classical or quantum advantage claim.'})
        print(model['model_id'],'retained',len(rows),flush=True)
    result={'schema':'qbp-boundary-constructor-family-v1','status':'passed' if all(r['arithmetic_status']=='passed' for r in rows) else 'contains rejected arithmetic cases',
        'financial_cases':36,'method_graphs':len(rows),'records':rows,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['boundary_value_oracle.py','boundary_value_certificate.py','boundary_derivative_certificate.py','results/boundary_derivative_norms.json','generalisation_models.py','results/generalisation_protocol.json']}}
    with gzip.open(ROOT/'results/boundary_constructor_assessment.json.gz','wt') as f:json.dump(result,f,separators=(',',':'));f.write('\n')
    print('Retained',len(rows),'graphs, maximum arithmetic money',max(float(F(r['monetary_value_error_upper'])) for r in rows),flush=True)

if __name__=='__main__':main()
