"""Fair finite opponent-inherited precision pool with independent budgets.

I retain the historical selected prefixes and every prefix actually counted
in their original shortlists. All conditional methods receive the union of
interval, generic full and cap-aware full prefixes, including the matched
interval diagnostic. Each transfer is rebuilt and certified for its target
method. Rejected transfers remain visible. This is not global optimisation.
"""
from decimal import Decimal,localcontext
from pathlib import Path
import hashlib
import json
import sys
import time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'checks'))
import generalisation as gen
from resource_study import propagate
from audit_generalisation_certificates import recurrence,distribution
from classical_reference_generalisation import exact_model_checks

D=Decimal


def prefix(row):
    return tuple(row[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits'))


def case_prefixes(rows,matched=None):
    conditional={};direct={}
    def add(target,key,source):target.setdefault(key,[]).append(source)
    for row in rows:
        dest=direct if row['method']=='direct' else conditional
        add(dest,prefix(row),f"{row['method']}:historical_selected")
        for i,candidate in enumerate(row['exact_counted_candidates']):
            add(dest,(candidate['q'],*candidate['function_address_bits']),f"{row['method']}:historical_counted_{i}")
    if matched is not None:add(conditional,prefix(matched),'interval:post_output_matched_prefix')
    # Direct can inherit the exp/q part of opponents without inventing a
    # log/probability function that direct pricing does not contain.
    for (q,e,_,_),provenance in conditional.items():
        for source in provenance:add(direct,(q,e,0,0),'exp/q_only_from_'+source)
    return conditional,direct


def build_candidate_graph(model,row):
    return gen.build_graph(model,row['method'],*[row[k] for k in
        ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],row['cap'])


def certify_pool(model,cap,epsilon,rows,matched=None):
    conditional,direct=case_prefixes(rows,matched);accepted=[];rejected=[]
    with localcontext() as ctx:
        ctx.prec=100
        for method in ('direct','interval','full_global','full_cap'):
            pool=direct if method=='direct' else conditional
            for pp,provenance in sorted(pool.items()):
                q,e,l,p=pp;record={'model_id':model['model_id'],'case_id':f"{model['model_id']}_cap{cap}",
                    'method':method,'cap':cap,'price_tolerance':epsilon,'q':q,
                    'exp_address_bits':e,'log_address_bits':l,'probability_address_bits':p,
                    'inherited_from':sorted(set(provenance))}
                dist=gen.distribution(model,method,cap,q,epsilon)
                if dist is None:
                    rejected.append({**record,'reason':'distribution budget infeasible'});continue
                record.update(dist)
                try:
                    graph=build_candidate_graph(model,record);cert=propagate(graph,cap)
                    independently,margin,domaincount=recurrence(graph,cap)
                    assert independently<=D(str(cert['arithmetic_price_error_upper']))+D('1e-30')
                    outer,inner=distribution(model,record)
                    assert outer<=D(str(record['outer_distribution_bound']))+D('1e-30')
                    assert inner<=D(str(record['inner_distribution_bound']))+D('1e-30')
                except (AssertionError,ValueError) as exc:
                    rejected.append({**record,'reason':'domain or independent certificate rejected','diagnostic':str(exc)});continue
                total=independently+outer+inner+D('1e-12')
                if independently>D(str(epsilon))*D('.1') or total>D(str(epsilon))*D('.5'):
                    rejected.append({**record,'reason':'target-method monetary budget exceeded',
                        'independent_arithmetic_upper_decimal':str(independently),
                        'independent_deterministic_bias_upper_decimal':str(total)});continue
                accepted.append({**record,'arithmetic_price_error_upper':cert['arithmetic_price_error_upper'],
                    'independent_arithmetic_upper_decimal':str(independently),
                    'independent_deterministic_bias_upper_decimal':str(total),
                    'minimum_function_domain_margin_diagnostic':margin,'function_domain_count':domaincount,
                    'reserved_statistical_price_allowance':epsilon*.5,
                    'selection_scope':'finite inherited pool; every target method is recertified; not globally optimal'})
    # Same graph and wider valid financial certificate imply generic-feasible
    # prefixes must also remain cap-aware-feasible. Check this constructive
    # containment before allowing a resource selector to choose a winner.
    gf={prefix(row) for row in accepted if row['method']=='full_global'}
    cf={prefix(row) for row in accepted if row['method']=='full_cap'}
    assert gf<=cf
    assert all(any(row['method']==method for row in accepted) for method in ('direct','interval','full_global','full_cap'))
    return accepted,rejected


def generate_candidates(data,matched_data=None):
    accepted=[];rejected=[];metrics=[exact_model_checks(model) for model in data['models']]
    for model in data['models']:
        caps=sorted({r['cap'] for r in data['configurations'] if r['model_id']==model['model_id']})
        for cap in caps:
            rows=[r for r in data['configurations'] if r['model_id']==model['model_id'] and r['cap']==cap]
            cid=rows[0]['case_id'];matched=None
            if matched_data is not None:matched=next((r for r in matched_data['configurations'] if r['case_id']==cid),None)
            aa,rr=certify_pool(model,cap,rows[0]['price_tolerance'],rows,matched)
            accepted+=aa;rejected+=rr
            print('Inherited pool certified',cid,len(aa),'accepted',len(rr),'rejected',flush=True)
    return accepted,rejected,metrics


def main():
    start=time.perf_counter();path=ROOT/'results/generalisation_resources.json'
    data=json.loads(path.read_text());mp=ROOT/'results/generalisation_matched_interval.json'
    matched=json.loads(mp.read_text()) if mp.exists() else None
    snapshot=json.loads((ROOT/'results/base_snapshot.json').read_text())
    before={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in snapshot['source_files']};assert before==snapshot['source_files']
    aa,rr,metrics=generate_candidates(data,matched)
    after={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in before};assert before==after
    out={'status':'independently certified finite opponent-inherited candidate pool',
        'historical_resource_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'historical_matched_sha256':hashlib.sha256(mp.read_bytes()).hexdigest() if mp.exists() else None,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'base_source_hashes':before,'model_checks':metrics,'models':data['models'],
        'feasible_candidates':aa,'rejected_candidates':rr,
        'summary':{'accepted':len(aa),'rejected':len(rr),
            'accepted_by_method':{m:sum(r['method']==m for r in aa) for m in ('direct','interval','full_global','full_cap')}},
        'scope':'union of historical selected and actually counted shortlist prefixes plus every matched interval prefix; all conditional opponents receive the same offered pool; direct receives exp/q transfers; finite-pool optimisation only',
        'elapsed_seconds':time.perf_counter()-start}
    output=ROOT/'results/inherited_precision_candidates.json';output.write_text(json.dumps(out,indent=2)+'\n');print(output)


if __name__=='__main__':main()
