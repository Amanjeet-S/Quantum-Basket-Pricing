"""Registered exact-factor-basis comparison of complete financial quantum plans.

Derived from the retained four-survivor producer without changing that source.
Identity and the analytically selected exact rational rotation are rebuilt on
one actual rounded Gaussian law. All failed scalar allocations are retained.
No price answer informs a quantum prior or factor-basis selection.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import gzip, hashlib, json, sys, time
import argparse
_BOOTSTRAP = argparse.ArgumentParser(add_help=False)
_BOOTSTRAP.add_argument('--repo', required=True)
_BOOTSTRAP.add_argument('--output', required=True)
_BOOTSTRAP.add_argument('--certificate', required=True)
_BOOTSTRAP_ARGS, _ = _BOOTSTRAP.parse_known_args()
sys.path.insert(0, str(Path(_BOOTSTRAP_ARGS.repo).resolve()))
from flint import arb, ctx
import higherdim_preparation as prep
import tensor_transfer as tt
import tensor_first_transfer as first
import boundary_value_oracle as oracle
import boundary_value_certificate as certifier
import boundary_precision_selection as retained
import correction_information_priors as informed
from boundary_grid_certificate import parameters, ball, serial, upper, endpoints, coefficient_certificate
from boundary_derivative_certificate import payoff_norms

ROOT = Path(_BOOTSTRAP_ARGS.repo).resolve()
RESULTS = Path(_BOOTSTRAP_ARGS.output).resolve()
RESULTS.mkdir(parents=True, exist_ok=True)
CERTIFICATE = Path(_BOOTSTRAP_ARGS.certificate).resolve()
PROTOCOL = RESULTS / 'factor_basis_quantum_protocol.json'
TARGETS = ('n2_lambda0.4_scale3', 'n4_lambda0.05_scale3/2',
           'n4_lambda0.2_scale3', 'n8_lambda0.2_scale3')
KINDS = retained.KINDS
ARCHITECTURES = ('exact7t', 'unitary8t')
EPS = F(1, 100)
ALPHA = F(1, 3600)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, default=str, indent=2) + '\n')

def source_hashes():
    return {p.name: sha(p) for p in sorted(ROOT.glob('*.py'))}

def input_hashes():
    names = ('dispersion_separation_protocol.json', 'generalisation_protocol.json',
             'boundary_derivative_norms.json', 'gaussian_derivative_norms.json',
             'qft_rotation_tapes.json', 'qft_rotation_certificates.json')
    return {n: sha(ROOT / 'results' / n) for n in names}

def selected_models():
    models = json.loads((ROOT / 'results' / 'dispersion_separation_protocol.json').read_text())['models']
    by_name = {m['model_id']: m for m in models}
    # Stored scale names use their exact rational representation.
    selected = []
    for assets, lam, scale in ((2, F(2,5), F(3)), (4,F(1,20),F(3,2)),
                              (4,F(1,5),F(3)), (8,F(1,5),F(3))):
        choices = [m for m in models if m['assets']==assets and F(m['lambda'])==lam
                   and F(m['dispersion_scale'])==scale]
        if len(choices)!=1: raise ValueError('Four-case target not unique')
        selected.append(choices[0])
    return selected

def register():
    if PROTOCOL.exists(): raise FileExistsError('Registration already exists; preserve it')
    protocol = {'schema':'qbp-factor-basis-complete-quantum-protocol-v1',
        'registered_before_new_execution': True,
        'registration_time_utc': datetime.now(timezone.utc).isoformat(),
        'research_record_date':'2026-10-10', 'models':selected_models(), 'cap':'1',
        'factor_basis_selection':{'screen_pool':'Identity and every single coordinate plane with cosine3/5 and sine+4/5 or-4/5, exact right rotations R times O',
            'screen_count_by_case':[1,7,7,43], 'screen_total':58,
            'select':'Minimum exact direct global R1 transfer upper; tie candidate_id lexicographically',
            'resource_pool':'Identity plus selected candidate per original case, duplicate identity removed',
            'full_cost_selection':'Minimum complete T among resource-pool records; all four methods and both architectures receive each selected basis',
            'not_full_screen_cost_search':True, 'classical_opponents_receive_same_transformed_model_and_preparation':True},
        'methods':list(KINDS), 'architectures':list(ARCHITECTURES),
        'financial_law':{'L':'6','b_y':10,'q':28,'Beta_bits':8,
                         'precision_bits':512,'MGF_series_order':120},
        'scalar_precision_baseline':12,
        'repair_policy':{'eligible_method':'boundary_r2',
            'trigger':'Certified all-12 deterministic budget exceeds 1/100 or arithmetic guard fails',
            'repair':'All function addresses at 14 bits; q28, by10, L6, Beta8 remain unchanged',
            'finite_search_scope':'No greedy precision search; every baseline failure retained'},
        'monetary_radius':str(EPS), 'alpha_per_price':str(ALPHA),
        'selected_four_price_union_failure_upper':str(4*ALPHA),
        'confidence_scope':'Alternative resource plans are deterministic hypothetical plans; only four final selected executions are covered by this union bound. No probability is consumed by offline architecture, method or basis selection.',
        'R1_financial_information':'Retain previous R2+U Pprime transport alongside direct global bounded-payoff derivative transport; use their deterministic minimum and recompute coefficient-informed probability information.',
        'classical_price_answer_excluded_from_inputs':True,
        'source_sha256':sha(__file__), 'dependency_source_sha256':source_hashes(),
        'certificate_producer_source_sha256':sha(CERTIFICATE.parent.parent/'factor_basis_certificate.py'),
        'input_sha256':input_hashes(),
        'data_origin':'Exact rational synthetic GBM inputs and synthetic oracle tables; no LSEG-derived data.',
        'scope':'Constructed ideal-logical arithmetic, finite-confidence and complete QPE resource plans. No hardware execution or physical-runtime comparison.'}
    write_json(PROTOCOL, protocol)
    print('Registered', [m['model_id'] for m in protocol['models']], flush=True)

def direct_first_transfer(model, transport, norms):
    # A_j means the full asset loading, not a centred loading.
    loadings = tuple(tuple(F(x) for x in row) for row in model['loadings'])
    _,_,_,_,rg,_ = parameters(model,False)
    p=len(rg); M1=F(norms['M1_upper']); M2=F(norms['M2_upper'])
    m1=upper(transport.rho1(0)); m2=upper(transport.rho2mass)
    other=upper(transport.mass**(p-1)); rows=[]; linear=quadratic=F(0)
    for j in range(p):
        A=max(abs(row[j]) for row in loadings); r=rg[j]
        D1=M1*(A+abs(r)); D2=(M1+M2)*(A*A+r*r)
        l=other*m1*D1; q=other*m2*D2; linear+=l; quadratic+=q
        rows.append({'coordinate':j,'full_asset_loading_max':str(A),'geometric_loading':str(r),
            'D1_upper':str(D1),'D2_upper':str(D2),'linear_error_upper':str(l),'quadratic_error_upper':str(q)})
    tail=upper(1-transport.mass**p) # cap H=1
    return {'model_id':model['model_id'],'payoff':'capped','cap':'1','order':1,
        'linear_transport_upper':str(linear),'quadratic_transport_upper':str(quadratic),
        'tail_upper':str(tail),'financial_transfer_upper':str(linear+quadratic+tail),
        'mass_interval':serial(transport.mass**p), 'coordinates':rows,
        'rho1_mass_upper':str(m1),'rho2_mass_upper':str(m2),'other_coordinate_mass_upper':str(other),
        'scope':'Bounded capped-call first residual; exact full-loading derivative envelopes with raw unnormalised hybrid masses. No price answer.'}

def build(model, kind, domain, table_bits, encoding='centred'):
    return oracle.build_graph(model, F(1),q=28,b_y=10,beta_bits=8,table_bits=table_bits,
        outer_bound=6,U_star=F(domain['U_star']),kind=kind,
        retained_bound=max(abs(F(v)) for v in domain['coordinate_bounds']),output_encoding=encoding)

class Archives:
    def __init__(self):
        self.path=RESULTS/'factor_basis_quantum_scalar_ordinates.jsonl.gz'
        self.raw=self.path.open('wb'); self.stream=gzip.GzipFile(fileobj=self.raw, mode='wb',mtime=0)
        self.seen=set(); self.entries=0
    def retain(self, graph, tables, certificates):
        for name, certificate in certificates.items():
            checksum=digest(certificate)
            key=oracle.table_key(name,graph.nodes[name])
            if checksum in self.seen: continue
            self.seen.add(checksum); self.entries+=1
            item={'certificate_sha256':checksum,'certificate':certificate,
                  'packed_table':list(tables[key]),'packed_table_sha256':digest(tables[key]),
                  'graph_word_width':graph.width,'q':graph.q,
                  'scope':'Synthetic one-dimensional ordinates and actual packed integer table, never a financial-data download.'}
            self.stream.write(json.dumps(item, default=str,separators=(',',':')).encode()+b'\n')
    def close(self):
        self.stream.close();self.raw.close()

def expanded_models(protocol, certificate_data):
    """Validate the frozen analytic choice and exact covariance, then form pairs."""
    originals={m['model_id']:m for m in protocol['models']}
    if certificate_data['source_sha256']!=protocol['certificate_producer_source_sha256']:
        raise ValueError('Registered analytic screen producer changed')
    for key,value in certificate_data['dependency_source_sha256'].items():
        if sha(ROOT/key)!=value: raise ValueError(('Analytic dependency changed',key))
    records=certificate_data['records']
    if len(records)!=58: raise ValueError('Full analytic screen not retained')
    chosen=certificate_data['selected_models']
    for mid,original in originals.items():
        p=original['outer_dimension'];screen=[r for r in records if r['original_model_id']==mid]
        expected={'identity'}|{f'plane_{a}_{b}_{word}' for a in range(p) for b in range(a+1,p) for word in ('minus','plus')}
        if {r['candidate_id'] for r in screen}!=expected or len(screen)!=len(expected):
            raise ValueError('Analytic screen differs from registered finite pool')
        winner=min(screen,key=lambda r:(F(r['direct_global_transfer']['financial_transfer_upper']),r['candidate_id']))
        actual=next(x for x in chosen if x['original_model_id']==mid)
        if actual!=winner: raise ValueError('Analytic factor choice is not exact registered minimum')
        R=tuple(tuple(map(F,row)) for row in original['loadings'])
        for row in screen:
            descriptor=row['candidate'];cid=row['candidate_id']
            O=[[F(int(a==b)) for b in range(p)] for a in range(p)]
            if cid!='identity':
                a,b=descriptor['plane'];c=F(descriptor['cosine']);s=F(descriptor['sine'])
                expected_id=f"plane_{a}_{b}_{'plus' if s>0 else 'minus'}"
                if not 0<=a<b<p or c!=F(3,5) or abs(s)!=F(4,5) or expected_id!=cid:
                    raise ValueError('Unregistered rotation descriptor')
                O[a][a]=O[b][b]=c;O[a][b]=-s;O[b][a]=s
            elif descriptor['plane'] is not None or F(descriptor['cosine'])!=1 or F(descriptor['sine'])!=0:
                raise ValueError('Invalid identity descriptor')
            expected_R=[[str(sum(R[i][k]*O[k][j] for k in range(p))) for j in range(p)] for i in range(len(R))]
            if row['rotated_model']['loadings']!=expected_R:
                raise ValueError('Analytic screen loading is not exact registered rotation')
    if set(x['original_model_id'] for x in chosen)!=set(originals):
        raise ValueError('Analytic result does not cover exactly the registered four cases')
    pairs=[]
    for selected in chosen:
        original=originals[selected['original_model_id']]
        pairs.append({'original_model_id':original['model_id'],'candidate_id':'identity','model':original})
        if selected['candidate_id']=='identity': continue
        model=selected['rotated_model']
        if model['model_id']==original['model_id']: raise ValueError('Rotated label collides with identity')
        R=tuple(tuple(map(F,row)) for row in original['loadings'])
        Rp=tuple(tuple(map(F,row)) for row in model['loadings'])
        if len(R)!=len(Rp) or any(len(a)!=len(b) for a,b in zip(R,Rp)):
            raise ValueError('Rotation changed financial dimensions')
        for i in range(len(R)):
            for j in range(len(R)):
                if sum(a*b for a,b in zip(R[i],R[j]))!=sum(a*b for a,b in zip(Rp[i],Rp[j])):
                    raise ValueError('Exact covariance not preserved')
        for key,value in original.items():
            if key not in ('model_id','loadings') and model.get(key)!=value:
                raise ValueError(('Rotation changed a financial field',key))
        pairs.append({'original_model_id':original['model_id'],'candidate_id':selected['candidate_id'],'model':model})
    if len(pairs)>7: raise ValueError('Frozen pair pool exceeded')
    return pairs

def paired_comparison(selected):
    output=[]
    for original in selected_models():
        rows=[x for x in selected if x['original_model_id']==original['model_id']]
        choices=[]
        for row in rows:
            for kind,value in row['methods'].items():
                if value: choices.append((value['selected']['result']['T_count'],row['factor_basis_candidate_id'],kind,value['selected']))
        if not choices: raise ArithmeticError('No viable quantum resource choice')
        choices.sort(key=lambda x:(x[0],x[1],x[2]))
        identity=next(x for x in rows if x['factor_basis_candidate_id']=='identity')
        ikind=identity['quantum_resource_winner']; ivalue=identity['methods'][ikind]['selected']
        cost, cid, kind, value=choices[0]
        output.append({'original_model_id':original['model_id'],'identity_winning_method':ikind,
            'identity_T_count':ivalue['result']['T_count'],'selected_candidate_id':cid,
            'selected_method':kind,'selected_T_count':cost,'relative_T_saving_exact':str(F(ivalue['result']['T_count']-cost,ivalue['result']['T_count'])),
            'selected_strategy':value,'basis_count':len(rows)})
    return output

def execute():
    if (RESULTS/'factor_basis_quantum_resources.json').exists() or (RESULTS/'factor_basis_quantum_scalar_ordinates.jsonl.gz').exists():
        raise FileExistsError('Preserve existing execution; use a fresh output directory')
    protocol=json.loads(PROTOCOL.read_text()); protocol_sha=sha(PROTOCOL)
    if protocol['source_sha256']!=sha(__file__) or protocol['dependency_source_sha256']!=source_hashes():
        raise RuntimeError('Registered source hashes changed')
    if protocol['input_sha256']!=input_hashes(): raise RuntimeError('Registered inputs changed')
    certificate_data = json.loads(CERTIFICATE.read_text())
    certificate_sha = sha(CERTIFICATE)
    model_pairs = expanded_models(protocol, certificate_data)
    start=time.perf_counter();ctx.prec=512
    transport=tt.TransportMGF(6,10,28); archive=Archives()
    records=[]; preparations=[]; block_cache={}; candidate_ix=0
    for pair in model_pairs:
        model=pair['model']
        mid=model['model_id']; p=model['outer_dimension']; print('Prepare',mid,flush=True)
        moments_record=prep.product_grid_moments(model,6,10,28,6,transport.grid)
        domain=prep.finite_domain(model,6,10,28,transport.grid)
        moments={r:F(moments_record['centred_U_moments'][str(r)]['upper_exact_fraction']) for r in (1,2,3,4,6)}
        gap=F(moments_record['correlated_underlying_gap_interval']['upper_exact_fraction'])
        coefficients=coefficient_certificate(model,F(1),'capped',1)
        continuous_moment=informed.continuous_second_moment(model)
        hybrid=tt.HybridMoments(model,transport).record(); tail=tt.tail_moment(model,6)
        r2=tt.financial_transfer(model,F(1),hybrid,tail); r2['mass_interval']=serial(transport.mass**p)
        extra=first.ExtraMoments(model,transport).record_extra(); first_tail=first.outside_U_tilt(model,6)
        old_first=first.transfer(model,F(1),r2,extra,first_tail); old_first['mass_interval']=serial(transport.mass**p)
        norms=payoff_norms(model['lambda'],F(1)); new_first=direct_first_transfer(model,transport,norms)
        best_first=dict(min((old_first,new_first),key=lambda x:F(x['financial_transfer_upper'])))
        best_first['selection']='Minimum of separately proved older and direct global transfer bounds'
        preparations.append({'original_model_id':pair['original_model_id'],'candidate_id':pair['candidate_id'],'model':model,'domain':domain,'finite_moments':moments_record,
            'coefficient_certificate':coefficients,'continuous_second_moment':continuous_moment,
            'R2_hybrid_geometry':hybrid,'R2_tail':tail,'R2_financial_transfer':r2,
            'R1_extra_hybrid_geometry':extra,'R1_tail':first_tail,
            'R1_older_financial_transfer':old_first,'R1_direct_global_financial_transfer':new_first,
            'R1_selected_financial_transfer':best_first,'payoff_norms':norms})
        table_cache={}
        for kind in KINDS:
            candidates=[12]
            for bits in candidates:
                label='all12_baseline' if bits==12 else 'all14_declared_boundary_repair'
                print('Certify',mid,kind,label,flush=True)
                graph=build(model,kind,domain,bits)
                tables,ordinate=oracle.make_tables(graph,table_cache,True)
                archive.retain(graph,tables,ordinate)
                certificate=certifier.certify_graph(graph,tables,ordinate)
                certificate['precision_search_certificate_sha256']=digest(certificate)
                financial=best_first if kind=='positive_r1' else r2
                budget=retained.budget(model,F(1),kind,financial,coefficients,certificate,moments) if certificate['status']=='passed' else None
                row={'original_model_id':pair['original_model_id'],'factor_basis_candidate_id':pair['candidate_id'],'candidate_index':candidate_ix,'model_id':mid,'payoff':'capped','cap':'1',
                    'kind':kind,'candidate_label':label,'scalar_bits':bits,
                    'node_address_bits':{n:node.bits for n,node in graph.nodes.items() if node.kind=='function'},
                    'base_graph_certificate':certificate,'complete_deterministic_budget':budget,
                    'profiles':{},'strategies':[],'selected_strategy':None}
                candidate_ix+=1
                rejected=budget is None or F(budget['total_deterministic_undiscounted_price_error_upper'])>=EPS
                if rejected:
                    row['baseline_failure']={'reason':'Arithmetic guard failure' if budget is None else 'Deterministic monetary budget exhausts tolerance',
                        'deterministic_error_upper':None if budget is None else budget['total_deterministic_undiscounted_price_error_upper']}
                    if bits==12 and kind=='boundary_r2': candidates.append(14)
                    records.append(row);print('Retained failure',mid,kind,row['baseline_failure'],flush=True);continue
                previous=retained.prior(kind,F(1),certificate,moments,gap);row['correction_prior']=previous
                prior_variants=[('previous_proved_prior',budget,previous)]
                if kind=='positive_r1':
                    stronger=informed.first_probability_interval(model,F(1),coefficients,best_first,certificate,previous,continuous_moment)
                    row['coefficient_informed_prior']=stronger
                    old_budget=retained.budget(model,F(1),kind,old_first,coefficients,certificate,moments)
                    old_informed=informed.first_probability_interval(model,F(1),coefficients,old_first,certificate,previous,continuous_moment)
                    row['older_financial_deterministic_budget']=old_budget;row['older_coefficient_informed_prior']=old_informed
                    prior_variants.extend([('coefficient_informed_direct_minimum',budget,stronger),
                        ('older_financial_previous_prior',old_budget,previous),
                        ('older_financial_coefficient_informed_prior',old_budget,old_informed)])
                positive=previous['native_positive']
                for architecture in ARCHITECTURES:
                    profiles={}
                    for encoding in (('centred',) if positive else ('centred','positive_part','negative_part')):
                        pg=graph if encoding=='centred' else build(model,kind,domain,bits,encoding)
                        pt,pc=(tables,ordinate) if encoding=='centred' else oracle.make_tables(pg,table_cache,True)
                        archive.retain(pg,pt,pc)
                        profiles[encoding]=oracle.counted_composition(pg,pt,block_cache,architecture)
                    row['profiles'][architecture]=profiles
                    for information,bud,prior in prior_variants:
                        for strategy in (('positive',) if positive else ('centred','split')):
                            result=(informed.informed_complete_plan(strategy,bud,certificate,prior,profiles)
                                if 'probability_lower' in prior else retained.complete(strategy,bud,certificate,prior,profiles))
                            row['strategies'].append({'architecture':architecture,'strategy':strategy,
                                'prior_information':information,'result':result})
                available=[s for s in row['strategies'] if s['result']['kind']!='rejected']
                row['selected_strategy']=min(available,key=lambda x:(x['result']['T_count'],x['architecture'],x['strategy'],x['prior_information'])) if available else None
                records.append(row)
                print('Costed',mid,kind,label,'T',None if row['selected_strategy'] is None else row['selected_strategy']['result']['T_count'],flush=True)
        table_cache.clear()
    archive.close()
    selected=[]
    for pair in model_pairs:
        model=pair['model']
        row={'original_model_id':pair['original_model_id'],'factor_basis_candidate_id':pair['candidate_id'],'model_id':model['model_id'],'methods':{}}
        for kind in KINDS:
            eligible=[r for r in records if r['model_id']==model['model_id'] and r['kind']==kind and r['selected_strategy']]
            winner=min(eligible,key=lambda r:(r['selected_strategy']['result']['T_count'],r['candidate_index'])) if eligible else None
            row['methods'][kind]=None if winner is None else {'candidate_index':winner['candidate_index'],
                'candidate_label':winner['candidate_label'],'selected':winner['selected_strategy']}
        available=[k for k,v in row['methods'].items() if v]
        row['quantum_resource_winner']=min(available,key=lambda k:row['methods'][k]['selected']['result']['T_count']) if available else None
        selected.append(row)
    if protocol['dependency_source_sha256']!=source_hashes() or protocol['input_sha256']!=input_hashes() or protocol_sha!=sha(PROTOCOL):
        raise RuntimeError('Registered sources or inputs changed during execution')
    result={'schema':'qbp-factor-basis-complete-quantum-resources-v1',
        'status':'completed with all failures retained','protocol_file':PROTOCOL.name,'protocol_sha256':protocol_sha,
        'records':records,'preparations':preparations,'selected_prices':selected,
        'certificate_file':CERTIFICATE.name,'certificate_sha256':certificate_sha,
        'resource_pool':model_pairs, 'paired_comparison':paired_comparison(selected),
        'transport_profile':transport.profile(),'actual_financial_node_integers':[c['actual_node_integer'] for c in transport.cells],
        'source_sha256':sha(__file__),'dependency_source_sha256':protocol['dependency_source_sha256'],
        'input_sha256':protocol['input_sha256'],'scalar_ordinate_archive':archive.path.name,
        'scalar_ordinate_archive_sha256':sha(archive.path),'unique_scalar_ordinate_certificates':archive.entries,
        'unique_emitted_block_profiles':len(block_cache),'monetary_radius':str(EPS),'alpha_per_price':str(ALPHA),
        'selected_four_price_union_failure_upper':str(4*ALPHA),
        'confidence_scope':'Alternative resource plans are deterministic hypothetical plans; only four final selected executions are covered by this union bound. No probability is consumed by offline architecture, method or basis selection.',
        'profiles_count':sum(len(p) for r in records for p in r['profiles'].values()),
        'confidence_plan_count':sum(len(s['result'].get('flags',())) for r in records for s in r['strategies']),
        'generation_seconds_diagnostic':time.perf_counter()-start,
        'data_origin':protocol['data_origin'],'scope':protocol['scope']}
    path=RESULTS/'factor_basis_quantum_resources.json';write_json(path,result)
    summary={k:v for k,v in result.items() if k not in ('records','preparations','actual_financial_node_integers')}
    summary['full_result_sha256']=sha(path);summary['full_result_file']=path.name
    summary['candidate_count']=len(records);summary['failed_baseline_count']=sum('baseline_failure' in r for r in records)
    summary['selected_method_targets']=sum(v is not None for r in selected for v in r['methods'].values())
    write_json(RESULTS/'factor_basis_quantum_resource_summary.json',summary)
    print('Complete',summary['candidate_count'],'candidates',summary['selected_method_targets'],'targets',
          result['profiles_count'],'profiles',result['confidence_plan_count'],'plans',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(parents=[_BOOTSTRAP])
    parser.add_argument('command',choices=('register','run'))
    args=parser.parse_args()
    register() if args.command=='register' else execute()
