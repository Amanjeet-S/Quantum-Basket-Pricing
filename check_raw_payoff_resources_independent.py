"""Independent complete-charge and actual winner-preparation review.
Does not import producer, producer's auditor or quantum confidence frontend.
Reuses established primitive emitters for actual chronological preparation.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from math import isqrt,comb
import gzip,hashlib,json,os,sys,time
sys.dont_write_bytecode=True
REPO=Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
STAGE=Path(os.environ['QBP_PRIVATE_STAGE'])/'raw_payoff_step'/'producer'
OUT=Path(__file__).parent/'independent_resource_review.json'
assert not OUT.exists()
sys.path.insert(0,str(REPO))
import boundary_value_oracle as emitter
import quantum_gates as gates
start=time.monotonic()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
digest=lambda x:hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
resources=json.load(gzip.open(STAGE/'results_tight/raw_payoff_resources.json.gz','rt'))
protocol=json.loads((STAGE/'results_tight/raw_payoff_protocol.json').read_text())
assert sha(STAGE/'raw_payoff_oracle_tight.py')==protocol['source_sha256']==resources['source_sha256']
assert sha(STAGE/'results_tight/raw_payoff_protocol.json')==resources['protocol_sha256']
for n,h in resources['dependencies_sha256'].items():assert sha(REPO/n)==h
assert len(resources['records'])==48
archive=json.load(gzip.open(STAGE/'results_tight/raw_payoff_scalar_ordinates.json.gz','rt'))
assert sha(STAGE/'results_tight/raw_payoff_scalar_ordinates.json.gz')==resources['scalar_archive_sha256']
# Certificate-to-table binding is derived from raw node integers, not table claims.
packed={}
for a in archive:
 c=a['certificate'];w=a['graph_word_width'];assert digest(c)==a['certificate_sha256']
 if c['schema']=='qbp-integral-quantile-v1':own=[n['offset_integer'] for n in c['nodes']]
 else:
  v=[n['integer'] for n in c['nodes']];own=[]
  for i,z in enumerate(v):
   d=v[i+1]-z if i+1<len(v) else 0
   assert 0<=z<1<<w and abs(d)<1<<w
   own.append(z | (abs(d)<<w) | ((d<0)<<(2*w)))
 assert own==a['packed_table'] and digest(own)==a['packed_table_sha256']
 key=(('quantile',c['b_y'],c['L']) if c['schema']=='qbp-integral-quantile-v1' else ('scalar',c['name'],c['lo'],c['span'],tuple(sorted(c['params'].items())),c['table_bits']))
 packed[key,w]=(tuple(own),c)

def own_graph(row):
 g=emitter.IntegralGraph(28,12)
 for j in range(row['model']['outer_dimension']):g.inp('iy'+str(j),10)
 g.inp('inner_normal_label',row['normal_bits'])
 for n in row['graph_nodes']:
  kind=n['kind'];meta=n['meta']
  if kind=='function':meta=emitter.ScalarSpec(meta['name'],F(meta['lo']),F(meta['span']),tuple(sorted((k,F(v)) for k,v in meta['parameters'].items())),tuple(map(F,meta['output_bounds'])))
  elif kind in ('constant','scale','clamp','clamp_signed','quantile'):meta=F(meta)
  elif meta=='None':meta=None
  assert g.op(kind,*n['args'],meta=meta,bits=n['bits'],bounds=n['ideal_bounds'])==n['name']
 g.semantic_nodes=row['graph_semantic_nodes'];g.configuration=dict(row['graph_configuration']);g.finish(g.semantic_nodes['probability'])
 assert g.configuration==row['graph_configuration'] and all(g.nodes[n['name']].output_bits==n['output_bits'] for n in row['graph_nodes'])
 tables={}
 for name,n in g.nodes.items():
  if n.kind=='quantile':key=('quantile',n.bits,str(n.meta))
  elif n.kind=='function':key=('scalar',n.meta.name,str(n.meta.lo),str(n.meta.span),tuple(sorted((k,str(v)) for k,v in n.meta.parameters)),n.bits)
  else:continue
  tables[emitter.table_key(name,n)]=packed[key,g.width][0]
 return g,tables

PI=F(355,113);POST=F(1,10**60)
def sqrt_cap(x):
 y=x*(1<<200);k=(y.numerator+y.denominator-1)//y.denominator;r=isqrt(k)
 return F(r+(r*r<k),1<<100)
def binom(n,p):return sum((F(comb(n,j))*p**j*(1-p)**(n-j) for j in range((n+1)//2)),F(0))
checks=0;bindings=[];summaries=[]
for r in resources['records']:
 assert not r.get('rejection') and len(r['strategies'])==2
 g,tables=own_graph(r);assert set(g.inputs)=={'iy'+str(j) for j in range(r['model']['outer_dimension'])}|{'inner_normal_label'}
 s=g.semantic_nodes;assert g.nodes[s['upper_endpoint']].args==(s['lower_endpoint'],s['logT'])
 assert g.nodes[s['PG_exp']].meta==g.nodes[s['PA_exp']].meta and g.nodes[s['PG_exp']].meta.output_bounds==(F(1),F(64,63))
 budget=r['complete_deterministic_budget'];cert=r['probability_certificate'];prior=r['coefficient_informed_prior']
 A=F(budget['external_factor_upper']);det=F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*F(cert['probability_error_upper'])+F(budget['external_factor_radius'])
 assert det==F(budget['total_deterministic_undiscounted_price_error_upper'])
 for st in r['strategies']:
  mode=st['architecture'];result=st['result'];p=r['profiles'][mode]['centred'];f=result['flags'][0];q=f['plan'];qc=f['complete_charge'];qft=f['fourier_profile'];n=sum(g.inputs.values())+29
  assert p['reflection_core_qubits']==n and p['uniform_H']==sum(g.inputs.values())+28 and p['endpoint_high_zero_bits']==1
  assert len(p['blocks'])==len(g.nodes) and all(a['node']==b for a,b in zip(p['blocks'],g.nodes))
  # Derive forward/clean/comparator totals from complete retained block list.
  fw=Counter();depth=five=lookups=0;retained=sum(g.inputs.values())+30;peak=retained
  for b in p['blocks']:
   fw.update(b['primitive_counts']);depth+=b['explicit_T_depth_serial_upper'];five+=b['conservative_5layer_T_depth_serial_upper'];lookups+=b['lookup_invocations'];peak=max(peak,retained+b['extra_qubits_upper']);retained+=b['output_bits']
  comparator=Counter(p['comparator_primitive_counts']);full=Counter({k:2*v for k,v in fw.items()})+comparator;full['H']+=p['uniform_H']
  assert full==Counter(p['primitive_counts']) and p['logical_T_count']==full['T']+full['Tdg']
  ccx=2*(n-1);reflection=Counter({'X':2*n,'CZ':2,'Z':1,'T':4*ccx,'Tdg':3*ccx,'H':2*ccx,'CX':6*ccx})
  assert reflection==Counter(qc['controlled_reflection']['gate_counts'])
  assert q['M']==1<<q['phase_register_bits'] and F(q['epsilon'])==(F(1,100)-det)/A and F(q['alpha'])==F(1,3600) and F(q['prior_upper'])==F(prior['probability_upper'])
  upper=F(prior['probability_upper']);variance=upper*(1-upper) if upper<=F(1,2) else F(1,4);assert variance==F(q['variance_upper'])
  M=q['M'];eps=2*PI*sqrt_cap(variance)/M+PI**2/M**2
  assert F(q['amplitude_error_upper'])==eps and F(q['postprocess_allowance'])==POST and eps+POST<=F(q['epsilon'])
  assert M==2 or 2*PI*sqrt_cap(variance)/(M//2)+PI**2/(M//2)**2+POST>F(q['epsilon'])
  success=8/(PI**2)-F(q['qft_norm_error_allowance']);rep=q['repetitions'];failure=binom(rep,success)
  assert rep%2 and failure==F(q['binomial_failure_upper'])<=F(1,3600) and (rep==1 or binom(rep-2,success)>F(1,3600))
  inverse=Counter()
  for op,v in full.items():inverse[{'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}.get(op,op)]+=v
  Q=full+inverse+reflection;shot=full+Counter({op:(M-1)*v for op,v in Q.items()})+Counter(qft['gate_counts']);shot['H']+=q['phase_register_bits']
  peak=max(p['logical_qubits_upper'],2*n)+q['phase_register_bits']+qft.get('extra_ancillas',0)
  assert qc['complete_oracle_T']==full['T']+full['Tdg'] and qc['controlled_Q_T']==Q['T']+Q['Tdg'] and qc['peak_logical_qubits']==peak
  assert qc['metrics']['T']==rep*(shot['T']+shot['Tdg'])==result['T_count']
  assert qc['metrics']['measurements']==rep*peak and qc['metrics']['reset_X']==rep*peak and qc['complete_reinitialisations']==rep and qc['controlled_Q_calls']==rep*(M-1)
  assert det+A*(eps+POST)==F(result['total_price_error_upper'])<=F(1,100)
  assert F(result['failure_probability_upper'])==failure
  checks+=1
 assert r['selected_strategy']==min(r['strategies'],key=lambda x:(x['result']['T_count'],x['architecture']))

numeric=json.loads((REPO/'results/strengthened_r1_resources.json').read_text());threshold=json.loads((REPO/'results/shared_threshold_resources.json').read_text())
for supplied in resources['selected_comparisons']:
 mid=supplied['original_model_id'];cohort=[r for r in resources['records'] if r['original_model_id']==mid]
 winner=min(cohort,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['payoff_exp_bits']))
 assert winner['selected_strategy']==supplied['selected_raw_payoff_strategy'] and len(cohort)==12
 raw=winner['selected_strategy']['result']['T_count'];num=next(c['selected_numerical_strategy']['result']['T_count'] for c in numeric['selected_comparisons'] if c['original_model_id']==mid);thr=next(c['selected_strategy']['result']['T_count'] for c in threshold['rows'] if c['original_model_id']==mid)
 strongest=min(raw,num);saving=F(strongest-thr,strongest);summaries.append({'model_id':mid,'raw_T':raw,'numerical_T':num,'shared_threshold_T':thr,'strongest_registered_comparator_T':strongest,'saving_exact':str(saving),'saving_per_cent_diagnostic':float(100*saving),'normal_bits':winner['normal_bits'],'exp_bits':winner['payoff_exp_bits'],'exp_upper_boundary':winner['payoff_exp_bits']==max(protocol['narrow_exp_address_bits'])})
 # Fresh literal whole preparation for both architectures, using only retained node tables.
 g,t=own_graph(winner)
 for mode in ('exact7t','unitary8t'):
  circuit,ins,compare,flag=emitter.compile_flag(g,t,mode,store_gates=False);m=circuit.resources();counts=Counter(m['clifford_T_expanded_gate_counts']);counts['H']+=len(circuit.uniform_preparation_gates);p=winner['profiles'][mode]['centred']
  assert counts==Counter(p['primitive_counts'])
  core=tuple(w for reg in ins.values() for w in reg)+tuple(compare)+(flag,);assert core==circuit.reflection_core and len(core)==p['reflection_core_qubits']
  excluded=flag-1;assert excluded not in core and len(compare)==28 and m['wires']<=p['logical_qubits_upper']
  work=tuple(w for w in range(m['wires']) if w not in core and w!=excluded)[:len(core)-1];assert len(work)==len(core)-1
  reflect=gates.expand_toffoli(gates.grover([],core,flag,work,control=m['wires']))
  f=next(st['result']['flags'][0] for st in winner['strategies'] if st['architecture']==mode)
  assert Counter(x[0] for x in reflect)==Counter(f['complete_charge']['controlled_reflection']['gate_counts'])
  bindings.append({'model_id':mid,'architecture':mode,'flag':flag,'excluded_high_zero':excluded,'reflection_core_qubits':len(core),'chronological_emitted_peak':m['wires'],'serial_profile_peak_upper':p['logical_qubits_upper'],'full_preparation_counts_equal':True,'reflection_counts_equal':True})
  print('Whole preparation verified',mid,mode,flush=True)
OUT.write_text(json.dumps({'schema':'qbp-independent-raw-resource-review-v1','status':'passed','resource_sha256':sha(STAGE/'results_tight/raw_payoff_resources.json.gz'),'protocol_sha256':sha(STAGE/'results_tight/raw_payoff_protocol.json'),'review_source_sha256':sha(__file__),'configurations':48,'independent_complete_charge_and_QPE_checks':checks,'whole_literal_preparations':bindings,'exact_comparator_minima':summaries,'elapsed_seconds_diagnostic':time.monotonic()-start,'limitations':['Reuses established primitive emitter, not a new low-level implementation.','Does not independently prove finite inner-normal mathematics, financial transfer or novelty.','Winning exp bits10 is registered upper boundary in every case, so no global precision optimality.','Counted logical upper bounds and canonical QPE only; no hardware outcomes or quantum price estimates.']},indent=2)+'\n')
print('Independent resource review passed',checks,'complete plans and',len(bindings),'whole preparations',flush=True)
