"""Whole literal preparation checks for smoothing winners that changed.
No producer or producer reviewer imported; graphs reconstructed from retained DAG.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,gzip,hashlib,json,sys
sys.dont_write_bytecode=True
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
digest=lambda x:hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def main(a):
 repo=Path(a.repo);raw=Path(a.raw_results);amended=Path(a.amended_results);out=Path(a.output);assert not out.exists();sys.path.insert(0,str(repo))
 import boundary_value_oracle as e
 import quantum_gates as qg
 data=json.load(gzip.open(raw/'raw_payoff_resources.json.gz','rt'));new=json.loads((amended/'raw_payoff_smoothing_resources.json').read_text());assert new['base_resources_sha256']==sha(raw/'raw_payoff_resources.json.gz')
 archives=json.load(gzip.open(raw/'raw_payoff_scalar_ordinates.json.gz','rt'));index={}
 for a in archives:
  c=a['certificate'];w=a['graph_word_width']
  if c['schema']=='qbp-integral-quantile-v1':key=('quantile',c['b_y'],c['L']);own=tuple(n['offset_integer'] for n in c['nodes'])
  else:
   key=('scalar',c['name'],c['lo'],c['span'],tuple(sorted(c['params'].items())),c['table_bits']);vs=[n['integer'] for n in c['nodes']];own=tuple(v | (abs(vs[j+1]-v if j+1<len(vs) else 0)<<w) | (((vs[j+1]-v if j+1<len(vs) else 0)<0)<<(2*w)) for j,v in enumerate(vs))
  assert own==tuple(a['packed_table']);index[key,w]=own
 bindings=[]
 for cmp in new['selected_comparisons']:
  before=next(c for c in data['selected_comparisons'] if c['original_model_id']==cmp['model_id'])
  if (before['selected_normal_bits'],before['selected_payoff_exp_bits'])==(cmp['selected_normal_bits'],cmp['selected_exp_bits']):continue
  row=next(r for r in data['records'] if r['original_model_id']==cmp['model_id'] and r['normal_bits']==cmp['selected_normal_bits'] and r['payoff_exp_bits']==cmp['selected_exp_bits'])
  g=e.IntegralGraph(28,12)
  for j in range(row['model']['outer_dimension']):g.inp('iy'+str(j),10)
  g.inp('inner_normal_label',row['normal_bits'])
  for n in row['graph_nodes']:
   kind=n['kind'];meta=n['meta']
   if kind=='function':meta=e.ScalarSpec(meta['name'],F(meta['lo']),F(meta['span']),tuple(sorted((k,F(v)) for k,v in meta['parameters'].items())),tuple(map(F,meta['output_bounds'])))
   elif kind in ('constant','scale','clamp','clamp_signed','quantile'):meta=F(meta)
   elif meta=='None':meta=None
   assert g.op(kind,*n['args'],meta=meta,bits=n['bits'],bounds=n['ideal_bounds'])==n['name']
  g.semantic_nodes=row['graph_semantic_nodes'];g.configuration=dict(row['graph_configuration']);g.finish(g.semantic_nodes['probability']);assert g.configuration==row['graph_configuration']
  tables={}
  for name,n in g.nodes.items():
   if n.kind=='quantile':key=('quantile',n.bits,str(n.meta))
   elif n.kind=='function':key=('scalar',n.meta.name,str(n.meta.lo),str(n.meta.span),tuple(sorted((k,str(v)) for k,v in n.meta.parameters)),n.bits)
   else:continue
   tables[e.table_key(name,n)]=index[key,g.width]
  for mode in ('exact7t','unitary8t'):
   c,ins,comparison,flag=e.compile_flag(g,tables,mode,store_gates=False);m=c.resources();p=row['profiles'][mode]['centred'];counts=Counter(m['clifford_T_expanded_gate_counts']);counts['H']+=len(c.uniform_preparation_gates);assert counts==Counter(p['primitive_counts'])
   core=tuple(w for reg in ins.values() for w in reg)+tuple(comparison)+(flag,);assert core==c.reflection_core and len(core)==p['reflection_core_qubits']
   highzero=flag-1;assert highzero not in core and m['wires']<=p['logical_qubits_upper'];work=tuple(w for w in range(m['wires']) if w not in core and w!=highzero)[:len(core)-1];assert len(work)==len(core)-1
   reflection=Counter(x[0] for x in qg.expand_toffoli(qg.grover([],core,flag,work,control=m['wires'])));ccx=2*(len(core)-1);assert reflection==Counter({'X':2*len(core),'CZ':2,'Z':1,'T':4*ccx,'Tdg':3*ccx,'H':2*ccx,'CX':6*ccx})
   bindings.append({'model_id':cmp['model_id'],'normal_bits':row['normal_bits'],'exp_bits':row['payoff_exp_bits'],'architecture':mode,'profile_sha256':digest(p),'flag':flag,'excluded_high_zero':highzero,'reflection_core_qubits':len(core),'full_literal_emitted_peak':m['wires'],'serial_profile_peak_upper':p['logical_qubits_upper'],'complete_primitive_counts_verified':True,'actual_core_and_controlled_reflection_verified':True})
   print('Changed smoothing winner verified',cmp['model_id'],mode,flush=True)
 out.write_text(json.dumps({'schema':'qbp-smoothing-changed-winner-preparation-audit-v1','status':'passed','source_sha256':sha(__file__),'raw_resources_sha256':sha(raw/'raw_payoff_resources.json.gz'),'amendment_resources_sha256':sha(amended/'raw_payoff_smoothing_resources.json'),'bindings':bindings,'scope':'Fresh whole literal preparations, actual flag/core, excluded high-zero and controlled reflections for amended winners absent from original whole-preparation binding set.'},indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--raw-results',required=True);p.add_argument('--amended-results',required=True);p.add_argument('--output',required=True);main(p.parse_args())
