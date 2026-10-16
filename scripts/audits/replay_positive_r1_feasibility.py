"""Portable replay of the unchanged historical finite feasibility source.

The original source and record are preserved. Only runtime path bindings are
changed. Fresh protocol/result files live in a separate supplied output path;
their numerical rows must reproduce the original rows exactly.
"""
from pathlib import Path
from importlib.machinery import SourceFileLoader
from importlib.util import spec_from_loader, module_from_spec
from fractions import Fraction
import argparse, hashlib, json, sys

sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True)
    p.add_argument('--factor-resources',required=True);p.add_argument('--historical-source',required=True)
    p.add_argument('--original-protocol',required=True);p.add_argument('--original-result',required=True)
    p.add_argument('--output',required=True);a=p.parse_args()
    repo=Path(a.repo).resolve();factor=Path(a.factor_resources).resolve();src=Path(a.historical_source).resolve()
    prot=Path(a.original_protocol).resolve();old=Path(a.original_result).resolve();out=Path(a.output).resolve()
    protocol=json.loads(prot.read_text());original=json.loads(old.read_text())
    assert original['protocol_sha256']==sha(prot)
    expected={Path(k).name:v for k,v in protocol['inputs_sha256'].items()}
    assert sha(src)==original['source_sha256']==expected['positive_r1_feasibility.py']
    assert sha(factor)==expected['factor_basis_quantum_resources.json']
    for name in ('taylor_pricing.py','continuous_quadrature.py','qpe_estimator.py','certified_scalar.py'):
        assert sha(repo/name)==expected[name]
    assert sha(repo/'results/gaussian_derivative_norms.json')==expected['gaussian_derivative_norms.json']
    sys.path.insert(0,str(repo));loader=SourceFileLoader('retained_feasibility_source',str(src))
    spec=spec_from_loader(loader.name,loader);m=module_from_spec(spec);loader.exec_module(m)
    out.mkdir(parents=True,exist_ok=True)
    m.ROOT=repo;m.STAGE=out;m.RESOURCE=factor;m.NORMS=repo/'results/gaussian_derivative_norms.json'
    m.PROTOCOL=out/'feasibility_replay_protocol.json';m.RESULT=out/'feasibility_replay_result.json'
    m.register();m.execute();new=json.loads(m.RESULT.read_text())
    assert new['rows']==original['rows']
    assert new['vertices_rebuilt']==original['vertices_rebuilt']
    record={'schema':'qbp-positive-r1-feasibility-portable-replay-v1','status':'passed',
      'audit_source_sha256':sha(__file__),'historical_source_sha256':sha(src),
      'original_protocol_sha256':sha(prot),'original_result_sha256':sha(old),
      'replay_protocol_sha256':sha(m.PROTOCOL),'replay_result_sha256':sha(m.RESULT),
      'numerical_rows_sha256':digest(new['rows']),'pairs_replayed':len(new['rows']),
      'vertices_rebuilt':new['vertices_rebuilt'],'unchanged_historical_source':True,
      'only_runtime_path_bindings_changed':True,'numerical_rows_identical':True}
    dest=out/'positive_r1_feasibility_portable_audit.json'
    if dest.exists():raise FileExistsError('Preserve existing audit')
    dest.write_text(json.dumps(record,indent=2)+'\n');print('Portable exact numerical replay passed')

if __name__=='__main__':main()
