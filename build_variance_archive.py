"""Preserve the full resource comparison as a deterministic portable gzip.

The driver calls the unchanged retained planner and adds a concrete reusable
shot reset protocol. This changes zero primary T costs, never selection or
confidence. Expanded JSON is a local optional copy, excluded from Git.
"""
from pathlib import Path
import gzip,hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reset_charge(candidate,L,q):
    if 'T' not in candidate:return
    charges=candidate['charges']
    if candidate['kind']=='dyadic':
        z=x=0
        for b in candidate['buckets']:
            live=L+(b['index']+1 if b['index']<q else 0)
            z+=b['maximum_shots']*live;x+=b['maximum_shots']*(live+1)
    else:
        live=L+candidate['c'];z=candidate['shots']*live;x=candidate['shots']*(live+1)
    charges['reinitialisation_additional_Z_measurements_upper']=z
    charges['reinitialisation_conditional_X_upper']=x
    charges['all_measurements_including_final_and_reinitialisation_upper']=charges.get('measurements',0)+candidate['shots']+z
    candidate['reinitialisation_scope']='measure each retained latent/comparison wire in Z after the final flag measurement; conditional X on each retained wire including flag; private work already clean'

def main():
    archive=ROOT/'results/variance_oracle_comparison.json.gz'
    manifest=ROOT/'results/variance_oracle_archive.json'
    if '--expand' in sys.argv:
        raw=gzip.decompress(archive.read_bytes());meta=json.loads(manifest.read_text())
        assert hashlib.sha256(raw).hexdigest()==meta['uncompressed_sha256']
        out=archive.with_suffix('')
        if out.exists() and out.read_bytes()!=raw:raise RuntimeError('Refusing to replace a different local expanded result')
        out.write_bytes(raw);print(out);return
    original_write=Path.write_text
    def write(self,data,*args,**kwargs):
        if self.resolve()!=ROOT/'results/variance_oracle_comparison.json':return original_write(self,data,*args,**kwargs)
        result=json.loads(data)
        source=json.loads((ROOT/'results/architecture_resource_candidates.json').read_text())
        lookup={(r['candidate_id'],r['architecture']):r for r in source['configurations']}
        for row in result['configurations']:
            L=lookup[(row['candidate_id'],row['architecture'])]['resources']['input_qubits'];q=row['q']
            for record in row['controllers'].values():
                for name,candidates in record.items():
                    for c in candidates if isinstance(candidates,list) else (candidates,):reset_charge(c,L,q)
        result['source_hashes']['build_variance_archive.py']=digest(Path(__file__))
        result['selection_protocol']['reinitialisation']='event counts plus explicit per-wire Z measurement and conditional X upper bounds; no T change'
        raw=(json.dumps(result,indent=2)+'\n').encode()
        packed=gzip.compress(raw,compresslevel=9,mtime=0);archive.write_bytes(packed)
        meta={'status':'complete comparison retained without loss as deterministic gzip','uncompressed_bytes':len(raw),
              'compressed_bytes':len(packed),'uncompressed_sha256':hashlib.sha256(raw).hexdigest(),
              'compressed_sha256':hashlib.sha256(packed).hexdigest(),'driver_sha256':digest(Path(__file__)),
              'rows':len(result['configurations']),'summary':result['summary'],
              'expanded_json':'optional local copy; ignored by Git; no private observations'}
        original_write(manifest,json.dumps(meta,indent=2)+'\n')
        print(json.dumps({'archive':str(archive),'uncompressed_bytes':len(raw),'compressed_bytes':len(packed)},indent=2),flush=True)
        return len(data)
    Path.write_text=write
    import variance_oracle_comparison
    sys.argv=['variance_oracle_comparison.py']+([] if '--compute' in sys.argv else ['--merge'])
    variance_oracle_comparison.main()
if __name__=='__main__':main()
