"""Reproduce the frozen four-case quantum experiment in an isolated directory.

Additional project producers do not enter the registered source set. The
original executed producer, dependency hashes and original evidence are
preserved. No private financial data are read or copied.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True,
        help='A new or empty directory for fresh registration and results')
    parser.add_argument('--prepare-only',action='store_true',
        help='Validate and stage the exact recorded inputs without executing')
    args=parser.parse_args();out=args.output.resolve()
    if out.exists() and any(out.iterdir()):raise SystemExit('Output directory must be new or empty')
    protocol=json.loads((ROOT/'results/four_survivor_quantum_protocol.json').read_text())
    sources=protocol['dependency_source_sha256'];inputs=protocol['input_sha256']
    for name,checksum in sources.items():
        if Path(name).name!=name or sha(ROOT/name)!=checksum:
            raise SystemExit('Frozen source failed verification: '+name)
    for name,checksum in inputs.items():
        if Path(name).name!=name or sha(ROOT/'results'/name)!=checksum:
            raise SystemExit('Frozen synthetic input failed verification: '+name)
    out.mkdir(parents=True,exist_ok=True);(out/'results').mkdir()
    for name in sources:shutil.copyfile(ROOT/name,out/name)
    for name in inputs:shutil.copyfile(ROOT/'results'/name,out/'results'/name)
    manifest={'schema':'qbp-isolated-four-survivor-reproduction-v1',
        'original_registration_sha256':sha(ROOT/'results/four_survivor_quantum_protocol.json'),
        'source_sha256':sources,'input_sha256':inputs,
        'wrapper_sha256':sha(__file__),'output_directory':str(out),
        'scope':'Exact isolated source/input copy. Fresh execution generates a new dated protocol and new diagnostic timings; original retained evidence remains unchanged.'}
    (out/'reproduction_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Verified and staged',len(sources),'sources and',len(inputs),'synthetic inputs at',out)
    if args.prepare_only:return
    producer=out/'four_survivor_quantum_resources.py'
    subprocess.run([sys.executable,str(producer),'--register'],cwd=out,check=True)
    subprocess.run([sys.executable,str(producer),'--run'],cwd=out,check=True)

if __name__=='__main__':main()
