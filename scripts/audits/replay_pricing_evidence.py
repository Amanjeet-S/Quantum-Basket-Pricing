"""Replay the checked pricing evidence with explicit isolated output paths.

The historical reviewers use a phase-directory layout. This helper copies
the installed sources/evidence into that layout before invoking them, so
even auditors that regenerate supplements never modify repository records.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    repo=a.repo.resolve();out=a.output.resolve()
    if out.exists():raise FileExistsError('Choose a fresh audit output directory')
    out.mkdir(parents=True);step=out/'phase_inputs';mapping={}
    def source(alias,name):mapping[alias]=repo/name
    def result(alias,name):source(alias,'results/'+name)
    source('certificate/factor_basis_certificate.py','factor_basis_certificate.py')
    source('resources/factor_basis_quantum_resources.py','factor_basis_quantum_resources.py')
    source('resources/audit_factor_basis_quantum_resources.py','checks/audit_factor_basis_quantum_resources.py')
    source('review/audit_factor_basis.py','checks/audit_factor_basis.py')
    source('adversarial/quantum/shared_threshold_oracle.py','shared_threshold_oracle.py')
    source('adversarial/quantum/strengthened_r1_baseline.py','strengthened_r1_baseline.py')
    source('plot_pricing_precision.py','figures/plot_pricing_precision.py')
    source('plot_dispersion_separation_v2.py','figures/plot_dispersion_separation_v2.py')
    for n in ('factor_basis_protocol.json','factor_basis_certificates.json','factor_basis_summary.json','factor_basis_independent_audit.json'):
        result('certificate/results/'+n,n)
    for n in ('factor_basis_quantum_protocol.json','factor_basis_quantum_resources.json',
              'factor_basis_quantum_scalar_ordinates.jsonl.gz','factor_basis_quantum_independent_audit.json',
              'factor_basis_quantum_table_width_reconstruction.json.gz'):
        result('resources/results/'+n,n)
    result('review/factor_basis_certificate_audit.json','factor_basis_certificate_audit.json')
    for n in ('shared_threshold_protocol.json','shared_threshold_resources.json','shared_threshold_scalar_ordinates.json.gz'):
        result('adversarial/quantum/threshold_results/'+n,n)
    for n in ('strengthened_r1_protocol.json','strengthened_r1_resources.json','strengthened_r1_scalar_ordinates.json.gz'):
        result('adversarial/quantum/strengthened_baseline_results/'+n,n)
    for prefix in ('dispersion_separation','pricing_precision'):
        for suffix in ('.png','.svg','_figure.tex','_coordinates.json','_plot_data.json','_render_record.json'):
            source('figure/'+prefix+suffix,'figures/'+prefix+suffix)
    before={str(src.relative_to(repo)):sha(src) for src in mapping.values()}
    for alias,src in mapping.items():
        dst=step/alias;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    env=dict(os.environ);env['PYTHONPATH']=str(repo);env['PYTHONDONTWRITEBYTECODE']='1'
    def run(name,*arguments):
        command=[sys.executable,str(repo/'checks'/name),*map(str,arguments)]
        log=out/(name+'.log')
        with log.open('w') as f:
            r=subprocess.run(command,cwd=repo,env=env,stdout=f,stderr=subprocess.STDOUT)
        if r.returncode:
            print(log.read_text()[-8000:]);raise RuntimeError('Audit failed: '+name)
        print('Passed '+name,flush=True)
    cr=step/'certificate/results';rr=step/'resources/results';rv=step/'review'
    sr=rv/'shared_threshold';sr.mkdir(parents=True)
    run('audit_factor_basis_certificate.py','--repo',repo,'--results',cr,'--output',cr/'factor_basis_independent_audit.json')
    run('audit_factor_basis_certificates.py','--repo',repo,'--certificate-dir',step/'certificate','--output',rv/'factor_basis_certificate_audit.json')
    run('audit_factor_basis_quantum_resources.py','--repo',repo,'--results',rr,
        '--producer',repo/'factor_basis_quantum_resources.py','--certificate',cr/'factor_basis_certificates.json')
    run('audit_factor_basis_financial.py','--step',step,'--output',rv/'factor_basis_financial_audit.json')
    run('audit_shared_threshold_predicate.py','--prototype',repo/'shared_threshold_oracle.py',
        '--repo',repo,'--factor-resources',repo/'results/factor_basis_quantum_resources.json',
        '--output',sr/'shared_threshold_predicate_audit.json')
    run('audit_shared_threshold_financial.py','--step',step,'--repo',repo,'--output',sr/'shared_threshold_financial_audit.json')
    run('audit_shared_threshold_depth.py','--step',step,'--repo',repo,'--output',sr/'shared_threshold_depth_audit.json')
    run('audit_strengthened_r1_baseline.py','--step',step,'--repo',repo,'--output',sr/'strengthened_r1_baseline_audit.json')
    run('audit_pricing_precision_figure.py','--step',step,'--repo',repo,'--output',sr/'pricing_precision_figure_portable_audit.json')
    run('replay_positive_r1_feasibility.py','--repo',repo,
        '--factor-resources',repo/'results/factor_basis_quantum_resources.json',
        '--historical-source',repo/'history/positive_r1_feasibility/positive_r1_feasibility.py.txt',
        '--original-protocol',repo/'results/positive_r1_feasibility_protocol.json',
        '--original-result',repo/'results/positive_r1_feasibility.json','--output',out/'feasibility_replay')
    for name,h in before.items():assert sha(repo/name)==h,name
    records={str(x.relative_to(out)):sha(x) for x in out.rglob('*audit.json') if x.is_file()}
    summary={'schema':'qbp-installed-pricing-evidence-replay-v1','status':'passed',
        'repository_inputs_unchanged':True,'checked_installed_input_sha256':before,
        'audit_records':records,'replay_source_sha256':sha(__file__),
        'scope':'Fresh isolated replay of mathematical, financial, coherent predicate, complete charge, strengthened rival, figure mapping and historical feasibility checks. No repository evidence is overwritten.'}
    (out/'installed_replay_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('All installed evidence replays passed; repository inputs unchanged.',flush=True)

if __name__=='__main__':main()
