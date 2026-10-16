"""Export portable audited metadata while retaining the full cell evidence."""
from pathlib import Path
import gzip,hashlib,json

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/'results/boundary_grid_independent_audit.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    full=ROOT/'results/boundary_financial_grid_transfer.json.gz';audit=json.loads(AUDIT.read_text())
    assert audit['status']=='passed' and audit['certificate_sha256']==sha(full)
    data=json.load(gzip.open(full,'rt'));records=[{k:v for k,v in row.items() if k!='cells'} for row in data['records']]
    result={k:v for k,v in data.items() if k!='records'};result['records']=records
    result['full_cell_certificate']={'file':full.name,'sha256':sha(full)}
    result['independent_audit']={'file':'boundary_grid_independent_audit.json','sha256':sha(AUDIT)}
    result['summary_export_source_sha256']=sha(Path(__file__))
    result['classical_setup_ledger']={'shared_quantile_cells':1024,'quantile_inverse_evaluations':2047,
        'nodal_residual_derivative_evaluations':sum(r['cell_count'] for r in records),
        'segment_residual_derivative_interval_evaluations':sum(r['derivative_subinterval_count'] for r in records),
        'geometric_coefficient_certificates':len(data['coefficients']),
        'scope':'Finite mathematical construction counts; no conversion of Arb work to T gates or seconds. The independently repeated audit is review work.'}
    (ROOT/'results/boundary_financial_grid_transfer_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'results/boundary_grid_independent_audit.json').write_bytes(AUDIT.read_bytes())
    assert (ROOT/'results/original_q28_quantile_nodes.json').exists(), 'Retained original quantile evidence must be included in the package'
    print(json.dumps({'summary':str(ROOT/'results/boundary_financial_grid_transfer_summary.json'),'records':len(records),'complete_source_sha256':data['source_sha256']}))
if __name__=='__main__':main()
