"""Complete, reproducible synthetic continuous-price certificates.

Run ``--freeze`` first.  No sample words are requested until the immutable
protocol exists and its source/input hashes agree.  The old interrupted
experiment is not replayed or used as a numerical input.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
from array import array
import argparse, hashlib, json, math, os, platform, sys, time
from flint import arb, ctx
import flint
from boundary_grid_certificate import (ball, endpoints, serial, interval,
    gaussian_cells, abs_upper, parameters, Phi)
from boundary_derivative_certificate import payoff_norms

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
PROTOCOL=RESULTS/'four_case_classical_protocol.json'
CASE_IDS=('n2_lambda0.4_scale3','n4_lambda0.05_scale3over2',
          'n4_lambda0.2_scale3','n8_lambda0.2_scale3')
INPUTS=('results/dispersion_separation_protocol.json',
        'results/dispersion_separation_results.json',
        'results/boundary_derivative_norms.json',
        'results/gaussian_derivative_norms.json',
        'boundary_grid_certificate.py','boundary_derivative_certificate.py')
N=100000
EPS=F(1,150)
Q_SCORE=40

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')
def fraction_record(lo,hi):
    return {'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)}

def freeze():
    source=sha(Path(__file__))
    p={'schema':'qbp-four-case-classical-fixed-protocol-v1',
       'created_utc':datetime.now(timezone.utc).isoformat(),
       'status':'Frozen before fresh operating-system random words are requested',
       'previous_interrupted_run':'Lost temporary scores and traces are not numerical inputs; all samples and enclosures will be generated afresh.',
       'case_ids':CASE_IDS,'cap':'1','draws_per_case':N,
       'grid':{'L':'6','node_bits':10,'node_rounding_bits':28},
       'statistical_epsilon':str(EPS),'per_case_failure_allowance':'1/3600',
       'family_failure_allowance':'1/900',
       'confidence_premise':'Captured operating-system random uint16 words are treated as independent uniformly distributed words. This probabilistic premise is explicit and is not physically proved by retaining a realised trace.',
       'sampling_rule':'Each coordinate label is the low ten bits of one unsigned little-endian uint16 word; 100000 fixed draws per case, without stopping based on scores.',
       'score_precision_bits':128,'outward_score_bits':Q_SCORE,
       'transfer_precision_bits':512,
       'method':'Direct analytically conditioned capped-call payoff on the exact dyadic CDF midpoint grid; global derivative norms, positive Gaussian cell transfer and bounded-payoff full tails; n2 also has complete 1024-node exact interval summation before sampling.',
       'reporting':'Undiscounted tolerance 1/100, then deterministic exp(-3/100) discount; every case retained. Timings are diagnostic, not quantum runtime comparisons.',
       'data':'Exact stipulated synthetic risk-neutral models only; no LSEG or other licensed observations, outputs or private data.',
       'producer_sha256':source,'input_sha256':{name:sha(ROOT/name) for name in INPUTS}}
    write_json(PROTOCOL,p)
    print('protocol frozen',sha(PROTOCOL),flush=True)

def validate_protocol():
    p=json.loads(PROTOCOL.read_text())
    if p['producer_sha256']!=sha(Path(__file__)):raise RuntimeError('Producer differs from frozen protocol')
    if p['input_sha256']!={name:sha(ROOT/name) for name in INPUTS}:raise RuntimeError('Input differs from frozen protocol')
    if list(p['case_ids'])!=list(CASE_IDS) or p['draws_per_case']!=N:raise RuntimeError('Protocol parameters differ')
    return p

def clip_score(value):
    lo,hi=endpoints(value);lo=max(F(0),lo);hi=min(F(1),hi)
    if lo>hi:raise ArithmeticError('Analytic score contradicts [0,1] payoff premise')
    scale=1<<Q_SCORE
    a=lo*scale;b=hi*scale
    lower=a.numerator//a.denominator
    upper=b.numerator//b.denominator+int(b.denominator!=1)
    if not 0<=lower<=upper<=scale:raise ArithmeticError('Invalid outward q40 score')
    return lower,upper

class Scores:
    def __init__(self,model,cells):
        ctx.prec=128
        self.w=tuple(ball(F(w)) for w in model['weights'])
        self.lam=ball(F(model['lambda']));self.hlog=ball(F(101,100)).log()
        self.drifts=tuple(ball(F(mu)).exp() for mu in model['log_drifts'])
        nodes=tuple(F(c['actual_node']) for c in cells)
        self.table=[]
        for row in model['loadings']:
            self.table.append([[ball(F(r)*c).exp() for c in nodes] for r in row])
        self.evaluations=0
    def __call__(self,labels):
        ctx.prec=128;total=arb(0)
        for w,drift,table in zip(self.w,self.drifts,self.table):
            term=w*drift
            for j,label in enumerate(labels):term*=table[j][label]
            total+=term
        m=total.log();x=m/self.lam;xh=(m-self.hlog)/self.lam
        growth=100*(m+self.lam*self.lam/2).exp()
        value=growth*(Phi(x+self.lam)-Phi(xh+self.lam))-100*(Phi(x)-Phi(xh))+Phi(xh)
        self.evaluations+=1
        return clip_score(value)

def save_u64(path,values):
    a=array('Q',values)
    if a.itemsize!=8:raise RuntimeError('Platform does not support uint64 array')
    if sys.byteorder!='little':a.byteswap()
    if path.exists():raise FileExistsError(path)
    path.write_bytes(a.tobytes())
    return {'file':path.name,'bytes':path.stat().st_size,'sha256':sha(path),
            'format':'unsigned uint64 little endian','count':len(values)}

def transfer(model,mass,cells):
    ctx.prec=512;p=model['outer_dimension'];lam=F(model['lambda'])
    norm=payoff_norms(lam,F(1));M1=F(norm['M1_upper']);M2=F(norm['M2_upper'])
    b1=sum(abs_upper(c['_first']) for c in cells)
    b2=sum((c['_second']/2 for c in cells),arb(0))
    b2upper=endpoints(b2)[1];a_upper=endpoints(mass)[1]
    A=[max(abs(F(row[j])) for row in model['loadings']) for j in range(p)]
    rows=[{'coordinate':j,'maximum_full_asset_loading':str(v),
           'first_derivative_upper':str(M1*v),
           'second_derivative_upper':str((M1+M2)*v*v)} for j,v in enumerate(A)]
    beta=a_upper**(p-1)*sum(M1*v*b1+(M1+M2)*v*v*b2upper for v in A)
    Fball=mass**p;Fl,Fh=endpoints(Fball)
    return {'global_M1_upper':str(M1),'global_M2_upper':str(M2),
      'rho1_mass_upper':str(b1),'rho2_mass_interval':serial(b2),
      'coordinate_derivative_bounds':rows,'financial_grid_transfer_upper':str(beta),
      'retained_mass_interval':fraction_record(Fl,Fh),'bounded_tail_upper':str(1-Fl)}

def continuous_price(lo,hi,t):
    Fl=F(t['retained_mass_interval']['lower_exact_fraction']);Fh=F(t['retained_mass_interval']['upper_exact_fraction'])
    beta=F(t['financial_grid_transfer_upper'])
    a=max(F(0),Fl*lo-beta);b=min(F(1),Fh*hi+beta+1-Fl)
    ctx.prec=512;discount=ball(F(-3,100)).exp();dl,dh=endpoints(discount)
    return {'undiscounted_continuous_price_interval':fraction_record(a,b),
      'undiscounted_radius_upper':str((b-a)/2),'certifies_undiscounted_cent':(b-a)/2<=F(1,100),
      'discount_interval':serial(discount),
      'discounted_continuous_price_interval':fraction_record(dl*a,dh*b),
      'discounted_radius_upper':str((dh*b-dl*a)/2)}

def main():
    start=time.perf_counter();p=validate_protocol();protocol_sha=sha(PROTOCOL)
    models=json.loads((RESULTS/'dispersion_separation_protocol.json').read_text())['models']
    survivors=json.loads((RESULTS/'dispersion_separation_results.json').read_text())['records']
    selected=[next(m for m in models if m['model_id']==name) for name in CASE_IDS]
    for m in selected:
        parameters(m,False)
        r=next(r for r in survivors if r['model_id']==m['model_id'] and r['cap']=='1')
        if not r['separates_after_degree5']:raise ValueError('Selected model is not a retained separator')
    grid_start=time.perf_counter();ctx.prec=512;mass,cells=gaussian_cells(F(6),10,28,512)
    grid_seconds=time.perf_counter()-grid_start
    b1=sum(abs_upper(c['_first']) for c in cells);b2=sum((c['_second']/2 for c in cells),arb(0))
    rows=[]
    for m in selected:
        begin=time.perf_counter();t=transfer(m,mass,cells);score=Scores(m,cells)
        preparation=time.perf_counter()-begin;det=None
        if m['outer_dimension']==1:
            dstart=time.perf_counter();ss=[score((j,)) for j in range(1024)]
            dlo=[v[0] for v in ss];dhi=[v[1] for v in ss]
            lo=F(sum(dlo),1024*(1<<Q_SCORE));hi=F(sum(dhi),1024*(1<<Q_SCORE))
            det={'finite_grid_mean_interval':fraction_record(lo,hi),
              'scores_lower':save_u64(RESULTS/(m['model_id']+'_exactsum_lower_u64.bin'),dlo),
              'scores_upper':save_u64(RESULTS/(m['model_id']+'_exactsum_upper_u64.bin'),dhi),
              'evaluations':1024,'seconds_diagnostic':time.perf_counter()-dstart,
              **continuous_price(lo,hi,t)}
        sample_start=time.perf_counter();raw=os.urandom(N*m['outer_dimension']*2)
        raw_path=RESULTS/(m['model_id']+'_fresh_random_u16.bin')
        if raw_path.exists():raise FileExistsError(raw_path)
        raw_path.write_bytes(raw);words=array('H');words.frombytes(raw)
        if words.itemsize!=2:raise RuntimeError('No 16-bit array')
        if sys.byteorder!='little':words.byteswap()
        lower=[];upper=[];dimension=m['outer_dimension']
        for i in range(N):
            labels=tuple(words[i*dimension+j]&1023 for j in range(dimension))
            lo,hi=score(labels);lower.append(lo);upper.append(hi)
        sample_mean_lo=F(sum(lower),N*(1<<Q_SCORE));sample_mean_hi=F(sum(upper),N*(1<<Q_SCORE))
        finite_lo=max(F(0),sample_mean_lo-EPS);finite_hi=min(F(1),sample_mean_hi+EPS)
        ctx.prec=512;failure=2*ball(-2*N*EPS*EPS).exp()
        if endpoints(failure)[1]>F(1,3600):raise ArithmeticError('Hoeffding failure allowance violated')
        mc={'draws':N,'statistical_epsilon':str(EPS),'hoeffding_failure_upper':str(endpoints(failure)[1]),
          'sample_mean_interval':fraction_record(sample_mean_lo,sample_mean_hi),
          'finite_grid_mean_confidence_interval':fraction_record(finite_lo,finite_hi),
          'random_words':{'file':raw_path.name,'bytes':len(raw),'sha256':sha(raw_path),
                          'format':'unsigned uint16 little endian','coordinate_count':N*dimension},
          'scores_lower':save_u64(RESULTS/(m['model_id']+'_fresh_lower_u64.bin'),lower),
          'scores_upper':save_u64(RESULTS/(m['model_id']+'_fresh_upper_u64.bin'),upper),
          'maximum_outward_score_width':str(F(max(b-a for a,b in zip(lower,upper)),1<<Q_SCORE)),
          'evaluations':N,'seconds_diagnostic':time.perf_counter()-sample_start,
          **continuous_price(finite_lo,finite_hi,t)}
        row={'model_id':m['model_id'],'assets':m['assets'],'outer_dimension':dimension,
             'lambda':m['lambda'],'dispersion_scale':m['dispersion_scale'],'cap':'1',
             'transfer_certificate':t,'table_preparation_seconds_diagnostic':preparation,
             'complete_exactsum':det,'complete_monte_carlo':mc,
             'primary_complete_method':'exact finite-grid interval summation' if det else 'fixed-sample Hoeffding finite-grid mean confidence interval',
             'total_score_evaluations':score.evaluations,
             'joint_case_seconds_diagnostic':time.perf_counter()-begin}
        rows.append(row)
        print('completed',m['model_id'],'MC radius',float(F(mc['undiscounted_radius_upper'])),
              'exact radius',None if det is None else float(F(det['undiscounted_radius_upper'])),flush=True)
    validate_protocol()
    if sha(PROTOCOL)!=protocol_sha:raise RuntimeError('Protocol changed during execution')
    public_cells=[{k:v for k,v in c.items() if not k.startswith('_')} for c in cells]
    result={'schema':'qbp-four-case-complete-classical-pricing-v1','status':'complete',
      'protocol_sha256':protocol_sha,'producer_sha256':sha(Path(__file__)),
      'input_sha256':p['input_sha256'],'records':rows,'gaussian_grid':{'mass_interval':serial(mass),
        'rho1_mass_upper':str(b1),'rho2_mass_interval':serial(b2),'actual_cells':public_cells,
        'construction_seconds_diagnostic':grid_seconds},
      'all_four_primary_prices_certify_cent':all((r['complete_exactsum'] or r['complete_monte_carlo'])['certifies_undiscounted_cent'] for r in rows),
      'all_four_monte_carlo_prices_certify_cent':all(r['complete_monte_carlo']['certifies_undiscounted_cent'] for r in rows),
      'family_hoeffding_failure_allowance':'1/900','family_confidence_lower':'899/900',
      'interpretation':'Complete continuous prices under the stated independent uniform sampling premise; exact n2 sum is deterministic. These fresh outcomes supersede no retained trace and do not establish quantum advantage or empirical validation.',
      'generation_seconds_diagnostic':time.perf_counter()-start,
      'python_version':platform.python_version(),'python_flint_version':flint.__version__}
    write_json(RESULTS/'four_case_classical_pricing.json',result)
    compact={k:v for k,v in result.items() if k!='gaussian_grid'}
    compact['full_result_sha256']=sha(RESULTS/'four_case_classical_pricing.json')
    compact['full_result_file']='four_case_classical_pricing.json'
    write_json(RESULTS/'four_case_classical_summary.json',compact)
    print('all four primary prices complete',result['all_four_primary_prices_certify_cent'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true')
    if parser.parse_args().freeze:freeze()
    else:main()
