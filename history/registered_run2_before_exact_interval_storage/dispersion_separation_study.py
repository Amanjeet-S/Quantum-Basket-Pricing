"""Controlled separation stress: residual transfer versus finite Taylor prices.

The pre-outcome protocol extends established frozen GBM inputs by exact
Gaussian dispersion scaling. This compares two different certified outputs.
It neither establishes a lower bound against every classical method nor
executes or demonstrates a quantum advantage.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial,prod
import hashlib,json,time,platform
import flint
from flint import arb,ctx
import tensor_transfer as tt
from boundary_grid_certificate import parameters,ball,interval,serial,endpoints,abs_upper,phi,coefficient_certificate
from boundary_derivative_certificate import payoff_norms
from taylor_pricing import compositions,dot
ROOT=Path(__file__).resolve().parent
PRECISION=512
PRICE_TARGET=F(1,100)
TRANSFER_TARGET=F(1,200)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pos(v,name):return tt.positive(v,name)

class ContinuousGeometry:
    def __init__(self,model):
        self.w,self.d,self.s,self.mg,self.rg,self.lam=parameters(model,False)
        self.terms=[];self.weighted=[];self.raw=[];self.moments=[];self.phase_work=[]
        self.extend(4)
    def extend(self,maximum):
        ctx.prec=PRECISION;start=time.perf_counter();first=len(self.terms);exps=0
        for order in range(first,maximum+1):
            row=[];weighted=[]
            for alpha in compositions(order,len(self.w)):
                coefficient=F(factorial(order),prod(factorial(a) for a in alpha))*prod(w**a for w,a in zip(self.w,alpha))
                v=tuple(sum(F(a)*s[j] for a,s in zip(alpha,self.s)) for j in range(len(self.rg)))
                constant=sum(F(a)*d for a,d in zip(alpha,self.d));exponent=constant+dot(v,v)/2;shift=dot(self.rg,v)
                factor=ball(coefficient)*ball(exponent).exp();exps+=1
                row.append((alpha,coefficient,exponent,shift));weighted.append((alpha,factor,shift))
            self.terms.append(row);self.weighted.append(weighted);self.raw.append(sum((f for _,f,_ in weighted),arb(0)))
            value=sum(((-1)**(order-j)*comb(order,j)*self.raw[j] for j in range(order+1)),arb(0))
            self.moments.append(pos(value,'continuous U^'+str(order)))
        work={'first_new_moment_degree':first,'maximum_moment_degree':maximum,'new_Gaussian_geometry_exponential_evaluations':exps,
          'new_exact_monomial_geometries':exps,'seconds_diagnostic':time.perf_counter()-start}
        self.phase_work.append(work);return work
    def record(self,asset_count,scale):
        return {'assets':asset_count,'dispersion_scale':scale,'raw_T_moment_intervals':{str(k):serial(v) for k,v in enumerate(self.raw)},
          'continuous_U_moment_intervals':{str(k):serial(v) for k,v in enumerate(self.moments)},
          'monomial_counts_by_degree':list(map(len,self.terms)),'phase_work':self.phase_work}

class BoundaryKernels:
    def __init__(self,m,sigma,cap):
        self.m=F(m);self.sigma=sigma;self.cap=cap;self.cache={};self.polys={2:[arb(1)]}
        self.counts={'distinct_tilted_boundary_points':0,'Gaussian_density_evaluations':0,'falling_boundary_derivative_evaluations':0,'polynomial_coefficient_recurrence_updates':0}
    def ensure_poly(self,maximum):
        for r in range(max(self.polys),maximum):
            if r+1 in self.polys:continue
            p=self.polys[r];new=[arb(0)]*(len(p)+1)
            for j,c in enumerate(p):
                new[j]-=r*c;new[j+1]-=c/self.sigma
                if j:new[j-1]+=j*c/self.sigma
                self.counts['polynomial_coefficient_recurrence_updates']+=3 if j else 2
            self.polys[r+1]=new
    def polynomial(self,p,z):
        value=arb(0)
        for coefficient in reversed(p):value=value*z+coefficient
        return value
    def values(self,shift,maximum):
        shift=F(shift);self.ensure_poly(maximum)
        if shift not in self.cache:
            z=ball(self.m+shift)/self.sigma;f=phi(z);self.counts['Gaussian_density_evaluations']+=1
            if self.cap is None:z2=f2=None
            else:
                ell=(1+ball(self.cap)/100).log();z2=(ball(self.m+shift)-ell)/self.sigma;f2=phi(z2);self.counts['Gaussian_density_evaluations']+=1
            self.cache[shift]={'z':z,'f':f,'z2':z2,'f2':f2,'values':{}};self.counts['distinct_tilted_boundary_points']+=1
        c=self.cache[shift]
        for r in range(2,maximum+1):
            if r in c['values']:continue
            value=100*self.polynomial(self.polys[r],c['z'])*c['f']/self.sigma
            if self.cap is not None:value-=(100+ball(self.cap))*self.polynomial(self.polys[r],c['z2'])*c['f2']/self.sigma
            c['values'][r]=value;self.counts['falling_boundary_derivative_evaluations']+=1
        return c['values']

class ContinuousCase:
    def __init__(self,model,cap,geometry):
        self.model=model;self.cap=cap;self.geometry=geometry;self.w,self.d,self.s,self.mg,self.rg,self.lam=parameters(model,False)
        sigma=ball(self.lam*self.lam+dot(self.rg,self.rg)).sqrt();self.kernels=BoundaryKernels(self.mg,sigma,cap)
        self.base=coefficient_certificate(model,cap,'vanilla' if cap is None else 'capped',1)
        self.coefficients=[interval(F(self.base[n]['lower_exact_fraction']),F(self.base[n]['upper_exact_fraction'])) for n in ('c0','c1')]
        self.norms=payoff_norms(self.lam,cap);self.orders=[];self.phase_work=[]
        self.EA=100*ball((F(model['risk_free_rate'])-F(model['dividend_yield']))*F(model['maturity'])).exp()
    def extend(self,maximum):
        start=time.perf_counter();previous=len(self.coefficients)-1;before=dict(self.kernels.counts);term_visits=0
        for k in range(max(2,previous+1),maximum+1):
            coefficient=arb(0)
            for j in range(k+1):
                row=arb(0)
                for _,factor,shift in self.geometry.weighted[j]:
                    row+=factor*self.kernels.values(shift,maximum)[k];term_visits+=1
                coefficient+=(-1)**(k-j)*comb(k,j)*row
            self.coefficients.append(coefficient/factorial(k))
        partial=sum(self.coefficients[:previous+1],arb(0)) if previous else self.coefficients[0]
        first=max(1,previous+1)
        if not self.orders:partial=self.coefficients[0];first=1
        for k in range(first,maximum+1):
            partial+=self.coefficients[k];M=F(self.norms['C_r_upper'][str(k+1)])
            radius=M*endpoints(self.geometry.moments[k+1])[1]/factorial(k+1)
            lo,hi=endpoints(partial);lo=max(F(0),lo-(F(0) if self.cap is None and k==1 else radius));hi+=radius
            if self.cap is not None:hi=min(F(self.cap),hi)
            if lo>hi:raise ArithmeticError('Continuous price enclosure empty')
            self.orders.append({'order':k,'coefficient_interval':serial(self.coefficients[k]),'partial_sum_interval':serial(partial),
              'global_C_r_norm_upper':str(M),'remainder_U_moment_interval':serial(self.geometry.moments[k+1]),'remainder_upper':str(radius),
              'continuous_price_interval':serial(interval(lo,hi)),'undiscounted_price_radius_upper':str((hi-lo)/2),
              'certifies_undiscounted_cent':(hi-lo)/2<=PRICE_TARGET})
        work={'starting_coefficient_degree':previous,'maximum_coefficient_degree':maximum,'weighted_kernel_term_visits':term_visits,
          'incremental_boundary_operation_counts':{k:self.kernels.counts[k]-before[k] for k in before},
          'seconds_diagnostic':time.perf_counter()-start}
        self.phase_work.append(work);return work
    def record(self):
        # Every retained interval contains the same continuous price. Exact
        # rational intersection is valid even when its constituent Taylor
        # bounds are individually wider than the requested radius.
        lo=F(0);hi=F(self.cap) if self.cap is not None else endpoints(self.EA)[1]
        intrinsic={'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)};progress=[]
        for row in self.orders:
            lo=max(lo,F(row['continuous_price_interval']['lower_exact_fraction']))
            hi=min(hi,F(row['continuous_price_interval']['upper_exact_fraction']))
            if lo>hi:raise ArithmeticError('Valid continuous price interval intersection is empty')
            progress.append({'through_degree':row['order'],'price_interval':{'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)},
              'price_radius_upper':str((hi-lo)/2),'certifies_undiscounted_cent':(hi-lo)/2<=PRICE_TARGET})
        return {'model_id':self.model['model_id'],'payoff':'vanilla' if self.cap is None else 'capped','cap':None if self.cap is None else str(self.cap),
          'c0_c1_certificate':self.base,'orders':list(self.orders),'first_successful_order':next((r['order'] for r in self.orders if r['certifies_undiscounted_cent']),None),
          'best_price_radius_upper':str(min(F(r['undiscounted_price_radius_upper']) for r in self.orders)),
          'intrinsic_price_interval':intrinsic,'EA_interval':serial(self.EA),'intrinsic_forward_exponential_evaluations':1,
          'intersection_progress':progress,'intersection_price_interval':progress[-1]['price_interval'],
          'intersection_radius_upper':progress[-1]['price_radius_upper'],'intersection_certifies_cent':progress[-1]['certifies_undiscounted_cent'],
          'intersection_first_successful_order':next((r['through_degree'] for r in progress if r['certifies_undiscounted_cent']),None),
          'phase_work':list(self.phase_work),'boundary_operation_counts':dict(self.kernels.counts)}


def main():
    ctx.prec=PRECISION;start=time.perf_counter();protocol_path=ROOT/'results/dispersion_separation_protocol.json';protocol_sha=sha(protocol_path)
    p=json.loads(protocol_path.read_text());models=p['models'];source_files=('dispersion_separation_study.py','tensor_transfer.py','boundary_grid_certificate.py','boundary_derivative_certificate.py','higherdim_preparation.py','taylor_pricing.py')
    source_before={n:sha(ROOT/n) for n in source_files};norm_before={n:sha(ROOT/'results'/n) for n in ('boundary_derivative_norms.json','gaussian_derivative_norms.json')}
    trans_start=time.perf_counter();trans=tt.TransportMGF(6,10,28);transport_seconds=time.perf_counter()-trans_start
    raw_transport=[]
    for c in trans.cells:
        raw_transport.append({'cell':c['cell'],'actual_node':c['actual_node'],'actual_node_integer':c['actual_node_integer'],
          'left':c['left'],'right':c['right'],'mass':c['mass'],'first_moment':c['first_moment'],'second_moment':c['second_moment'],
          'left_cdf':c['left_cdf'],'right_cdf':c['right_cdf'],'rho1_positive_atom_weight':str(abs_upper(c['_first']))})
    groups={};geometry_records=[];tails=[];records=[];continuous_geometry={};timings=[];extended_keys=set()
    for n in (2,4,8):
        for scale in p['scales']:
            ms=[m for m in models if m['assets']==n and m['dispersion_scale']==scale];m0=ms[0];base=parameters(m0,False)
            for m in ms:
                q=parameters(m,False)
                if q[:3]!=base[:3] or q[4]!=base[4]:raise ArithmeticError('Scaled geometry reuse premise failed')
            group_start=time.perf_counter();hybrid=tt.HybridMoments(m0,trans).record();tail=tt.tail_moment(m0)
            geometry_records.append({'assets':n,'dispersion_scale':scale,**hybrid});tails.append({'assets':n,'dispersion_scale':scale,**tail})
            g=ContinuousGeometry(m0);continuous_geometry[(n,scale)]=g
            for model in ms:
                for cap in (F(1),F(10),F(25),None):
                    begin=time.perf_counter();transfer=tt.financial_transfer(model,cap,hybrid,tail)
                    transfer['dispersion_scale']=scale;transfer['original_model_id']=model['original_model_id']
                    case=ContinuousCase(model,cap,g);case.extend(3);degree3=case.record();degree3['dispersion_scale']=scale;degree3['original_model_id']=model['original_model_id']
                    separator=F(transfer['financial_transfer_upper'])<TRANSFER_TARGET and not degree3['intersection_certifies_cent']
                    later=None
                    if separator:
                        if (n,scale) not in extended_keys:g.extend(6);extended_keys.add((n,scale))
                        case.extend(5);later=case.record()
                    records.append({'model_id':model['model_id'],'original_model_id':model['original_model_id'],'assets':n,'dispersion_scale':scale,
                      'lambda':model['lambda'],'payoff':transfer['payoff'],'cap':transfer['cap'],'aggregate_R2_transport':transfer,
                      'degree3_continuous_price':degree3,'selected_for_degree5':separator,'degree5_continuous_price':later,
                      'separates_after_degree3':separator,'separates_after_degree5':separator and not later['intersection_certifies_cent'],
                      'joint_case_evaluation_seconds_diagnostic':time.perf_counter()-begin})
            elapsed=time.perf_counter()-group_start;timings.append({'assets':n,'dispersion_scale':scale,'joint_hybrid_tail_continuous_geometry_and12cases_seconds_diagnostic':elapsed})
            print('completed',n,scale,'seconds',elapsed,flush=True)
    # Scale one must reproduce previous mathematical inputs and outputs.
    original=json.loads((ROOT/'results/tensor_transfer.json').read_text())['records'];agreements=[]
    prior=json.loads((ROOT/'results/original_degree3_price_intervals.json').read_text())['records']
    for r in records:
        if r['dispersion_scale']!='1':continue
        old=next(v for v in original if v['model_id']==r['original_model_id'] and v['cap']==r['cap'])
        error=abs(F(old['financial_transfer_upper'])-F(r['aggregate_R2_transport']['financial_transfer_upper']))
        cp=next(v for v in prior if v['model_id']==r['original_model_id'] and v['cap']==r['cap'])
        discrepancies=[]
        for neworder,oldorder in zip(r['degree3_continuous_price']['orders'],cp['orders']):
            oldradius=F(oldorder['undiscounted_radius_upper']);newradius=F(neworder['undiscounted_price_radius_upper'])
            discrepancies.append(abs(oldradius-newradius))
        ok=error<F(1,10**60) and max(discrepancies)<F(1,10**60)
        if not ok:raise ArithmeticError('Scale-one agreement with original36 failed')
        agreements.append({'model_id':r['model_id'],'cap':r['cap'],'transfer_upper_absolute_difference':str(error),
          'degree3_radius_max_absolute_difference':str(max(discrepancies)),'outward_tolerance':'1e-60','passed':ok})
    if sha(protocol_path)!=protocol_sha or source_before!={n:sha(ROOT/n) for n in source_files} or norm_before!={n:sha(ROOT/'results'/n) for n in norm_before}:raise RuntimeError('A registered input or producer changed during execution')
    summaries=[]
    for scale in p['scales']:
        rows=[r for r in records if r['dispersion_scale']==scale]
        summaries.append({'dispersion_scale':scale,'cases':len(rows),'transport_passes':sum(r['aggregate_R2_transport']['status']=='passed' for r in rows),
          'degree3_price_passes':sum(r['degree3_continuous_price']['intersection_certifies_cent'] for r in rows),
          'degree3_individual_interval_price_passes':sum(r['degree3_continuous_price']['first_successful_order'] is not None for r in rows),'degree3_separators':sum(r['separates_after_degree3'] for r in rows),'degree5_evaluated':sum(r['selected_for_degree5'] for r in rows),
          'degree5_rescued':sum(r['selected_for_degree5'] and r['degree5_continuous_price']['intersection_certifies_cent'] for r in rows),
          'remaining_degree5_separators':sum(r['separates_after_degree5'] for r in rows)})
    result={'schema':'qbp-controlled-dispersion-separation-v1','status':'complete; every allocation/price failure and targeted extension retained',
      'protocol_sha256':protocol_sha,'source_sha256':source_before,'norm_input_sha256':norm_before,'records':records,'summary_by_scale':summaries,
      'transport_actual_inputs':raw_transport,'transport_profile':trans.profile(),'shared_transport_construction_seconds_diagnostic':transport_seconds,
      'hybrid_geometry_records':geometry_records,'tail_records':tails,
      'continuous_geometry_records':[g.record(n,scale) for (n,scale),g in continuous_geometry.items()],
      'joint_group_timings':timings,'scale1_agreement_records':agreements,'generation_seconds_diagnostic':time.perf_counter()-start,
      'python_version':platform.python_version(),'python_flint_version':flint.__version__,
      'timing_scope':'Shared transport setup separate; group times combinehybridmoments/fulltails/continuousmoments/12payoffs and admitteddegree5 extensions. Classical phase counters/times separateinsideeachcase andgeometry. No cold-vs-warm or cross-output speedratio claimed.',
      'claim_boundary':'An aggregate actual-grid residual transfer guarantee maypass while a declared finite Taylor price comparator fails. This is separation of outputs, not a quantum speed-up, market validation, generic algorithm priority or lower bound against classical pricing.'}
    output=ROOT/'results/dispersion_separation_results.json';output.write_text(json.dumps(result,indent=2,default=str)+'\n')
    compact={k:v for k,v in result.items() if k not in ('records','transport_actual_inputs','hybrid_geometry_records','tail_records','continuous_geometry_records')}
    compact['records']=[{'model_id':r['model_id'],'original_model_id':r['original_model_id'],'assets':r['assets'],'dispersion_scale':r['dispersion_scale'],'lambda':r['lambda'],'cap':r['cap'],
      'transfer_upper':r['aggregate_R2_transport']['financial_transfer_upper'],'transport_pass':r['aggregate_R2_transport']['status']=='passed',
      'degree3_best_radius_upper':r['degree3_continuous_price']['best_price_radius_upper'],'degree3_first_successful_order':r['degree3_continuous_price']['first_successful_order'],
      'degree3_intersection_radius_upper':r['degree3_continuous_price']['intersection_radius_upper'],'degree3_intersection_first_successful_order':r['degree3_continuous_price']['intersection_first_successful_order'],
      'degree5_best_radius_upper':None if not r['selected_for_degree5'] else r['degree5_continuous_price']['best_price_radius_upper'],
      'degree5_first_successful_order':None if not r['selected_for_degree5'] else r['degree5_continuous_price']['first_successful_order'],
      'degree5_intersection_radius_upper':None if not r['selected_for_degree5'] else r['degree5_continuous_price']['intersection_radius_upper'],'degree5_intersection_first_successful_order':None if not r['selected_for_degree5'] else r['degree5_continuous_price']['intersection_first_successful_order'],
      'separates_after_degree3':r['separates_after_degree3'],'separates_after_degree5':r['separates_after_degree5']} for r in records]
    compact['full_result_sha256']=sha(output);compact['full_result_file']=output.name
    (ROOT/'results/dispersion_separation_summary.json').write_text(json.dumps(compact,indent=2)+'\n')
    print(json.dumps(summaries,indent=2),flush=True)
if __name__=='__main__':main()
