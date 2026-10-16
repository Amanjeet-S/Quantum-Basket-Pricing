"""Matched positive-measure block certificates for an actual rounded CDF grid.

A block refines Cauchy's inequality under positive first- or second-order
transport measures. This is established measure/quadrature mathematics. The
experiment tests the price-independent residual envelope, not priority.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import argparse,gzip,hashlib,json,time
from flint import arb,ctx
import tensor_transfer as tt
from boundary_grid_certificate import ball,endpoints,serial,upper,abs_upper,phi,Phi,parameters
from boundary_derivative_certificate import payoff_norms

ROOT=Path(__file__).resolve().parent
PRECISION=512
PARTITIONS=(1,8,32,128,1024)
ORDER=120

def positive(x,what):return tt.positive(x,what)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def normal_moments(a,b,n):
    """Outward integrals of z**k phi(z), k=0,...,n."""
    pa,pb=phi(a),phi(b);out=[Phi(b)-Phi(a),pa-pb]
    ap=bp=arb(1)
    for k in range(2,n+1):
        ap*=a;bp*=b;out.append(ap*pa-bp*pb+(k-1)*out[k-2])
    return out[:n+1]

class BlockTransport:
    def __init__(self,L=6,b_y=10,q=28,order=ORDER):
        ctx.prec=PRECISION;start=time.perf_counter()
        self.base=tt.TransportMGF(L,b_y,q)
        self.L=F(L);self.N=len(self.base.cells);self.order=order
        self.cells=self.base.cells;self.cache={};self.fallbacks=[]
        self.counters={'normal_moment_recurrences':0,'prefix_coefficient_additions':0,
          'block_polynomial_horner_steps':0,'block_continuous_erf_evaluations':0,
          'block_continuous_exponential_evaluations':0,'grouped_geometry_coefficients':0,
          'block_geometry_dot_multiplications':0,'geometry_factor_exponential_evaluations':0}
        # Prefixes share coefficients across all nested contiguous partitions.
        self.prefix={k:[[arb(0)] for _ in range(order+1)] for k in ('rounded','signed','rho1','rho0')}
        self.second_mass=[arb(0)];self.cell_metadata=[]
        cellmass=self.base.mass/self.N
        for i,c in enumerate(self.cells):
            y=F(c['actual_node']);yb=ball(y);a,b=c['_left'],c['_right']
            pa,pb=c['_left_cdf'],c['_right_cdf'];pc=Phi(yb)
            left=normal_moments(a,yb,order+1);right=normal_moments(yb,b,order+1)
            self.counters['normal_moment_recurrences']+=2*order
            twice=2*pc-pa-pb;first_upper=abs_upper(c['_first'])
            coef={k:[] for k in self.prefix};power=arb(1)
            for r in range(order+1):
                coef['rounded'].append(cellmass*power)
                coef['signed'].append(c['_first']*power)
                coef['rho1'].append(ball(first_upper)*power)
                # Integral of z**r under positive FTC displacement density.
                coef['rho0'].append((yb**(r+1)*twice+right[r+1]-left[r+1])/factorial(r+1))
                power*=yb/(r+1)
                for name in self.prefix:
                    self.prefix[name][r].append(self.prefix[name][r][-1]+coef[name][-1])
                    self.counters['prefix_coefficient_additions']+=1
            self.second_mass.append(self.second_mass[-1]+c['_second']/2)
            rho0mass=positive(coef['rho0'][0],'FTC cell mass')
            self.cell_metadata.append({'cell':i,'actual_node':str(y),
              'first_moment_interval':serial(c['_first']),'rho1_positive_atom_weight':str(first_upper),
              'rho0_mass_interval':serial(rho0mass),'rho2_mass_interval':serial(c['_second']/2)})
        self.elapsed=time.perf_counter()-start

    def blocks(self,B):
        if B<1 or self.N%B:raise ValueError('Partition must divide actual cell count')
        width=self.N//B
        return [(i*width,(i+1)*width) for i in range(B)]

    def coefficients(self,name,start,end):
        return [p[end]-p[start] for p in self.prefix[name]]

    def polynomial(self,name,rate,start,end):
        cs=self.coefficients(name,start,end);v=ball(rate);value=arb(0)
        for c in reversed(cs):value=value*v+c
        self.counters['block_polynomial_horner_steps']+=len(cs)
        if name=='signed':
            variation=self.prefix['rho1'][0][end]-self.prefix['rho1'][0][start]
        else:variation=cs[0]
        argument=ball(abs(rate)*self.L)
        rem=positive(variation,'series total variation')*argument.exp()*argument**(self.order+1)/factorial(self.order+1)
        value+=arb(0,ball(upper(rem)).upper())
        return value

    def continuous_block(self,rate,start,end):
        v=ball(rate);a=self.cells[start]['_left'];b=self.cells[end-1]['_right']
        self.counters['block_continuous_erf_evaluations']+=2
        self.counters['block_continuous_exponential_evaluations']+=1
        return (v*v/2).exp()*(Phi(b-v)-Phi(a-v))

    def vector(self,kind,rate,B):
        rate=F(rate);key=(kind,rate,B)
        if key not in self.cache:
            out=[]
            for start,end in self.blocks(B):
                if kind in (0,1):
                    name='rho0' if kind==0 else 'rho1'
                    value=positive(self.polynomial(name,rate,start,end),name+' block MGF')
                elif kind==2:
                    mass=positive(self.second_mass[end]-self.second_mass[start],'Peano block mass')
                    if not rate:value=mass
                    elif abs(rate)<F(1,10**20):
                        lo,hi=endpoints(mass);argument=ball(abs(rate)*self.L)
                        value=tt.interval(endpoints(ball(lo)*(-argument).exp())[0],upper(ball(hi)*argument.exp()))
                        self.fallbacks.append({'kind':2,'rate':str(rate),'B':B,'block_start':start,'block_end':end})
                    else:
                        v=ball(rate)
                        value=positive((self.continuous_block(rate,start,end)-self.polynomial('rounded',rate,start,end)-v*self.polynomial('signed',rate,start,end))/(v*v),'Peano block MGF')
                else:raise ValueError('Transport kind must be 0,1,2')
                out.append(value)
            self.cache[key]=tuple(out)
        return self.cache[key]

    def mass_record(self,B):
        out=[]
        for k,(start,end) in enumerate(self.blocks(B)):
            out.append({'block':k,'cell_start':start,'cell_end_exclusive':end,
              'rho0_mass_interval':serial(self.vector(0,F(0),B)[k]),
              'rho1_mass_interval':serial(self.vector(1,F(0),B)[k]),
              'rho2_mass_interval':serial(self.vector(2,F(0),B)[k])})
        return out

    def profile(self):
        return {'L':str(self.L),'b_y':self.base.b_y,'q':self.base.q,'N':self.N,
          'precision_bits':PRECISION,'series_order':self.order,'partitions':list(PARTITIONS),
          'node_certificate':self.base.grid.node_certificate,
          'retained_mass_interval':serial(self.base.mass),'construction_seconds_diagnostic':self.elapsed,
          'operation_counts':self.counters,'distinct_block_MGF_queries':len(self.cache),
          'actual_block_MGF_evaluations':sum(len(v) for v in self.cache.values()),
          'tiny_rate_support_fallbacks':self.fallbacks,
          'construction_scope':'Original rounded nodes and masses. Prefix factorial moments and grouped scalar rate queries; no Cartesian grid, payoff price or derivative samples.'}

class GroupedGeometry:
    def __init__(self,model,transport):
        self.model=model;self.transport=transport
        self.w,self.d,self.s,self.mg,self.rg,self.lam=parameters(model,False)
        self.p=len(self.rg);self.group_cache={};self.alpha_cache={}

    def alpha(self,alpha):
        if alpha not in self.alpha_cache:
            constant=sum(F(a)*d for a,d in zip(alpha,self.d))
            rates=tuple(sum(F(a)*si[k] for a,si in zip(alpha,self.s)) for k in range(self.p))
            self.alpha_cache[alpha]=(ball(constant).exp(),rates)
            self.transport.counters['geometry_factor_exponential_evaluations']+=1
        return self.alpha_cache[alpha]

    def grouped(self,poly,j):
        key=(tuple(sorted(poly.items())),j)
        if key not in self.group_cache:
            groups={}
            for alpha,c in poly.items():
                value,rates=self.alpha(alpha);value=ball(c)*value
                for k,rate in enumerate(rates):
                    if k<j:value*=self.transport.base.continuous(rate)
                    elif k>j:value*=self.transport.base.rounded(rate)
                r=rates[j];groups[r]=groups.get(r,arb(0))+value
            self.transport.counters['grouped_geometry_coefficients']+=len(groups)
            self.group_cache[key]=groups
        return self.group_cache[key]

    def moments(self,poly,j,kind,B,name):
        groups=self.grouped(poly,j);out=[arb(0) for _ in range(B)]
        for rate,coefficient in groups.items():
            vec=self.transport.vector(kind,rate,B)
            for k,value in enumerate(vec):out[k]+=coefficient*value
        self.transport.counters['block_geometry_dot_multiplications']+=len(groups)*B
        return [positive(v,name+' block '+str(k)) for k,v in enumerate(out)]

    def record(self,B):
        rows=[]
        for j in range(self.p):
            polys=tt.geometry_polynomials(self.model,j)
            blocks=[{'block':k,'geometric_loading':str(self.rg[j]),'moments':{str(kind):{} for kind in (0,1,2)}} for k in range(B)]
            for kind in (0,1,2):
                # rho0 only needs the first derivative's positive moments.
                names=('U2','J2','U2_J2','U4') if kind==0 else tuple(polys)
                for name in names:
                    values=self.moments(polys[name],j,kind,B,name)
                    for k,v in enumerate(values):blocks[k]['moments'][str(kind)][name]=serial(v)
            rows.append({'coordinate':j,'blocks':blocks})
        return {'model_id':self.model['model_id'],'B':B,'coordinates':rows,
           'distinct_geometry_alphas':len(self.alpha_cache),'grouped_rate_counts':[len(v) for v in self.group_cache.values()]}

def first_only(norms,row,kind=0):
    """Exactly the existing first-derivative envelope under another measure."""
    m={k:tt.B(v) for k,v in row['moments'][str(kind)].items()}
    M2=F(norms['M2_upper']);M3=F(norms['C_r_upper']['3'])
    D1=M2*F(norms['scalar_boundary_density_first_upper']) if norms['payoff']=='vanilla' else 2*M2*F(norms['scalar_boundary_h_first_upper'])
    r=abs(F(row['geometric_loading']));uj=upper((m['U2']*m['J2']).sqrt())
    v1=upper((m['U2']*m['U2_J2']).sqrt());v2=upper((m['U4']*m['J2']).sqrt())
    c={'U_Tj':M2*uj,'r_U2':r*D1*upper(m['U2'])/2,'U2_Tj':M3*min(v1,v2)/6}
    return {'error_upper':str(sum(c.values(),F(0))),'components_upper':{k:str(v) for k,v in c.items()},
       'cauchy_variant_a_upper':str(v1),'cauchy_variant_b_upper':str(v2)}

def transfer(model,cap,geometry,tail):
    norms=payoff_norms(F(model['lambda']),cap);sums={k:F(0) for k in (0,1,2)};rows=[]
    for row in geometry['coordinates']:
        sub=[];coord={k:F(0) for k in (0,1,2)}
        for block in row['blocks']:
            envelopes,n=tt.derivative_envelope(norms,block)
            ftc=first_only(norms,block,0)
            coord[0]+=F(ftc['error_upper'])
            for kind in (1,2):coord[kind]+=F(envelopes[str(kind)]['error_upper'])
            sub.append({'block':block['block'],'first_order_FTC':ftc,'centred_first_atoms':envelopes['1'],'second_Peano':envelopes['2']})
        for k in sums:sums[k]+=coord[k]
        rows.append({'coordinate':row['coordinate'],'component_sums_upper':{str(k):str(v) for k,v in coord.items()},'blocks':sub})
    tailup=F(norms['M2_upper'])*F(tail['U2_outside_interval']['upper_exact_fraction'])/2
    return {'model_id':model['model_id'],'cap':None if cap is None else str(cap),'payoff':norms['payoff'],'order':2,'B':geometry['B'],
      'FTC_transport_upper':str(sums[0]),'centred_first_transport_upper':str(sums[1]),
      'second_Peano_transport_upper':str(sums[2]),'tail_upper':str(tailup),
      'FTC_financial_transfer_upper':str(sums[0]+tailup),
      'second_order_financial_transfer_upper':str(sums[1]+sums[2]+tailup),'coordinates':rows,
      'norms':n,'target_allocation':'1/200',
      'FTC_status':'passed' if sums[0]+tailup<F(1,200) else 'retained allocation failure',
      'second_order_status':'passed' if sums[1]+sums[2]+tailup<F(1,200) else 'retained allocation failure'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'results/block_certificate_comparison.json.gz')
    ap.add_argument('--summary',type=Path,default=ROOT/'results/block_certificate_summary.json')
    ap.add_argument('--partitions',default=','.join(map(str,PARTITIONS)));args=ap.parse_args()
    parts=tuple(map(int,args.partitions.split(',')))
    if not parts or sorted(set(parts))!=list(parts) or any(b not in PARTITIONS for b in parts):raise ValueError('Declared nested partitions required')
    start=time.perf_counter();ctx.prec=PRECISION;trans=BlockTransport();models=json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']
    records=[];geometries=[];timings=[];tails=[]
    for n in (2,4,8):
        ms=[m for m in models if m['assets']==n];engine=GroupedGeometry(ms[0],trans);tail=tt.tail_moment(ms[0]);tails.append({'assets':n,**tail})
        p0=parameters(ms[0],False)
        for m in ms:
            pm=parameters(m,False)
            if pm[:3]!=p0[:3] or pm[4]!=p0[4]:raise ArithmeticError('Geometry reuse premise failed')
        for B in parts:
            begun=time.perf_counter();geometry=engine.record(B);geometries.append(geometry)
            for model in ms:
                for cap in (F(1),F(10),F(25),None):
                    record=transfer(model,cap,geometry,tail)
                    previous=[r for r in records if r['model_id']==model['model_id'] and r['cap']==record['cap']]
                    for kind,raw in (('FTC','FTC_financial_transfer_upper'),('second_order','second_order_financial_transfer_upper')):
                        parent=min((F(r[kind+'_hierarchy_upper']) for r in previous),default=F(record[raw]))
                        chosen=min(parent,F(record[raw]));record[kind+'_hierarchy_upper']=str(chosen)
                        record[kind+'_hierarchy_status']='passed' if chosen<F(1,200) else 'retained allocation failure'
                    records.append(record)
            elapsed=time.perf_counter()-begun;timings.append({'assets':n,'B':B,'geometry_and_12_payoffs_seconds_diagnostic':elapsed})
            print('completed assets',n,'B',B,'seconds',elapsed,flush=True)
    source_files=('block_certificate_comparison.py','tensor_transfer.py','higherdim_preparation.py','boundary_grid_certificate.py','boundary_derivative_certificate.py')
    full={'schema':'qbp-matched-block-certificate-v1','status':'complete; all adverse allowances retained',
      'records':records,'geometry_records':geometries,'tail_records':tails,'transport_profile':trans.profile(),
      'cell_masses':trans.cell_metadata,'block_masses':{str(B):trans.mass_record(B) for B in parts},
      'timings':timings,'generation_seconds_diagnostic':time.perf_counter()-start,
      'source_sha256':{f:digest(ROOT/f) for f in source_files},'model_protocol_sha256':digest(ROOT/'results/generalisation_protocol.json'),
      'matched_inputs':'36 frozen GBM/payoff cases; identical original rounded equal-CDF nodes, L6, by10, q28, 512-bit outward arithmetic, global boundary norms, geometry envelopes and tails',
      'scope':'R2 only. FTC first order and centred second order with nested positive-measure Cauchy blocks. Paid shared preparation. No tensor-grid enumeration or unknown price input. No priority or quantum advantage claim.'}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(args.out,'wt') as f:json.dump(full,f,separators=(',',':'))
    summary={k:v for k,v in full.items() if k not in ('geometry_records','records','cell_masses','block_masses')}
    summary['records']=[{k:v for k,v in r.items() if k!='coordinates'} for r in records]
    summary['archive_sha256']=digest(args.out);summary['archive_file']=args.out.name
    args.summary.write_text(json.dumps(summary,indent=2)+'\n')
    print('saved',args.out,'elapsed',time.perf_counter()-start,flush=True)

if __name__=='__main__':main()
