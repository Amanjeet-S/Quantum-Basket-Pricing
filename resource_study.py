"""Matched-bias-budget study of three complete finite evaluator constructions."""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from encoding_resources import (Graph,build_graph,make_tables,counted_composition,
                               check_small,decimal_cdf,PI,INV_SQRT_2PI,
                               function_span,capped_kernel_decimal)

D=Decimal
ROOT=Path(__file__).resolve().parent
ALLOWANCE=D('1e-38')
KERNEL_ALLOWANCE=D('1e-34')


@dataclass
class Bound:
    lo: Decimal
    hi: Decimal
    error: Decimal


def assert_ordered_cdf_difference(graph, first, second):
    """Prove the retained ordered-CDF premise from exact integer structure.

    Nonlinear/rounded nodes are independent atoms. Only exact addition,
    subtraction and integer scaling are expanded. The accepted difference
    must be a nonnegative constant smaller than the CDF input span; the
    propagator separately checks both arguments are in [0,16).
    """
    def linear(name):
        node=graph.nodes[name]
        if node.kind=='constant':
            value=node.meta
            if value >= 1 << (graph.width-1):value-=1 << graph.width
            return {'constant':value}
        if node.kind in ('add','sub'):
            left=linear(node.args[0]);right=linear(node.args[1])
            sign=1 if node.kind=='add' else -1
            for key,value in right.items():left[key]=left.get(key,0)+sign*value
            return {key:value for key,value in left.items() if value}
        if node.kind=='scale' and node.meta[1]==0:
            return {key:value*node.meta[0] for key,value in linear(node.args[0]).items()}
        return {name:1}
    a,b=graph.nodes[first],graph.nodes[second]
    if not (a.kind==b.kind=='function' and a.meta==b.meta=='cdf'):
        raise AssertionError(('unproved ordered CDF subtraction',first,second))
    difference=linear(a.args[0]);right=linear(b.args[0])
    for key,value in right.items():difference[key]=difference.get(key,0)-value
    difference={key:value for key,value in difference.items() if value}
    delta=difference.pop('constant',0)
    if difference or not 0<=delta<16*(1 << graph.q):
        raise AssertionError(('unproved ordered CDF subtraction',first,second,delta,difference))
    return delta


def propagate(graph,cap):
    """Directed conservative analytical bounds using 70-digit decimal centres.

    A fixed 1e-38 outward allowance per non-rational bound operation is
    enormously larger than its 70-digit correctly-rounded centre error.
    It is included in, rather than substituted for, interpolation error.
    """
    q=graph.q
    unit=D(1)/D(1 << q)
    vals={}
    domains=[]
    with localcontext() as ctx:
        ctx.prec=70
        for name,node in graph.nodes.items():
            a=[vals[v] for v in node.args if v in vals]
            if node.kind=='constant':
                real=graph.ideal_constants[name]
                v=Bound(real,real,unit/2+ALLOWANCE)
            elif node.kind=='quantile':
                # For b<=18, Phi(-4.75)<1/(2*2^b), checked separately.
                v=Bound(D('.25'),D('9.75'),unit)
            elif node.kind=='uniform':
                n=D(1 << node.bits)
                v=Bound(1/(2*n),1-1/(2*n),D(0))
            elif node.kind in ('add','sub','positive_sub'):
                if node.kind=='add':lo,hi=a[0].lo+a[1].lo,a[0].hi+a[1].hi
                else:lo,hi=a[0].lo-a[1].hi,a[0].hi-a[1].lo
                if node.kind=='positive_sub':
                    # Both operands are CDFs at ordered arguments t0>=t1.
                    # The monotone integer interpolator retains that order.
                    assert_ordered_cdf_difference(graph,*node.args)
                    lo=max(D(0),lo);hi=min(D(1),hi)
                v=Bound(lo,hi,a[0].error+a[1].error+ALLOWANCE)
            elif node.kind=='scale':
                k=D(node.meta[2]);integer,drop,_=node.meta
                ck=unit/2 if drop else D(0)
                if a[0].lo < 0 or k < 0:raise AssertionError(('negative constant product',name))
                v=Bound(a[0].lo*k,a[0].hi*k,
                        k*a[0].error+a[0].hi*ck+a[0].error*ck+(unit if drop else D(0))+ALLOWANCE)
            elif node.kind=='shift':
                div=D(1 << node.meta)
                v=Bound(a[0].lo/div,a[0].hi/div,a[0].error/div+unit)
            elif node.kind=='mul':
                if a[0].lo < 0 or a[1].lo < 0:
                    raise AssertionError(('unsigned multiply negative mathematical input',name,a))
                lo=a[0].lo*a[1].lo;hi=a[0].hi*a[1].hi
                error=a[0].hi*a[1].error+a[1].hi*a[0].error+a[0].error*a[1].error+unit
                v=Bound(lo,hi,error+ALLOWANCE)
            elif node.kind=='clamp':
                ceiling=D(node.meta)/D(1 << q)
                # Ideal comparison cap differs by <=.5 unit from its exact
                # requested cap, recorded by the graph's rounded threshold.
                # Using rounded cap consistently changes price by <=.5 unit.
                lo=min(max(a[0].lo,D(0)),max(D(0),ceiling-unit/2))
                hi=min(max(a[0].hi,D(0)),ceiling+unit/2)
                v=Bound(lo,hi,a[0].error+unit/2)
            elif node.kind=='function':
                aa=a[0]
                generator_allowance=ALLOWANCE
                span=D(function_span(node.meta))
                spacing=span/D(1 << node.bits)
                if aa.lo-aa.error < 0 or aa.hi+aa.error >= span:
                    raise AssertionError(('function address out of domain',name,node.meta,aa))
                if node.meta=='exp':
                    lo=(aa.lo-2).exp();hi=(aa.hi-2).exp()
                    maximum=(min(span,aa.hi+aa.error+spacing)-2).exp()+ALLOWANCE
                    curvature=maximum;lip=maximum
                elif node.meta=='log':
                    lo=(aa.lo+D('.25')).ln()+2;hi=(aa.hi+D('.25')).ln()+2
                    minimum=max(D('.25'),aa.lo-aa.error-spacing+D('.25'))
                    curvature=1/minimum**2;lip=1/(aa.lo-aa.error+D('.25'))
                elif node.meta.startswith('capped_'):
                    generator_allowance=KERNEL_ALLOWANCE
                    kernel_cap=int(node.meta.split('_')[1])
                    lo=capped_kernel_decimal(aa.lo-1,kernel_cap)
                    hi=capped_kernel_decimal(aa.hi-1,kernel_cap)
                    curvature=25/(2*PI*D(1).exp()).sqrt()+ALLOWANCE
                    lip=5*INV_SQRT_2PI+ALLOWANCE
                elif node.meta=='scaled_cdf':
                    lo=decimal_cdf(5*(aa.lo-1));hi=decimal_cdf(5*(aa.hi-1))
                    curvature=25/(2*PI*D(1).exp()).sqrt()+ALLOWANCE
                    lip=5*INV_SQRT_2PI+ALLOWANCE
                elif node.meta.startswith('strike_log_'):
                    h=D(int(node.meta.split('_')[-1]))/100
                    lo=1+(1+h*aa.lo).ln();hi=1+(1+h*aa.hi).ln()
                    curvature=h*h+ALLOWANCE;lip=h+ALLOWANCE
                else:
                    lo=decimal_cdf(aa.lo-8);hi=decimal_cdf(aa.hi-8)
                    curvature=1/(2*PI*D(1).exp()).sqrt()+ALLOWANCE
                    lip=INV_SQRT_2PI+ALLOWANCE
                interpolation=curvature*spacing*spacing/8
                error=lip*aa.error+interpolation+D('1.5')*unit+generator_allowance
                v=Bound(lo-generator_allowance,hi+generator_allowance,error)
                domains.append({'node':name,'kernel':node.meta,'argument_min':float(aa.lo),
                    'argument_max':float(aa.hi),'argument_error':float(aa.error),
                    'interpolation_error':float(interpolation),'output_error':float(error),
                    'generator_absolute_allowance':float(generator_allowance)})
            else:raise ValueError(node.kind)
            # All mathematical signed values and their perturbations stay in
            # the signed safe range; modular intermediate wraps are harmless.
            safe_magnitude=max(abs(v.lo),abs(v.hi))+v.error
            if node.kind=='clamp':
                safe_magnitude=D(node.meta)/D(1 << q)
            elif node.kind=='scale' and graph.nodes[node.args[0]].kind=='clamp':
                # The integer multiplier has a classical coefficient, which
                # needs no fixed-point quantum constant register. Its input
                # is physically clipped, so even a loose bias certificate
                # cannot make the implemented product overflow.
                ceiling=D(graph.nodes[node.args[0]].meta)/D(1 << q)
                safe_magnitude=ceiling*D(node.meta[0])/D(1 << node.meta[1])
            if safe_magnitude >= 32:
                raise AssertionError(('register overflow certificate',name,v))
            vals[name]=v
        bound=vals[graph.output].error*D(cap)
        return {'arithmetic_price_error_upper':math.nextafter(float(bound),math.inf),
                'domain_checks':domains,'register_absolute_limit':32,
                'quantile_coordinate_error_upper':float(unit),
                'endpoint_enclosure_allowance':float(ALLOWANCE),
                'ordered_cdf_subtractions_structurally_verified':sum(
                    n.kind=='positive_sub' for n in graph.nodes.values())}


def coupling_bound(bits):
    with localcontext() as ctx:
        ctx.prec=70
        # clip quantile function atR=5: midpoint step error<=R/N;
        # E(|Y|-R)+<=2phi(R)/(R^2+1) by a Mills lower bound.
        phi=INV_SQRT_2PI*D('-12.5').exp()
        return D(5)/D(1 << bits)+2*phi/26+ALLOWANCE


def choose_distribution(kind,cap,q,price_allowance):
    with localcontext() as ctx:
        ctx.prec=70
        lr=D('.25')*(100+D(cap))
        lz=D('.4')*(100+D(cap))
        choices=[]
        bu=q-1 if kind in ('interval','interval_kernel') else 0
        for by in range(6,19):
            outer=min(lr*coupling_bound(by),D('1.5')*D(cap)/D(1 << by))
            if kind=='direct':
                for bz in range(6,19):
                    inner=min(lz*coupling_bound(bz),D(cap)/D(1 << bz))
                    if outer+inner <= price_allowance:
                        choices.append((2**by+2**bz,by,bz,bu,outer,inner))
            else:
                inner=D(cap)/D(1 << bu) if kind in ('interval','interval_kernel') else D(0)
                if outer+inner <= price_allowance:
                    choices.append((2**by,by,0,bu,outer,inner))
        if not choices:return None
        _,by,bz,bu,outer,inner=min(choices)
        return {'b_y':by,'b_z':bz,'b_u':bu,'outer_coupling_or_BV_price_bound':math.nextafter(float(outer),math.inf),
                'inner_coupling_BV_or_strike_price_bound':math.nextafter(float(inner),math.inf),
                'distribution_price_bound':math.nextafter(float(outer+inner),math.inf)}


def gate_proxy(kind,q,e,l,c):
    # Search ranking only. Actual results are counted by emitting every block.
    counts={'direct':(3,0,0,5),'interval':(2,2,2,8),'full':(3,1,8,12),
            'interval_kernel':(2,2,2,6),'full_kernel':(2,1,2,5)}[kind]
    ne,nl,nc,nm=counts
    return (ne*(2**e)+nl*(2**l)+nc*(2**c))*28+nm*28*(q+6)**2


def tune(kind,cap,epsilon):
    """Independent finite precision search; no pilot-variance tuning."""
    price_arithmetic=D(str(epsilon))*D('.1')
    price_distribution=D(str(epsilon))*D('.4')-D('1e-12')
    records=[]
    for q in range(18,31):
        dist=choose_distribution(kind,cap,q,price_distribution)
        if dist is None:continue
        params=[4,4,4]
        if kind=='direct':params=[4,0,0]
        if kind=='full':params=[4,4,4]
        def evaluate(pp):
            graph=build_graph(kind,q,dist['b_y'],dist['b_z'],dist['b_u'],*pp,cap)
            cert=propagate(graph,cap)
            return cert['arithmetic_price_error_upper'],graph,cert
        # A domain-reduced kernel can reject a coarse prefix enclosure even
        # when a finer prefix is valid. Refine that prefix before the ordinary
        # search; do not silently discard a feasible specialised comparator.
        for _ in range(14):
            try:
                current,graph,cert=evaluate(params)
                break
            except AssertionError as exc:
                if not (kind.endswith('_kernel') and exc.args and
                        exc.args[0][0]=='function address out of domain'):
                    raise
                for i in (0,1):params[i]+=1
                if max(params[:2])>min(17,q-1):raise
        else:raise AssertionError('prefix refinement did not certify kernel domain')
        for _ in range(45):
            if current <= float(price_arithmetic):break
            candidates=[]
            for i in range(3):
                if params[i]==0 or params[i]>=min(17,q-1):continue
                pp=params.copy();pp[i]+=1
                error,gg,cc=evaluate(pp)
                improvement=current-error
                extra=gate_proxy(kind,q,*pp)-gate_proxy(kind,q,*params)
                candidates.append((improvement/max(extra,1),error,pp,gg,cc))
            if not candidates:break
            _,current,params,graph,cert=max(candidates,key=lambda x:x[0])
        if current <= float(price_arithmetic):
            records.append((gate_proxy(kind,q,*params),q,params,dist,cert))
    if not records:raise ValueError(('no feasible precision',kind,cap,epsilon))
    _,q,params,dist,cert=min(records,key=lambda x:x[0])
    graph=build_graph(kind,q,dist['b_y'],dist['b_z'],dist['b_u'],*params,cap)
    return graph, {'method':kind,'cap':cap,'price_tolerance':epsilon,
                   'discounted_error_conservative':True,'q':q,'exp_address_bits':params[0],
                   'log_address_bits':params[1],'cdf_address_bits':params[2],
                   'probability_address_bits':params[2],
                   'probability_kernel':('cap-specific conditional payoff' if kind=='full_kernel' else
                                         'domain-reduced conditional CDF' if kind=='interval_kernel' else 'normal CDF' if kind!='direct' else 'none'),
                   **dist,**cert,'bias_allowance_fraction':.5,
                   'reserved_statistical_price_allowance':epsilon*.5,
                   'bias_bound_including_correction':math.nextafter(dist['distribution_price_bound']+cert['arithmetic_price_error_upper']+1e-12,math.inf),
                   'precision_selection':'independent greedy table refinement across q=18..30, with prefix-domain refinement for specialised kernels; analytic gate proxy ranks feasible configurations; not globally optimal',
                   'feasible_configurations':len(records)}


def geometric_correction(cap):
    """Offline continuous capped geometric mean, with separate enclosure."""
    with localcontext() as ctx:
        ctx.prec=70
        variance=D('.0425');sd=variance.sqrt();mu=D('-.0225')
        h=D(cap)/100
        t0=mu/sd;t1=(mu-(1+h).ln())/sd
        value=100*((mu+variance/2).exp()*(decimal_cdf(t0+sd)-decimal_cdf(t1+sd))
                   -(decimal_cdf(t0)-decimal_cdf(t1))+h*decimal_cdf(t1))
        return {'cap':cap,'undiscounted_geometric_correction_decimal':str(value),
                'certified_error_allowance':1e-12,
                'method':'70-digit capped lognormal interval formula; normal-series absolute allowance1e-38 per call, overall retained allowance1e-12'}


def check_decimal_certificates():
    # Independent exact values at0 and symmetry, erf-based diagnostics only.
    max_difference=0
    for x in [D(j)/4 for j in range(-32,33)]:
        v=decimal_cdf(x)
        assert abs(v+decimal_cdf(-x)-1)<D('1e-38')
        diagnostic=.5*math.erfc(-float(x)/math.sqrt(2))
        max_difference=max(max_difference,abs(float(v)-diagnostic))
    tail=decimal_cdf(D('-4.75'))+ALLOWANCE
    assert tail < D(1)/D(2*(1 << 18))
    assert decimal_cdf(D('-5'))+ALLOWANCE < D(1)/D(2*(1 << 18))
    return {'cdf_grid_checks':65,'max_libm_difference_diagnostic':max_difference,
            'quantile_domain_bound':4.75,'maximum_gaussian_bits':18,
            'cdf_at_negative_4_75_upper':float(tail),
            'cdf_series_remainder_allowance':float(ALLOWANCE)}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=ROOT/'results'/'conditioning_resources.json')
    p.add_argument('--tune-only',action='store_true')
    p.add_argument('--certify-quantiles',action='store_true')
    p.add_argument('--small-checks',action='store_true')
    p.add_argument('--caps',type=int,nargs='+',default=[1,10,25])
    p.add_argument('--tolerances',type=float,nargs='+',default=[.1,.01])
    p.add_argument('--methods',nargs='+',choices=['direct','interval','full','interval_kernel','full_kernel'],
                   default=['direct','interval','full','interval_kernel','full_kernel'])
    args=p.parse_args()
    start=time.perf_counter()
    results=[];tables={};block_cache={};audits=[]
    checks=check_decimal_certificates()
    for cap in args.caps:
        for epsilon in args.tolerances:
            for kind in args.methods:
                graph,row=tune(kind,cap,epsilon)
                if not args.tune_only:
                    setup_start=time.perf_counter()
                    previous_keys=set(tables)
                    audits+=make_tables(graph,tables,args.certify_quantiles)
                    row['incremental_classical_table_setup']={
                        'new_table_keys':[list(k) for k in sorted(set(tables)-previous_keys)],
                        'measured_seconds':time.perf_counter()-setup_start,
                        'scope':'cached classical constant generation and optional quantile enclosure; order-dependent diagnostic, not a quantum execution cost'}
                    row['resources']=counted_composition(graph,tables,block_cache)
                results.append(row)
                print(f'Completed {kind} cap={cap} epsilon={epsilon}: q={graph.q}',flush=True)
    result={'status':'constructed evaluator upper bounds; not end-to-end quantum speedup',
            'model':{'assets':2,'currency':'GBP','spot':100,'strike':100,'lambda':.2,'residual_loadings':[-.1,.2],
                     'log_drifts_relative_spot':[-.015,-.03],'risk_free_rate':.03,'dividend_yield':.02,
                     'covariance':[[.05,.02],[.02,.08]],'maturity':1,'weights':[.5,.5]},
            'data':'fully stipulated public synthetic instance',
            'python_version':__import__('sys').version,
            'table_manifest':[{'kernel':k[0],'address_bits':k[1],'fraction_bits':k[2],'entries':len(v),'sha256':hashlib.sha256(json.dumps(v,separators=(',',':')).encode()).hexdigest()} for k,v in sorted(tables.items())],
            'quantum_encoding':'uniform input and q-bit comparison registers, coherent reversible Boolean arithmetic, flag predicate r<threshold, exact body inverse',
            'cost_boundary':'one complete price residual flag evaluator; includes input quantile QROM, function QROM, fixed-point arithmetic, flag and uncomputation; excludes mean-estimation reflections, controlled calls, fault tolerance and physical runtime',
            'primitive_checks':checks,'geometric_corrections':[geometric_correction(cap) for cap in args.caps],'quantile_tables':audits,'configurations':results,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'kernel_script_sha256':hashlib.sha256((ROOT/'encoding_resources.py').read_bytes()).hexdigest(),
            'gate_library_sha256':hashlib.sha256(Path(__import__('reversible').__file__).read_bytes()).hexdigest(),
            'reflection_source_sha256':hashlib.sha256((ROOT/'reflections.py').read_bytes()).hexdigest(),
            'elapsed_classical_workflow_seconds':time.perf_counter()-start}
    if args.small_checks:result['complete_toy_gate_replay']=check_small()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
