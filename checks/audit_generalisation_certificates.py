"""Independent numerical check of uniform certificate premises and budgets.

This supplements, and does not substitute for, the mathematical proofs of
linear interpolation, Mills tails, cap-kernel derivatives and total variation.
The independent recurrence uses 100 decimal digits and preserves source hashes.
"""
from decimal import Decimal,localcontext
from pathlib import Path
import hashlib
import json
import math
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'checks'))
import generalisation as gen
from audit_generalisation import endpoint_enclosure,cdf_enclosure,PHI0
from classical_reference_generalisation import exact_model_checks
D=Decimal;PAD=D('1e-38')


def recurrence(graph,cap):
    q=graph.q;u=D(1)/D(1<<q);bounds={};domains=[];minmargin=D(100)
    with localcontext() as ctx:
        ctx.prec=100
        for name,node in graph.nodes.items():
            a=[bounds[x] for x in node.args if x in bounds]
            if node.kind=='constant':
                lo=hi=graph.ideal_constants[name];err=u/2+PAD
            elif node.kind=='quantile':
                lo=D('.25');hi=D('9.75');err=u
            elif node.kind=='uniform':
                lo=1/D(2<<node.bits);hi=1-lo;err=D(0)
            elif node.kind in ('add','sub'):
                if node.kind=='add':lo=a[0][0]+a[1][0];hi=a[0][1]+a[1][1]
                else:lo=a[0][0]-a[1][1];hi=a[0][1]-a[1][0]
                err=a[0][2]+a[1][2]+PAD
            elif node.kind=='scale':
                integer,drop,original=node.meta;k=D(original)
                assert a[0][0]>=0 and k>=0
                ck=u/2 if drop else D(0)
                lo=k*a[0][0];hi=k*a[0][1]
                err=k*a[0][2]+a[0][1]*ck+a[0][2]*ck+(u if drop else D(0))+PAD
            elif node.kind=='clamp':
                ceiling=D(node.meta)*u
                lo=min(max(a[0][0],D(0)),max(D(0),ceiling-u/2));hi=min(max(a[0][1],D(0)),ceiling+u/2)
                err=a[0][2]+u/2
            elif node.kind=='function':
                kind=node.meta;alo,ahi,ae=a[0]
                if kind in ('general_exp','general_log'):sp=D(8)
                elif kind.startswith('general_strike:'):sp=D(1)
                else:sp=D(4)
                spacing=sp/D(1<<node.bits)
                assert alo-ae>=0 and ahi+ae<sp,(name,kind,alo,ahi,ae)
                minmargin=min(minmargin,alo-ae,sp-ahi-ae)
                low,lowerr=endpoint_enclosure(kind,alo);high,higherr=endpoint_enclosure(kind,ahi)
                if kind=='general_exp':
                    lip=curv=(min(sp,ahi+ae+spacing)-4).exp()+PAD;generator=PAD
                elif kind=='general_log':
                    lip=1/(alo-ae+D('.0625'));curv=1/max(D('.0625'),alo-ae-spacing+D('.0625'))**2;generator=PAD
                else:
                    parts=kind.split(':')
                    if parts[0]=='general_strike':h=D(parts[1])/100;lip=h+PAD;curv=h*h+PAD;generator=PAD
                    else:
                        lam=D(parts[1]);lip=PHI0/lam;curv=PHI0/D(1).exp().sqrt()/(lam*lam)
                        iscap=parts[0]=='general_cap'
                        if iscap and graph.kernel_bound_mode=='cap_aware':
                            h=D(parts[2])/100;lip=min(lip,(1+h)/h);curv=min(curv,(1+h)*PHI0/(h*lam))
                        lip+=PAD;curv+=PAD;generator=D('2e-12') if iscap else D('6.4e-16')
                lo=low-generator;hi=high+generator
                err=lip*ae+curv*spacing*spacing/8+D('1.5')*u+generator
                domains.append((name,kind,float(alo),float(ahi),float(ae)))
            else:raise ValueError(node.kind)
            magnitude=max(abs(lo),abs(hi))+err
            if node.kind=='clamp':magnitude=D(node.meta)*u
            elif node.kind=='scale' and graph.nodes[node.args[0]].kind=='clamp':
                magnitude=D(graph.nodes[node.args[0]].meta)*u*D(node.meta[0])/D(1<<node.meta[1])
            assert magnitude<D(64),(name,magnitude)
            bounds[name]=(lo,hi,err)
        return bounds[graph.output][2]*cap,float(minmargin),len(domains)


def distribution(model,row):
    cap=D(row['cap']);lam=D(model['lambda']);h=cap/100;d=model['outer_dimension']
    r=[[D(x) for x in rr] for rr in model['loadings']];w=[D(x) for x in model['weights']]
    by,bz,bu=row['b_y'],row['b_z'],row['b_u']
    with localcontext() as ctx:
        ctx.prec=100
        rg=[sum(w[i]*r[i][j] for i in range(len(w))) for j in range(d)]
        plip=min(PHI0/lam,(1+h)/h)+PAD
        financial_lip=min(100+cap,cap*plip)
        def coupling(b):return D(5)/D(1<<b)+2*PHI0*D('-12.5').exp()/26+PAD
        outer=sum(min(D('1.5')*cap/D(1<<by),financial_lip*(max(abs(rr[j]) for rr in r)+abs(rg[j]))*coupling(by)) for j in range(d))
        if row['method']=='direct':inner=min(cap/D(1<<bz),2*lam*(100+cap)*coupling(bz))
        elif row['method']=='interval':inner=cap/D(1<<bu)
        else:inner=D(0)
        return outer,inner


def main():
    path=ROOT/'results/generalisation_resources.json';data=json.loads(path.read_text());sources=data['source_hashes']
    before={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources};assert before==sources
    models=data['models'];records=[]
    with localcontext() as ctx:
        ctx.prec=100
        tail=PHI0*D('-32').exp()/8
        tail_amplification=(2*(D(2)+D('.4')**2/2).exp()+2+D('.25'))/D('.01')+1
        assert tail<D('6.4e-16') and tail*tail_amplification<D('2e-12')
        saturation={'normal_tail_at_8_upper':str(tail),'worst_cap_and_projection_amplification_upper':str(tail_amplification),
                    'cap_generator_tail_upper':str(tail*tail_amplification),'charged_cdf_allowance':'6.4e-16','charged_cap_allowance':'2e-12'}
        # The mathematical quantile domain premise covers every b<=18.
        phi475,e475=cdf_enclosure(D('-4.75'));phi5,e5=cdf_enclosure(D('-5'))
        assert phi475+e475<D(1)/D(2<<18)
        assert phi5+e5<D(1)/D(2<<18)
        for row in data['configurations']:
            model=next(m for m in models if m['model_id']==row['model_id'])
            exact_model_checks(model)
            g=gen.build_graph(model,row['method'],*[row[k] for k in ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],row['cap'])
            arithmetic,margin,domaincount=recurrence(g,row['cap']);outer,inner=distribution(model,row)
            # Float JSON bounds were directed upward. A tiny decimal-pad
            # mismatch is permitted only well below the 1e-12 offline reserve.
            assert arithmetic<=D(str(row['arithmetic_price_error_upper']))+D('1e-30')
            assert outer<=D(str(row['outer_distribution_bound']))+D('1e-30')
            assert inner<=D(str(row['inner_distribution_bound']))+D('1e-30')
            total=arithmetic+outer+inner+D('1e-12')
            assert total<=D('.005') and row['reserved_statistical_price_allowance']==.005
            records.append({'case_id':row['case_id'],'method':row['method'],
                'independent_arithmetic_upper_decimal':str(arithmetic),'independent_outer_bias_upper_decimal':str(outer),
                'independent_inner_bias_upper_decimal':str(inner),'independent_deterministic_bias_upper_decimal':str(total),
                'function_domain_count':domaincount,'minimum_domain_margin_diagnostic':margin,'signed_safe_limit':64})
        matched_records=[]
        matched_path=ROOT/'results/generalisation_matched_interval.json'
        if matched_path.exists():
            matched=json.loads(matched_path.read_text())
            for entry in matched['configurations']:
                interval=next(row for row in data['configurations'] if row['case_id']==entry['case_id'] and row['method']=='interval')
                model=next(m for m in models if m['model_id']==interval['model_id'])
                for method in ('full_global','full_cap'):
                    g=gen.build_graph(model,method,entry['q'],entry['b_y'],0,0,
                        entry['exp_address_bits'],entry['log_address_bits'],entry['probability_address_bits'],interval['cap'])
                    arithmetic,margin,domaincount=recurrence(g,interval['cap'])
                    assert arithmetic<=D(str(entry['results'][method]['arithmetic_price_error_upper']))+D('1e-30')
                    assert arithmetic<=D('.001')
                    matched_records.append({'case_id':entry['case_id'],'method':method,
                        'independent_arithmetic_upper_decimal':str(arithmetic),
                        'minimum_domain_margin_diagnostic':margin,'signed_safe_limit':64})
        corrections=[]
        for entry in data['geometric_corrections']:
            model=next(m for m in models if m['model_id']==entry['model_id']);cap=entry['cap']
            w=list(map(D,model['weights']));mu=sum(a*D(b) for a,b in zip(w,model['log_drifts']))
            sd=(D(model['lambda'])**2+D(model['geometric_outer_variance'])).sqrt()
            value,error=endpoint_enclosure(f'general_cap:{sd}:{cap}',mu+2)
            value*=cap;error=error*cap+D('1e-42')
            rootval=D(entry['undiscounted_value_decimal']);assert abs(value-rootval)+error<=D('1e-12')
            corrections.append({'model_id':entry['model_id'],'cap':cap,'value_decimal':str(value),
                'absolute_decimal_reference_enclosure':str(error),'retained_offline_constant_within_1e-12':True})
    after={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources};assert before==after
    out={'status':'independent uniform recurrence, distribution budgets and offline constants passed',
         'configuration_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_hashes':before,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'saturation_certificate':saturation,
         'quantile_domain':{'largest_b':18,'Phi_minus4.75_upper':str(phi475+e475),'minimum_grid_tail_probability':str(D(1)/D(2<<18))},
         'configurations':records,'matched_prefix_certificates':matched_records,
         'matched_prefix_file_sha256':hashlib.sha256(matched_path.read_bytes()).hexdigest() if matched_path.exists() else None,
         'offline_geometric_constants':corrections,
         'qualification':'independent recurrence and premise check, combined with separately stated analytical interpolation/BV proofs; no empirical or end-to-end quantum confidence claim'}
    output=ROOT/'results/generalisation_independent_certificates.json';output.write_text(json.dumps(out,indent=2)+'\n');print(output)


if __name__=='__main__':main()
