"""Independent exact search, monetary prior, QPE and complete-charge replay."""
from collections import Counter
from dataclasses import asdict
from fractions import Fraction as F
from pathlib import Path
from math import comb,isqrt
import argparse,gzip,hashlib,json,sys,time
from flint import arb,ctx

KINDS=('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1')
PI=F(355,113);POST=F(1,10**60);EPS=F(1,100);ALPHA=F(1,3600)
def archive_path(p):
    """Resolve a canonical evidence archive or a checked lossless shard set."""
    p=Path(p)
    if p.exists():return p
    if p.parent.name!='results':raise FileNotFoundError(p)
    root=p.parent.parent
    sys.path.insert(0,str(root))
    from artifact_archive import resolve_archive
    return resolve_archive(root,str(p.relative_to(root)))
def sha(p):return hashlib.sha256(archive_path(p).read_bytes()).hexdigest()
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()
def digest(v):return hashlib.sha256(canonical(v)).hexdigest()
def read(p):
    p=archive_path(p);return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def up(d):return F(d['upper_exact_fraction'])
def low(d):return F(d['lower_exact_fraction'])
def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator

def endpoint(x,upper):
    v=x.upper() if upper else x.lower();m,e=map(int,v.man_exp());return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def signature(c):
    if c['schema']=='qbp-integral-quantile-v1':return ('quantile',F(c['L']),c['b_y'],c['q'])
    return ('scalar',c['name'],F(c['lo']),F(c['span']),tuple(sorted(c['params'].items())),c['q'],F(c['output_offset']))
def exp2(x):
    x=F(x);assert x>0 and not x.numerator&(x.numerator-1) and not x.denominator&(x.denominator-1)
    return x.numerator.bit_length()-x.denominator.bit_length()
def exact_binomial(n,s):return sum((F(comb(n,j))*(1-s)**j*s**(n-j) for j in range((n+1)//2,n+1)),F(0))
def sqrt_upper(x):
    scaled=x*(1<<200);integer=(scaled.numerator+scaled.denominator-1)//scaled.denominator;r=isqrt(integer)
    if r*r<integer:r+=1
    return F(r,1<<100)


def independent_qft(bits,rotations,bounds):
    tape=[];phase=F(0);norm=F(0)
    for j in range(bits//2):k=bits-1-j;tape.extend([('CX',j,k),('CX',k,j),('CX',j,k)])
    for j in range(bits-1,-1,-1):
        for k in range(bits-1,j,-1):
            d=k-j
            if d>=32:norm+=PI/(1<<d);continue
            key=d+1;minus=rotations[f'-1:{key}'];plus=rotations[f'1:{key}'];r=minus['chronological_gates'];s=plus['chronological_gates']
            local=[(op,k) for op in r]+[(op,j) for op in r]+[('CX',k,j)]+[(op,j) for op in s]+[('CX',k,j)]
            if not r and not s:local=[]
            tape.extend(local);phase+=-F(1,1<<(d+2))+2*F(minus['phase_pi'])+F(plus['phase_pi']);norm+=2*bounds[f'-1:{key}']+bounds[f'1:{key}']
        tape.append(('H',j))
    return tape,phase,norm


def check_plan(plan,target,alpha,prior,qft,rotations,bounds):
    assert F(plan['epsilon'])==target and F(plan['alpha'])==alpha and F(plan['prior_upper'])==prior
    var=prior*(1-prior) if prior<=F(1,2) else F(1,4);assert F(plan['variance_upper'])==var
    assert plan['repetitions']>0 and plan['repetitions']%2==1 and plan['M']==1<<plan['phase_register_bits']
    M=plan['M'];root=sqrt_upper(var);err=2*PI*root/M+PI*PI/(M*M)
    assert F(plan['amplitude_error_upper'])==err and F(plan['postprocess_allowance'])==POST
    assert err+POST<=target
    if M>2:assert 2*PI*root/(M//2)+PI*PI/F((M//2)**2)+POST>target
    tape,phase,norm=independent_qft(plan['phase_register_bits'],rotations,bounds);counts=Counter(g[0] for g in tape)
    assert counts==Counter(qft['gate_counts']) and qft['T_count']==counts['T']+counts['Tdg'] and F(qft['norm_error_upper'])==norm
    assert F(qft['global_phase_correction_pi'])==phase and qft['chronological_tape_sha256']==hashlib.sha256(json.dumps(tuple(tape),separators=(',',':')).encode()).hexdigest()
    scale=1<<80;z=norm*scale;rounded=F((z.numerator+z.denominator-1)//z.denominator,scale);assert F(plan['qft_norm_error_allowance'])==rounded
    success=8/(PI*PI)-rounded;assert success>F(1,2) and F(plan['good_probability_lower'])==success
    failure=exact_binomial(plan['repetitions'],success);assert F(plan['binomial_failure_upper'])==failure<=alpha
    if plan['repetitions']>1:assert exact_binomial(plan['repetitions']-2,success)>alpha
    return counts


def charge_check(charge,profile,qft,plan):
    ap=Counter(profile['primitive_counts']);n=profile['reflection_core_qubits'];ccx=2*(n-1)
    reflection=Counter({'X':2*n,'CZ':2,'Z':1,'T':4*ccx,'Tdg':3*ccx,'H':2*ccx,'CX':6*ccx})
    r=charge['controlled_reflection'];assert Counter(r['gate_counts'])==reflection and r['T']==7*ccx
    assert r['logical_qubits']==2*n and r['explicit_T_depth_serial_upper']==4*ccx and r['conservative_5layer_T_depth_serial_upper']==5*ccx
    assert r['expanded_primitive_operations']==sum(reflection.values()) and r['global_minus_Z_charged'] and r['controlled_good_CZ_charged']
    at=ap['T']+ap['Tdg'];adj=Counter()
    for op,count in ap.items():adj[{'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}.get(op,op)]+=count
    qp=ap+adj+reflection;qmT=qp['T']+qp['Tdg'];fp=Counter(qft['gate_counts']);bits=plan['phase_register_bits'];M=plan['M'];reps=plan['repetitions']
    widths=max(profile['logical_qubits_upper'],2*n)+bits+qft.get('extra_ancillas',0)
    per_counts=ap+Counter({k:(M-1)*v for k,v in qp.items()})+fp;per_counts['H']+=bits
    per={'T':at+(M-1)*qmT+qft['T_count'],'measurements':widths,'reset_X':widths,'feedback_CZ':0,
        'expanded_primitive_operations':sum(per_counts.values())+2*widths,
        'uniform_H':profile['uniform_H']*(2*M-1)+bits,
        'explicit_T_depth_serial_upper':profile['explicit_T_depth_serial_upper']*(2*M-1)+(M-1)*4*ccx+qft['T_count'],
        'conservative_5layer_T_depth_serial_upper':profile['conservative_5layer_T_depth_serial_upper']*(2*M-1)+(M-1)*5*ccx+qft['T_count']}
    per.update({f'gate_{op}':value for op,value in per_counts.items()})
    expected=Counter(charge['per_repetition_metrics']);assert Counter(per)==expected,('Per-shot charge',per,dict(expected))
    assert Counter(charge['metrics'])==Counter({k:reps*v for k,v in per.items()})
    assert charge['peak_logical_qubits']==widths and charge['complete_oracle_T']==at and charge['controlled_Q_T']==qmT
    assert charge['controlled_Q_calls']==reps*(M-1) and charge['phase_register_measurements']==reps*bits and charge['complete_reinitialisations']==reps


def fresh_scalar(name,x,params):
    p={k:F(v) for k,v in params.items()};lam=ball(p.get('lambda',F(1)));K=ball(p.get('K',F(100)))
    def pdf(v):return (-v*v/2).exp()/(2*arb.pi()).sqrt()
    def cdf(v):return (1+(v/arb(2).sqrt()).erf())/2
    if name=='exp':return x.exp()
    if name=='log':return x.log()
    if name=='inv_square':return 1/(x*x)
    if name=='square':return x*x
    if name=='boundary_density':return K*pdf(x/lam)/(lam*ball(p['M2']))
    if name=='boundary_h':
        H=ball(p['cap']);c=1+H/K;a=c.log();return (1+K*(pdf(x/lam)-c*pdf((x-a)/lam))/(lam*ball(p['M2'])))/2
    forward=K*(x+lam*lam/2).exp()*cdf(x/lam+lam)
    if name=='vanilla_call_prime':return forward
    if name=='vanilla_call':return forward-K*cdf(x/lam)
    if name in ('cap','cap_prime'):
        H=ball(p['cap']);a=(1+H/K).log();xh=(x-a)/lam
        value=K*(x+lam*lam/2).exp()*(cdf(x/lam+lam)-cdf(xh+lam))
        if name=='cap':value-=K*(cdf(x/lam)-cdf(xh));value+=H*cdf(xh)
        return value/ball(p['cap']) if p.get('normalised',F(0)) else value
    if name in ('h','h2'):
        value=1+(-1/lam-x/(lam*lam))*pdf(x)/ball(p['D2']);return value/2 if name=='h' else value
    raise ValueError('Unknown independent scalar '+name)


def independent_scalar_bounds(c,norm_rows,hermites):
    name=c['name'];p={k:F(v) for k,v in c['params'].items()};lo=F(c['lo']);hi=lo+F(c['span']);rlo,rhi=map(F,c['exact_function_range']);first=F(c['first_derivative_upper']);second=F(c['second_derivative_upper']);ctx.prec=768
    def contains(value):
        vl,vu=endpoint(value,False),endpoint(value,True);assert rlo<=vl<=vu<=rhi
    if name=='exp':
        contains(ball(lo).exp());contains(ball(hi).exp());assert first==second==rhi
    elif name=='log':
        assert lo>0;contains(ball(lo).log());contains(ball(hi).log());assert first==1/lo and second==1/(lo*lo)
    elif name=='inv_square':assert lo>0 and (rlo,rhi)==(1/(hi*hi),1/(lo*lo)) and first==2/lo**3 and second==6/lo**4
    elif name=='square':
        assert (rlo,rhi)==(F(0) if lo<=0<=hi else min(lo*lo,hi*hi),max(lo*lo,hi*hi));assert first==2*max(abs(lo),abs(hi)) and second==2
    elif name in ('boundary_h','boundary_density','vanilla_call','vanilla_call_prime','cap','cap_prime'):
        lam=p['lambda'];K=p.get('K',F(100));assert K==100 and lam>0;vanilla=name in ('boundary_density','vanilla_call','vanilla_call_prime');cap=None if vanilla else str(p['cap']);n=norm_rows[lam,cap];M2=F(n['M2_upper']);M3=F(n['C_r_upper']['3']);M4=F(n['C_r_upper']['4'])
        if name in ('boundary_h','boundary_density'):
            denom=p['M2'];assert denom>=M2 and (rlo,rhi)==(0,1)
            if name=='boundary_density':expected=(K*hermites[1]/(lam**2*denom),K*hermites[2]/(lam**3*denom))
            else:
                H=p['cap'];expected=(min((2*K+H)*hermites[1]/lam**2,M3+2*M2)/(2*denom),min((2*K+H)*hermites[2]/lam**3,M4+5*M3+4*M2)/(2*denom))
            assert (first,second)==expected
        elif vanilla:
            assert rlo==0 and rhi>=endpoint(ball(K)*ball(hi+lam*lam/2).exp(),True)
            assert (first,second)==((rhi,rhi+M2) if name=='vanilla_call' else (rhi+M2,rhi+3*M2+M3))
        else:
            H=p['cap'];M1=F(n['M1_upper']);scale=1/H if p.get('normalised',F(0)) else F(1);oldfirst=min(K+H,H*hermites[0]/lam);oldsecond=H*hermites[1]/lam**2
            if name=='cap':expected=(min(oldfirst,M1)*scale,min(oldsecond,M1+M2)*scale);assert (rlo,rhi)==(0,H*scale)
            else:expected=(min(oldsecond,M1+M2)*scale,min(H*hermites[2]/lam**3,M1+3*M2+M3)*scale);assert (rlo,rhi)==(0,M1*scale)
            assert (first,second)==expected
    elif name in ('h','h2'):
        lam=p['lambda'];D=p['D2'];factor=F(1,2) if name=='h' else F(1);assert (rlo,rhi)==(0,2*factor);assert (first,second)==(factor*(hermites[1]/lam+hermites[2]/lam**2)/D,factor*(hermites[2]/lam+hermites[3]/lam**2)/D)
    else:raise ValueError(name)

def exact_optimisation_replay(graph,tables,oracle):
    """Independent masked constant induction and backwards liveness replay."""
    import copy
    unit=1<<graph.q;mask=(1<<graph.width)-1
    def rnd(x):
        t=F(x)+F(1,2);return t.numerator//t.denominator
    def scale(x):
        x=F(x);return (x.numerator,x.denominator.bit_length()-1) if not x.denominator&(x.denominator-1) else (rnd(x*unit),graph.q)
    def description(g):
        return {'q':g.q,'width':g.width,'inputs':g.inputs,'output':g.output,'configuration':g.configuration,'nodes':[{'name':name,'kind':n.kind,'args':list(n.args),'meta':n.meta.as_dict() if isinstance(n.meta,oracle.ScalarSpec) else str(n.meta),'bits':n.bits,'output_bits':n.output_bits} for name,n in g.nodes.items()]}
    def const(n,args,name):
        k=n.kind
        if k=='constant':v=rnd(n.meta*unit)
        elif k=='add':v=sum(args)
        elif k in ('sub','signed_sub'):v=args[0]-args[1]
        elif k=='positive_sub':v=max(0,args[0]-args[1])
        elif k=='mul':v=args[0]*args[1]>>graph.q
        elif k=='scale':a,b=scale(n.meta);v=args[0]*a>>b
        elif k=='uniform':v=(2*args[0]+1)*(unit>>(n.meta+1))
        elif k=='minimum':v=min(args)
        elif k=='quantile':v=tables[oracle.table_key(name,n)][args[0]]
        elif k=='function':
            f=graph.q+exp2(n.meta.span)-n.bits;index=args[0]>>f;rem=args[0]&((1<<f)-1);word=tables[oracle.table_key(name,n)][index];base=word&mask;mag=(word>>graph.width)&mask;drop=mag*rem>>f;v=base-drop if word>>(2*graph.width)&1 else base+drop
        elif k in ('clamp','clamp_signed'):
            v=args[0]
            if k=='clamp_signed' and v>=1<<(graph.width-1):v-=1<<graph.width
            v=max(0,min(v,rnd(n.meta*unit)))
        else:raise ValueError(k)
        return v&((1<<n.output_bits)-1)
    before=description(graph);g=copy.deepcopy(graph);known={};rewrites=[];narrow=[];preceding=set(graph.inputs)
    for name,n in graph.nodes.items():
        assert set(n.args)<=preceding;preceding.add(name);reason=None
        if n.kind=='constant' or all(a in known for a in n.args):v=const(n,[known[a] for a in n.args],name);reason='all_arguments_are_exact_masked_constants'
        elif n.kind=='scale' and scale(n.meta)[0]==0:v=0;reason='exact_implemented_scale_integer_is_zero'
        elif n.kind=='mul' and any(known.get(a)==0 for a in n.args):v=0;reason='exact_masked_zero_multiplicand'
        elif n.kind in ('sub','signed_sub','positive_sub') and n.args[0]==n.args[1]:v=0;reason='identical_integer_arguments_subtract_to_zero'
        elif n.kind in ('clamp','clamp_signed') and rnd(n.meta*unit)==0:v=0;reason='clamp_cap_integer_is_zero'
        else:continue
        known[name]=v
        if n.output_bits!=graph.width:narrow.append(name);continue
        if n.kind=='constant':continue
        g.nodes[name]=oracle.IntegralNode('constant',(),meta=F(v,unit),bits=0,output_bits=n.output_bits)
        rewrites.append({'node':name,'original_kind':n.kind,'original_args':list(n.args),'original_output_bits':n.output_bits,'integer_value':v,'exact_dyadic_constant':str(F(v,unit)),'reason':reason})
    live=set();stack=[g.output]
    while stack:
        name=stack.pop()
        if name in g.inputs or name in live:continue
        assert name in g.nodes;live.add(name);stack.extend(g.nodes[name].args)
    removed=[name for name in graph.nodes if name not in live];g.nodes={name:n for name,n in g.nodes.items() if name in live}
    required={oracle.table_key(name,n) for name,n in g.nodes.items() if n.kind in ('quantile','function')}
    expected={'source_graph_sha256':digest(before),'source_tables_sha256':digest({str(k):list(v) for k,v in sorted(tables.items(),key=lambda x:str(x[0]))}),'rewritten_graph_sha256':digest(description(g)),'source_node_count':len(graph.nodes),'retained_node_count':len(g.nodes),'constant_rewrites':rewrites,'dead_nodes_removed':removed,'known_narrow_nodes_left_unchanged':narrow,'original_input_registers_preserved':True,'word_width_preserved':True,'fractional_width_preserved':True,'output_node_preserved':True,'reflection_core_qubits':sum(g.inputs.values())+g.q+1,'source_table_count':len(tables),'retained_table_count':len(required)}
    assert description(graph)==before
    return g,expected

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));p.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/higherdim_precision_independent_audit.json'));args=p.parse_args();root=Path(args.root);start=time.perf_counter();sys.path.insert(0,str(root))
    import boundary_value_oracle as oracle
    import boundary_value_certificate as certifier
    print('Loading frozen comparison and scalar evidence',flush=True)
    data=read(root/'results/higherdim_precision_comparison.json.gz');assert data['source_sha256']==sha(root/'higherdim_precision_selection.py')
    for f,h in data['dependency_source_hashes'].items():assert sha(root/f)==h
    for f,h in data['input_sha256'].items():assert sha(root/'results'/f)==h
    assert F(data['monetary_radius'])==EPS and F(data['alpha_per_selected_price'])==ALPHA and F(data['union_over_original36_failure_upper'])==36*ALPHA
    path=archive_path(root/'results'/data['scalar_ordinate_archive']);assert sha(path)==data['scalar_ordinate_archive_sha256'];archive={};nodes=0;headers=set();archive_rows=[];anchors={};fresh_ordinates=0
    with gzip.open(path,'rt') as z:
        for line in z:
            row=json.loads(line);c=row['certificate'];assert row['sha256']==digest(c);header=digest({k:v for k,v in c.items() if k!='nodes'});assert header not in headers;headers.add(header);archive_rows.append(c)
            sig=signature(c)
            if c['schema']=='qbp-integral-quantile-v1' or c['table_bits']==12:
                assert sig not in anchors;anchors[sig]=c
    import boundary_grid_certificate as bg
    scalar_norm_rows={(F(r['lambda']),r['cap']):r for r in read(root/'results/boundary_derivative_norms.json')['records']}
    scalar_hermites={r['unit_derivative_order']:F(r['norm_upper']) for r in read(root/'results/gaussian_derivative_norms.json')['unit_norms'] if r['passed']}
    print('Verifying fresh 768-bit scalar values and exact subsets',flush=True)
    for c in archive_rows:
        sig=signature(c);anchor=anchors[sig]
        if c['schema']=='qbp-integral-quantile-v1':
            assert c==anchor and F(c['L'])==6 and c['b_y']==10 and c['q']==28;key=sig
            ctx.prec=768;L=ball(6);cdf=lambda x:(1+(x/arb(2).sqrt()).erf())/2;c0=cdf(-L);mass=cdf(L)-c0
            for j,n in enumerate(c['nodes']):
                value=arb(2).sqrt()*(2*(c0+mass*ball(F(2*j+1,2048)))-1).erfinv();vl,vu=endpoint(value,False),endpoint(value,True)
                assert F(n['quantile_lower'])<=vl<=vu<=F(n['quantile_upper'])
                actual=F(n['offset_integer'],1<<28)-6;assert actual==F(n['actual_coordinate'])
                assert max(abs(actual-vl),abs(actual-vu))<=F(n['node_error_upper']);fresh_ordinates+=1
        else:
            assert 6<=c['table_bits']<=12 and c['q']==28 and c['table_bits']<=anchor['table_bits']
            assert c['nodes']==anchor['nodes'][::1<<(12-c['table_bits'])]
            for field in ('source_sha256','norm_source_sha256','exact_function_range','first_derivative_upper','second_derivative_upper','precision_bits'):assert c[field]==anchor[field]
            assert c['source_sha256']==sha(root/'boundary_value_certificate.py') and c['norm_source_sha256']==sha(root/'results/boundary_derivative_norms.json')
            if c is anchor:
                independent_scalar_bounds(c,scalar_norm_rows,scalar_hermites)
                ctx.prec=768
                for n in c['nodes']:
                    value=fresh_scalar(c['name'],ball(F(n['input'])),c['params'])-ball(F(c['output_offset']));vl,vu=endpoint(value,False),endpoint(value,True)
                    assert F(n['lower'])<=vl<=vu<=F(n['upper']),('Scalar ordinate',c['name'],n['input']);fresh_ordinates+=1
            unit=1<<28;largest=max(max(abs(F(n['integer'],unit)-F(n['lower'])),abs(F(n['integer'],unit)-F(n['upper']))) for n in c['nodes'])
            rem=F(c['second_derivative_upper'])*(F(c['span'])/(1<<c['table_bits']))**2/8
            assert largest==F(c['nodal_error_upper']) and rem==F(c['linear_interpolation_error_upper']) and F(c['output_error_upper'])==largest+rem+F(1,unit)
            assert c['fractional_remainder_bits']==28+exp2(c['span'])-c['table_bits'];key=(*sig,c['table_bits'])
        assert key not in archive;archive[key]=c;nodes+=len(c['nodes'])
    assert len(headers)==data['unique_scalar_certificates']
    print(f'Scalar audit passed: {len(archive)} certificates, {fresh_ordinates} fresh ordinates',flush=True)
    models={m['model_id']:m for m in read(root/'results/generalisation_protocol.json')['models']}
    preparation={r['model_id']:r for r in read(root/'results/higherdim_preparation.json')['records']};assert len(preparation)==6
    domains={mid:r['domain'] for mid,r in preparation.items()};moments={mid:r['moments'] for mid,r in preparation.items()}
    coef={(mid,r['payoff'],r['cap']):r['coefficients'] for mid,prep in preparation.items() for r in prep['cases']};assert len(coef)==24
    trans={};mass_iv=read(root/'results/higherdim_preparation.json')['shared_MGF_profile']['node_certificate']['mass_interval'];ctx.prec=512
    masses={p:bg.serial(bg.interval(low(mass_iv),up(mass_iv))**p) for p in (3,7)}
    for filename in ('tensor_transfer.json','tensor_first_transfer.json'):
        for r in read(root/'results'/filename)['records']:
            if r['model_id'] not in preparation:continue
            cap=None if r['cap'] is None else int(F(r['cap']));tr=dict(r);tr['mass_interval']=masses[models[r['model_id']]['outer_dimension']];trans[r['model_id'],r['payoff'],cap,r['order']]=tr
    assert len(trans)==48
    norms={(F(r['lambda']),r['cap']):r for r in read(root/'results/boundary_derivative_norms.json')['records']}
    original=anchors['quantile',F(6),10,28];gapcache={};literal_mgfs={};packedcache={};evalcache={};programs=0;profiles=0;flagplans=0;zero=0;rejected=0;optimised_profiles=0;original_profiles=0
    rotationdata=read(root/'results/qft_rotation_tapes.json')['rotations'];rotationbounds={r['key']:F(r['operator_norm_upper_rational']) for r in read(root/'results/qft_rotation_certificates.json')['certificates']}

    def gap(mid):
        if mid not in gapcache:
            ctx.prec=768;m=models[mid];w=list(map(F,m['weights']));mu=list(map(F,m['log_drifts']));load=[list(map(F,a)) for a in m['loadings']];p=m['outer_dimension'];mg=sum(a*b for a,b in zip(w,mu));rg=[sum(wi*row[j] for wi,row in zip(w,load)) for j in range(p)];lam=F(m['lambda']);d=[v-mg for v in mu];s=[[v-r for v,r in zip(row,rg)] for row in load]
            def empirical(v):
                if v not in literal_mgfs:literal_mgfs[v]=sum((ball(v*F(n['actual_coordinate'])).exp() for n in original['nodes']),arb(0))/1024
                return literal_mgfs[v]
            def product(rates):
                value=arb(1)
                for rate in rates:value*=empirical(rate)
                return value
            total=sum((ball(wi)*ball(di).exp()*product([x+y for x,y in zip(rg,si)]) for wi,di,si in zip(w,d,s)),arb(0))-product(rg)
            g=100*ball(mg+lam*lam/2).exp()*total;gapcache[mid]=endpoint(g,False),endpoint(g,True)
        return gapcache[mid]

    def make(mid,payoff,cap,kind,bits,encoding='centred'):
        d=domains[mid];g=oracle.build_graph(models[mid],cap,q=28,b_y=10,beta_bits=8,table_bits=12,outer_bound=6,U_star=F(d['U_star']),kind=kind,retained_bound=max(abs(F(x)) for x in d['coordinate_bounds']),output_encoding=encoding)
        funcs={n:v.bits for n,v in g.nodes.items() if v.kind=='function'};assert set(funcs)==set(bits)
        for name,b in bits.items():assert type(b) is int and 6<=b<=12;g.nodes[name].bits=b
        g.configuration['function_address_bits']=dict(bits);assert g.configuration['payoff']==payoff
        expected_inputs={f'iy{j}':10 for j in range(models[mid]['outer_dimension'])}
        if kind=='boundary_r2':expected_inputs.update({'iz0':8,'iz1':8})
        assert g.inputs==expected_inputs
        assert g.configuration['strike_bits']==0 and g.configuration['q']==28
        tables={};cc={}
        for name,n in g.nodes.items():
            if n.kind=='quantile':key=('quantile',F(n.meta),n.bits,g.q);c=archive[key];tk=('quantile',n.bits,str(n.meta))
            elif n.kind=='function':
                s=n.meta;key=('scalar',s.name,s.lo,s.span,tuple(sorted((k,str(v)) for k,v in s.parameters)),g.q,F(0),n.bits);c=archive[key];tk=('scalar',json.dumps(s.as_dict(),sort_keys=True),n.bits)
            else:continue
            pk=(key,g.width)
            if pk not in packedcache:
                if n.kind=='quantile':words=tuple(r['offset_integer'] for r in c['nodes'])
                else:
                    values=[r['integer'] for r in c['nodes']];words=[]
                    for j,v in enumerate(values):
                        delta=values[j+1]-v if j+1<len(values) else 0;assert min(v,abs(delta))>=0 and max(v,abs(delta))<1<<g.width
                        words.append(v|(abs(delta)<<g.width)|(int(delta<0)<<(2*g.width)))
                    words=tuple(words)
                packedcache[pk]=words
            tables[tk]=packedcache[pk];cc[name]=c
        return g,tables,cc

    def money(key,c):
        mid,payoff,cap,kind=key;cfg=c['configuration'];B=F(cfg['monetary_signal_normaliser']);decode=F(cfg['monetary_probability_decode_factor']);eta=F(c['probability_error_upper'])
        tr=trans[mid,payoff,cap,1 if kind=='positive_r1' else 2];head=coef[mid,payoff,cap]['c0' if kind=='positive_r1' else 'coefficient'];mass=tr['mass_interval'];ml,mh=low(mass),up(mass)
        mu={r:up(moments[mid]['centred_U_moments'][str(r)]) for r in (2,3,4,6)};mu[1]=up(moments[mid]['centred_U_moments']['1']);nr=norms[F(models[mid]['lambda']),None if cap is None else str(cap)]
        M2=F(nr['M2_upper']);M3=F(nr['C_r_upper']['3']);M4=F(nr['C_r_upper']['4']);M6=F(nr['C_r_upper']['6']);aux=F(0)
        if kind=='boundary_r2':aux=(M3*mu[3]/12+M4*mu[4]/48)/(1<<16)
        elif kind=='stable_gauss_r2':
            constants=cfg['gauss_constants'];zn=sum(F(constants['w_'+s]['ideal_upper'])*F(constants['z_'+s]['error_upper']) for s in ('minus','plus'));wn=sum(F(constants['w_'+s]['error_upper']) for s in ('minus','plus'));aux=M6*mu[6]/14400+M3*mu[3]*zn/2+M2*mu[2]*wn/2
        total=(up(head)-low(head))/2+F(tr['financial_transfer_upper'])+mh*(aux+decode*B*eta)+(mh-ml)*B/2
        return total,aux,mu

    def evaluate(key,bits,full=False):
        nonlocal programs
        ek=(*key,tuple(sorted(bits.items())))
        if ek in evalcache and not full:return evalcache[ek]
        g,t,cc=make(*key,bits);c=certifier.certify_graph(g,t,cc);h=digest(c);det=money(key,c)[0] if c['status']=='passed' else None
        short=(h,c['status'],c['failures'],det);evalcache[ek]=short;programs+=1
        return (g,t,c,det) if full else short

    def profile_check(g,t,p):
        nonlocal profiles
        profiles+=1;assert p['configuration']==g.configuration and p['source_sha256']==sha(root/'boundary_value_oracle.py') and p['retained_emitter_sha256']==sha(root/'integral_value_oracle.py')
        assert p['comparison_bits']==28 and p['integer_fractional_bits']==28 and p['integer_word_width']==g.width and p['maximum_probability_integer']==1<<28
        assert p['reflection_core_qubits']==sum(g.inputs.values())+29 and p['uniform_H']==sum(g.inputs.values())+28
        assert [r['node'] for r in p['blocks']]==list(g.nodes);forward=Counter();depth=F(0);five=F(0);lookups=0;retained=sum(g.inputs.values())+30;peak=retained
        for r,(name,n) in zip(p['blocks'],g.nodes.items()):
            distinct=tuple(dict.fromkeys(n.args));alias=[distinct.index(a) for a in n.args];widths=[g.inputs[a] if a in g.inputs else g.nodes[a].output_bits for a in n.args]
            assert r['kind']==n.kind and r['argument_widths']==widths and r['argument_alias_pattern']==alias and r['output_bits']==n.output_bits
            assert r['function']==(n.meta.as_dict() if n.kind=='function' else None)
            table=t.get(oracle.table_key(name,n)) if n.kind in ('function','quantile') else None;th=hashlib.sha256(json.dumps(table,separators=(',',':')).encode()).hexdigest() if table is not None else None;assert r['table_sha256']==th
            counts=Counter(r['primitive_counts']);assert not set(counts)-{'X','H','CX','CZ','Z','S','Sdg','T','Tdg'} and all(type(v) is int and v>=0 for v in counts.values())
            assert r['logical_T_count']==counts['T']+counts['Tdg'];forward.update(counts);depth+=r['explicit_T_depth_serial_upper'];five+=r['conservative_5layer_T_depth_serial_upper'];lookups+=r['lookup_invocations'];peak=max(peak,retained+r['extra_qubits_upper']);retained+=n.output_bits
        assert forward==Counter(p['primitive_counts_forward_body']);complete=Counter({k:2*v for k,v in forward.items()});complete.update(p['comparator_primitive_counts']);complete['H']+=p['uniform_H'];assert complete==Counter(p['primitive_counts']) and p['logical_T_count']==complete['T']+complete['Tdg']
        comp=oracle.old.vo.ValueCircuit();a,b=comp.alloc(29),comp.alloc(29);flag=comp.alloc(1)[0];comp.compare_into(a,b,flag);cm=comp.resources()
        assert Counter(cm['clifford_T_expanded_gate_counts'])==Counter(p['comparator_primitive_counts'])
        assert p['explicit_T_depth_serial_upper']==2*depth+cm['explicit_7T_sequence_T_depth'] and p['conservative_5layer_T_depth_serial_upper']==2*five+cm['conservative_5layer_CCX_T_depth_upper']
        assert p['lookup_invocations_forward_and_cleanup']==2*lookups;peak=max(peak,retained+cm['wires']-59);assert p['logical_qubits_upper']==peak

    print('Replaying candidate profiles, integer simplification and complete confidence accounting',flush=True)
    records=data['records'];groups={};bases={};expected={(*key,kind) for key in coef for kind in KINDS}
    for ix,r in enumerate(records):
        if ix and ix%50==0:print(f'Checked {ix}/{len(records)} candidate profiles',flush=True)
        key=r['model_id'],r['payoff'],r['cap'],r['kind'];assert key in expected;groups.setdefault(key,[]).append((ix,r));bits=r['node_address_bits'];g,t,c,det=evaluate(key,bits,True)
        saved=dict(r['base_graph_certificate']);h=saved.pop('precision_search_certificate_sha256');assert h==digest(saved) and c==saved
        if 'unchanged_all12_baseline' in r['allocation_limits']:assert all(b==12 for b in bits.values()) and key not in bases;bases[key]=r
        for allocation in r['allocation_limits']:
            if allocation!='unchanged_all12_baseline':assert F(allocation) in (F(1,1000),F(1,400),F(1,200)) and det<=F(allocation)
        if det is None:assert r['complete_deterministic_budget'] is None and not r['profiles'] and not r['strategies'];continue
        bud=r['complete_deterministic_budget'];total,aux,mu=money(key,c);assert total==F(bud['total_deterministic_undiscounted_price_error_upper']) and aux==F(bud['auxiliary_components']['bias_upper'])
        B=F(c['configuration']['monetary_signal_normaliser']);eta=F(c['probability_error_upper']);U=F(c['configuration']['U_star']);positive=c['configuration']['output_encoding']=='positive';pr=r['correction_prior']
        assert positive==pr['native_positive'] and F(bud['monetary_probability_decode_factor'])==(1 if positive else 2)
        if key[3]=='positive_r1':
            gl,gh=gap(key[0]);gapupper=F(pr['underlying_gap_upper']);assert gapupper>=gh and gapupper-gl<F(1,10**90)
            means=[gapupper]
            if key[2] is not None:means.extend([F(key[2]),F(c['configuration']['M1'])*mu[1]])
            mean=min(means);pupper=min(F(1),mean/B+eta);assert F(pr['ideal_mean_upper'])==mean
        else:pupper=min(F(1),mu[2]/(U*U)+(1 if positive else 2)*eta)
        assert F(pr['probability_upper'])==pupper
        assert set(r['profiles'])=={'original','exact_integer_simplification'}
        for implementation,implementations in r['profiles'].items():
            assert set(implementations)=={'exact7t','unitary8t'}
            for architecture,encoded in implementations.items():
                assert set(encoded)==({'centred'} if positive else {'centred','positive_part','negative_part'})
                for encoding,p in encoded.items():
                    pg,pt,_=make(*key,bits,encoding);assert p['architecture']==architecture and p['implementation']==implementation
                    assert p['integer_optimisation_source_sha256']==sha(root/'boundary_integer_optimisation.py')
                    if implementation=='exact_integer_simplification':
                        pg,expected_optimisation=exact_optimisation_replay(pg,pt,oracle)
                        for name,value in expected_optimisation.items():assert p['integer_optimisation'][name]==value
                        assert p['integer_optimisation']['source_sha256']==sha(root/'boundary_integer_optimisation.py');optimised_profiles+=1
                    else:assert p['integer_optimisation']['applied'] is False;original_profiles+=1
                    profile_check(pg,pt,p)
        for strategy in r['strategies']:
            result=strategy['result'];remaining=EPS-det;C=F(bud['external_factor_upper'])*B;radius=C*pupper/(2 if positive else 1)
            if remaining<=0:assert result['kind']=='rejected';rejected+=1;continue
            if radius<=remaining:
                assert result['kind']=='deterministic_correction_prior' and result['T_count']==0 and result['quantum_shots']==0 and F(result['total_price_error_upper'])==det+radius;zero+=1;continue
            assert result['kind']=='complete_quantum_resource_plan';encoding=strategy['strategy'];labels=['centred'] if positive or encoding=='centred' else ['positive_part','negative_part'];target=remaining/C if positive else remaining/(2*C);priors=[pupper]*len(labels) if positive or encoding=='split' else [F(1)];alphas=[ALPHA/len(labels)]*len(labels);decode=1 if positive or encoding=='split' else 2
            assert [f['flag'] for f in result['flags']]==labels;stat=F(0);failure=F(0);T=0;shots=0;peak=0
            for label,u,alpha,flag in zip(labels,priors,alphas,result['flags']):
                profile=r['profiles'][strategy['implementation']][strategy['architecture']][label];assert flag['financial_profile_sha256']==digest(profile) and F(flag['probability_tolerance'])==target and F(flag['prior_upper'])==u
                check_plan(flag['plan'],target,alpha,u,flag['fourier_profile'],rotationdata,rotationbounds);charge_check(flag['complete_charge'],profile,flag['fourier_profile'],flag['plan']);flagplans+=1
                stat+=F(flag['plan']['amplitude_error_upper'])+F(flag['plan']['postprocess_allowance']);failure+=F(flag['plan']['binomial_failure_upper']);T+=flag['complete_charge']['metrics']['T'];shots+=flag['plan']['repetitions'];peak=max(peak,flag['complete_charge']['peak_logical_qubits'])
            stat*=C*decode;assert F(result['statistical_price_error_upper'])==stat and F(result['total_price_error_upper'])==det+stat<=EPS and F(result['failure_probability_upper'])==failure<=ALPHA
            assert result['T_count']==T and result['quantum_shots']==shots and result['logical_qubits_upper']==peak
        valid=[s for s in r['strategies'] if s['result']['kind']!='rejected'];assert r['selected_strategy']==(min(valid,key=lambda s:(s['result']['T_count'],s['implementation'],s['architecture'],s['strategy'])) if valid else None)
    assert set(groups)==set(bases)==expected
    print('Replaying every precision-search proposal and final full-cost selection',flush=True)
    trials={};proposal_count=0;chosen_count=0;skipped=0
    for row in data['search_trials']:
        key=row['model_id'],row['payoff'],row['cap'],row['kind'];trials.setdefault((*key,F(row['allocation_limit'])),[]).append(row)
    for key in expected:
        initial=dict(bases[key]['node_address_bits']);initial_det=F(bases[key]['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper']) if bases[key]['complete_deterministic_budget'] else None
        for limit in (F(1,1000),F(1,400),F(1,200)):
            rows=trials.get((*key,limit),[])
            if initial_det is None or initial_det>limit:
                assert len(rows)==1 and rows[0]['status']=='starting_program_outside_this_allocation';skipped+=1;assert all(str(limit) not in r['allocation_limits'] for _,r in groups[key]);continue
            bits=dict(initial);previous=initial_det;offset=0;step=0
            while True:
                names=[n for n,b in bits.items() if b>6];batch=rows[offset:offset+len(names)];assert len(batch)==len(names);choices=[]
                for n,row in zip(names,batch):
                    proposed=dict(bits);proposed[n]-=1;h,status,failures,det=evaluate(key,proposed);accepted=status=='passed' and det is not None and det<=limit
                    assert row['node']==n and row['from_bits']==bits[n] and row['to_bits']==proposed[n] and row['step']==step and row['current_width_map']==bits and row['proposed_width_map']==proposed
                    assert row['graph_certificate_sha256']==h and row['graph_failures']==failures and row['accepted_by_certificate']==accepted
                    assert (F(row['deterministic_error_upper']) if row['deterministic_error_upper'] is not None else None)==det
                    if accepted:
                        score=F(1<<(bits[n]-1))/max(F(1,1<<28),det-previous);assert F(row['greedy_score'])==score;choices.append((score,n,proposed,det))
                    else:assert row['greedy_score'] is None and not row['chosen']
                    proposal_count+=1
                offset+=len(names)
                if not choices:break
                choices.sort(key=lambda x:(-x[0],x[1]));_,winner,bits,previous=choices[0];assert [r['node'] for r in batch if r['chosen']]==[winner];step+=1;chosen_count+=1
            assert offset==len(rows);final=[r for _,r in groups[key] if str(limit) in r['allocation_limits']];assert len(final)==1 and final[0]['node_address_bits']==bits
    assert len(data['selected_prices'])==len(coef)==24;seen_cases=set()
    for row in data['selected_prices']:
        key=row['model_id'],row['payoff'],row['cap'];assert key in coef and key not in seen_cases;seen_cases.add(key)
        for kind in KINDS:
            valid=[(ix,r) for ix,r in groups[(*key,kind)] if r['selected_strategy']];winner=min(valid,key=lambda x:(x[1]['selected_strategy']['result']['T_count'],x[0])) if valid else None
            sel=row['methods'][kind];assert (sel is None)==(winner is None)
            if winner:assert sel['candidate_index']==winner[0] and sel['selected']==winner[1]['selected_strategy']
        available=[k for k in KINDS if row['methods'][k]];winner=min(available,key=lambda k:row['methods'][k]['selected']['result']['T_count']) if available else None;assert row['quantum_resource_winner']==winner
    result={'status':'passed','audit_source_sha256':sha(Path(__file__)),'comparison_sha256':sha(root/'results/higherdim_precision_comparison.json.gz'),'scalar_archive_sha256':sha(path),
        'retained_candidates':len(records),'method_groups':len(groups),'profiles_checked':profiles,'candidate_profiles_freshly_emitted_by_this_audit':False,'driver_unique_counted_blocks':data['unique_counted_blocks'],'scalar_certificates_checked':len(archive),'scalar_ordinates_linked_to_checked_base':nodes,'fresh_768bit_ordinates_checked':fresh_ordinates,'exact_optimisation_profiles_replayed':optimised_profiles,'original_profiles_checked':original_profiles,
        'unique_program_certificates_reconstructed':len(evalcache),'greedy_proposals_checked':proposal_count,'starting_infeasible_allocations_checked':skipped,'chosen_greedy_steps_checked':chosen_count,
        'complete_QPE_flag_plans_checked':flagplans,'deterministic_prior_strategy_records_checked':zero,'rejected_strategy_records_checked':rejected,'selected_price_cases_checked':24,
        'scope':'Fresh768bit base scalar and quantile values; exact retained scalar subsets, preserved input law, full graph guards, declared greedy tie/choice replay, deterministic monetary and auxiliary envelopes, correlated native positive priors, independent binary QPE ceilings/binomial tails, complete controlled reflection/QFT/observation/reset accounting and full-cost selection. Norm and literal primitive-body audits are retained separately.',
        'elapsed_seconds_diagnostic':time.perf_counter()-start}
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
