"""Independent positive transport and block-Cauchy audit.

The producer is imported only after the caller has frozen its source. Literal
Gaussian cell sums are evaluated at 768 bits. Positive interval quadrature is
an additional check, not a replacement for the exact transform enclosure.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib, json, random, sys
sys.dont_write_bytecode = True

from flint import arb, ctx


def root_module(root):
    sys.path.insert(0, str(root))
    return importlib.import_module('boundary_grid_certificate')


def literal_rho0(cell, rate, bg):
    a, b, c = cell['_left'], cell['_right'], bg.ball(F(cell['actual_node']))
    if not rate:
        return c*(2*bg.Phi(c)-bg.Phi(a)-bg.Phi(b))+2*bg.phi(c)-bg.phi(a)-bg.phi(b)
    v = bg.ball(rate)
    left = (v*v/2).exp()*(bg.Phi(c-v)-bg.Phi(a-v))
    right = (v*v/2).exp()*(bg.Phi(b-v)-bg.Phi(c-v))
    return ((v*c).exp()*(2*bg.Phi(c)-bg.Phi(a)-bg.Phi(b))+right-left)/v


def literal_rho2(cell, rate, bg):
    if not rate:
        return cell['_second']/2
    a,b,c = cell['_left'],cell['_right'],bg.ball(F(cell['actual_node']))
    v=bg.ball(rate)
    cont=(v*v/2).exp()*(bg.Phi(b-v)-bg.Phi(a-v))
    return (cont-bg.interval(F(cell['mass']['lower_exact_fraction']),F(cell['mass']['upper_exact_fraction']))*(v*c).exp()
            -v*cell['_first']*(v*c).exp())/(v*v)


def positive_quadrature(cell, rate, kind, bg, pieces=128):
    """Validated rectangles after mapping each segment to [0,1]."""
    a,b,c=cell['_left'],cell['_right'],bg.ball(F(cell['actual_node']))
    result=arb(0)
    for left,right,endpoint,side in ((a,c,a,'left'),(c,b,b,'right')):
        length=right-left
        for j in range(pieces):
            t=bg.interval(F(j,pieces),F(j+1,pieces))
            z=left+length*t
            if kind==0:
                weight=bg.Phi(z)-bg.Phi(endpoint) if side=='left' else bg.Phi(endpoint)-bg.Phi(z)
            else:
                weight=bg.phi(z)+z*bg.Phi(z)-bg.phi(endpoint)-z*bg.Phi(endpoint)
            lo,hi=bg.endpoints(weight)
            if hi<0:raise AssertionError('A positive transport density has a negative upper endpoint')
            weight=bg.interval(max(F(0),lo),hi)
            result+=length/pieces*weight*(bg.ball(rate)*z).exp()
    return result


def cauchy_audit(bg):
    rng=random.Random(28061024)
    checks=0;strict=0
    for case in range(120):
        weights=[F(rng.randrange(1,31),17) for _ in range(16)]
        us=[F(rng.randrange(0,31),11) for _ in range(16)]
        js=[F(rng.randrange(-30,31),13) for _ in range(16)]
        if case==0:us=[F(0)]*16
        if case==1:js=[F(0)]*16
        values=[]
        for B in (1,2,4,8,16):
            first=arb(0);second=arb(0);aa=arb(0);bb=arb(0)
            for b in range(B):
                inds=range(b*16//B,(b+1)*16//B)
                def mom(fn):return sum((w*fn(u,j) for w,u,j in ((weights[i],us[i],js[i]) for i in inds)),F(0))
                U2,J2,U4,U2J2=[mom(fn) for fn in (lambda u,j:u*u,lambda u,j:j*j,lambda u,j:u**4,lambda u,j:u*u*j*j)]
                a=(bg.ball(U2)*bg.ball(U2J2)).sqrt()
                bnd=(bg.ball(U4)*bg.ball(J2)).sqrt()
                first+=(bg.ball(U2)*bg.ball(J2)).sqrt()
                aa+=a;bb+=bnd
                second+=bg.ball(min(bg.upper(a),bg.upper(bnd)))
            values.append((B,first,second,aa,bb))
        for old,new in zip(values,values[1:]):
            # Exact inequalities are verified by independently enclosed
            # algebraic square roots. Exact degenerate equality is permitted.
            for oi,ni in ((1,1),(2,2)):
                ol,oh=bg.endpoints(old[oi]);nl,nh=bg.endpoints(new[ni])
                if nh>oh and nl>oh:raise AssertionError('Cauchy block refinement grew')
                if nh<ol:strict+=1
                checks+=1
        global_a,global_b=values[0][3],values[0][4]
        for B,first,second,aa,bb in values:
            if bg.endpoints(second)[0]>min(bg.upper(global_a),bg.upper(global_b)):
                raise AssertionError('Minimum of Cauchy alternatives failed refinement')
            checks+=1
    return {'exact_positive_atom_families':120,'hierarchy_checks':checks,'strict_refinements':strict,
            'includes_zero_U_and_zero_derivative':True,'random_seed':28061024}


def production_transform_audit(root,cells,rates,bg,producer_name):
    producer=importlib.import_module(producer_name)
    source=root/(producer_name+'.py');initial=hashlib.sha256(source.read_bytes()).hexdigest()
    transport=producer.BlockTransport()
    ctx.prec=768
    partitions=(1,8,128,1024)
    checks=0;additive_checks=0;records=[]
    for rate in rates:
        prefixes={kind:[arb(0)] for kind in (0,1,2)}
        for i,c in enumerate(cells):
            prefixes[0].append(prefixes[0][-1]+literal_rho0(c,rate,bg))
            # The actual rational dominating atom weights are producer inputs,
            # not silently changed to higher-precision abs(mu1) values.
            weight=F(transport.cell_metadata[i]['rho1_positive_atom_weight'])
            node=F(c['actual_node'])
            prefixes[1].append(prefixes[1][-1]+bg.ball(weight)*(bg.ball(rate*node)).exp())
            prefixes[2].append(prefixes[2][-1]+literal_rho2(c,rate,bg))
        for B in partitions:
            for kind in (0,1,2):
                values=transport.vector(kind,rate,B)
                references=[]
                for block,(start,end) in enumerate(transport.blocks(B)):
                    reference=prefixes[kind][end]-prefixes[kind][start]
                    al,ah=bg.endpoints(values[block]);rl,rh=bg.endpoints(reference)
                    if not al<=rl<=rh<=ah:
                        raise AssertionError(('Producer block transform missed independent literal sum',B,kind,rate,block))
                    references.append(reference);checks+=1
                total=sum(values,arb(0));literal_total=prefixes[kind][-1]
                tl,th=bg.endpoints(total);ll,lh=bg.endpoints(literal_total)
                if not tl<=ll<=lh<=th:raise AssertionError('Block additivity lost transform target')
                additive_checks+=1
                records.append({'B':B,'rate':str(rate),'measure':kind,
                                'block_checks':len(values),'sum_interval':bg.serial(total),
                                'literal_768bit_total':bg.serial(literal_total),'passed':True})
    if hashlib.sha256(source.read_bytes()).hexdigest()!=initial:
        raise AssertionError('Producer source changed during independent audit')
    return {'status':'passed','production_block_transform_checks':checks,
            'block_additivity_checks':additive_checks,'partitions':list(partitions),
            'producer_source_sha256':initial,'records':records}


def production_hierarchy_audit(root,summary_path):
    data=json.loads(summary_path.read_text())
    records=data['records'];models=json.loads((root/'results/generalisation_protocol.json').read_text())['models']
    partitions=tuple(data['transport_profile']['partitions'])
    expected={(m['model_id'],cap,B) for m in models for cap in ('1','10','25',None) for B in partitions}
    keys=[(r['model_id'],r['cap'],r['B']) for r in records]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:raise AssertionError('Hierarchy cohort incomplete or duplicated')
    source_checks=0
    for name,digest in data['source_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError(('Source hash differs from hierarchy evidence',name))
        source_checks+=1
    arithmetic_checks=0;refinement_checks=0;raw_increases=0
    for model in models:
        for cap in ('1','10','25',None):
            cohort=sorted((r for r in records if r['model_id']==model['model_id'] and r['cap']==cap),key=lambda r:r['B'])
            tail=F(cohort[0]['tail_upper']);parents={};previous_raw={}
            for r in cohort:
                if F(r['tail_upper'])!=tail:raise AssertionError('Tail changes across matched hierarchy')
                if F(r['FTC_financial_transfer_upper'])!=F(r['FTC_transport_upper'])+tail:
                    raise AssertionError('FTC full financial allowance omits or duplicates tail')
                if F(r['second_order_financial_transfer_upper'])!=F(r['centred_first_transport_upper'])+F(r['second_Peano_transport_upper'])+tail:
                    raise AssertionError('Peano allowance omits a retained component or duplicates tail')
                arithmetic_checks+=2
                for kind,field in (('FTC','FTC_financial_transfer_upper'),('second_order','second_order_financial_transfer_upper')):
                    raw=F(r[field]);hierarchy=F(r[kind+'_hierarchy_upper'])
                    if hierarchy!=min(raw,parents.get(kind,raw)):raise AssertionError('Parent/current minimum not retained exactly')
                    if kind in parents and hierarchy>parents[kind]:raise AssertionError('Numerical hierarchy is not monotone')
                    if kind in previous_raw and raw>previous_raw[kind]:raw_increases+=1
                    if r[kind+'_hierarchy_status']!=('passed' if hierarchy<F(1,200) else 'retained allocation failure'):
                        raise AssertionError('Admitted/failure status has changed allocation')
                    parents[kind]=hierarchy;previous_raw[kind]=raw;refinement_checks+=1
    original=json.loads((root/'results/tensor_transfer.json').read_text())['records']
    original_map={(r['model_id'],r['cap']):r for r in original}
    aggregate=[];maximum_difference=F(0);maximum_relative=F(0)
    for r in records:
        if r['B']!=1:continue
        old=original_map[(r['model_id'],r['cap'])]
        amount=F(r['second_order_financial_transfer_upper']);old_amount=F(old['financial_transfer_upper'])
        difference=abs(amount-old_amount);relative=difference/old_amount if old_amount else difference
        maximum_difference=max(maximum_difference,difference);maximum_relative=max(maximum_relative,relative)
        if (amount<F(1,200))!=(old_amount<F(1,200)):
            raise AssertionError('One-block aggregate has changed the original allocation decision')
        aggregate.append({'model_id':r['model_id'],'cap':r['cap'],'absolute_upper_difference':str(difference),
                          'relative_upper_difference':str(relative),'same_allocation_decision':True})
    return {'status':'passed','cohort_records':len(records),'source_hash_checks':source_checks,
            'full_tail_arithmetic_checks':arithmetic_checks,'numerical_hierarchy_checks':refinement_checks,
            'raw_interval_upper_increases_retained':raw_increases,
            'aggregate_comparison':{'records':aggregate,'maximum_absolute_upper_difference':str(maximum_difference),
                                    'maximum_relative_upper_difference':str(maximum_relative),
                                    'scope':'Same mathematical aggregate closure; literal upper bounds may differ under outward arithmetic schedules. Numerical difference is reported, not treated as exact equality.'},
            'summary_sha256':hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            'scope':'Complete frozen cohort, common full tails, all failure statuses and exact parent/current minimum hierarchy. This audit does not infer that an endpoint bound is globally optimal.'}


def independent_geometry_audit(root,bg,producer_name):
    """Literal small-grid factorisation with independently formed polynomials."""
    producer=importlib.import_module(producer_name)
    protocol=json.loads((root/'results/generalisation_protocol.json').read_text())['models']
    transport=producer.BlockTransport(b_y=4)
    actual=[]
    for assets in (2,4,8):
        model=next(m for m in protocol if m['assets']==assets)
        ctx.prec=512
        actual.append((model,producer.GroupedGeometry(model,transport).record(8)))
    ctx.prec=1024
    mass,cells=bg.gaussian_cells(F(6),4,28,1024)
    cache={};rounded_cache={};cont_cache={}
    def transform(kind,rate,start,end):
        key=(kind,rate,start,end)
        if key not in cache:
            fn=literal_rho0 if kind==0 else literal_rho2
            if kind==1:
                cache[key]=sum((bg.ball(F(transport.cell_metadata[i]['rho1_positive_atom_weight']))*(bg.ball(rate*F(cells[i]['actual_node']))).exp() for i in range(start,end)),arb(0))
            else:cache[key]=sum((fn(cells[i],rate,bg) for i in range(start,end)),arb(0))
        return cache[key]
    def continuous(rate):
        if rate not in cont_cache:
            v=bg.ball(rate);cont_cache[rate]=(v*v/2).exp()*(bg.Phi(bg.ball(6)-v)-bg.Phi(bg.ball(-6)-v))
        return cont_cache[rate]
    def rounded(rate):
        if rate not in rounded_cache:
            rounded_cache[rate]=mass/16*sum(((bg.ball(rate*F(c['actual_node']))).exp() for c in cells),arb(0))
        return rounded_cache[rate]
    def add(a,b):
        out=dict(a)
        for k,v in b.items():out[k]=out.get(k,F(0))+v
        return {k:v for k,v in out.items() if v}
    def mul(a,b):
        out={}
        for x,a0 in a.items():
            for y,b0 in b.items():
                key=tuple(i+j for i,j in zip(x,y));out[key]=out.get(key,F(0))+a0*b0
        return {k:v for k,v in out.items() if v}
    checks=0;models=[]
    for model,evidence in actual:
        n=model['assets'];p=model['outer_dimension'];w=tuple(map(F,model['weights']))
        raw_s=tuple(tuple(map(F,x)) for x in model['loadings']);mu=tuple(map(F,model['log_drifts']))
        mg=sum(a*b for a,b in zip(w,mu));rg=tuple(sum(w[i]*raw_s[i][j] for i in range(n)) for j in range(p))
        d=tuple(x-mg for x in mu);s=tuple(tuple(x-y for x,y in zip(row,rg)) for row in raw_s)
        monomial_cache={}
        for row in evidence['coordinates']:
            j=row['coordinate'];zero=(0,)*n;T={};J={};JJ={}
            for i in range(n):
                a=tuple(int(k==i) for k in range(n));T[a]=w[i];J[a]=w[i]*s[i][j];JJ[a]=w[i]*s[i][j]**2
            U=add(T,{zero:F(-1)});U2=mul(U,U);J2=mul(J,J)
            polys={'U2':U2,'J2':J2,'U_JJ':mul(U,JJ),'U_J2':mul(U,J2),'U2_JJ':mul(U2,JJ),'U2_J2':mul(U2,J2),'U4':mul(U2,U2)}
            for block in row['blocks']:
                start,end=transport.blocks(8)[block['block']]
                for kind in (0,1,2):
                    for name,rec in block['moments'][str(kind)].items():
                        reference=arb(0)
                        for alpha,coefficient in polys[name].items():
                            key=(alpha,j,kind,start,end)
                            if key not in monomial_cache:
                                rates=tuple(sum(F(a)*si[k] for a,si in zip(alpha,s)) for k in range(p))
                                value=(bg.ball(sum(F(a)*di for a,di in zip(alpha,d)))).exp()
                                for k,rate in enumerate(rates):
                                    value*=transform(kind,rate,start,end) if k==j else (continuous(rate) if k<j else rounded(rate))
                                monomial_cache[key]=value
                            reference+=bg.ball(coefficient)*monomial_cache[key]
                        rl,rh=bg.endpoints(reference);rl=max(F(0),rl)
                        al,ah=F(rec['lower_exact_fraction']),F(rec['upper_exact_fraction'])
                        if not al<=rl<=rh<=ah:
                            raise AssertionError(('Signed-polynomial block moment missed literal factorisation',n,j,kind,name,block['block']))
                        checks+=1
        models.append({'assets':n,'outer_dimension':p,'blocks':8,'nodes_per_coordinate':16,'passed':True})
    return {'status':'passed','literal_hybrid_moment_checks':checks,'reference_precision_bits':1024,
            'models':models,'distinct_literal_transport_queries':len(cache),
            'scope':'Independent centred exponential-polynomial construction and literal scalar transforms on a 16-node audit grid. Tests implementation for dimensions 1,3,7, all admitted moment kinds and transport measures; not an additional financial result.'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--producer',default='block_certificate_comparison')
    ap.add_argument('--skip-producer',action='store_true')
    ap.add_argument('--hierarchy-summary',type=Path)
    ap.add_argument('--out',type=Path)
    args=ap.parse_args();audit_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    bg=root_module(args.root);ctx.prec=768
    mass,cells=bg.gaussian_cells(F(6),10,28,768)
    rates=[F(-2),F(-1,2),F(-1,20),F(-1,1000),F(-1,10**40),F(0),F(1,10**40),F(1,1000),F(1,20),F(1,2),F(2)]
    selected=[0,1,127,255,511,512,768,1022,1023]
    rows=[]
    for ci in selected:
        for rate in rates:
            for kind,fn in ((0,literal_rho0),(2,literal_rho2)):
                exact=fn(cells[ci],rate,bg)
                quad=positive_quadrature(cells[ci],rate,kind,bg)
                el,eh=bg.endpoints(exact);ql,qh=bg.endpoints(quad)
                if el<0 or not ql<=el<=eh<=qh:
                    raise AssertionError(('Positive quadrature missed the literal transform',ci,rate,kind))
                rows.append({'cell':ci,'rate':str(rate),'measure':kind,'literal_interval':bg.serial(exact),
                             'positive_interval_quadrature':bg.serial(quad),'passed':True})
    production=None if args.skip_producer else production_transform_audit(args.root,cells,rates,bg,args.producer)
    hierarchy=production_hierarchy_audit(args.root,args.hierarchy_summary) if args.hierarchy_summary else None
    geometry=None if args.skip_producer else independent_geometry_audit(args.root,bg,args.producer)
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=audit_source_sha256:
        raise AssertionError('Independent audit source changed during execution')
    result={'schema':'qbp-independent-block-transport-review-v1','status':'passed',
            'precision_bits':768,'literal_vs_positive_quadrature_checks':len(rows),
            'cauchy':cauchy_audit(bg),'records':rows,'production_transform_audit':production,'production_hierarchy_audit':hierarchy,'independent_geometry_audit':geometry,
            'scope':'Exact rho0/rho2 cell transforms checked against independent positive interval quadrature; positive-atom block closures checked under refinement. Production block transforms are compared with literal 768-bit Gaussian cell sums after source freeze.',
            'source_sha256':audit_source_sha256}
    out=args.out or args.root/'results/block_transport_math_audit.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print({k:result[k] for k in ('status','literal_vs_positive_quadrature_checks')},result['cauchy'])


if __name__=='__main__':main()
