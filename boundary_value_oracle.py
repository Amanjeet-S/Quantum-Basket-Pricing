"""Complete coherent graphs for the conditional cash-Gamma remainder.

The new graphs share the retained integer compiler and clean comparator.
No joint latent/payoff table is used. Each scalar table is one-dimensional.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math
import integral_value_oracle as old
from boundary_derivative_certificate import payoff_norms

IntegralGraph=old.IntegralGraph
ScalarSpec=old.ScalarSpec
IntegralNode=old.IntegralNode
table_key=old.table_key
compile_value=old.compile_value
compile_flag=old.compile_flag
physical_sequence=old.physical_sequence
emit_node=old.emit_node


def build_graph(model,cap=None,q=28,b_y=10,beta_bits=8,table_bits=12,outer_bound=6,
                U_star=None,M2=None,M1=None,kind='boundary_r2',payoff=None,
                headroom_bits=None,distribution='quantile',retained_bound=None,
                output_encoding='centred'):
    payoff=('vanilla' if cap is None else 'capped') if payoff is None else payoff
    if payoff not in ('capped','vanilla'):raise ValueError('Unknown payoff')
    if kind not in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):raise ValueError('Unknown correction graph')
    if distribution not in ('quantile','uniform_box'):raise ValueError('Unknown input law')
    if output_encoding not in ('centred','positive_part','negative_part'):raise ValueError('Unknown output encoding')
    if payoff=='vanilla' and output_encoding!='centred':raise ValueError('Vanilla remainder has its native positive encoding')
    if kind=='positive_r1' and output_encoding!='centred':raise ValueError('First remainder has its native positive encoding')
    cap=None if payoff=='vanilla' else F(cap);L=F(outer_bound);M=F(retained_bound) if retained_bound is not None else L
    if L<=0 or M<=0 or M>L or cap is not None and cap<=0:raise ValueError('Positive valid financial domain required')
    if distribution=='uniform_box' and old.power_span(L)!=L:raise ValueError('Uniform box needs a dyadic half-span')
    if U_star is None or F(U_star)<0:raise ValueError('Certified finite U_star is required')
    U_star=F(U_star);w,d,s,mg,rg,lam=old.model_parameters(model);norm=payoff_norms(lam,cap)
    M2=F(norm['M2_upper']) if M2 is None else F(M2)
    if M2<F(norm['M2_upper']):raise ValueError('M2 must dominate the retained payoff certificate')
    mhi=mg+M*sum(map(abs,rg))
    if payoff=='capped':
        M1=F(norm['M1_upper']) if M1 is None else F(M1)
        if M1<F(norm['M1_upper']):raise ValueError('M1 must dominate the retained capped certificate')
        B1=min(cap,M1*U_star)
    else:
        from boundary_value_certificate import scalar_bounds
        finite_M1=scalar_bounds('vanilla_call',mhi-1,mhi,{'lambda':lam,'K':100})['first']
        M1=finite_M1 if M1 is None else F(M1)
        if M1<finite_M1:raise ValueError('Finite vanilla first slope is not certified')
        B1=M1*U_star
    B2=M2*U_star**2/2;normaliser=B1 if kind=='positive_r1' else B2
    positive_output=payoff=='vanilla' or kind=='positive_r1'
    g=IntegralGraph(q,table_bits,headroom_bits)
    g.configuration={'model_id':model.get('model_id'),'kind':kind,'payoff':payoff,
        'cap':None if cap is None else str(cap),'lambda':str(lam),'K':'100','L':str(L),
        'retained_coordinate_bound':str(M),'U_star':str(U_star),'M1':str(M1),'M2':str(M2),
        'B2':str(B2),'positive_R1_normaliser':str(B1),'monetary_signal_normaliser':str(normaliser),
        'distribution':distribution,'b_y':b_y,'beta_bits':beta_bits if kind=='boundary_r2' else 0,
        'strike_bits':0,'table_bits_requested':table_bits,'q':q,
        'output_encoding':'positive' if positive_output else output_encoding,
        'monetary_probability_decode_factor':1 if positive_output or output_encoding!='centred' else 2,
        'correction_order':1 if kind=='positive_r1' else 2,'degenerate_zero_correction':U_star==0,
        'financial_target':'Unchanged classical risk-neutral arithmetic basket, exact analytic Gaussian c0 or c0+c1 retained separately.',
        'payoff_norm_sha256':hashlib.sha256((Path(__file__).resolve().parent/'results/boundary_derivative_norms.json').read_bytes()).hexdigest()}
    if U_star==0:
        g.semantic_nodes={}
        return g.finish(g.const(F(1,2) if not positive_output and output_encoding=='centred' else 0))
    coordinates=[];Mc=L
    if distribution=='quantile' and M<L:Mc=min(L,F(math.ceil(M*(1<<q))+2,1<<q))
    for j in range(len(rg)):
        raw=g.inp(f'iy{j}',b_y)
        y=g.op('quantile',raw,meta=L,bits=b_y,bounds=(L-M,L+M)) if distribution=='quantile' else g.uniform(raw,b_y,2*L)
        if Mc<L:y=g.sub(y,g.const(L-Mc))
        coordinates.append(y)
    exps=[]
    for drift,exposure in zip(d,s):
        encoded,B=old.shifted_affine(g,coordinates,exposure,drift,Mc,Mc)
        exps.append(g.scalar('exp',encoded,-B,2*B,output_bounds=(0,old.conservative_elementary_upper('exp',-B,2*B))))
    T=g.const(0)
    for exponential,weight in zip(exps,w):T=g.add(T,g.scale(exponential,weight))
    U=g.clamp(g.positive_sub(T,g.const(1)),U_star)
    mG,BG=old.shifted_affine(g,coordinates,rg,mg,Mc,Mc)
    g.semantic_nodes={'T':T,'U':U,'mG_offset':mG}
    span=old.power_span(max(U_star,F(1,1<<q)))
    params={'lambda':lam,'K':100}
    if cap is not None:params['cap']=cap
    if kind in ('boundary_r2','stable_gauss_r2'):
        ratio=g.clamp(g.scale(U,1/U_star),1);ratio2=g.mul(ratio,ratio)
        g.semantic_nodes.update({'U_ratio':ratio,'U_ratio_square':ratio2})
        def at(z,label):
            ZU=g.mul(z,U)
            logt=g.scalar('log',ZU,1,span,output_bounds=(0,span));g.bounds[logt]=(0,U_star)
            invt=g.scalar('inv_square',ZU,1,span,output_bounds=(0,1))
            BT=old.power_span(max(BG,abs(mg)+Mc*sum(map(abs,rg))+U_star+F(8,1<<q)))
            m_t=g.add(g.add(mG,g.const(BT-BG)),logt)
            # Correlated shifted argument lies inside the declared padded span.
            g.bounds[m_t]=(mg-Mc*sum(map(abs,rg))+BT,mg+Mc*sum(map(abs,rg))+BT+U_star)
            p=dict(params,M2=M2)
            h=g.scalar('boundary_density' if positive_output else 'boundary_h',m_t,-BT,2*BT,p,(0,1))
            amplitude=g.clamp(g.mul(ratio2,invt),1)
            probability=g.mul(amplitude,h) if positive_output else g.add(g.scale(g.sub(g.const(1),amplitude),F(1,2)),g.mul(amplitude,h))
            for key,value in {'Z':z,'ZU':ZU,'logt':logt,'invt_square':invt,'m_t_offset':m_t,'kernel':h,'amplitude_envelope':amplitude,'probability':probability}.items():g.semantic_nodes[label+key]=value
            return probability
        if kind=='boundary_r2':
            az=g.inp('iz0',beta_bits);bz=g.inp('iz1',beta_bits)
            Z=g.uniform(g.minimum(az,bz,beta_bits),beta_bits)
            probability=at(Z,'beta_')
        else:
            from boundary_value_certificate import gauss_constants
            constants=gauss_constants(q);g.configuration['gauss_constants']=constants
            zm,zp=(g.const(F(constants[k]['actual'])) for k in ('z_minus','z_plus'))
            pm,pp=at(zm,'minus_'),at(zp,'plus_')
            wm,wp=(F(constants[k]['actual']) for k in ('w_minus','w_plus'))
            probability=g.add(g.scale(pm,wm),g.scale(pp,wp))
            g.semantic_nodes.update({'minus_probability':pm,'plus_probability':pp})
    else:
        logT=g.scalar('log',U,1,span,output_bounds=(0,span));g.bounds[logT]=(0,U_star)
        BA=old.power_span(max(BG,abs(mg)+Mc*sum(map(abs,rg))+U_star+F(8,1<<q)))
        mA=g.add(g.add(mG,g.const(BA-BG)),logT)
        names=('vanilla_call','vanilla_call_prime') if payoff=='vanilla' else ('cap','cap_prime')
        if payoff=='vanilla':
            from boundary_value_certificate import scalar_bounds
            pa_upper=scalar_bounds(names[0],-BA,BA,params)['range'][1]
            pg_upper=scalar_bounds(names[0],-BG,BG,params)['range'][1]
            prime_upper=scalar_bounds(names[1],-BG,BG,params)['range'][1]
        else:pa_upper=pg_upper=cap;prime_upper=M1
        PA=g.scalar(names[0],mA,-BA,2*BA,params,(0,pa_upper))
        PG=g.scalar(names[0],mG,-BG,2*BG,params,(0,pg_upper))
        if kind=='raw_direct_r2':
            prime=g.scalar(names[1],mG,-BG,2*BG,params,(0,prime_upper))
            negative=g.scale(g.add(PG,g.mul(U,prime)),1/(B2 if positive_output else 2*B2))
            positive=g.scale(PA,1/(B2 if positive_output else 2*B2))
            if not positive_output:positive=g.add(g.const(F(1,2)),positive)
            probability=g.clamp(g.sub(positive,negative,final_signed=True),1,signed=True)
            g.semantic_nodes['PG_prime']=prime
        else:
            probability=g.clamp(g.sub(g.scale(PA,1/B1),g.scale(PG,1/B1),final_signed=True),1,signed=True)
        g.semantic_nodes.update({'logT':logT,'mA_offset':mA,'PA':PA,'PG':PG})
    if distribution=='uniform_box':
        squares=[g.scalar('square',y,-L,2*L,output_bounds=(0,L*L)) for y in coordinates]
        sumsq=g.const(0)
        for square in squares:sumsq=g.add(sumsq,square)
        BD=old.power_span(max(F(1),len(coordinates)*L*L/2))
        argument=g.sub(g.const(BD),g.scale(sumsq,F(1,2)))
        likelihood=g.scalar('exp',argument,-BD,BD,output_bounds=(0,1))
        probability=g.mul(likelihood,probability) if positive_output else g.add(g.scale(g.sub(g.const(1),likelihood),F(1,2)),g.mul(likelihood,probability))
        g.semantic_nodes['density_relative_weight']=likelihood
        g.configuration['uniform_box_monetary_factor']='(2L/sqrt(2pi))^p'
        g.configuration['uniform_box_dimension']=len(coordinates)
    probability=g.clamp(probability,1);g.semantic_nodes['centred_probability']=probability
    if not positive_output:
        if output_encoding=='positive_part':probability=g.positive_sub(g.scale(probability,2),g.const(1))
        elif output_encoding=='negative_part':probability=g.positive_sub(g.const(1),g.scale(probability,2))
    return g.finish(g.clamp(probability,1))


def make_tables(graph,cache=None,return_certificates=False):
    from boundary_value_certificate import scalar_certificate,quantile_certificate
    cache={} if cache is None else cache;tables={};certificates={}
    for name,node in graph.nodes.items():
        if node.kind not in ('quantile','function'):continue
        key=table_key(name,node);full=(key,graph.q,graph.width)
        if full not in cache:
            if node.kind=='quantile':
                certificate=quantile_certificate(node.meta,node.bits,graph.q)
                words=[n['offset_integer'] for n in certificate['nodes']]
            else:
                spec=node.meta;certificate=scalar_certificate(spec.name,spec.lo,spec.span,dict(spec.parameters),graph.q,node.bits)
                values=[n['integer'] for n in certificate['nodes']];words=[]
                for j,value in enumerate(values):
                    delta=values[j+1]-value if j+1<len(values) else 0
                    if min(value,abs(delta))<0 or max(value,abs(delta))>=1<<graph.width:raise ArithmeticError('Scalar table overflow')
                    words.append(value|(abs(delta)<<graph.width)|(int(delta<0)<<(2*graph.width)))
            cache[full]=(tuple(words),certificate)
        tables[key],certificates[name]=cache[full]
    graph.table_certificates=certificates
    return (tables,certificates) if return_certificates else tables


def counted_composition(graph,tables,cache=None,mode='exact7t'):
    result=old.counted_composition(graph,tables,cache,mode)
    result['schema']='qbp-boundary-full-financial-oracle-v1'
    result['retained_emitter_sha256']=result['source_sha256']
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result
