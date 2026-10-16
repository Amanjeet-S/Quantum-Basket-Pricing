"""Full finite financial graphs for integral and direct Taylor corrections.

All financial scalar kernels are one-dimensional tables. Interpolation stores
an unsigned magnitude and a Boolean decreasing flag; it never multiplies a
modular negative integer. Complete clean financial preparations are emitted.
No precision or competitive-cost claim is inferred from the fractional width.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_CEILING
from pathlib import Path
from statistics import NormalDist
from collections import Counter
import hashlib,json,math
import value_oracle as vo
import architecture as arch
import quantum_gates as qg
import reversible as base


@dataclass(frozen=True)
class ScalarSpec:
    name:str
    lo:F
    span:F
    parameters:tuple
    output_bounds:tuple

    def as_dict(self):
        return {'name':self.name,'lo':str(self.lo),'span':str(self.span),
                'parameters':{k:str(v) for k,v in self.parameters},
                'output_bounds':list(map(str,self.output_bounds))}


@dataclass
class IntegralNode:
    kind:str
    args:tuple
    meta:object=None
    bits:int=0
    output_bits:int=0


def rounded(x):
    x=F(x)+F(1,2);return x.numerator//x.denominator


def power_span(x):
    x=F(x)
    if x<=0:raise ValueError('Positive span required')
    value=F(1)
    while value<x:value*=2
    while value/2>=x:value/=2
    return value


def span_log(span):
    span=F(span)
    if span<=0 or span.numerator&(span.numerator-1) or span.denominator&(span.denominator-1):raise ValueError('Scalar span must be a power of two')
    return span.numerator.bit_length()-span.denominator.bit_length()


def conservative_elementary_upper(name,lo,span):
    with localcontext() as ctx:
        ctx.prec=100
        x=Decimal((lo+span).numerator)/Decimal((lo+span).denominator)
        if name=='exp':value=x.exp()
        else:value=x.ln()
        # An integer guard exceeds decimal rounding in the declared domains.
        return F(int(value.to_integral_value(rounding=ROUND_CEILING))+1)


class IntegralGraph:
    def __init__(self,q,table_bits,headroom_bits=None):
        if not isinstance(q,int) or q<2 or not isinstance(table_bits,int) or table_bits<1:raise ValueError('Positive binary precisions required')
        self.q=q;self.table_bits=table_bits;self.width=0;self.requested_headroom=headroom_bits
        self.nodes={};self.inputs={};self.bounds={};self.output='';self.constant_names={};self.configuration={}

    def inp(self,name,bits):
        self.inputs[name]=bits;self.bounds[name]=(F(0),F((1<<bits)-1));return name

    def op(self,kind,*args,meta=None,bits=0,output_bits=0,bounds=None):
        name=f'n{len(self.nodes)}';self.nodes[name]=IntegralNode(kind,tuple(args),meta,bits,output_bits)
        self.bounds[name]=tuple(map(F,bounds));return name

    def const(self,x):
        x=F(x)
        if x<0:raise ValueError('Only unsigned constants are allowed')
        if x in self.constant_names:return self.constant_names[x]
        name=self.op('constant',meta=x,bounds=(x,x));self.constant_names[x]=name;return name

    def add(self,a,b):
        aa,bb=self.bounds[a],self.bounds[b]
        return self.op('add',a,b,bounds=(aa[0]+bb[0],aa[1]+bb[1]))

    def sub(self,a,b,final_signed=False):
        aa,bb=self.bounds[a],self.bounds[b];lower,upper=aa[0]-bb[1],aa[1]-bb[0]
        if lower<0 and not final_signed:raise ValueError(('Unsigned subtraction lacks a nonnegative domain',a,b,lower))
        return self.op('signed_sub' if final_signed else 'sub',a,b,bounds=(lower,upper))

    def positive_sub(self,a,b):
        aa,bb=self.bounds[a],self.bounds[b]
        return self.op('positive_sub',a,b,bounds=(max(F(0),aa[0]-bb[1]),max(F(0),aa[1]-bb[0])))

    def mul(self,a,b):
        aa,bb=self.bounds[a],self.bounds[b]
        if aa[0]<0 or bb[0]<0:raise ValueError('Signed multiplication is forbidden')
        return self.op('mul',a,b,bounds=(aa[0]*bb[0],aa[1]*bb[1]))

    def scale(self,a,x):
        x=F(x);aa=self.bounds[a]
        if x<0 or aa[0]<0:raise ValueError('Signed scaling is forbidden')
        return self.op('scale',a,meta=x,bounds=(aa[0]*x,aa[1]*x))

    def clamp(self,a,cap=1,signed=False):
        aa=self.bounds[a];cap=F(cap)
        return self.op('clamp_signed' if signed else 'clamp',a,meta=cap,
                       bounds=(max(F(0),min(cap,aa[0])),max(F(0),min(cap,aa[1]))))

    def scalar(self,name,a,lo,span,params=None,output_bounds=None):
        lo,span=F(lo),F(span);log=span_log(span)
        if output_bounds is None:raise ValueError('A scalar output bound is required')
        clipped=self.clamp(a,span)
        spec=ScalarSpec(name,lo,span,tuple(sorted((k,F(v)) for k,v in (params or {}).items())),tuple(map(F,output_bounds)))
        bits=min(self.table_bits,self.q+log)
        if bits<1:raise ValueError('Insufficient fractional resolution for scalar span')
        return self.op('function',clipped,meta=spec,bits=bits,bounds=output_bounds)

    def uniform(self,a,bits,span=1):
        if self.q<bits+1:raise ValueError('Uniform midpoint needs q>=bits+1')
        value=self.op('uniform',a,meta=bits,bounds=(F(1,1<<(bits+1)),1-F(1,1<<(bits+1))))
        return self.scale(value,span) if F(span)!=1 else value

    def minimum(self,a,b,bits):
        return self.op('minimum',a,b,meta=bits,output_bits=bits,bounds=(0,(1<<bits)-1))

    def finish(self,output):
        self.output=output
        largest=max([F(1)]+[abs(x) for name,node in self.nodes.items() if not node.output_bits for x in self.bounds[name]])
        # Interpolation candidates and signed endpoints also need safe headroom.
        h=max(3,power_span(8*largest+8).numerator.bit_length())
        if self.requested_headroom is not None:
            if self.requested_headroom<h:raise ValueError(('Insufficient integer headroom',self.requested_headroom,h))
            h=self.requested_headroom
        self.width=self.q+h
        for node in self.nodes.values():
            if not node.output_bits:node.output_bits=self.width
        self.configuration['headroom_bits']=h
        return self

    def evaluate(self,inputs,tables,trace=False):
        unit=1<<self.q;mask=(1<<self.width)-1;values=dict(inputs)
        for name,node in self.nodes.items():
            args=[values[a] for a in node.args]
            if node.kind=='constant':v=rounded(node.meta*unit)
            elif node.kind=='add':v=args[0]+args[1]
            elif node.kind in ('sub','signed_sub'):v=args[0]-args[1]
            elif node.kind=='positive_sub':v=max(0,args[0]-args[1])
            elif node.kind=='mul':v=args[0]*args[1]>>self.q
            elif node.kind=='scale':
                integer,drop=scale_integer(node.meta,self.q);v=args[0]*integer>>drop
            elif node.kind=='uniform':v=(2*args[0]+1)*(unit>>(node.meta+1))
            elif node.kind=='minimum':v=min(args)
            elif node.kind=='quantile':v=tables[table_key(name,node)][args[0]]
            elif node.kind=='function':
                f=self.q+span_log(node.meta.span)-node.bits
                index,rem=args[0]>>f,args[0]&((1<<f)-1)
                word=tables[table_key(name,node)][index]
                baseword=word&mask;mag=(word>>self.width)&mask;decreasing=(word>>(2*self.width))&1
                adjustment=mag*rem>>f
                v=baseword-adjustment if decreasing else baseword+adjustment
            elif node.kind in ('clamp','clamp_signed'):
                v=args[0]
                if node.kind=='clamp_signed' and v>=(1<<(self.width-1)):v-=1<<self.width
                v=max(0,min(v,rounded(node.meta*unit)))
            else:raise ValueError(node.kind)
            values[name]=v&((1<<node.output_bits)-1)
        return values if trace else values[self.output]


def scale_integer(x,q):
    x=F(x)
    if x.denominator&(x.denominator-1)==0:return x.numerator,x.denominator.bit_length()-1
    return rounded(x*(1<<q)),q


def table_key(name,node):
    if node.kind=='quantile':return ('quantile',node.bits,str(node.meta))
    return ('scalar',json.dumps(node.meta.as_dict(),sort_keys=True),node.bits)


def clean_selected_difference(c,a,b,positive=False):
    if not positive:return c.unsigned_subtract(a,b)
    out=c.alloc(len(a));snap=c.snapshot();difference=c.unsigned_subtract(a,b);flag=c.unsigned_less_than(b,a)[0];end=c.gate_count
    for x,y in zip(difference,out):c.toffoli(flag,x,y)
    c.uncompute(snap,end_gate_count=end);return out


def clean_clamp(c,a,cap,signed=False):
    out=c.alloc(len(a));snap=c.snapshot();upper=c.alloc(len(a))
    for k,w in enumerate(upper):
        if cap>>k&1:c.x(w)
    below=c.unsigned_less_than(a,upper)[0];positive=c.alloc(1)[0];c.x(positive)
    if signed:c.cnot(a[-1],positive)
    take_value=c.alloc(1)[0];c.toffoli(positive,below,take_value)
    take_cap=c.alloc(1)[0];c.x(below);c.toffoli(positive,below,take_cap);c.x(below)
    end=c.gate_count
    for k,w in enumerate(out):
        c.toffoli(take_value,a[k],w)
        if cap>>k&1:c.cnot(take_cap,w)
    c.uncompute(snap,end_gate_count=end);return out


def emit_node(c,node,args,q,width,table=None):
    kind=node.kind
    if kind=='constant':
        out=c.alloc(width);integer=rounded(node.meta*(1<<q))
        for k,w in enumerate(out):
            if integer>>k&1:c.x(w)
        return out
    if kind=='add':return c.unsigned_add(*args)
    if kind in ('sub','signed_sub','positive_sub'):return clean_selected_difference(c,*args,positive=kind=='positive_sub')
    if kind=='mul':return c.multiply(*args,output_bits=width,drop_lsb=q)
    if kind=='scale':
        integer,drop=scale_integer(node.meta,q);return c.multiply_constant(args[0],integer,output_bits=width,drop_lsb=drop)
    if kind=='uniform':
        out=c.alloc(width);start=q-node.meta;c.x(out[start-1])
        for k,w in enumerate(args[0]):c.cnot(w,out[start+k])
        return out
    if kind=='minimum':
        out=c.alloc(node.output_bits);snap=c.snapshot();less=c.unsigned_less_than(*args)[0];result=c.alloc(node.output_bits)
        for a,b,o in zip(args[0],args[1],result):
            difference=c.alloc(1)[0];c.cnot(a,difference);c.cnot(b,difference);c.cnot(b,o);c.toffoli(less,difference,o)
        end=c.gate_count
        for x,y in zip(result,out):c.cnot(x,y)
        c.uncompute(snap,end_gate_count=end);return out
    if kind=='quantile':return c.lookup(args[0],table,width)
    if kind=='function':
        fractional=q+span_log(node.meta.span)-node.bits
        address=args[0][fractional:fractional+node.bits+1]
        out=c.alloc(width);snap=c.snapshot();packed=c.lookup(address,table,2*width+1)
        value,mag,decreasing=packed[:width],packed[width:2*width],packed[-1]
        directions={int(word>>(2*width)&1) for word in table if (word>>width)&((1<<width)-1)}
        adjustment=c.multiply(mag,args[0][:fractional],output_bits=width,drop_lsb=fractional) if fractional and directions else c.alloc(width)
        if not directions:result=value
        elif directions=={0}:result=c.unsigned_add(value,adjustment)
        elif directions=={1}:result=c.unsigned_subtract(value,adjustment)
        else:
            plus,minus=c.alloc(width),c.alloc(width)
            for x,p,n in zip(adjustment,plus,minus):c.toffoli(decreasing,x,n)
            c.x(decreasing)
            for x,p in zip(adjustment,plus):c.toffoli(decreasing,x,p)
            c.x(decreasing)
            result=c.unsigned_subtract(c.unsigned_add(value,plus),minus)
        end=c.gate_count
        for x,y in zip(result,out):c.cnot(x,y)
        c.uncompute(snap,end_gate_count=end);return out
    if kind in ('clamp','clamp_signed'):return clean_clamp(c,args[0],rounded(node.meta*(1<<q)),kind=='clamp_signed')
    raise ValueError(kind)


def model_parameters(model):
    w=tuple(map(F,model['weights']));rows=tuple(tuple(map(F,r)) for r in model['loadings']);mu=tuple(map(F,model['log_drifts']));lam=F(model['lambda'])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    mg=dot(w,mu);rg=tuple(sum(wi*r[j] for wi,r in zip(w,rows)) for j in range(len(rows[0])))
    d=tuple(a-mg for a in mu);s=tuple(tuple(a-b for a,b in zip(row,rg)) for row in rows)
    if sum(w)!=1 or min(w)<=0 or lam<=0:raise ValueError('Positive model weights and common volatility required')
    if dot(w,d)!=0 or any(sum(wi*row[j] for wi,row in zip(w,s))!=0 for j in range(len(rg))):raise ValueError('Uncentred geometry')
    return w,d,s,mg,rg,lam


def shifted_affine(g,coordinates,coefficients,constant,L,bound):
    # Offset is chosen from the actual interval, with a rounding guard.
    padding=F(2+len(coefficients)*(2*L+1),1<<g.q)
    B=power_span(abs(constant)+bound*sum(abs(x) for x in coefficients)+padding)
    positive=g.const(constant+B-L*sum(coefficients));negative=g.const(0)
    for value,coefficient in zip(coordinates,coefficients):
        if coefficient>0:positive=g.add(positive,g.scale(value,coefficient))
        elif coefficient<0:negative=g.add(negative,g.scale(value,-coefficient))
    encoded=g.sub(positive,negative)
    # The tighter correlated affine range is exact and exposes the semantics.
    low=constant+B-bound*sum(abs(x) for x in coefficients);high=constant+B+bound*sum(abs(x) for x in coefficients)
    g.bounds[encoded]=(low,high)
    return encoded,B


def build_graph(model,cap,q=28,b_y=10,beta_bits=8,strike_bits=8,table_bits=12,
                outer_bound=6,U_star=None,D2=None,D1=None,kind='integral_r2',
                headroom_bits=None,distribution='quantile',retained_bound=None,output_encoding='centred'):
    if kind not in ('integral_r2','direct_r2','positive_r1'):raise ValueError('Unknown financial correction')
    if distribution not in ('quantile','uniform_box'):raise ValueError('Unknown factor preparation')
    if output_encoding not in ('centred','positive_part','negative_part') or kind=='positive_r1' and output_encoding!='centred':raise ValueError('Unsupported correction encoding')
    cap,L=F(cap),F(outer_bound);M=F(retained_bound) if retained_bound is not None else L
    if cap<=0 or L<=0 or M<=0 or M>L:raise ValueError('Positive cap and valid retained coordinate bound required')
    if distribution=='uniform_box' and power_span(L)!=L:raise ValueError('Uniform-box prototype requires power-of-two L')
    if U_star is None or D2 is None or F(U_star)<0 or F(D2)<=0:raise ValueError('Certified U_star and D2 must be supplied')
    U_star,D2=F(U_star),F(D2);w,d,s,mg,rg,lam=model_parameters(model)
    if D1 is None:
        # D1=phi(0)/lambda, conservatively bounded without an irrational constant.
        D1=F(1,2)/lam
    D1=F(D1)
    g=IntegralGraph(q,table_bits,headroom_bits)
    B2=cap*D2*U_star*U_star/2;C1=min(cap,cap*D1*U_star)
    normaliser=C1 if kind=='positive_r1' else B2
    g.configuration={'model_id':model.get('model_id'),'kind':kind,'cap':str(cap),'lambda':str(lam),
        'L':str(L),'retained_coordinate_bound':str(M),'U_star':str(U_star),'D2':str(D2),'D1':str(D1),
        'B2':str(B2),'positive_R1_normaliser':str(C1),'monetary_signal_normaliser':str(normaliser),
        'distribution':distribution,'b_y':b_y,'beta_bits':beta_bits,'strike_bits':strike_bits,
        'table_bits_requested':table_bits,'q':q,'output_encoding':output_encoding,'degenerate_zero_correction':U_star==0,
        'financial_target':'Classical risk-neutral capped basket, continuous coefficients retained externally; box/grid/tail and value error require independent certificates.'}
    if U_star==0:return g.finish(g.const(F(1,2) if kind!='positive_r1' and output_encoding=='centred' else 0))
    coordinates=[]
    centred_offset=L
    if distribution=='quantile' and M<L:
        centred_offset=min(L,F(math.ceil(M*(1<<q))+2,1<<q))
    for j in range(len(rg)):
        raw=g.inp(f'iy{j}',b_y)
        if distribution=='quantile':
            y=g.op('quantile',raw,meta=L,bits=b_y,bounds=(L-M,L+M))
        else:y=g.uniform(raw,b_y,2*L)
        if centred_offset<L:y=g.sub(y,g.const(L-centred_offset))
        coordinates.append(y)
    exponential=[]
    for drift,exposure in zip(d,s):
        encoded,B=shifted_affine(g,coordinates,exposure,drift,centred_offset,centred_offset)
        upper=conservative_elementary_upper('exp',-B,2*B)
        exponential.append(g.scalar('exp',encoded,-B,2*B,output_bounds=(0,upper)))
    T=g.const(0)
    for value,weight in zip(exponential,w):T=g.add(T,g.scale(value,weight))
    U=g.clamp(g.positive_sub(T,g.const(1)),U_star)
    g.semantic_nodes={'T':T,'U':U}
    mG,BG=shifted_affine(g,coordinates,rg,mg,centred_offset,centred_offset);g.semantic_nodes['mG_offset']=mG
    if kind=='integral_r2':
        a=g.inp('iz0',beta_bits);b=g.inp('iz1',beta_bits)
        Z=g.uniform(g.minimum(a,b,beta_bits),beta_bits);S=g.uniform(g.inp('istrike',strike_bits),strike_bits)
        ZU=g.mul(Z,U);span=power_span(max(U_star,F(1,1<<q)))
        logt=g.scalar('log',ZU,1,span,output_bounds=(0,span));g.bounds[logt]=(0,U_star)
        invt=g.scalar('inv_square',ZU,1,span,output_bounds=(0,1))
        strike_input=g.scale(S,cap/100);strike_span=power_span(max(cap/100,F(1,1<<q)))
        strike=g.scalar('log',strike_input,1,strike_span,output_bounds=(0,strike_span))
        log_bound=g.bounds[logt][1];guard=F(8,1<<q)
        xiB=power_span((BG+log_bound+strike_span+guard)/lam)
        positive=g.add(g.add(mG,logt),g.const(lam*xiB-BG))
        encoded=g.scale(g.sub(positive,strike),1/lam)
        g.bounds[encoded]=(0,2*xiB)
        h=g.scalar('h',encoded,-xiB,2*xiB,{'lambda':lam,'D2':D2},(0,1))
        ratio=g.clamp(g.scale(U,1/U_star),1)
        amplitude=g.clamp(g.mul(g.mul(ratio,ratio),invt),1)
        probability=g.add(g.scale(g.sub(g.const(1),amplitude),F(1,2)),g.mul(amplitude,h))
        g.semantic_nodes.update({'Z':Z,'S':S,'ZU':ZU,'logt':logt,'invt_square':invt,'strike_log':strike,'xi_offset':encoded,'h':h,'amplitude_envelope':amplitude})
    else:
        span=power_span(max(U_star,F(1,1<<q)))
        logT=g.scalar('log',U,1,span,output_bounds=(0,span));g.bounds[logT]=(0,U_star)
        BA=power_span(max(BG,abs(mg)+centred_offset*sum(abs(x) for x in rg)+U_star+F(8,1<<q)))
        mA=g.add(g.add(mG,g.const(BA-BG)),logT)
        params={'lambda':lam,'cap':cap,'K':100}
        PA=g.scalar('cap',mA,-BA,2*BA,params,(0,cap))
        PG=g.scalar('cap',mG,-BG,2*BG,params,(0,cap))
        if kind=='direct_r2':
            prime=g.scalar('cap_prime',mG,-BG,2*BG,params,(0,100+cap))
            positive=g.add(g.const(F(1,2)),g.scale(PA,1/(2*B2)))
            negative=g.scale(g.add(PG,g.mul(U,prime)),1/(2*B2))
            probability=g.clamp(g.sub(positive,negative,final_signed=True),1,signed=True)
            g.semantic_nodes['PG_prime']=prime
        else:
            probability=g.clamp(g.sub(g.scale(PA,1/C1),g.scale(PG,1/C1),final_signed=True),1,signed=True)
        g.semantic_nodes.update({'logT':logT,'mA_offset':mA,'PA':PA,'PG':PG})
    if distribution=='uniform_box':
        squares=[]
        for coordinate in coordinates:squares.append(g.scalar('square',coordinate,-L,2*L,output_bounds=(0,L*L)))
        square_sum=g.const(0)
        for square in squares:square_sum=g.add(square_sum,square)
        BD=power_span(max(F(1),len(coordinates)*L*L/2))
        argument=g.sub(g.const(BD),g.scale(square_sum,F(1,2)))
        likelihood=g.scalar('exp',argument,-BD,BD,output_bounds=(0,1))
        if kind=='positive_r1':probability=g.mul(probability,likelihood)
        else:
            # Multiply the signed observable, rather than the probability offset.
            probability=g.add(g.scale(g.sub(g.const(1),likelihood),F(1,2)),g.mul(likelihood,probability))
        g.semantic_nodes['density_relative_weight']=likelihood
        g.configuration['uniform_box_monetary_factor']='(2L/sqrt(2pi))^p'
        g.configuration['uniform_box_dimension']=len(coordinates)
    probability=g.clamp(probability,1);g.semantic_nodes['centred_probability']=probability
    if output_encoding=='positive_part':probability=g.positive_sub(g.scale(probability,2),g.const(1))
    elif output_encoding=='negative_part':probability=g.positive_sub(g.const(1),g.scale(probability,2))
    return g.finish(g.clamp(probability,1))


def make_tables(graph,cache=None,return_certificates=False):
    from integral_value_certificate import scalar_certificate,quantile_certificate
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
                    if value<0 or value>=(1<<graph.width) or abs(delta)>=(1<<graph.width):raise ArithmeticError('Scalar table output overflow')
                    words.append(value|(abs(delta)<<graph.width)|(int(delta<0)<<(2*graph.width)))
            cache[full]=(tuple(words),certificate)
        tables[key],certificates[name]=cache[full]
    graph.table_certificates=certificates
    return (tables,certificates) if return_certificates else tables


def emit_body(c,graph,tables,inputs):
    registers=dict(inputs);blocks=[]
    for name,node in graph.nodes.items():
        start=c.gate_count
        registers[name]=emit_node(c,node,[registers[a] for a in node.args],graph.q,graph.width,tables.get(table_key(name,node)) if node.kind in ('quantile','function') else None)
        blocks.append({'node':name,'kind':node.kind,'args':node.args,'register':registers[name],'logical_gate_start':start,'logical_gate_end':c.gate_count})
    c.financial_blocks=blocks;c.node_registers=registers
    return registers[graph.output]


def compile_value(graph,tables,mode='exact7t',store_gates=True):
    if mode not in ('exact7t','unitary8t'):raise ValueError('Coherent lookup architecture required')
    c=vo.ValueCircuit(store_gates=store_gates,mode=mode);inputs={name:c.alloc(bits,name) for name,bits in graph.inputs.items()}
    out=c.alloc(graph.q+1,'normalised probability XOR output');snap=c.snapshot();value=emit_body(c,graph,tables,inputs);end=c.gate_count
    for a,b in zip(value[:graph.q+1],out):c.cnot(a,b)
    c.uncompute(snap,end_gate_count=end);c.uniform_preparation_gates=[]
    return c,inputs,out


def compile_flag(graph,tables,mode='exact7t',store_gates=True):
    if mode not in ('exact7t','unitary8t'):raise ValueError('Coherent lookup architecture required')
    c=vo.ValueCircuit(store_gates=store_gates,mode=mode);inputs={name:c.alloc(bits,name) for name,bits in graph.inputs.items()}
    comparison=c.alloc(graph.q,'uniform comparator');zero=c.alloc(1,'comparator high zero');flag=c.alloc(1,'success flag')[0]
    snap=c.snapshot();value=emit_body(c,graph,tables,inputs);end=c.gate_count
    c.compare_into(comparison+zero,value[:graph.q+1],flag);c.uncompute(snap,end_gate_count=end)
    c.uniform_preparation_gates=[('H',w) for reg in inputs.values() for w in reg]+[('H',w) for w in comparison]
    c.reflection_core=tuple(w for reg in inputs.values() for w in reg)+comparison+(flag,)
    return c,inputs,comparison,flag


def physical_sequence(c,inverse=False,include_uniform=True):
    gates=list(getattr(c,'uniform_preparation_gates',())) if include_uniform else []
    aliases={'TD':'Tdg','SD':'Sdg'}
    for op,a,b,target in c.gates:
        if op==arch.MACRO:
            for kind,wires,condition in c._macros[a].operations():
                if condition is not None:raise ValueError('Measurement instrument is not supported here')
                gates.append((aliases.get(kind,kind),*wires))
        elif op==base.X:gates.append(('X',target))
        elif op==base.CNOT:gates.append(('CX',a,target))
        elif op==base.TOFFOLI:gates.append(('CCX',a,b,target))
    gates=qg.expand_toffoli(gates)
    return qg.adjoint(gates) if inverse else gates


def counted_composition(graph,tables,cache=None,mode='exact7t'):
    """Actual emitted blocks, fresh clean inverse and complete flag charge.

    Depths are valid sums of literal expanded block depths, not inherited
    three-layer Toffoli assignments. All original inputs remain retained.
    """
    if mode not in ('exact7t','unitary8t'):raise ValueError('Coherent architecture required')
    cache={} if cache is None else cache;forward=Counter();raw_forward=Counter();blocks=[]
    retained=sum(graph.inputs.values())+graph.q+2;peak=retained;explicit_depth=conservative_depth=all_depth=0;lookups=0
    for name,node in graph.nodes.items():
        widths=tuple(graph.inputs[a] if a in graph.inputs else graph.nodes[a].output_bits for a in node.args)
        distinct=tuple(dict.fromkeys(node.args));aliases=tuple(distinct.index(a) for a in node.args)
        table=tables.get(table_key(name,node)) if node.kind in ('function','quantile') else None
        meta=node.meta.as_dict() if isinstance(node.meta,ScalarSpec) else str(node.meta)
        key=(mode if table is not None else 'shared_exact',node.kind,json.dumps(meta,sort_keys=True),node.bits,node.output_bits,graph.q,graph.width,widths,aliases)
        if table is not None:key+=(hashlib.sha256(json.dumps(table,separators=(',',':')).encode()).hexdigest(),)
        if key not in cache:
            c=vo.ValueCircuit(store_gates=False,mode=mode);distinct_registers=[c.alloc(widths[node.args.index(a)]) for a in distinct];args=[distinct_registers[k] for k in aliases];incoming=sum(map(len,distinct_registers))
            out=emit_node(c,node,args,graph.q,graph.width,table);metrics=c.resources()
            assert len(out)==node.output_bits and metrics['active_wires_at_end']-incoming==len(out)
            cache[key]=(metrics,metrics['wires']-incoming,len(out))
        metrics,extra,live=cache[key];counts=metrics['clifford_T_expanded_gate_counts']
        assert not set(counts)-{'X','H','CX','CZ','Z','S','Sdg','T','Tdg'}
        forward.update(counts);raw_forward.update(metrics['primitive_gate_counts'])
        peak=max(peak,retained+extra);retained+=live
        explicit_depth+=metrics['explicit_7T_sequence_T_depth'];conservative_depth+=metrics['conservative_5layer_CCX_T_depth_upper']
        all_depth+=metrics['explicit_Clifford_T_all_gate_dependency_depth'];lookups+=metrics['lookup_macro_invocations']
        blocks.append({'node':name,'kind':node.kind,'function':node.meta.as_dict() if isinstance(node.meta,ScalarSpec) else None,
                       'argument_widths':list(widths),'argument_alias_pattern':list(aliases),'output_bits':node.output_bits,'primitive_counts':dict(counts),
                       'logical_T_count':metrics['logical_T_count'],'explicit_T_depth_serial_upper':metrics['explicit_7T_sequence_T_depth'],
                       'conservative_5layer_T_depth_serial_upper':metrics['conservative_5layer_CCX_T_depth_upper'],
                       'lookup_invocations':metrics['lookup_macro_invocations'],'extra_qubits_upper':extra,
                       'table_sha256':key[-1] if table is not None else None})
    comp=vo.ValueCircuit(store_gates=False,mode=mode);a,b=comp.alloc(graph.q+1),comp.alloc(graph.q+1);flag=comp.alloc(1)[0]
    comp.compare_into(a,b,flag);cm=comp.resources()
    uniform=sum(graph.inputs.values())+graph.q
    complete=Counter({g:2*n for g,n in forward.items()});complete.update(cm['clifford_T_expanded_gate_counts']);complete['H']+=uniform
    peak=max(peak,retained+cm['wires']-2*(graph.q+1)-1)
    source_path=Path(__file__)
    keys={table_key(name,node) for name,node in graph.nodes.items() if node.kind in ('quantile','function')}
    return {'schema':'qbp-full-financial-oracle-v1','architecture':mode,'configuration':graph.configuration,
      'primitive_counts':dict(complete),'primitive_counts_forward_body':dict(forward),'pre_expansion_forward_counts':dict(raw_forward),
      'logical_T_count':complete['T']+complete['Tdg'],'reflection_core_qubits':sum(graph.inputs.values())+graph.q+1,
      'uniform_H':uniform,'logical_qubits_upper':peak,'active_qubits_after_clean_prep':sum(graph.inputs.values())+graph.q+2,
      'explicit_T_depth_serial_upper':2*explicit_depth+cm['explicit_7T_sequence_T_depth'],
      'conservative_5layer_T_depth_serial_upper':2*conservative_depth+cm['conservative_5layer_CCX_T_depth_upper'],
      'explicit_all_gate_depth_serial_upper':1+2*all_depth+cm['explicit_Clifford_T_all_gate_dependency_depth'],
      'lookup_invocations_forward_and_cleanup':2*lookups,'blocks':blocks,
      'comparator_primitive_counts':cm['clifford_T_expanded_gate_counts'],
      'maximum_probability_integer':1<<graph.q,'comparison_bits':graph.q,'endpoint_high_zero_bits':1,
      'integer_fractional_bits':graph.q,'integer_word_width':graph.width,
      'scalar_table_entries':sum(len(tables[k]) for k in keys),
      'scalar_table_program_bits':sum(len(tables[k])*(graph.width if k[0]=='quantile' else 2*graph.width+1) for k in keys),
      'source_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),
      'scope':'Fully emitted financial preparation and clean comparator, including uniform Hadamards. Continuous-target price accuracy and complete mean estimation are separate certified obligations.'}
