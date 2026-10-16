"""Executable finite fixed-point encodings and conservative logical resources.

Shared kernel library, distinct arithmetic precisions and table sizes. Full
counted circuits are upper bounds for this constructor, not optimal circuits.
The quantum mean estimator is deliberately outside these evaluator counts.
"""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reversible import Circuit
from reflections import resources as reflection_resources

ROOT = Path(__file__).resolve().parent
D = Decimal
NORMAL = NormalDist()


def decimal_pi():
    with localcontext() as ctx:
        ctx.prec = 70
        def atan(x):
            total, term, n = x, x, 0
            while abs(term) > D('1e-68'):
                n += 1
                term *= -x*x
                total += term / (2*n+1)
            return total
        return +(16*atan(D(1)/5)-4*atan(D(1)/239))

PI = decimal_pi()
with localcontext() as ctx:
    ctx.prec = 70
    INV_SQRT_2PI = 1/(2*PI).sqrt()


def decimal_cdf(x):
    """Normal integral series, with an explicit generous absolute enclosure.

    For |x|<=8 the first omitted alternating term is <1e-50, after
    the terms become decreasing.  At 70 decimal digits the total absolute
    rounding allowance 1e-38 exceeds the operation-count bound (see report).
    This routine returns the centre; certificate comparisons use +/-1e-38.
    """
    x = D(x)
    if abs(x) > 8:
        raise ValueError('CDF series domain exceeded')
    with localcontext() as ctx:
        ctx.prec = 70
        term, total, n = x, x, 0
        while n < int(x*x/2)+2 or abs(term) > D('1e-50'):
            n += 1
            term *= -x*x*D(2*n-1)/(2*n*(2*n+1))
            total += term
        return +(D('.5')+INV_SQRT_2PI*total)


def rounded_scaled(value, q):
    with localcontext() as ctx:
        ctx.prec = 70
        return int((value * (1 << q)).to_integral_value())


def function_span(kind):
    if kind.startswith('general_'):
        from general_kernels import function_span as general_span
        return general_span(kind)
    if kind.startswith('strike_log_'): return 1
    if kind.startswith('capped_') or kind=='scaled_cdf': return 2
    return 4 if kind in ('exp','log') else 16


def capped_kernel_decimal(m, cap):
    """Known scalar conditional capped payoff, normalised by the cap.

    This is not a financial-input lookup: basket log bases are computed by
    the graph. The cap-specific table's whole programme and gates are costed.
    On m in [-1,1], all formula CDF arguments remain inside [-8,8].
    """
    with localcontext() as ctx:
        ctx.prec=70
        h=D(cap)/100
        t0=5*m; t1=t0-5*(1+h).ln()
        p0,p1,a0,a1=[decimal_cdf(t) for t in (t0,t1,t0+D('.2'),t1+D('.2'))]
        return +((m+D('.02')).exp()*(a0-a1)-(p0-p1)+h*p1)/h


def function_table(kind, bits, q):
    """Packed endpoint and endpoint difference, with dyadic input domains."""
    n = 1 << bits
    with localcontext() as ctx:
        ctx.prec = 70
        values = []
        for j in range(n+1):
            if kind == 'exp':
                x = D(-2)+D(4)*j/n
                value = x.exp()
            elif kind == 'log':
                x = D('.25')+D(4)*j/n
                value = x.ln()+2
            elif kind == 'cdf':
                x = D(-8)+D(16)*j/n
                value = decimal_cdf(x)
            elif kind.startswith('capped_'):
                x = D(-1)+D(2)*j/n
                value = capped_kernel_decimal(x, int(kind.split('_')[1]))
            elif kind == 'scaled_cdf':
                x = D(-1)+D(2)*j/n
                value = decimal_cdf(5*x)
            elif kind.startswith('strike_log_'):
                h=D(int(kind.split('_')[-1]))/100
                value=1+(1+h*D(j)/n).ln()
            else:
                raise ValueError(kind)
            values.append(rounded_scaled(value, q))
    width = q+6
    assert all(0 <= v < 1 << width for v in values)
    assert all(a <= b for a,b in zip(values, values[1:]))
    return [values[j] | ((values[j+1]-values[j]) << width) for j in range(n)]


def quantile_table(bits, q, certify=False):
    """Uniform midpoint quantile register Y+5; never a joint payoff table.

    Double precision supplies a seed only. The optional exhaustive check
    encloses the true inverse by the stored value +/- one fixed-point unit.
    """
    n = 1 << bits
    values, largest_residual = [], D(0)
    for j in range(n):
        p = (j+.5)/n
        if j < n//2:
            positive = n-j-1
            value = 10*(1 << q)-round((NORMAL.inv_cdf((positive+.5)/n)+5)*(1 << q))
        else:
            value = round((NORMAL.inv_cdf(p)+5)*(1 << q))
        if certify and j >= n//2:
            with localcontext() as ctx:
                ctx.prec = 70
                centre = D(value)/D(1 << q)-5
                exact_p = (D(j)+D('.5'))/n
                lo = decimal_cdf(centre-D(1)/D(1 << q))+D('1e-38')
                hi = decimal_cdf(centre+D(1)/D(1 << q))-D('1e-38')
                if not lo <= exact_p <= hi:
                    raise AssertionError(('inverse-normal enclosure', bits, q, j))
                largest_residual = max(largest_residual, abs(lo-exact_p),abs(hi-exact_p))
        values.append(value)
    # Rounded binary seeds are exactly symmetric, including midpoint ties.
    assert all(values[j]+values[-j-1] == 10*(1 << q) for j in range(n))
    return values, float(largest_residual)


@dataclass
class Node:
    kind: str
    args: tuple[str, ...]
    meta: object = None
    bits: int = 0


class Graph:
    def __init__(self, q):
        self.q, self.width = q, q+6
        self.nodes: dict[str, Node] = {}
        self.inputs: dict[str, int] = {}
        self.output = ''
        self.constants = {}
        self.ideal_constants = {}

    def inp(self, name, width):
        self.inputs[name] = width
        return name

    def const(self, value):
        value = D(str(value))
        integer = rounded_scaled(value, self.q) % (1 << self.width)
        if integer in self.constants:
            return self.constants[integer]
        name = f'c{len(self.constants)}'
        self.constants[integer] = name
        self.nodes[name] = Node('constant', (), integer)
        self.ideal_constants[name] = value
        return name

    def op(self, kind, *args, meta=None, bits=0):
        name = f'n{len(self.nodes)}'
        self.nodes[name] = Node(kind, args, meta, bits)
        return name

    def add(self,a,b): return self.op('add',a,b)
    def sub(self,a,b): return self.op('sub',a,b)
    def positive_sub(self,a,b): return self.op('positive_sub',a,b)
    def mul(self,a,b): return self.op('mul',a,b)
    def scale(self,a,k):
        ideal=D(str(k))
        if ideal == int(ideal):
            integer,drop=int(ideal),0
        else:
            integer,drop=rounded_scaled(ideal,self.q),self.q
        return self.op('scale',a,meta=(integer,drop,str(ideal)))
    def shift(self,a,bits): return self.op('shift',a,meta=bits)
    def clamp(self,a,cap=1): return self.op('clamp',a,meta=rounded_scaled(D(str(cap)),self.q))
    def function(self,kind,a,bits): return self.op('function',a,meta=kind,bits=bits)

    def affine(self,y,z,which):
        # Exact real model: r1=-.1, r2=.2, lambda=.2; y,z are shifted by5.
        if which == 0:
            out = self.sub(self.const('2.485'), self.scale(y,'.1'))
        elif which == 1:
            out = self.add(self.const('.97'), self.scale(y,'.2'))
        else:
            out = self.add(self.const('1.7275'), self.scale(y,'.05'))
        if z:
            out = self.add(self.sub(out,self.const('1')),self.scale(z,'.2'))
        return out

    def payoff(self,a,cap):
        return self.clamp(self.sub(a,self.const(1)),cap)

    def evaluate(self, inputs, tables):
        vals = dict(inputs)
        mask, unit = (1 << self.width)-1, 1 << self.q
        for name,node in self.nodes.items():
            a = [vals[v] for v in node.args]
            if node.kind == 'constant': v=node.meta
            elif node.kind == 'add': v=a[0]+a[1]
            elif node.kind in ('sub','positive_sub'): v=a[0]-a[1]
            elif node.kind == 'mul': v=(a[0]*a[1]) >> self.q
            elif node.kind == 'scale': v=(a[0]*node.meta[0]) >> node.meta[1]
            elif node.kind == 'shift': v=a[0] >> node.meta
            elif node.kind == 'quantile': v=tables[(node.kind,node.bits,self.q)][a[0]]
            elif node.kind == 'uniform':
                v=(2*a[0]+1)*(unit >> (node.bits+1))
            elif node.kind == 'function':
                span = function_span(node.meta)
                fractional = self.q+int(math.log2(span))-node.bits
                index, rem = a[0] >> fractional, a[0] & ((1 << fractional)-1)
                packed=tables[(node.meta,node.bits,self.q)][index]
                base,delta=packed & mask,packed >> self.width
                v=base+((delta*rem) >> fractional)
            elif node.kind == 'clamp':
                signed=a[0] if a[0] < 1 << (self.width-1) else a[0]-(1 << self.width)
                v=min(max(signed,0),node.meta)
            else: raise ValueError(node.kind)
            vals[name]=v & mask
        return vals[self.output]


def build_graph(kind,q,b_y,b_z,b_u,exp_bits,log_bits,cdf_bits,cap=10):
    if q < max(b_u+1,exp_bits,log_bits,cdf_bits):
        raise ValueError('Insufficient fixed-point fractional bits')
    g=Graph(q)
    iy=g.inp('iy',b_y)
    y=g.op('quantile',iy,bits=b_y)
    h=D(cap)/100
    z=None
    if kind=='direct':
        iz=g.inp('iz',b_z)
        z=g.op('quantile',iz,bits=b_z)
        prices=[g.function('exp',g.affine(y,z,j),exp_bits) for j in (0,1)]
        ba=g.shift(g.add(*prices),1)
        bg=g.function('exp',g.affine(y,z,2),exp_bits)
        res=g.clamp(g.sub(g.payoff(ba,h),g.payoff(bg,h)),h)
        g.output=g.clamp(g.scale(res,D(1)/h),1)
        return g
    logs=[g.affine(y,None,j) for j in (0,1)]
    prices=[g.function('exp',v,exp_bits) for v in logs]
    ba=g.shift(g.add(*prices),1)
    ma=g.function('log',g.sub(ba,g.const('.25')),log_bits)
    mg=g.shift(g.add(*logs),1)
    if kind=='interval_kernel':
        iu=g.inp('iu',b_u)
        u=g.op('uniform',iu,bits=b_u)
        mk=g.function(f'strike_log_{cap}',u,log_bits)
        probs=[g.function('scaled_cdf',g.sub(m,mk),cdf_bits) for m in (ma,mg)]
        g.output=g.clamp(g.sub(*probs),1)
    elif kind=='interval':
        iu=g.inp('iu',b_u)
        u=g.op('uniform',iu,bits=b_u)
        k=g.add(g.const(1),g.scale(u,h))
        mk=g.function('log',g.sub(k,g.const('.25')),log_bits)
        probs=[]
        for m in (ma,mg):
            t=g.scale(g.sub(g.add(m,g.const('1.6')),mk),5)
            probs.append(g.function('cdf',t,cdf_bits))
        g.output=g.clamp(g.sub(*probs),1)
    elif kind=='full_kernel':
        # A specialised one-dimensional known function replaces the eight
        # CDF calls and products. This computes the exact same full conditional
        # expectation before rounding, without a direct input-payoff table.
        probs=[g.function(f'capped_{cap}',g.sub(m,g.const(1)),cdf_bits)
               for m in (ma,mg)]
        g.output=g.clamp(g.sub(*probs),1)
    elif kind=='full':
        with localcontext() as ctx:
            ctx.prec=70
            delta=((1+h).ln())*5
            growth=D('.02').exp()
        outputs=[]
        for m,b in ((ma,ba),(mg,g.function('exp',mg,exp_bits))):
            t0=g.scale(g.sub(m,g.const('.4')),5) #8+5*log(B)
            t1=g.sub(t0,g.const(delta))
            p0=g.function('cdf',t0,cdf_bits)
            p1=g.function('cdf',t1,cdf_bits)
            a0=g.function('cdf',g.add(t0,g.const('.2')),cdf_bits)
            a1=g.function('cdf',g.add(t1,g.const('.2')),cdf_bits)
            term1=g.mul(g.scale(b,growth),g.positive_sub(a0,a1))
            term2=g.positive_sub(p0,p1)
            term3=g.scale(p1,h)
            outputs.append(g.clamp(g.add(g.sub(term1,term2),term3),h))
        res=g.clamp(g.sub(*outputs),h)
        g.output=g.clamp(g.scale(res,D(1)/h),1)
    else: raise ValueError(kind)
    return g


def make_tables(graph,cache,certify_quantiles=False):
    audits=[]
    for node in graph.nodes.values():
        if node.kind=='quantile':
            key=('quantile',node.bits,graph.q)
            if key not in cache:
                cache[key],residual=quantile_table(node.bits,graph.q,certify_quantiles)
                audits.append({'bits':node.bits,'q':graph.q,'entries':1 << node.bits,
                    'exhaustive_inverse_enclosure':certify_quantiles,'largest_cdf_residual':residual})
        elif node.kind=='function':
            key=(node.meta,node.bits,graph.q)
            if key not in cache: cache[key]=function_table(*key)
    return audits


def emit_node(c,node,args,q,width,table):
    if node.kind=='constant':
        out=c.alloc(width)
        for i,w in enumerate(out):
            if (node.meta >> i)&1:c.x(w)
        return out
    if node.kind=='add': return c.unsigned_add(*args)
    if node.kind in ('sub','positive_sub'): return c.unsigned_subtract(*args)
    if node.kind=='mul': return c.multiply(*args,output_bits=width,drop_lsb=q)
    if node.kind=='scale': return c.multiply_constant(args[0],node.meta[0],output_bits=width,drop_lsb=node.meta[1])
    if node.kind=='shift':
        # Copy logical shifted values, including explicitly reset top bits.
        out=c.alloc(width)
        for a,b in zip(args[0][node.meta:],out):c.cnot(a,b)
        return out
    if node.kind=='quantile': return c.lookup(args[0],table,width)
    if node.kind=='uniform':
        out=c.alloc(width)
        start=q-node.bits
        c.x(out[start-1])
        for i,w in enumerate(args[0]):c.cnot(w,out[start+i])
        return out
    if node.kind=='function':
        span=function_span(node.meta)
        fractional=q+int(math.log2(span))-node.bits
        address=args[0][fractional:fractional+node.bits]
        packed=c.lookup(address,table,2*width)
        base,delta=packed[:width],packed[width:]
        product=c.multiply(delta,args[0][:fractional],output_bits=width,drop_lsb=fractional)
        return c.unsigned_add(base,product)
    if node.kind=='clamp':
        upper=c.alloc(width)
        for i,w in enumerate(upper):
            if(node.meta >> i)&1:c.x(w)
        flag=c.unsigned_less_than(args[0],upper)[0]
        out=c.alloc(width)
        # Signed lower clipping; overflow-safe because all real intermediates
        # are strictly smaller than half the unsigned register range.
        positive=c.alloc(1)[0]; c.x(positive);c.cnot(args[0][-1],positive)
        take_value=c.alloc(1)[0];c.toffoli(positive,flag,take_value)
        take_cap=c.alloc(1)[0];c.x(flag);c.toffoli(positive,flag,take_cap);c.x(flag)
        for i,w in enumerate(out):
            c.toffoli(take_value,args[0][i],w)
            if(node.meta >> i)&1:c.cnot(take_cap,w)
        return out
    raise ValueError(node.kind)


def compile_small(graph,tables):
    c=Circuit()
    inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    # Uniform comparator has q+1 bits so it represents values below1 exactly.
    comparison=c.alloc(graph.q,'comparison integer')
    zero=c.alloc(1,'comparison high zero')
    flag=c.alloc(1,'good flag')[0]
    snap=c.snapshot()
    regs=dict(inputs)
    for name,node in graph.nodes.items():
        key=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        regs[name]=emit_node(c,node,[regs[a] for a in node.args],graph.q,graph.width,tables.get(key))
    body_end=c.gate_count
    c.compare_into(comparison+zero,regs[graph.output][:graph.q+1],flag)
    c.uncompute(snap,end_gate_count=body_end)
    return c,inputs,comparison,flag


def counted_composition(graph,tables,cache):
    """Emit each instantiated block, then compose a conservative serial schedule.

    Large complete tapes are unnecessary: every block expands actual gates.
    The resource total contains every forward block and its exact inverse,
    plus the clean output comparator. Depth is an executable serial upper
    bound, not the global ASAP depth of a materialised complete circuit.
    """
    counts={'X':0,'CNOT':0,'Toffoli':0}
    retained=sum(graph.inputs.values())+graph.q+2
    peak=retained
    serial_depth=0
    serial_toffoli_depth=0
    blocks=[]
    for name,node in graph.nodes.items():
        keytable=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        table=tables.get(keytable)
        key=(node.kind,node.meta,node.bits,graph.q,tuple(graph.inputs.get(a,graph.width) for a in node.args))
        if node.kind in ('function','quantile'):key+=(hashlib.sha256(str(table).encode()).hexdigest(),)
        if key not in cache:
            c=Circuit(store_gates=False)
            args=[c.alloc(graph.inputs.get(a,graph.width)) for a in node.args]
            start_inputs=sum(map(len,args))
            out=emit_node(c,node,args,graph.q,graph.width,table)
            metrics=c.resources()
            cache[key]=(metrics,c.wire_count-start_inputs,metrics['active_wires_at_end']-start_inputs)
        metrics,extra,live=cache[key]
        for gate in counts:counts[gate]+=metrics[{'X':'x','CNOT':'cnot','Toffoli':'toffoli'}[gate]]
        peak=max(peak,retained+extra)
        # Keep output; primitive work is zero or can be recovered by its
        # local inverse. Non-clean function/clamp blocks keep all their
        # intermediate outputs until the complete body is reversed.
        retained+=live
        serial_depth+=metrics['all_gate_dependency_depth']
        serial_toffoli_depth+=metrics['toffoli_layer_depth']
        blocks.append({'node':name,'kind':node.kind,'function':node.meta if node.kind=='function' else None,
                       'address_bits':node.bits,'counts':{k:metrics[v] for k,v in [('X','x'),('CNOT','cnot'),('Toffoli','toffoli')]},'wires':metrics['wires']})
    # Comparison is clean; its inputs are pre-existing and flag is reserved.
    comp=Circuit(store_gates=False)
    a,b=comp.alloc(graph.q+1),comp.alloc(graph.q+1)
    f=comp.alloc(1)[0]
    comp.compare_into(a,b,f)
    cm=comp.resources()
    for gate in counts:counts[gate]=2*counts[gate]+cm[{'X':'x','CNOT':'cnot','Toffoli':'toffoli'}[gate]]
    peak=max(peak,retained+comp.wire_count-(2*(graph.q+1)+1))
    td=2*serial_toffoli_depth+cm['toffoli_layer_depth']
    coherent_inputs=sum(graph.inputs.values())+graph.q
    reflection_wires=coherent_inputs+1
    reflection=reflection_resources(reflection_wires)
    reflection_toffolis=reflection["gate_counts"]["CCX"]
    table_keys={('quantile',n.bits,graph.q) if n.kind=='quantile' else (n.meta,n.bits,graph.q)
                for n in graph.nodes.values() if n.kind in ('quantile','function')}
    table_storage_bits=sum(len(tables[k])*(graph.width if k[0]=='quantile' else 2*graph.width) for k in table_keys)
    return {'expanded_gate_counts':counts,'logical_T_count_exact_7T_Toffoli':7*counts['Toffoli'],
            'logical_T_depth_serial_3_layers_per_Toffoli':3*td,
            'serial_all_gate_depth':2*serial_depth+cm['all_gate_dependency_depth'],
            'logical_qubits_serial_upper_bound':peak,
            'input_qubits':sum(graph.inputs.values()),'comparison_qubits':graph.q,
            'includes_arithmetic_inverse':True,'includes_finite_uniform_preparation_H_gates':sum(graph.inputs.values())+graph.q,
            'includes_mean_estimator':False,
            'classical_table_entries':sum(len(tables[k]) for k in table_keys),
            'classical_table_program_bits':table_storage_bits,
            'state_preparation_depth_including_uniform_H_layer_upper':2*serial_depth+cm['all_gate_dependency_depth']+1,
            'grover_iterate_T_count_exact_7T_Toffoli':14*counts['Toffoli']+7*reflection_toffolis,
            'grover_iterate_T_depth_serial_upper':6*td+reflection['T_depth_3layer_upper'],
            'grover_iterate_logical_qubits_upper':max(peak,reflection['logical_qubits']),
            'grover_zero_reflection':reflection,
            'grover_scope':'two clean flag evaluators, two uniform Hadamard layers, all-zero reflection on input/comparison/flag invariant subspace and good-flag Z; no mean-estimator repetition/confidence schedule',
            'schedule':'sum of fully emitted blocks and exact inverse; serial upper bound',
            'blocks':blocks}


def ideal_value(kind,iy,iz,iu,by,bz,bu,cap):
    y=NORMAL.inv_cdf((iy+.5)/(1 << by))
    ma1,ma2=-.015-.1*y,-.03+.2*y
    ma=math.log((math.exp(ma1)+math.exp(ma2))/2)
    mg=(ma1+ma2)/2
    cdf=lambda x:.5*math.erfc(-x/math.sqrt(2))
    h=cap/100
    if kind=='direct':
        z=NORMAL.inv_cdf((iz+.5)/(1 << bz))
        f=lambda x:min(max(math.exp(x+.2*z)-1,0),h)
        return min(max((f(ma)-f(mg))/h,0),1)
    if kind in ('interval','interval_kernel'):
        u=(iu+.5)/(1 << bu)
        mk=math.log(1+h*u)
        return min(max(cdf(5*(ma-mk))-cdf(5*(mg-mk)),0),1)
    delta=5*math.log(1+h)
    def g(m):
        t0=5*m;t1=t0-delta
        return math.exp(m+.02)*(cdf(t0+.2)-cdf(t1+.2))-(cdf(t0)-cdf(t1))+h*cdf(t1)
    return min(max((g(ma)-g(mg))/h,0),1)


def check_small():
    results=[]; tables={}
    # Coarse precisions intentionally separate from the monetary configurations.
    for kind in ('direct','interval','full','interval_kernel','full_kernel'):
        graph=build_graph(kind,8,2,2,2,3,3,4,10)
        make_tables(graph,tables)
        c,inputs,comparison,flag=compile_small(graph,tables)
        max_error=0
        # Exhaustive model inputs and boundary comparator levels around output.
        from itertools import product
        tested=0
        for indices in product(*[range(1 << b) for b in graph.inputs.values()]):
            mapping=dict(zip(graph.inputs,indices))
            threshold=graph.evaluate(mapping,tables)
            max_error=max(max_error,abs(threshold/(1 << graph.q)-ideal_value(kind,mapping['iy'],mapping.get('iz',0),mapping.get('iu',0),2,2,2,10)))
            for level in sorted(set([0,min(threshold,(1 << graph.q)-1),max(0,threshold-1),(1 << graph.q)-1])):
                values={inputs[n]:v for n,v in mapping.items()};values[comparison]=level
                state=c.run(values)
                assert state[flag]==int(level<threshold),(kind,mapping,level,threshold)
                for n,v in mapping.items():assert c.value(state,inputs[n])==v
                assert c.value(state,comparison)==level
                reserved={flag,*comparison,*[w for reg in inputs.values() for w in reg]}
                assert all(v==0 for i,v in enumerate(state) if i not in reserved)
                tested+=1
        results.append({'kind':kind,'basis_replay_checks':tested,'q':graph.q,'gaussian_bits':2,
                        'max_probability_error_diagnostic':max_error,'gate_level_resources':c.resources(),
                        'all_work_ancillas_reset':True,'all_inputs_preserved':True})
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check-small',action='store_true')
    args=p.parse_args()
    if args.check_small:
        print(json.dumps(check_small(),indent=2))
