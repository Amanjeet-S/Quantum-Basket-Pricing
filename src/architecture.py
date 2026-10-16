"""Expanded phase-correct QROM architectures, not abstract four-T Toffolis.

The measured implementation uses an explicit four-T clean-target AND
isometry and corrected X-basis erasure. Each erasure branch is proportional
to the same logical operation. Lookup inversion invokes a fresh complete
instrument; measurements are never reversed. The unitary alternative uses
the literal inverse of the compute sequence at eight T per selector pair.

The Boolean outer tape records an atomic logical XOR lookup together with
its fully expandable primitive plan. This permits existing Boolean arithmetic
to use the exact logical unitary without pretending that a measurement has
a Boolean inverse. All lookup work must start and end in zero.
"""
from collections import Counter
from dataclasses import dataclass
import hashlib
import math
from typing import Sequence
import reversible as base
import encoding_resources as er

MACRO=20
MODES=('measured4t','unitary8t','exact7t')


def and_compute(a,b,t):
    # On a clean target, the RCCX-like four-T sequence followed by S†
    # gives |a,b,0> -> |a,b,ab> with no relative phase. This is not a
    # four-T full Toffoli on an arbitrary target.
    return [('H',(t,),None),('T',(t,),None),('CX',(b,t),None),
            ('TD',(t,),None),('CX',(a,t),None),('T',(t,),None),
            ('CX',(b,t),None),('TD',(t,),None),('H',(t,),None),('SD',(t,),None)]


def and_inverse(a,b,t):
    inverse={'H':'H','T':'TD','TD':'T','CX':'CX','SD':'S'}
    return [(inverse[op],wires,None) for op,wires,_ in reversed(and_compute(a,b,t))]


@dataclass(frozen=True)
class LookupPlan:
    address:tuple
    output:tuple
    work:tuple
    table:tuple
    mode:str
    enable:int|None=None

    def operations(self):
        n=len(self.address);out=self.output;words=self.table
        if self.mode not in MODES:raise ValueError(self.mode)
        if n==0:
            value=words[0] if words else 0
            for j,w in enumerate(out):
                if value>>j&1:yield ('X',(w,),None) if self.enable is None else ('CX',(self.enable,w),None)
            return
        if len(self.work)!=n+(self.enable is None):raise ValueError('QROM selector work width mismatch')
        nonzero=[0]
        for value in words:nonzero.append(nonzero[-1]+(value!=0))
        def nonempty(lo,hi):
            lo=min(lo,len(words));hi=min(hi,len(words))
            return nonzero[hi]!=nonzero[lo]
        if not nonempty(0,1<<n):return
        root=self.work[0] if self.enable is None else self.enable
        selector=self.work[1:] if self.enable is None else self.work;measurement=0
        if self.enable is None:yield 'X',(root,),None
        def visit(prefix,bit,lo,hi):
            nonlocal measurement
            if not nonempty(lo,hi):return
            if bit<0:
                value=words[lo]
                for j,w in enumerate(out):
                    if value>>j&1:yield 'CX',(prefix,w),None
                return
            child=selector[bit];control=self.address[bit];middle=(lo+hi)//2
            if self.mode=='exact7t':yield 'CCX',(prefix,control,child),None
            else:yield from and_compute(prefix,control,child)
            yield from visit(child,bit-1,middle,hi)
            yield 'CX',(prefix,child),None
            yield from visit(child,bit-1,lo,middle)
            yield 'CX',(prefix,child),None
            if self.mode=='exact7t':yield 'CCX',(prefix,control,child),None
            elif self.mode=='unitary8t':yield from and_inverse(prefix,control,child)
            else:
                label=measurement;measurement+=1
                yield 'H',(child,),None
                yield 'MZ',(child,),label
                yield 'CZ_IF',(prefix,control),label
                yield 'X_IF',(child,),label
        yield from visit(root,n-1,0,1<<n)
        if self.enable is None:yield 'X',(root,),None

    def summary(self):
        counts=Counter(op for op,_,_ in self.operations())
        return {'primitive_gate_counts':dict(counts),'T_count':counts['T']+counts['TD']+7*counts['CCX'],
            'measurement_count':counts['MZ'],'classical_feedback_CZ_count':counts['CZ_IF'],
            'classical_reset_X_count':counts['X_IF'],
            'selector_AND_computations':counts['SD'],
            'selector_unitary_erases':counts['S'],
            'physical_operations':sum(counts.values()),
            'logical_inverse':'fresh identical XOR instrument' if self.mode=='measured4t' else 'clean-subspace XOR lookup is its own logical inverse'}


def lookup_plan(address,output,table,mode='measured4t',work=None,enable=None):
    address=tuple(address);output=tuple(output);table=tuple(table)
    if mode not in MODES or len(table)>1<<len(address):raise ValueError('Invalid QROM mode or table size')
    if any(not isinstance(v,int) or v<0 or v>=1<<len(output) for v in table):raise ValueError('QROM word outside output width')
    if work is None:
        start=max(address+output+(() if enable is None else (enable,)),default=-1)+1
        work=tuple(range(start,start+len(address)+(enable is None))) if address else ()
    work=tuple(work)
    wires=address+output+work+(() if enable is None else (enable,))
    if len(set(wires))!=len(wires):raise ValueError('Address, dirty output and clean work must be disjoint')
    return LookupPlan(address,output,work,table,mode,enable)


class ArchitectureCircuit(base.Circuit):
    """Boolean composition with explicit, phase-audited QROM instruments."""
    def __init__(self,store_gates=True,mode='measured4t'):
        if mode not in MODES:raise ValueError(mode)
        super().__init__(store_gates);self.mode=mode;self._macros=[]
        self._physical=Counter();self._tdepth=[];self._measurement_depth=[]
        self._lookup_measurements=0

    def _ensure_depths(self):
        while len(self._tdepth)<self.wire_count:self._tdepth.append(0);self._measurement_depth.append(0)

    def _physical_operation(self,op,wires,condition,classical):
        self._ensure_depths();self._physical[op]+=1
        if not set(wires)<=self._active:raise ValueError('Physical lookup uses inactive wires')
        layer=1+max((self._depth[w] for w in wires),default=0)
        td=max((self._tdepth[w] for w in wires),default=0)
        md=max((self._measurement_depth[w] for w in wires),default=0)
        if condition is not None and op!='MZ':
            cl,ct,cm=classical[condition];layer=max(layer,cl+1);td=max(td,ct);md=max(md,cm)
        if op in ('T','TD'):td+=1
        elif op=='CCX':td+=3
        if op=='MZ':
            md+=1;classical[condition]=(layer,td,md);self._lookup_measurements+=1
        for w in wires:self._depth[w]=layer;self._tdepth[w]=td;self._measurement_depth[w]=md

    def _emit(self,op,c1,c2,target):
        if op!=MACRO:
            super()._emit(op,c1,c2,target);self._ensure_depths()
            touched=(target,) if op==base.X else ((c1,target) if op==base.CNOT else (c1,c2,target))
            td=max(self._tdepth[w] for w in touched)+(3 if op==base.TOFFOLI else 0)
            md=max(self._measurement_depth[w] for w in touched)
            for w in touched:self._tdepth[w]=td;self._measurement_depth[w]=md
            return
        plan=self._macros[c1]
        for w in plan.address+plan.output+plan.work+(() if plan.enable is None else (plan.enable,)):
            if w not in self._active:
                self._active.add(w);self._free.remove(w)
        self._max_active=max(self._max_active,len(self._active))
        self._counts[MACRO]+=1
        if self.store_gates or self._snapshots:self._tape.extend((MACRO,c1,-1,-1))
        classical={}
        # This iterates every actual gate and conditional operation. The
        # resource model does not substitute a four-T integer multiplier.
        for gate,wires,condition in plan.operations():self._physical_operation(gate,wires,condition,classical)

    def lookup(self,address,table,out_bits,enable=None):
        out=self.alloc(out_bits,'QROM output')
        if not any(table):return out
        work=self.alloc(len(address)+(enable is None),'QROM clean selectors') if address else ()
        plan=lookup_plan(address,out,table,self.mode,work,enable)
        index=len(self._macros);self._macros.append(plan)
        self._emit(MACRO,index,-1,-1)
        self.release(work)
        return out

    def _reverse_segment(self,start,end):
        for i in range(end-4,start-1,-4):
            op,c1,c2,target=self._tape[i:i+4]
            if op==MACRO:
                plan=self._macros[c1];touched=plan.address+plan.output+plan.work+(() if plan.enable is None else (plan.enable,))
            else:touched=(target,) if op==base.X else ((c1,target) if op==base.CNOT else (c1,c2,target))
            for w in touched:
                if w not in self._active:self._active.add(w);self._free.remove(w)
            self._max_active=max(self._max_active,len(self._active))
            # XOR lookup's inverse is a fresh valid instrument. In
            # measured4t mode this is explicitly not measurement reversal.
            self._emit(op,c1,c2,target)

    def run(self,initial=None):
        if not self.store_gates:raise RuntimeError('Basis replay requires a stored logical tape')
        state=[0]*self.wire_count;assigned=set()
        for reg,value in (initial or {}).items():
            if value<0 or value>=1<<len(reg) or assigned.intersection(reg):raise ValueError('Invalid initial registers')
            assigned.update(reg)
            for i,w in enumerate(reg):state[w]=value>>i&1
        for op,c1,c2,target in self.gates:
            if op==MACRO:
                plan=self._macros[c1]
                assert all(state[w]==0 for w in plan.work),'QROM private work was not clean'
                index=self.value(state,plan.address);word=plan.table[index] if index<len(plan.table) else 0
                if plan.enable is not None and not state[plan.enable]:word=0
                for j,w in enumerate(plan.output):state[w]^=word>>j&1
            elif op==base.X or op==base.CNOT and state[c1] or op==base.TOFFOLI and state[c1] and state[c2]:state[target]^=1
        return state

    def resources(self):
        self._ensure_depths();p=self._physical.copy()
        p['X']+=self._counts[base.X];p['CX']+=self._counts[base.CNOT];p['CCX']+=self._counts[base.TOFFOLI]
        td=p['T']+p['TD']+7*p['CCX']
        internal=(self._physical['CCX']//2 if self.mode=='exact7t' else self._physical['SD'])
        return {'x':p['X'],'cnot':p['CX'],'toffoli':p['CCX'],'primitive_gate_counts':dict(p),
            'logical_T_count':td,'expanded_gates':sum(p.values()),'logical_tape_tokens':self.gate_count,
            'wires':self.wire_count,'peak_active_wires':self._max_active,'active_wires_at_end':len(self._active),
            'all_gate_dependency_depth':max(self._depth,default=0),
            'logical_T_layer_depth_upper':max(self._tdepth,default=0),
            'measurement_layer_depth_upper':max(self._measurement_depth,default=0),
            'measurement_count':p['MZ'],'classical_feedback_CZ_count':p['CZ_IF'],'classical_reset_X_count':p['X_IF'],
            'lookup_macro_invocations':self._counts[MACRO],'selector_AND_computations':p['SD'],
            'lookup_selector_internal_nodes':internal,
            'arithmetic_exact_Toffoli_count':self._counts[base.TOFFOLI],
            'selector_T_count':{'exact7t':14,'unitary8t':8,'measured4t':4}[self.mode]*internal,
            'architecture':self.mode,'counting':'all primitive lookup gates emitted; conditional Clifford operations charged at worst-case execution'}


def compile_graph(graph,tables,mode='measured4t'):
    c=ArchitectureCircuit(mode=mode);inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    comparison=c.alloc(graph.q,'comparison');zero=c.alloc(1,'comparison high zero');flag=c.alloc(1,'good flag')[0]
    snap=c.snapshot();regs=dict(inputs)
    for name,node in graph.nodes.items():
        key=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        regs[name]=er.emit_node(c,node,[regs[a] for a in node.args],graph.q,graph.width,tables.get(key))
    end=c.gate_count;c.compare_into(comparison+zero,regs[graph.output][:graph.q+1],flag);c.uncompute(snap,end_gate_count=end)
    return c,inputs,comparison,flag


def counted_composition(graph,tables,cache,mode='measured4t'):
    counts=Counter();retained=sum(graph.inputs.values())+graph.q+2;peak=retained
    depth=tdepth=mdepth=0;internal=arithmetic=0;blocks=[]
    for name,node in graph.nodes.items():
        keytable=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        table=tables.get(keytable)
        # Pure arithmetic blocks have identical emitted gates in all modes;
        # only lookup-containing blocks require an architecture-specific key.
        cache_mode=mode if node.kind in ('quantile','function') else 'shared_exact_arithmetic'
        key=(cache_mode,node.kind,node.meta,node.bits,graph.q,graph.width,tuple(graph.inputs.get(a,graph.width) for a in node.args))
        if table is not None:key+=(hashlib.sha256(str(table).encode()).hexdigest(),)
        if key not in cache:
            c=ArchitectureCircuit(store_gates=False,mode=mode);args=[c.alloc(graph.inputs.get(a,graph.width)) for a in node.args]
            incoming=sum(map(len,args));er.emit_node(c,node,args,graph.q,graph.width,table)
            metrics=c.resources();cache[key]=(metrics,c.wire_count-incoming,metrics['active_wires_at_end']-incoming)
        metrics,extra,live=cache[key];counts.update(metrics['primitive_gate_counts'])
        peak=max(peak,retained+extra);retained+=live
        depth+=metrics['all_gate_dependency_depth'];tdepth+=metrics['logical_T_layer_depth_upper'];mdepth+=metrics['measurement_layer_depth_upper']
        internal+=metrics['lookup_selector_internal_nodes'];arithmetic+=metrics['arithmetic_exact_Toffoli_count']
        blocks.append({'node':name,'kind':node.kind,'function':node.meta if node.kind=='function' else None,
            'address_bits':node.bits,'primitive_gate_counts':metrics['primitive_gate_counts'],
            'logical_T_count':metrics['logical_T_count'],'measurements':metrics['measurement_count']})
    comp=ArchitectureCircuit(store_gates=False,mode=mode);a,b=comp.alloc(graph.q+1),comp.alloc(graph.q+1)
    flag=comp.alloc(1)[0];comp.compare_into(a,b,flag);cm=comp.resources()
    counts=Counter({k:2*v for k,v in counts.items()});counts.update(cm['primitive_gate_counts'])
    peak=max(peak,retained+comp.wire_count-(2*(graph.q+1)+1))
    tcount=counts['T']+counts['TD']+7*counts['CCX'];reflection=er.reflection_resources(sum(graph.inputs.values())+graph.q+1)
    internal*=2;arithmetic=2*arithmetic+cm['arithmetic_exact_Toffoli_count']
    assert tcount=={'exact7t':14,'unitary8t':8,'measured4t':4}[mode]*internal+7*arithmetic
    keys={('quantile',n.bits,graph.q) if n.kind=='quantile' else (n.meta,n.bits,graph.q)
          for n in graph.nodes.values() if n.kind in ('quantile','function')}
    return {'architecture':mode,'expanded_gate_counts':dict(counts),'logical_T_count':tcount,
        'logical_T_depth_serial_upper':2*tdepth+cm['logical_T_layer_depth_upper'],
        'serial_all_gate_depth_upper':2*depth+cm['all_gate_dependency_depth'],
        'measurement_count':counts['MZ'],'measurement_depth_serial_upper':2*mdepth+cm['measurement_layer_depth_upper'],
        'classical_feedback_CZ_count':counts['CZ_IF'],'classical_reset_X_count':counts['X_IF'],
        'lookup_selector_internal_nodes':internal,'arithmetic_exact_Toffoli_count':arithmetic,
        'logical_qubits_serial_upper_bound':peak,'input_qubits':sum(graph.inputs.values()),'comparison_qubits':graph.q,
        'includes_complete_arithmetic_and_lookup_logical_inverse':True,
        'lookup_inverse':'fresh instrument with independent outcomes' if mode=='measured4t' else 'clean-subspace logical XOR inverse',
        'includes_finite_uniform_preparation_H_gates':sum(graph.inputs.values())+graph.q,
        'includes_mean_estimator':False,'grover_iterate_T_count':2*tcount+reflection['T_count_exact_7T'],
        'grover_iterate_measurement_count':2*counts['MZ'],
        'grover_iterate_T_depth_serial_upper':2*(2*tdepth+cm['logical_T_layer_depth_upper'])+reflection['T_depth_3layer_upper'],
        'grover_iterate_measurement_depth_serial_upper':2*(2*mdepth+cm['measurement_layer_depth_upper']),
        'grover_iterate_logical_qubits_upper':max(peak,reflection['logical_qubits']),
        'grover_zero_reflection':reflection,
        'classical_table_entries':sum(len(tables[k]) for k in keys),
        'classical_table_program_bits':sum(len(tables[k])*(graph.width if k[0]=='quantile' else 2*graph.width) for k in keys),
        'schedule':'sum of emitted forward blocks, valid logical inverses and comparator; measurement/feedforward counted separately',
        'blocks':blocks}
