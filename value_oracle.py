"""Exact finite-value and dyadic-bucket quantum preparations.

The financial arithmetic is the existing checked integer graph. The logical
value oracle is an XOR into a q+1-bit dirty output. The bucket flag uses only
exact Boolean gates and uniform Hadamards; no amplitude rotation is assumed.
All lookup macros have fully emitted phase-correct implementations. Their
logical inverse is a fresh corrected instrument in measured4t mode.

The finite dyadic identity is established mathematics, not a novelty claim.
Production counters emit and cache each complete forward block and the
predicate, then compose their clean logical inverses without retaining a
multi-million-gate full production tape.
"""
from collections import Counter
import hashlib
import json
import architecture as arch
import encoding_resources as er
import quantum_gates as qg
import reversible as base


class ValueCircuit(arch.ArchitectureCircuit):
    """ArchitectureCircuit with literal Clifford+T dependency accounting.

    The historical three-layer assigned CCX metric is inherited. A second
    schedule expands every CCX into the actual seven-T sequence in
    quantum_gates.expand_toffoli. A conservative five-layer serial CCX upper
    bound is also retained explicitly; neither replaces a hidden primitive.
    """
    def __init__(self, store_gates=True, mode='measured4t'):
        super().__init__(store_gates=store_gates, mode=mode)
        self._expanded_tdepth=[];self._expanded_depth=[];self._five_tdepth=[]
        self._expanded_records={};self._five_records={};self._expanded_counts=Counter()

    def _expanded_ensure(self):
        while len(self._expanded_tdepth)<self.wire_count:
            self._expanded_tdepth.append(0);self._expanded_depth.append(0);self._five_tdepth.append(0)

    def _track_atom(self, op, wires, condition=None):
        self._expanded_ensure()
        if op=='CCX':
            gates=qg.expand_toffoli([('CCX',*wires)])
        else:
            gates=[({'TD':'Tdg','SD':'Sdg'}.get(op,op),*wires)]
        for gate in gates:
            kind=gate[0];ww=gate[1:]
            td=max((self._expanded_tdepth[w] for w in ww),default=0)
            gd=1+max((self._expanded_depth[w] for w in ww),default=0)
            if condition is not None and op!='MZ':
                recorded_gd,recorded_td=self._expanded_records[condition]
                gd=max(gd,recorded_gd+1);td=max(td,recorded_td)
            if kind in ('T','Tdg'):td+=1
            for w in ww:self._expanded_tdepth[w]=td;self._expanded_depth[w]=gd
            self._expanded_counts[kind]+=1
            if op=='MZ':self._expanded_records[condition]=(gd,td)
        upper=max((self._five_tdepth[w] for w in wires),default=0)
        if condition is not None and op!='MZ':
            upper=max(upper,self._five_records[condition])
        upper += 5 if op=='CCX' else 1 if op in ('T','TD') else 0
        for w in wires:self._five_tdepth[w]=upper
        if op=='MZ':self._five_records[condition]=upper

    def _physical_operation(self,op,wires,condition,classical):
        super()._physical_operation(op,wires,condition,classical)
        self._track_atom(op,wires,condition)

    def _emit(self,op,c1,c2,target):
        super()._emit(op,c1,c2,target)
        if op!=arch.MACRO:
            wires=(target,) if op==base.X else ((c1,target) if op==base.CNOT else (c1,c2,target))
            self._track_atom({base.X:'X',base.CNOT:'CX',base.TOFFOLI:'CCX'}[op],wires)

    def resources(self):
        metrics=super().resources();self._expanded_ensure()
        assert metrics['logical_T_count']==self._expanded_counts['T']+self._expanded_counts['Tdg']
        assert max(self._five_tdepth,default=0)>=max(self._expanded_tdepth,default=0)
        metrics.update({'clifford_T_expanded_gate_counts':dict(self._expanded_counts),
            'explicit_7T_sequence_T_depth':max(self._expanded_tdepth,default=0),
            'explicit_Clifford_T_all_gate_dependency_depth':max(self._expanded_depth,default=0),
            'conservative_5layer_CCX_T_depth_upper':max(self._five_tdepth,default=0),
            'explicit_CCX_sequence_T_depth_from_clean_start':4,
            'CCX_sequence_source':'quantum_gates.expand_toffoli, literal emitted order',
            'arbitrary_rotation_count':0})
        return metrics


def _validate_value_graph(graph):
    if graph.q<1 or graph.width<graph.q+1 or graph.output not in graph.nodes:
        raise ValueError('Value graph must have a declared q+1-bit normalised output')


def _body(c,graph,tables,inputs):
    regs=dict(inputs)
    for name,node in graph.nodes.items():
        key=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        regs[name]=er.emit_node(c,node,[regs[a] for a in node.args],graph.q,graph.width,tables.get(key))
    return regs[graph.output][:graph.q+1]


def value_into(c,graph,tables,inputs,output):
    """Emit |latent>|z>|0_work> -> |latent>|z XOR a(latent)>|0_work>.

    output must be allocated before the body snapshot, preserving arbitrary
    dirty logical targets under Bennett compute/copy/fresh-inverse. The caller
    supplies the previously certified range 0 <= a <= 2**q.
    """
    _validate_value_graph(graph)
    if len(output)!=graph.q+1 or set(output)&set(sum((tuple(v) for v in inputs.values()),())):
        raise ValueError('Value output has wrong width or overlaps inputs')
    snap=c.snapshot();value=_body(c,graph,tables,inputs);end=c.gate_count
    for a,z in zip(value,output):c.cnot(a,z)
    c.uncompute(snap,end_gate_count=end)


def compile_value(graph,tables,mode='measured4t'):
    _validate_value_graph(graph)
    c=ValueCircuit(mode=mode);inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    output=c.alloc(graph.q+1,'clean numerical value XOR output')
    value_into(c,graph,tables,inputs,output)
    return c,inputs,output


def _conjunction(c,controls):
    controls=tuple(controls)
    if not controls:raise ValueError('Conjunction requires a control')
    if len(controls)==1:return controls[0]
    work=c.alloc(len(controls)-1,'bucket conjunction work')
    c.toffoli(controls[0],controls[1],work[0])
    for j in range(2,len(controls)):c.toffoli(work[j-2],controls[j],work[j-1])
    return work[-1]


def bucket_predicate_into(c,value,comparison,bucket,flag):
    """XOR a dyadic bucket's good predicate into flag and reset private work.

    For b<q: [2**b <= a < 2**(b+1)] [U < a], U in [0,2**(b+1)).
    For b=q: a_q, equal to [a=2**q] on the certified domain 0<=a<=2**q,
    with no comparison register. The construction is a Boolean permutation
    on all basis states; the top-bucket financial meaning uses that domain.
    """
    value=tuple(value);comparison=tuple(comparison);q=len(value)-1
    if q<1 or not 0<=bucket<=q or len(comparison)!=(bucket+1 if bucket<q else 0):
        raise ValueError('Invalid dyadic bucket or comparison width')
    if len(set(value+comparison+(flag,)))!=len(value)+len(comparison)+1:
        raise ValueError('Value, comparison and flag must be disjoint')
    if bucket==q:
        c.cnot(value[q],flag)
        return
    snap=c.snapshot()
    if bucket<q:
        for wire in value[bucket+1:]:c.x(wire)
        member=_conjunction(c,value[bucket:])
        smaller=c._ripple(comparison,value[:bucket+1],True)[1]
        end=c.gate_count;c.toffoli(member,smaller,flag)
    c.uncompute(snap,end_gate_count=end)


def compile_bucket(graph,tables,bucket,mode='measured4t'):
    _validate_value_graph(graph)
    if not 0<=bucket<=graph.q:raise ValueError('Bucket is outside 0,...,q')
    c=ValueCircuit(mode=mode);inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    comparison=c.alloc(bucket+1 if bucket<graph.q else 0,'bucket uniform comparator integer')
    flag=c.alloc(1,'dyadic bucket good flag')[0]
    snap=c.snapshot();value=_body(c,graph,tables,inputs);end=c.gate_count
    bucket_predicate_into(c,value,comparison,bucket,flag)
    c.uncompute(snap,end_gate_count=end)
    return c,inputs,comparison,flag


def range_predicate_into(c,value,comparison,range_bits,flag):
    """Mark U<a with c bits, given the separately proved a<=2**c.

    For c>0 the comparison input is zero-extended to c+1 bits so that the
    endpoint a=2**c is represented without wrap. For c=0, a in {0,1}, a
    single CNOT produces the exact Bernoulli flag. The promise is used only
    to interpret the Born mean; the emitted operation is a Boolean unitary.
    """
    value=tuple(value);comparison=tuple(comparison)
    if not 0<=range_bits<len(value) or len(comparison)!=range_bits:
        raise ValueError('Invalid certified comparison range')
    if len(set(value+comparison+(flag,)))!=len(value)+len(comparison)+1:
        raise ValueError('Range predicate registers overlap')
    if range_bits==0:
        c.cnot(value[0],flag)
        return
    zero=c.alloc(1,'range comparison high zero')
    c.compare_into(comparison+zero,value[:range_bits+1],flag)
    c.release(zero)


def compile_range(graph,tables,range_bits,mode='measured4t'):
    """Clean fused single-flag preparation with a certified smaller range."""
    _validate_value_graph(graph)
    if not 0<=range_bits<=graph.q:raise ValueError('Invalid range exponent')
    c=ValueCircuit(mode=mode);inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    comparison=c.alloc(range_bits,'certified range uniform comparison')
    flag=c.alloc(1,'certified range good flag')[0]
    snap=c.snapshot();value=_body(c,graph,tables,inputs);end=c.gate_count
    range_predicate_into(c,value,comparison,range_bits,flag)
    c.uncompute(snap,end_gate_count=end)
    return c,inputs,comparison,flag


def clipped_range_predicate_into(c,value,comparison,range_bits,flag):
    """Mark U_c<a using the full value, so the mean clips automatically.

    No smaller-range promise is needed. The comparison register is padded
    with clean zero bits to the full q+1 value width, including when c=0.
    Thus the Born flag mean is E[min(a,2**c)]/2**c exactly on a>=0.
    """
    value=tuple(value);comparison=tuple(comparison)
    if not 0<=range_bits<len(value) or len(comparison)!=range_bits:
        raise ValueError('Invalid clipped comparison range')
    if len(set(value+comparison+(flag,)))!=len(value)+len(comparison)+1:
        raise ValueError('Clipped range predicate registers overlap')
    zeros=c.alloc(len(value)-range_bits,'clipped comparison high zeros')
    c.compare_into(comparison+zeros,value,flag)
    c.release(zeros)


def compile_clipped_range(graph,tables,range_bits,mode='measured4t'):
    """Clean fused truncated-range flag, with its clipping bias explicit."""
    _validate_value_graph(graph)
    if not 0<=range_bits<=graph.q:raise ValueError('Invalid clipped range exponent')
    c=ValueCircuit(mode=mode);inputs={n:c.alloc(b,n) for n,b in graph.inputs.items()}
    comparison=c.alloc(range_bits,'clipped range uniform comparison')
    flag=c.alloc(1,'clipped range good flag')[0]
    snap=c.snapshot();value=_body(c,graph,tables,inputs);end=c.gate_count
    clipped_range_predicate_into(c,value,comparison,range_bits,flag)
    c.uncompute(snap,end_gate_count=end)
    return c,inputs,comparison,flag


def logical_gates(c):
    """Expand Boolean atoms and retain QROM instrument plans explicitly.

    Tuples for QROM contain the concrete LookupPlan, rather than claiming that
    measurements possess an operational reverse. Used by instrument tests.
    """
    for op,a,b,t in c.gates:
        if op==arch.MACRO:yield ('LOOKUP',c._macros[a])
        elif op==base.X:yield ('X',t)
        elif op==base.CNOT:yield ('CX',a,t)
        elif op==base.TOFFOLI:yield ('CCX',a,b,t)
        else:raise ValueError(op)


def _forward_blocks(graph,tables,cache,mode):
    """Fully emit every instantiated forward block once and cache metrics."""
    counts=Counter();expanded=Counter();retained=sum(graph.inputs.values());peak=retained
    depth=tdepth=mdepth=edepth=etdepth=fdepth=0;internal=arithmetic=0;blocks=[]
    for name,node in graph.nodes.items():
        keytable=('quantile',node.bits,graph.q) if node.kind=='quantile' else (node.meta,node.bits,graph.q)
        table=tables.get(keytable)
        cache_mode=mode if node.kind in ('quantile','function') else 'shared_exact_arithmetic'
        alias={a:j for j,a in enumerate(dict.fromkeys(node.args))}
        key=('value_v3',cache_mode,node.kind,node.meta,node.bits,graph.q,graph.width,
            tuple(graph.inputs.get(a,graph.width) for a in node.args),tuple(alias[a] for a in node.args))
        if table is not None:key+=(hashlib.sha256(str(table).encode()).hexdigest(),)
        if key not in cache:
            # Only an individual forward block is retained, never the whole
            # production financial graph. Its logical inverse is re-emitted
            # independently to avoid assuming equal directional depths.
            c=ValueCircuit(store_gates=True,mode=mode)
            unique_args={a:c.alloc(graph.inputs.get(a,graph.width)) for a in alias}
            args=[unique_args[a] for a in node.args]
            incoming=sum(map(len,unique_args.values()));er.emit_node(c,node,args,graph.q,graph.width,table)
            metrics=c.resources();inverse=ValueCircuit(store_gates=False,mode=mode)
            inverse.alloc(c.wire_count);inverse._macros=list(c._macros)
            for op,a,b,t in reversed(list(c.gates)):inverse._emit(op,a,b,t)
            inv=inverse.resources()
            assert metrics['primitive_gate_counts']==inv['primitive_gate_counts']
            assert metrics['clifford_T_expanded_gate_counts']==inv['clifford_T_expanded_gate_counts']
            cache[key]=(metrics,c.wire_count-incoming,metrics['active_wires_at_end']-incoming,inv)
        metrics,extra,live,inv=cache[key]
        counts.update(metrics['primitive_gate_counts']);expanded.update(metrics['clifford_T_expanded_gate_counts'])
        peak=max(peak,retained+extra);retained+=live
        depth+=metrics['all_gate_dependency_depth']+inv['all_gate_dependency_depth']
        tdepth+=metrics['logical_T_layer_depth_upper']+inv['logical_T_layer_depth_upper']
        mdepth+=metrics['measurement_layer_depth_upper']+inv['measurement_layer_depth_upper']
        edepth+=metrics['explicit_Clifford_T_all_gate_dependency_depth']+inv['explicit_Clifford_T_all_gate_dependency_depth']
        etdepth+=metrics['explicit_7T_sequence_T_depth']+inv['explicit_7T_sequence_T_depth']
        fdepth+=metrics['conservative_5layer_CCX_T_depth_upper']+inv['conservative_5layer_CCX_T_depth_upper']
        internal+=metrics['lookup_selector_internal_nodes'];arithmetic+=metrics['arithmetic_exact_Toffoli_count']
        blocks.append({'node':name,'kind':node.kind,'address_bits':node.bits,
            'logical_T_count':metrics['logical_T_count'],'measurements':metrics['measurement_count'],
            'primitive_gate_counts':metrics['primitive_gate_counts']})
    return {'counts':counts,'expanded':expanded,'retained':retained,'peak':peak,'depth':depth,'tdepth':tdepth,
        'mdepth':mdepth,'edepth':edepth,'etdepth':etdepth,'fdepth':fdepth,'internal':internal,'arithmetic':arithmetic,'blocks':blocks}


def _counts_twice(counts):return Counter({k:2*v for k,v in counts.items()})


def _reflection_profile(core_qubits):
    work=tuple(range(core_qubits,2*core_qubits-2))
    gates=qg.zero_phase(tuple(range(core_qubits)),work)
    old=qg.resources(gates,max(2*core_qubits-2,core_qubits))
    exact=qg.resources(qg.expand_toffoli(gates),old['logical_qubits'])
    return {'gate_counts':old['gate_counts'],'clifford_T_expanded_gate_counts':exact['gate_counts'],
        'T_count':old['T_count_exact_7T_Toffoli'],
        'historical_assigned_3layer_T_depth':old['T_depth_assigned_3layer_Toffoli_or_explicit_T_upper'],
        'explicit_sequence_T_depth':exact['T_depth_assigned_3layer_Toffoli_or_explicit_T_upper'],
        'conservative_5layer_CCX_T_depth_upper':5*old['gate_counts'].get('CCX',0),
        'explicit_all_gate_dependency_depth':exact['all_gate_dependency_depth'],
        'logical_qubits':old['logical_qubits'],'core_qubits':core_qubits,'arbitrary_rotations':0}


def _compose(graph,tables,cache,mode,bucket):
    _validate_value_graph(graph)
    if bucket is not None and not 0<=bucket<=graph.q:raise ValueError('Invalid bucket')
    f=_forward_blocks(graph,tables,cache,mode)
    extra=graph.q+1 if bucket is None else (bucket+1 if bucket<graph.q else 0)+1
    predicate=ValueCircuit(store_gates=False,mode=mode);value=predicate.alloc(graph.q+1)
    if bucket is None:
        target=predicate.alloc(graph.q+1)
        for a,z in zip(value,target):predicate.cnot(a,z)
        pm=predicate.resources()
    else:
        comparison=predicate.alloc(bucket+1 if bucket<graph.q else 0);flag=predicate.alloc(1)[0]
        bucket_predicate_into(predicate,value,comparison,bucket,flag);pm=predicate.resources()
    counts=_counts_twice(f['counts']);counts.update(pm['primitive_gate_counts'])
    expanded=_counts_twice(f['expanded']);expanded.update(pm['clifford_T_expanded_gate_counts'])
    peak=max(f['peak']+extra,f['retained']+extra+predicate.wire_count-(graph.q+1+extra))
    tcount=counts['T']+counts['TD']+7*counts['CCX']
    assert tcount==expanded['T']+expanded['Tdg']
    internal=2*f['internal'];arithmetic=2*f['arithmetic']+pm['arithmetic_exact_Toffoli_count']
    assert tcount=={'exact7t':14,'unitary8t':8,'measured4t':4}[mode]*internal+7*arithmetic
    keys={('quantile',n.bits,graph.q) if n.kind=='quantile' else (n.meta,n.bits,graph.q)
          for n in graph.nodes.values() if n.kind in ('quantile','function')}
    metrics={'architecture':mode,'oracle_kind':'numerical_value_XOR' if bucket is None else 'dyadic_bucket_flag',
        'q':graph.q,'bucket':bucket,'expanded_gate_counts':dict(counts),
        'clifford_T_expanded_gate_counts':dict(expanded),'logical_T_count':tcount,
        'logical_T_depth_serial_3layer_upper':f['tdepth']+pm['logical_T_layer_depth_upper'],
        'explicit_sequence_T_depth_serial_upper':f['etdepth']+pm['explicit_7T_sequence_T_depth'],
        'conservative_5layer_CCX_T_depth_serial_upper':f['fdepth']+pm['conservative_5layer_CCX_T_depth_upper'],
        'explicit_Clifford_T_all_gate_depth_serial_upper':f['edepth']+pm['explicit_Clifford_T_all_gate_dependency_depth'],
        'legacy_all_gate_depth_serial_upper':f['depth']+pm['all_gate_dependency_depth'],
        'measurement_count':counts['MZ'],'measurement_depth_serial_upper':f['mdepth']+pm['measurement_layer_depth_upper'],
        'classical_feedback_CZ_count':counts['CZ_IF'],'classical_reset_X_count':counts['X_IF'],
        'lookup_selector_internal_nodes':internal,'arithmetic_exact_Toffoli_count':arithmetic,
        'logical_qubits_serial_upper_bound':peak,'input_qubits':sum(graph.inputs.values()),
        'value_output_qubits':graph.q+1 if bucket is None else 0,
        'comparison_qubits':0 if bucket is None or bucket==graph.q else bucket+1,
        'flag_qubits':0 if bucket is None else 1,'arbitrary_rotation_count':0,
        'includes_complete_arithmetic_and_lookup_logical_inverse':True,
        'lookup_inverse':'fresh corrected instrument' if mode=='measured4t' else 'clean-work XOR inverse',
        'predicate':pm,'blocks':f['blocks'],
        'classical_table_entries':sum(len(tables[k]) for k in keys),
        'classical_table_program_bits':sum(len(tables[k])*(graph.width if k[0]=='quantile' else 2*graph.width) for k in keys),
        'counting':'fully emitted forward blocks plus fresh inverse; fully emitted clean predicate; serial depth upper bounds'}
    if bucket is not None:
        coherent=sum(graph.inputs.values())+metrics['comparison_qubits'];reflection=_reflection_profile(coherent+1)
        metrics.update({'weight_numerator':1 if bucket==graph.q else 1,
            'weight_power_of_two':0 if bucket==graph.q else bucket+1-graph.q,
            'uniform_H_count':coherent,'preparation_expanded_gate_counts':dict(expanded+Counter({'H':coherent})),
            'grover_zero_reflection':reflection,'grover_iterate_T_count':2*tcount+reflection['T_count'],
            'grover_iterate_measurement_count':2*counts['MZ'],
            'grover_iterate_classical_feedback_CZ_count':2*counts['CZ_IF'],
            'grover_iterate_classical_reset_X_count':2*counts['X_IF'],
            'grover_iterate_explicit_sequence_T_depth_serial_upper':2*metrics['explicit_sequence_T_depth_serial_upper']+reflection['explicit_sequence_T_depth'],
            'grover_iterate_conservative_5layer_T_depth_serial_upper':2*metrics['conservative_5layer_CCX_T_depth_serial_upper']+reflection['conservative_5layer_CCX_T_depth_upper'],
            'grover_iterate_logical_qubits_upper':max(peak,reflection['logical_qubits']),
            'grover_iterate_expanded_gate_counts':dict(_counts_twice(expanded)+Counter(reflection['clifford_T_expanded_gate_counts'])+Counter({'H':2*coherent,'Z':1})),
            'Born_probability':'E[(a/2**(b+1))*1{2**b<=a<2**(b+1)}]' if bucket<graph.q else 'P(a=2**q)'})
        # The fused compiler touches the known financial body twice. In a
        # numerical-value black-box formulation, V;predicate;V† touches that
        # body four times and copies the value twice. Record this alternative
        # rather than equating the two oracle models.
        black_counts=Counter({k:4*v for k,v in f['counts'].items()});black_counts.update(pm['primitive_gate_counts'])
        black_counts['CX']+=2*(graph.q+1)
        black_expanded=Counter({k:4*v for k,v in f['expanded'].items()});black_expanded.update(pm['clifford_T_expanded_gate_counts'])
        black_expanded['CX']+=2*(graph.q+1)
        bt=black_expanded['T']+black_expanded['Tdg']
        black_depth=2*f['etdepth']+pm['explicit_7T_sequence_T_depth']
        black_five=2*f['fdepth']+pm['conservative_5layer_CCX_T_depth_upper']
        black_peak=max(f['peak']+graph.q+1+extra,
            sum(graph.inputs.values())+predicate.wire_count)
        metrics['blackbox_value_oracle_alternative']={
            'contract':'V; clean bucket predicate; V†, with a q+1-bit numerical value output',
            'financial_body_invocations_per_A':4,'numerical_value_oracle_invocations_per_A':2,
            'expanded_gate_counts':dict(black_counts),'clifford_T_expanded_gate_counts':dict(black_expanded),
            'logical_T_count':bt,'measurement_count':black_counts['MZ'],
            'classical_feedback_CZ_count':black_counts['CZ_IF'],'classical_reset_X_count':black_counts['X_IF'],
            'explicit_sequence_T_depth_serial_upper':black_depth,
            'conservative_5layer_CCX_T_depth_serial_upper':black_five,
            'logical_qubits_serial_upper_bound':black_peak,'uniform_H_count':coherent,
            'grover_iterate_T_count':2*bt+reflection['T_count'],
            'grover_iterate_measurement_count':2*black_counts['MZ'],
            'grover_iterate_classical_feedback_CZ_count':2*black_counts['CZ_IF'],
            'grover_iterate_classical_reset_X_count':2*black_counts['X_IF'],
            'grover_iterate_explicit_sequence_T_depth_serial_upper':2*black_depth+reflection['explicit_sequence_T_depth'],
            'grover_iterate_conservative_5layer_T_depth_serial_upper':2*black_five+reflection['conservative_5layer_CCX_T_depth_upper'],
            'grover_iterate_logical_qubits_upper':max(black_peak,reflection['logical_qubits']),
            'grover_iterate_expanded_gate_counts':dict(_counts_twice(black_expanded)+Counter(reflection['clifford_T_expanded_gate_counts'])+Counter({'H':2*coherent,'Z':1})),
            'scope':'conservative composed black-box alternative; fused known-body implementation is the headline compiler'}
        metrics['financial_body_invocations_per_A']=2
        metrics['financial_body_invocations_per_Q']=4
    return metrics


def counted_value_composition(graph,tables,cache,mode='measured4t'):
    return _compose(graph,tables,cache,mode,None)


def counted_bucket_composition(graph,tables,cache,bucket,mode='measured4t'):
    return _compose(graph,tables,cache,mode,bucket)


def resource_api(profile,blackbox=False):
    """Additive A/Q cost interface for the confidence allocation controller."""
    if profile['oracle_kind'] not in ('dyadic_bucket_flag','certified_range_flag','clipped_range_flag'):raise ValueError('A flagged profile is required')
    p=profile['blackbox_value_oracle_alternative'] if blackbox else profile
    a_counts=Counter(p['clifford_T_expanded_gate_counts'])+Counter({'H':p['uniform_H_count']})
    q_counts=Counter(p['grover_iterate_expanded_gate_counts'])
    A={'T':p['logical_T_count'],'measurements':p['measurement_count'],
        'feedback_CZ':p['classical_feedback_CZ_count'],'reset_X':p['classical_reset_X_count'],
        'explicit_T_depth_serial_upper':p['explicit_sequence_T_depth_serial_upper'],
        'conservative_5layer_T_depth_serial_upper':p['conservative_5layer_CCX_T_depth_serial_upper'],
        'expanded_primitive_operations':sum(a_counts.values()),'uniform_H':p['uniform_H_count']}
    Q={'T':p['grover_iterate_T_count'],'measurements':p['grover_iterate_measurement_count'],
        'feedback_CZ':p['grover_iterate_classical_feedback_CZ_count'],
        'reset_X':p['grover_iterate_classical_reset_X_count'],
        'explicit_T_depth_serial_upper':p['grover_iterate_explicit_sequence_T_depth_serial_upper'],
        'conservative_5layer_T_depth_serial_upper':p['grover_iterate_conservative_5layer_T_depth_serial_upper'],
        'expanded_primitive_operations':sum(q_counts.values()),'uniform_H':2*p['uniform_H_count']}
    for kind,n in a_counts.items():A['gate_'+kind]=n
    for kind,n in q_counts.items():Q['gate_'+kind]=n
    return {'A':A,'Q':Q,'qubits':p['grover_iterate_logical_qubits_upper'],
        'model':'blackbox_numerical_value' if blackbox else 'fused_known_financial_body',
        'extra_flag_measurements_per_shot':1,
        'extra_complete_reinitialisations_per_shot':1,
        'counts_scope':'A includes uniform H and clean evaluator; Q includes two A, good Z and all-zero reflection; final shot measurement/reset separate'}


def _expand_count_dictionary(counts):
    """Exact primitive multiset of the declared literal CCX expansion."""
    result=Counter()
    for kind,n in counts.items():
        if kind=='CCX':result.update({'H':2*n,'CX':6*n,'T':4*n,'Tdg':3*n})
        else:result[{'TD':'Tdg','SD':'Sdg'}.get(kind,kind)]+=n
    return result


def _historical_body(resources,q,input_qubits=None,provenance_sha256=None):
    """Recover, and verify, an already emitted finite financial body."""
    mode=resources['architecture'];L=resources['input_qubits'] if input_qubits is None else input_qubits
    if mode not in arch.MODES or q<1 or L<1:raise ValueError('Unsupported historical financial profile')
    body=Counter()
    for block in resources['blocks']:body.update(block['primitive_gate_counts'])
    old=ValueCircuit(store_gates=False,mode=mode)
    a,b=old.alloc(q+1),old.alloc(q+1);flag=old.alloc(1)[0]
    old.compare_into(a,b,flag);om=old.resources()
    reconstructed=_counts_twice(body)+Counter(om['primitive_gate_counts'])
    if reconstructed!=Counter(resources['expanded_gate_counts']):
        raise ValueError('Retained forward blocks do not reproduce historical complete flag primitives')
    expected=body['T']+body['TD']+7*body['CCX']
    if 2*expected+om['logical_T_count']!=resources['logical_T_count']:
        raise ValueError('Historical T count does not match retained body/comparator')
    source=provenance_sha256 or hashlib.sha256(json.dumps(resources,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    # The old comparator and flag are allocated before the financial body.
    # Removing those q+2 reserved wires leaves a conservative body-space
    # bound, including any larger old comparator workspace in its peak.
    body_peak=resources['logical_qubits_serial_upper_bound']-(q+2)
    if body_peak<L:raise ValueError('Historical qubit bound is inconsistent')
    return {'counts':body,'expanded':_expand_count_dictionary(body),'T':expected,
        'qubits':body_peak,'input_qubits':L,'mode':mode,'source_sha256':source,
        'table_entries':resources.get('classical_table_entries',0),
        'table_program_bits':resources.get('classical_table_program_bits',0),
        'blocks':resources['blocks']}


def _profile_from_body(body,q,predicate,comparison_qubits,kind,coordinate,weight_power):
    """Exact additive counts and safe serial depths, without emitting B."""
    pm=predicate.resources();mode=body['mode'];L=body['input_qubits']
    counts=_counts_twice(body['counts']);counts.update(pm['primitive_gate_counts'])
    expanded=_counts_twice(body['expanded']);expanded.update(pm['clifford_T_expanded_gate_counts'])
    tcount=expanded['T']+expanded['Tdg']
    body_tphases=body['counts']['T']+body['counts']['TD'];body_ccx=body['counts']['CCX']
    four=2*(body_tphases+4*body_ccx)+pm['explicit_7T_sequence_T_depth']
    five=2*(body_tphases+5*body_ccx)+pm['conservative_5layer_CCX_T_depth_upper']
    three=2*(body_tphases+3*body_ccx)+pm['logical_T_layer_depth_upper']
    # Every body atom may be placed in serial order. This is intentionally
    # conservative and remains valid for a reverse Boolean tape or a fresh
    # corrected lookup instrument; no directional-depth identity is needed.
    expanded_depth=2*sum(body['expanded'].values())+pm['explicit_Clifford_T_all_gate_dependency_depth']
    retained=comparison_qubits+1
    private=predicate.wire_count-(q+1+retained)
    peak=body['qubits']+retained+private
    coherent=L+comparison_qubits;reflection=_reflection_profile(coherent+1)
    p={'architecture':mode,'oracle_kind':kind,'q':q,
        'bucket':coordinate if kind=='dyadic_bucket_flag' else None,
        'range_bits':coordinate if kind in ('certified_range_flag','clipped_range_flag') else None,
        'weight_numerator':1,'weight_power_of_two':weight_power,
        'expanded_gate_counts':dict(counts),'clifford_T_expanded_gate_counts':dict(expanded),
        'logical_T_count':tcount,'logical_T_depth_serial_3layer_upper':three,
        'explicit_sequence_T_depth_serial_upper':four,
        'conservative_5layer_CCX_T_depth_serial_upper':five,
        'explicit_Clifford_T_all_gate_depth_serial_upper':expanded_depth,
        'measurement_count':counts['MZ'],'measurement_depth_serial_upper':counts['MZ'],
        'classical_feedback_CZ_count':counts['CZ_IF'],'classical_reset_X_count':counts['X_IF'],
        'logical_qubits_serial_upper_bound':peak,'input_qubits':L,'comparison_qubits':comparison_qubits,
        'value_output_qubits':0,'flag_qubits':1,'arbitrary_rotation_count':0,
        'predicate':pm,'uniform_H_count':coherent,
        'preparation_expanded_gate_counts':dict(expanded+Counter({'H':coherent})),
        'grover_zero_reflection':reflection,
        'grover_iterate_T_count':2*tcount+reflection['T_count'],
        'grover_iterate_measurement_count':2*counts['MZ'],
        'grover_iterate_classical_feedback_CZ_count':2*counts['CZ_IF'],
        'grover_iterate_classical_reset_X_count':2*counts['X_IF'],
        'grover_iterate_explicit_sequence_T_depth_serial_upper':2*four+reflection['explicit_sequence_T_depth'],
        'grover_iterate_conservative_5layer_T_depth_serial_upper':2*five+reflection['conservative_5layer_CCX_T_depth_upper'],
        'grover_iterate_logical_qubits_upper':max(peak,reflection['logical_qubits']),
        'grover_iterate_expanded_gate_counts':dict(_counts_twice(expanded)+Counter(reflection['clifford_T_expanded_gate_counts'])+Counter({'H':2*coherent,'Z':1})),
        'includes_complete_arithmetic_and_lookup_logical_inverse':True,
        'lookup_inverse':'fresh corrected instrument' if mode=='measured4t' else 'clean-work XOR inverse',
        'financial_body_invocations_per_A':2,'financial_body_invocations_per_Q':4,
        'historical_emitted_resource_sha256':body['source_sha256'],
        'historical_complete_flag_reconstruction_verified':True,
        'classical_table_entries':body['table_entries'],'classical_table_program_bits':body['table_program_bits'],
        'blocks':body['blocks'],
        'counting':'exact multisets from previously emitted and hash-retained forward blocks; fresh emitted predicate/reflection; conservative serial depth and workspace bounds'}
    black_counts=Counter({k:4*v for k,v in body['counts'].items()});black_counts.update(pm['primitive_gate_counts'])
    black_counts['CX']+=2*(q+1)
    black_expanded=Counter({k:4*v for k,v in body['expanded'].items()});black_expanded.update(pm['clifford_T_expanded_gate_counts'])
    black_expanded['CX']+=2*(q+1)
    bt=black_expanded['T']+black_expanded['Tdg']
    black_four=4*(body_tphases+4*body_ccx)+pm['explicit_7T_sequence_T_depth']
    black_five=4*(body_tphases+5*body_ccx)+pm['conservative_5layer_CCX_T_depth_upper']
    black_peak=peak+(q+1)
    p['blackbox_value_oracle_alternative']={
        'contract':'V; clean flagged predicate; V†, with a q+1-bit held numerical value',
        'financial_body_invocations_per_A':4,'numerical_value_oracle_invocations_per_A':2,
        'expanded_gate_counts':dict(black_counts),'clifford_T_expanded_gate_counts':dict(black_expanded),
        'logical_T_count':bt,'measurement_count':black_counts['MZ'],
        'classical_feedback_CZ_count':black_counts['CZ_IF'],'classical_reset_X_count':black_counts['X_IF'],
        'explicit_sequence_T_depth_serial_upper':black_four,
        'conservative_5layer_CCX_T_depth_serial_upper':black_five,
        'logical_qubits_serial_upper_bound':black_peak,'uniform_H_count':coherent,
        'grover_iterate_T_count':2*bt+reflection['T_count'],
        'grover_iterate_measurement_count':2*black_counts['MZ'],
        'grover_iterate_classical_feedback_CZ_count':2*black_counts['CZ_IF'],
        'grover_iterate_classical_reset_X_count':2*black_counts['X_IF'],
        'grover_iterate_explicit_sequence_T_depth_serial_upper':2*black_four+reflection['explicit_sequence_T_depth'],
        'grover_iterate_conservative_5layer_T_depth_serial_upper':2*black_five+reflection['conservative_5layer_CCX_T_depth_upper'],
        'grover_iterate_logical_qubits_upper':max(black_peak,reflection['logical_qubits']),
        'grover_iterate_expanded_gate_counts':dict(_counts_twice(black_expanded)+Counter(reflection['clifford_T_expanded_gate_counts'])+Counter({'H':2*coherent,'Z':1})),
        'scope':'exact composed black-box counts; serial depth and workspace bounds'}
    return p


def prepare_existing_body(resources,q,input_qubits=None,provenance_sha256=None):
    """Verify a retained body once, then reuse it across bucket/range flags."""
    body=_historical_body(resources,q,input_qubits,provenance_sha256)
    body['q']=q
    return body


def profile_from_prepared_body(body,bucket):
    q=body['q']
    if not 0<=bucket<=q:raise ValueError('Invalid dyadic bucket')
    c=ValueCircuit(store_gates=False,mode=body['mode']);value=c.alloc(q+1)
    comparison=c.alloc(bucket+1 if bucket<q else 0);flag=c.alloc(1)[0]
    bucket_predicate_into(c,value,comparison,bucket,flag)
    p=_profile_from_body(body,q,c,len(comparison),'dyadic_bucket_flag',bucket,0 if bucket==q else bucket+1-q)
    p['Born_probability']='P(a=2**q)' if bucket==q else 'E[(a/2**(b+1))*1{2**b<=a<2**(b+1)}]'
    return p


def profile_from_existing_blocks(resources,q,bucket,input_qubits=None,provenance_sha256=None):
    """Reuse hash-retained emitted blocks; emit only new bucket/reflection."""
    if not 0<=bucket<=q:raise ValueError('Invalid dyadic bucket')
    return profile_from_prepared_body(prepare_existing_body(resources,q,input_qubits,provenance_sha256),bucket)


def profile_from_prepared_range(body,range_bits):
    q=body['q']
    if not 0<=range_bits<=q:raise ValueError('Invalid certified range exponent')
    c=ValueCircuit(store_gates=False,mode=body['mode']);value=c.alloc(q+1)
    comparison=c.alloc(range_bits);flag=c.alloc(1)[0]
    range_predicate_into(c,value,comparison,range_bits,flag)
    p=_profile_from_body(body,q,c,len(comparison),'certified_range_flag',range_bits,range_bits-q)
    p['Born_probability']='E[a/2**c]; requires the independently certified range a<=2**c'
    p['range_certificate_scope']='caller must retain the proof/certificate; this counter does not establish the range'
    return p


def profile_from_existing_range(resources,q,range_bits,input_qubits=None,provenance_sha256=None):
    """Fair single-flag opponent under a separately certified a<=2**c."""
    return profile_from_prepared_range(prepare_existing_body(resources,q,input_qubits,provenance_sha256),range_bits)


def profile_from_prepared_clipped_range(body,range_bits):
    q=body['q']
    if not 0<=range_bits<=q:raise ValueError('Invalid clipped range exponent')
    c=ValueCircuit(store_gates=False,mode=body['mode']);value=c.alloc(q+1)
    comparison=c.alloc(range_bits);flag=c.alloc(1)[0]
    clipped_range_predicate_into(c,value,comparison,range_bits,flag)
    p=_profile_from_body(body,q,c,len(comparison),'clipped_range_flag',range_bits,range_bits-q)
    p['Born_probability']='E[min(a,2**c)]/2**c; no smaller-range promise'
    p['reconstructed_finite_mean']='E[min(X,2**(c-q))]; unestimated excess-tail mean must be accounted for'
    p['range_certificate_scope']='valid for every finite a in 0,...,2**q; caller separately certifies clipping-tail allowance'
    return p


def profile_from_existing_clipped_range(resources,q,range_bits,input_qubits=None,provenance_sha256=None):
    """Truncated single-flag opponent using all q+1 value bits."""
    return profile_from_prepared_clipped_range(prepare_existing_body(resources,q,input_qubits,provenance_sha256),range_bits)
