"""Independent exact replay of the retained greedy scalar precision pool.

The audit reconstructs the finite input law, local scalar metadata, greedy
proposal order and certificate-conditioned choices. Its optional literal
emission check independently expands every distinct retained oracle block.
It does not infer a global optimum or a physical quantum price experiment.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse,gzip,hashlib,json,sys,time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()
def content_sha(value):return hashlib.sha256(canonical(value)).hexdigest()
def read(path):
    path=Path(path)
    return json.loads(gzip.open(path,'rt').read() if path.suffix=='.gz' else path.read_text())
def scalar_key(c):
    if c['schema']=='qbp-integral-quantile-v1':return ('quantile',F(c['L']),c['b_y'],c['q'])
    return ('scalar',c['name'],F(c['lo']),F(c['span']),tuple(sorted(c['params'].items())),c['q'],c['table_bits'],F(c['output_offset']))
def span_exp(x):
    x=F(x);assert x>0 and not x.numerator&(x.numerator-1) and not x.denominator&(x.denominator-1)
    return x.numerator.bit_length()-x.denominator.bit_length()
def strict_int(n,lo=0,hi=None):
    assert type(n) is int and n>=lo and (hi is None or n<=hi),('Invalid exact integer',n,lo,hi)
def exact_meta(c):
    assert c['q']==28
    if c['schema']=='qbp-integral-quantile-v1':
        assert c['b_y']==10 and F(c['L'])==6 and len(c['nodes'])==1024
        return len(c['nodes'])
    strict_int(c['table_bits'],6,12);N=1<<c['table_bits'];unit=1<<c['q']
    assert len(c['nodes'])==N+1
    assert c['fractional_remainder_bits']==c['q']+span_exp(c['span'])-c['table_bits']>=0
    lo,span=F(c['lo']),F(c['span']);largest=F(0);values=[]
    for j,node in enumerate(c['nodes']):
        assert F(node['input'])==lo+span*F(j,N)
        low,high=F(node['lower']),F(node['upper']);assert low<=high
        strict_int(node['integer']);value=node['integer'];values.append(value)
        centre=(low+high)/2;expected=max(0,(centre*unit+F(1,2)).numerator//(centre*unit+F(1,2)).denominator)
        assert value==expected
        error=max(abs(F(value,unit)-low),abs(F(value,unit)-high));assert error==F(node['error_upper']);largest=max(largest,error)
    assert largest==F(c['nodal_error_upper'])
    assert min(values)==c['value_min_integer'] and max(values)==c['value_max_integer']
    assert c['integer_width']==max(1,max(values).bit_length())
    assert F(c['first_derivative_upper'])>=0 and F(c['second_derivative_upper'])>=0
    interpolation=F(c['second_derivative_upper'])*(span/N)**2/8
    assert interpolation==F(c['linear_interpolation_error_upper'])
    assert F(c['interpolation_floor_error_upper'])==F(1,unit)
    assert largest+interpolation+F(1,unit)==F(c['output_error_upper'])
    for left,right in zip(values,values[1:]):
        if right<left:assert abs(right-left)<=left
    return len(values)


def audit(root,out,emit=False):
    start=time.perf_counter();root=Path(root);sys.path.insert(0,str(root))
    import integral_value_oracle as oracle
    import integral_value_certificate as certifier
    import integral_auxiliary_quadrature as aux
    dataset=root/'results/integral_precision_comparison.json.gz';data=read(dataset)
    assert data['schema']=='qbp-independent-precision-comparison-v1'
    assert F(data['monetary_radius'])==F(1,100) and F(data['alpha_per_selected_price'])==F(1,2700)
    assert data['source_sha256']==sha(root/'integral_precision_selection.py')
    for name,digest in data['dependency_source_hashes'].items():assert digest==sha(root/name)
    for name,digest in data['input_sha256'].items():assert digest==sha(root/'results'/name)
    archive_path=root/'results'/data['scalar_ordinate_archive'];assert sha(archive_path)==data['scalar_ordinate_archive_sha256']
    old_path=root/'results/integral_scalar_table_archive.json.gz'
    old_audit_path=Path('/private/tmp/qbp-integral-audit/results/integral_scalar_archive_audit.json')
    old_audit=read(old_audit_path);assert old_audit['archive_sha256']==sha(old_path)
    assert old_audit['reference_precision_bits']==416 and 'passed' in old_audit['status']
    old_archive=read(old_path)['records'];anchors={}
    for c in old_archive.values():
        k=scalar_key(c)
        if c['schema']=='qbp-integral-scalar-v1':
            assert c['table_bits']==12
            anchors[(k[0],*k[1:6],k[7])]=c
        else:anchors[k]=c
    archive={};digests=set();node_count=0;anchored_nodes=0
    with gzip.open(archive_path,'rt') as stream:
        for line in stream:
            row=json.loads(line);c=row['certificate'];assert row['sha256']==content_sha(c)
            assert row['sha256'] not in digests;digests.add(row['sha256']);key=scalar_key(c)
            if key in archive:assert content_sha(archive[key])==row['sha256']
            archive[key]=c;node_count+=exact_meta(c)
            if c['schema']=='qbp-integral-scalar-v1':
                anchor=anchors[(key[0],*key[1:6],key[7])]
                stride=1<<(12-c['table_bits'])
                assert c['nodes']==anchor['nodes'][::stride],('New ordinate is not an audited baseline subset',key)
                for field in ('exact_function_range','first_derivative_upper','second_derivative_upper','precision_bits'):
                    assert c[field]==anchor[field]
            else:
                assert c==anchors[key]
            anchored_nodes+=len(c['nodes'])
    assert len(digests)==data['unique_scalar_certificates']
    protocol=read(root/'results/generalisation_protocol.json');models={m['model_id']:m for m in protocol['models']}
    domains={d['model_id']:d for d in read(root/'results/integral_finite_input_domains.json')['domains']}
    finances=read(root/'results/integral_financial_grid_assessment.json')['records']
    finite_moments={m['model_id']:m for m in read(root/'results/integral_finite_grid_moments.json')['records']}
    cases={(f['model_id'],f['cap']):f for f in finances};kinds=('integral_r2','direct_r2','positive_r1')
    expected={(mid,cap,kind) for mid,cap in cases for kind in kinds}
    quantile=archive[('quantile',F(6),10,28)];mass=quantile['mass_interval'];factor=F(mass['lower_exact_fraction']),F(mass['upper_exact_fraction'])
    programs={};build_count=0;profile_count=0;emitted={};basis_lanes=0
    if emit:
        sys.path.insert(0,'/private/tmp/qbp-integral-audit/checks')
        from integral_audit_support import independent_gate_accounting,packed_replay,register_values

    def graph_for(mid,cap,kind,bits,encoding='centred'):
        model,domain=models[mid],domains[mid]
        g=oracle.build_graph(model,cap,q=28,b_y=10,beta_bits=8,strike_bits=8,table_bits=12,
            outer_bound=6,U_star=F(domain['U_star']),D2=certifier.joint_norm(2,F(model['lambda'])),
            D1=certifier.joint_norm(1,F(model['lambda'])),kind=kind,
            retained_bound=max(abs(F(x)) for x in domain['coordinate_bounds']),output_encoding=encoding)
        original={n:v.bits for n,v in g.nodes.items() if v.kind=='function'}
        assert set(bits)==set(original)
        for name,b in bits.items():strict_int(b,6,12);g.nodes[name].bits=b
        assert g.inputs==({'iy0':10,'iz0':8,'iz1':8,'istrike':8} if kind=='integral_r2' else {'iy0':10})
        assert all(n.bits==10 for n in g.nodes.values() if n.kind=='quantile')
        g.configuration['function_address_bits']=dict(bits)
        return g,original

    def tables_for(g):
        tables={};certificates={}
        for name,node in g.nodes.items():
            if node.kind=='quantile':
                c=archive[('quantile',F(node.meta),node.bits,g.q)];words=tuple(n['offset_integer'] for n in c['nodes'])
                key=('quantile',node.bits,str(node.meta))
            elif node.kind=='function':
                s=node.meta;keyc=('scalar',s.name,s.lo,s.span,tuple(sorted((k,str(v)) for k,v in s.parameters)),g.q,node.bits,F(0))
                assert keyc in archive,('Scalar certificate absent',name,keyc)
                c=archive[keyc];values=[n['integer'] for n in c['nodes']];words=[]
                assert c['fractional_remainder_bits']==g.q+span_exp(s.span)-node.bits
                for j,left in enumerate(values):
                    delta=values[j+1]-left if j+1<len(values) else 0
                    assert left<(1<<g.width) and abs(delta)<(1<<g.width)
                    words.append(left|(abs(delta)<<g.width)|(int(delta<0)<<(2*g.width)))
                words=tuple(words);key=('scalar',json.dumps(s.as_dict(),sort_keys=True),node.bits)
            else:continue
            tables[key]=words;certificates[name]=c
        return tables,certificates

    def evaluate(mid,cap,kind,bits):
        nonlocal build_count
        key=(mid,cap,kind,tuple(sorted(bits.items())))
        if key in programs:return programs[key]
        g,original=graph_for(mid,cap,kind,bits);tables,cc=tables_for(g)
        gc=certifier.certify_graph(g,tables,cc);build_count+=1
        if gc['status']!='passed':result=(g,tables,gc,None)
        else:
            transfer=F(cases[(mid,cap)]['kinds']['R1' if kind=='positive_r1' else 'R2']['financial_transfer_including_tail_and_quantile_rounding_upper'])
            moments={r:F(finite_moments[mid]['centred_U_moments'][str(r)]['upper_exact_fraction']) for r in (1,2,3)}
            budget=certifier.complete_bias_budget(models[mid],cap,kind,transfer,gc,factor,8,8,moments)
            if kind=='integral_r2':
                replacement=aux.auxiliary_bias(models[mid],cap,8,8,F(domains[mid]['U_star']),moments[2],moments[3])
                budget['total_deterministic_undiscounted_price_error_upper']=str(F(budget['total_deterministic_undiscounted_price_error_upper'])+F(budget['external_factor_upper'])*(F(replacement['total_auxiliary_bias_upper'])-F(budget['auxiliary_before_factor_upper'])))
            result=(g,tables,gc,budget)
        programs[key]=result;return result

    def profile_check(g,tables,p):
        nonlocal profile_count,basis_lanes
        profile_count+=1
        assert p['configuration']==g.configuration and p['source_sha256']==sha(root/'integral_value_oracle.py')
        assert p['integer_fractional_bits']==28 and p['integer_word_width']==g.width
        assert p['reflection_core_qubits']==sum(g.inputs.values())+29
        assert p['uniform_H']==sum(g.inputs.values())+28
        assert p['active_qubits_after_clean_prep']==sum(g.inputs.values())+30
        assert p['maximum_probability_integer']==1<<28 and p['comparison_bits']==28 and p['endpoint_high_zero_bits']==1
        assert [b['node'] for b in p['blocks']]==list(g.nodes)
        forward=Counter();retained=sum(g.inputs.values())+g.q+2;peak=retained;depth=0;five=0;lookups=0
        for record,(name,node) in zip(p['blocks'],g.nodes.items()):
            distinct=tuple(dict.fromkeys(node.args));aliases=tuple(distinct.index(x) for x in node.args)
            widths=tuple(g.inputs[a] if a in g.inputs else g.nodes[a].output_bits for a in node.args)
            assert record['kind']==node.kind and record['function']==(node.meta.as_dict() if isinstance(node.meta,oracle.ScalarSpec) else None)
            assert record['argument_widths']==list(widths) and record['argument_alias_pattern']==list(aliases) and record['output_bits']==node.output_bits
            table=tables.get(oracle.table_key(name,node)) if node.kind in ('quantile','function') else None
            th=hashlib.sha256(json.dumps(table,separators=(',',':')).encode()).hexdigest() if table is not None else None
            assert record['table_sha256']==th
            counts=Counter(record['primitive_counts']);assert not set(counts)-{'X','H','CX','CZ','Z','S','Sdg','T','Tdg'}
            for n in counts.values():strict_int(n)
            assert record['logical_T_count']==counts['T']+counts['Tdg']
            if emit:
                meta=node.meta.as_dict() if isinstance(node.meta,oracle.ScalarSpec) else str(node.meta)
                ek=(p['architecture'] if table is not None else 'shared',node.kind,json.dumps(meta,sort_keys=True),node.bits,node.output_bits,g.q,g.width,widths,aliases,th)
                if ek not in emitted:
                    c=oracle.vo.ValueCircuit(mode=p['architecture']);regs=[c.alloc(widths[node.args.index(a)]) for a in distinct];args=[regs[k] for k in aliases]
                    incoming=sum(map(len,regs));output=oracle.emit_node(c,node,args,g.q,g.width,table);metrics=c.resources();account=independent_gate_accounting(c)
                    literal=Counter()
                    for gate,n in account['expanded_literal_gate_counts'].items():literal[{'TD':'Tdg','SD':'Sdg'}.get(gate,gate)]+=n
                    er={'counts':literal,'T':account['logical_T_count'],'depth':account['actual_literal_T_dependency_depth'],
                        'five':metrics['conservative_5layer_CCX_T_depth_upper'],'lookups':metrics['lookup_macro_invocations'],'extra':c.wire_count-incoming,'output':len(output),'expanded_unitary_tape_sha256':account['expanded_unitary_tape_sha256']}
                    if node.kind=='function':
                        unit=1<<g.q;spanint=int(node.meta.span*unit);xs=[0,spanint,spanint//2,spanint//3,max(0,spanint-1),1,spanint//4,3*spanint//4]
                        assignments=[{regs[0]:x} for x in xs]
                        state,_=packed_replay(c,assignments);mask=(1<<g.width)-1;f=g.q+span_exp(node.meta.span)-node.bits;expectedvalues=[]
                        for x in xs:
                            word=table[x>>f];base=word&mask;delta=(word>>g.width)&mask;down=(word>>(2*g.width))&1;adj=delta*(x&((1<<f)-1))>>f
                            expectedvalues.append(base-adj if down else base+adj)
                        assert register_values(state,output,len(xs))==expectedvalues
                        assert register_values(state,regs[0],len(xs))==xs
                        keep=set(output+regs[0]);assert all(v==0 for w,v in enumerate(state) if w not in keep)
                        er['basis_inputs_integer']=xs;er['basis_outputs_integer']=expectedvalues;er['clean_scratch_verified']=True;er['argument_preserved']=True
                        basis_lanes+=len(xs)
                    emitted[ek]=er;del c
                er=emitted[ek]
                assert er['counts']==counts and er['T']==record['logical_T_count'] and er['depth']==record['explicit_T_depth_serial_upper']
                assert er['five']==record['conservative_5layer_T_depth_serial_upper'] and er['lookups']==record['lookup_invocations']
                assert er['extra']==record['extra_qubits_upper'] and er['output']==record['output_bits']
            forward.update(counts);peak=max(peak,retained+record['extra_qubits_upper']);retained+=node.output_bits
            depth+=record['explicit_T_depth_serial_upper'];five+=record['conservative_5layer_T_depth_serial_upper'];lookups+=record['lookup_invocations']
        assert forward==Counter(p['primitive_counts_forward_body'])
        complete=Counter({k:2*v for k,v in forward.items()});complete.update(p['comparator_primitive_counts']);complete['H']+=p['uniform_H']
        assert complete==Counter(p['primitive_counts']) and p['logical_T_count']==complete['T']+complete['Tdg']
        assert p['lookup_invocations_forward_and_cleanup']==2*lookups
        comp=oracle.vo.ValueCircuit(mode=p['architecture']);a,b=comp.alloc(g.q+1),comp.alloc(g.q+1);flag=comp.alloc(1)[0];comp.compare_into(a,b,flag);cm=comp.resources()
        assert Counter(cm['clifford_T_expanded_gate_counts'])==Counter(p['comparator_primitive_counts'])
        assert p['explicit_T_depth_serial_upper']==2*depth+cm['explicit_7T_sequence_T_depth']
        assert p['conservative_5layer_T_depth_serial_upper']==2*five+cm['conservative_5layer_CCX_T_depth_upper']
        peak=max(peak,retained+cm['wires']-2*(g.q+1)-1);assert p['logical_qubits_upper']==peak
        keys={oracle.table_key(n,v) for n,v in g.nodes.items() if v.kind in ('quantile','function')}
        assert p['scalar_table_entries']==sum(len(tables[k]) for k in keys)
        assert p['scalar_table_program_bits']==sum(len(tables[k])*(g.width if k[0]=='quantile' else 2*g.width+1) for k in keys)

    groups={};baseline={};records=data['records']
    for i,r in enumerate(records):
        key=r['model_id'],r['cap'],r['kind'];assert key in expected;groups.setdefault(key,[]).append((i,r))
        bits=r.get('node_address_bits',r['base_graph_certificate']['configuration']['function_address_bits']);g,tables,gc,budget=evaluate(*key,bits)
        saved_gc=dict(r['base_graph_certificate']);saved_hash=saved_gc.pop('precision_search_certificate_sha256')
        assert saved_hash==content_sha(saved_gc) and gc==saved_gc
        if budget is None:
            assert gc['status']=='rejected' and gc['failures'] and r['selected_strategy'] is None
            assert r['allocation_limits']==['unchanged_all12_baseline'] and all(b==12 for b in bits.values())
            assert 'complete_deterministic_budget' not in r and 'profiles' not in r
            assert key not in baseline;baseline[key]=r
            continue
        assert gc['status']=='passed' and not gc['failures']
        assert F(budget['total_deterministic_undiscounted_price_error_upper'])==F(r['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper'])
        assert len(set(r['allocation_limits']))==len(r['allocation_limits'])
        for label in r['allocation_limits']:
            if label=='unchanged_all12_baseline':assert all(b==12 for b in bits.values()) and key not in baseline;baseline[key]=r
            else:assert F(label) in (F(1,1000),F(1,400),F(1,200)) and F(budget['total_deterministic_undiscounted_price_error_upper'])<=F(label)
        for architecture,encodings in r['profiles'].items():
            assert architecture in ('exact7t','unitary8t')
            assert set(encodings)==({'centred'} if r['kind']=='positive_r1' else {'centred','positive_part','negative_part'})
            for encoding,p in encodings.items():
                pg,_=graph_for(*key,bits,encoding);pt,_=tables_for(pg);assert p['architecture']==architecture;profile_check(pg,pt,p)
        eligible=[s for s in r['strategies'] if s['result']['kind']!='rejected']
        assert r['selected_strategy']==(min(eligible,key=lambda s:(s['result']['T_count'],s['architecture'],s['strategy'])) if eligible else None)
    assert set(groups)==set(baseline)==expected
    assert all(len(v)<=4 for v in groups.values())
    trialgroups={};downsteps=0;infeasible=0;trialchecks=0
    for row in data['search_trials']:
        key=row['model_id'],row['cap'],row['kind'];assert key in expected
        limit=F(row['allocation_limit']);assert limit in (F(1,1000),F(1,400),F(1,200));trialgroups.setdefault((*key,limit),[]).append(row)
    for key in sorted(expected):
        base=baseline[key];initial=dict(base.get('node_address_bits',base['base_graph_certificate']['configuration']['function_address_bits']));initial_bias=F(base['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper']) if 'complete_deterministic_budget' in base else None
        for limit in (F(1,1000),F(1,400),F(1,200)):
            rows=trialgroups.get((*key,limit),[])
            if initial_bias is None or initial_bias>limit:
                assert len(rows)==1 and rows[0]['status']=='starting_program_outside_this_allocation';infeasible+=1
                assert all(str(limit) not in r['allocation_limits'] for _,r in groups[key]);continue
            bits=dict(initial);old_bias=initial_bias;offset=0;step=0
            while True:
                names=[name for name,b in bits.items() if b>6];batch=rows[offset:offset+len(names)]
                assert len(batch)==len(names),('Truncated greedy proposals',key,limit,step)
                choices=[]
                for name,row in zip(names,batch):
                    assert row['node']==name and row['from_bits']==bits[name] and row['to_bits']==bits[name]-1
                    proposed=dict(bits);proposed[name]-=1
                    g,t,cc,b=evaluate(*key,proposed);passed=cc['status']=='passed' and b is not None and F(b['total_deterministic_undiscounted_price_error_upper'])<=limit
                    assert row['accepted_by_certificate']==passed
                    if b is not None:assert F(row['deterministic_error_upper'])==F(b['total_deterministic_undiscounted_price_error_upper'])
                    assert row['current_width_map']==bits and row['proposed_width_map']==proposed and row['step']==step
                    assert row['graph_certificate_sha256']==content_sha(cc)
                    assert row['graph_failures']==cc['failures']
                    assert type(row['chosen']) is bool
                    if passed:
                        det=F(row['deterministic_error_upper']);increase=max(F(1,1<<28),det-old_bias);score=F(1<<(bits[name]-1))/increase
                        assert F(row['greedy_score'])==score
                        choices.append((score,name,proposed,det,row))
                    else:assert row['greedy_score'] is None and not row['chosen']
                    trialchecks+=1
                offset+=len(names)
                if not choices:break
                choices.sort(key=lambda c:(-c[0],c[1]));score,name,bits,old_bias,row=choices[0];downsteps+=1;step+=1
                assert [item['node'] for item in batch if item['chosen']]==[name]
            assert offset==len(rows),('Unconsumed trial rows',key,limit,len(rows)-offset)
            matches=[r for _,r in groups[key] if str(limit) in r['allocation_limits']];assert len(matches)==1 and matches[0]['node_address_bits']==bits
    assert len(data['selected_prices'])==len(cases)
    selected_cases=set()
    for row in data['selected_prices']:
        key=row['model_id'],row['cap'];assert key in cases and key not in selected_cases;selected_cases.add(key)
        assert set(row['methods'])==set(kinds)
        for kind in kinds:
            candidates=[(i,r) for i,r in groups[(*key,kind)] if r['selected_strategy'] is not None]
            winner=min(candidates,key=lambda p:p[1]['selected_strategy']['result']['T_count']) if candidates else None
            result=row['methods'][kind]
            assert (result is None)==(winner is None)
            if winner:assert result['candidate_index']==winner[0] and result['selected']==winner[1]['selected_strategy']
        choice=min((k for k in kinds if row['methods'][k]),key=lambda k:row['methods'][k]['selected']['result']['T_count'])
        assert row['quantum_only_method']==choice
    assert len(data['all27_scope'])==27
    assert {(r['model_id'],r['cap']) for r in data['all27_scope']}=={(m['model_id'],cap) for m in protocol['models'] for cap in (1,10,25)}
    for r in data['all27_scope']:
        if models[r['model_id']]['outer_dimension']>1:assert 'not certified' in r['status'] and 'not ranked' in r['status']
    result={'status':'passed','comparison_sha256':sha(dataset),'scalar_archive_sha256':sha(archive_path),'audit_source_sha256':sha(Path(__file__)),
        'audit_dependency_sha256':{'integral_audit_support.py':sha('/private/tmp/qbp-integral-audit/checks/integral_audit_support.py')} if emit else {},'finite_financial_cases':len(cases),'matched_method_groups':len(groups),'retained_candidates':len(records),'profiles_checked':profile_count,
        'unique_scalar_certificates':len(digests),'scalar_nodes_checked':node_count,'ordinates_linked_exactly_to_independently_audited_all12_baseline':anchored_nodes,'original_scalar_audit_sha256':sha(old_audit_path),'unique_reconstructed_programs':build_count,
        'greedy_proposals_replayed':trialchecks,'greedy_downsteps':downsteps,'starting_infeasible_allocations_retained':infeasible,
        'literal_emission_enabled':emit,'distinct_blocks_literally_expanded':len(emitted),'scalar_basis_lanes':basis_lanes,
        'deterministic_tie_rule':'Descending exact rational cell-saving/error proxy, then ascending node name; retained strategies use complete T, architecture and strategy, and candidate ties use stable retention order.',
        'scope':'Fixed q28, quantile10 and auxiliary8 input laws; archived scalar lattice/error metadata, exact ordinate links to the prior 416-bit independent all12 audit, and full-graph certificate replay; exact greedy search and retained full-cost choice. Literal emission and scalar basis replay are included only when enabled. No global precision optimum, physical experiment, theoretical originality or practical quantum advantage claim.',
        'emitted_block_records':[{'signature':signature,'signature_sha256':content_sha(signature),**record} for signature,record in sorted(emitted.items(),key=lambda pair:content_sha(pair[0]))],
        'elapsed_seconds_diagnostic':time.perf_counter()-start}
    Path(out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='emitted_block_records'}));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='/private/tmp/qbp-integral-value-step');p.add_argument('--out',default='/private/tmp/qbp-precision-selector-review/results/precision_selector_audit.json');p.add_argument('--emit',action='store_true')
    a=p.parse_args();audit(a.root,a.out,a.emit)
