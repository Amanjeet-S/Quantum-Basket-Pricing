"""A fully costed positive shared-threshold R1 preparation.

The layer-cake and common-random-number identity is established mathematics.
This source implements one project-specific finite financial preparation and
its conditional probability certificate. No theoretical novelty is claimed.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
from collections import Counter
from itertools import product
import argparse, copy, gzip, hashlib, json, sys, time

sys.dont_write_bytecode = True
ROOT = None
RESULTS = None
FACTOR = None

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v): return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def write(p,v):
    if p.exists(): raise FileExistsError('Preserve existing evidence')
    p.write_text(json.dumps(v,indent=2,default=str)+'\n')
def configure(repo, factor_resources, output):
    global ROOT, RESULTS, FACTOR, oracle, bodycert, vo, base, qg, retained, informed, arb, ctx, exact_endpoint
    ROOT=Path(repo).resolve();RESULTS=Path(output).resolve();FACTOR=Path(factor_resources).resolve()
    RESULTS.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(ROOT))
    import boundary_value_oracle as oracle
    import boundary_value_certificate as bodycert
    import value_oracle as vo
    import reversible as base
    import quantum_gates as qg
    import boundary_precision_selection as retained
    import correction_information_priors as informed
    from flint import arb, ctx
    from continuous_quadrature import exact_endpoint

DEPENDENCIES=('boundary_value_oracle.py','boundary_value_certificate.py',
 'integral_value_oracle.py','integral_value_certificate.py','boundary_derivative_certificate.py',
 'value_oracle.py','reversible.py','architecture.py','quantum_gates.py',
 'boundary_precision_selection.py','correction_information_priors.py','integral_pricing_budget.py',
 'qpe_estimator.py','qft_synthesis.py','certified_scalar.py','taylor_pricing.py','continuous_quadrature.py')
INPUTS=('gaussian_derivative_norms.json','boundary_derivative_norms.json',
 'qft_rotation_tapes.json','qft_rotation_certificates.json')
def signatures():
    return {'source_sha256':sha(__file__),'factor_resources_sha256':sha(FACTOR),
      'dependencies_sha256':{n:sha(ROOT/n) for n in DEPENDENCIES},
      'inputs_sha256':{n:sha(ROOT/'results'/n) for n in INPUTS}}

def register():
    d=json.loads(FACTOR.read_text());pairs=[]
    for row in d['paired_comparison']:
        pairs.append([row['original_model_id'],row['selected_candidate_id']])
    write(RESULTS/'shared_threshold_protocol.json',{'schema':'qbp-shared-threshold-protocol-v1',
      'registered_before_fresh_production_execution':True,
      'registration_utc':datetime.now(timezone.utc).isoformat(),'research_record_date':'2026-10-10',
      **signatures(),'selected_basis_policy':'Previously selected complete-resource winning basis, fixed before this new implementation.',
      'pairs':pairs,'cap':'1','outer_L':'6','outer_quantile_bits':10,'q':28,
      'outer_scalar_table_bits':12,'inner_normal_L':'6','inner_normal_bits':10,
      'strike_uniform_bits':4,'strike_log_table_bits':6,
      'architectures':['exact7t','unitary8t'],'no_precision_repairs':True,
      'representation':'H times indicator {mG < log(1+H S/K)-lambda N <= mA}; common S,N for both endpoints.',
      'strict_endpoints':'Lower strict, upper non-strict, implemented by lt(lower,threshold) AND NOT lt(upper,threshold).',
      'normal_decode':'Inner normal labels are a normalised approximation to the full Gaussian. Only the original financial Gaussian box mass is decoded externally.',
      'certificate':'Strike midpoint curvature + two-sided inner normal CDF error + Gaussian-density endpoint-arithmetic sandwich + unchanged outer transfer.',
      'confidence':'Fresh canonical QPE and exact binomial amplification, alpha=1/3600 per price, monetary radius=1/100.',
      'financial_price_answers_excluded':True,'hardware_execution':False,
      'data_origin':'Synthetic exact rational no-arbitrage models only; no licensed observations.',
      'prior_screen_known':'The already frozen mathematical feasibility screen predicts unchanged ideal QPE phase sizes. No emitted cost or endpoint-arithmetic result was inspected before this registration.'})

def build_graph(model,domain):
    old=oracle.build_graph(model,F(1),q=28,b_y=10,table_bits=12,outer_bound=6,
       U_star=F(domain['U_star']),retained_bound=max(abs(F(v)) for v in domain['coordinate_bounds']),kind='positive_r1')
    mgold=old.semantic_nodes['mG_offset'];logold=old.semantic_nodes['logT']
    BG=-old.nodes[old.semantic_nodes['PG']].meta.lo
    needed=set()
    def visit(n):
        if n in old.inputs or n in needed:return
        needed.add(n)
        for a in old.nodes[n].args:visit(a)
    visit(mgold);visit(logold)
    g=oracle.IntegralGraph(28,12);remap={n:g.inp(n,b) for n,b in old.inputs.items()}
    for name,node in old.nodes.items():
        if name not in needed:continue
        new=g.op(node.kind,*(remap[a] for a in node.args),meta=node.meta,bits=node.bits,
            output_bits=0,bounds=old.bounds[name]);remap[name]=new
        if node.kind=='constant':g.constant_names[node.meta]=new
    mG,logT=remap[mgold],remap[logold];lam=F(model['lambda']);L=F(6)
    J=oracle.old.power_span(max(BG,lam*L)+F(1,64)+F(16,1<<g.q))
    lower=g.add(mG,g.const(J-BG));upper=g.add(lower,logT)
    # Both endpoints share the same rounded lower register and logT>=0.
    normal_label=g.inp('inner_normal_label',10)
    nplusL=g.op('quantile',normal_label,meta=L,bits=10,bounds=(0,2*L))
    normal_scaled=g.scale(nplusL,lam)
    strike_label=g.inp('uniform_strike_label',4)
    strike_uniform=g.uniform(strike_label,4)
    strike_ratio=g.scale(strike_uniform,F(1,100))
    strike_span=oracle.old.power_span(F(1,100))
    strike=g.scalar('log',strike_ratio,1,strike_span,output_bounds=(0,strike_span))
    g.nodes[strike].bits=6
    threshold=g.sub(g.add(strike,g.const(J+lam*L)),normal_scaled)
    g.semantic_nodes={'lower_endpoint':lower,'upper_endpoint':upper,'threshold':threshold,
        'mG_offset':mG,'logT':logT,'normal_plus_L':nplusL,'strike_log':strike}
    g.configuration=copy.deepcopy(old.configuration)
    g.configuration.update({'pricing_representation':'shared_threshold_interval',
       'q':28,'kind':'positive_r1','output_encoding':'positive',
       'monetary_signal_normaliser':'1','positive_R1_normaliser':'1',
       'monetary_probability_decode_factor':1,'strike_bits':4,'inner_normal_bits':10,
       'inner_normal_L':'6','strike_log_table_bits':6,'common_endpoint_offset':str(J),
       'placeholder_output_scope':'Constant-zero numerical output is only an arithmetic-graph interface. It is not the binary good flag or its probability.'})
    return g.finish(g.const(0))

def emit_interval(c,lower,threshold,upper,flag):
    """XOR the exact interval predicate, preserving all argument registers."""
    snap=c.snapshot()
    lower_lt=c.unsigned_less_than(lower,threshold)[0]
    upper_lt=c.unsigned_less_than(upper,threshold)[0]
    end=c.gate_count
    c.x(upper_lt);c.toffoli(lower_lt,upper_lt,flag);c.x(upper_lt)
    c.uncompute(snap,end_gate_count=end)

def compile_interval(graph,tables,mode='exact7t',store_gates=True):
    c=vo.ValueCircuit(store_gates=store_gates,mode=mode)
    inputs={name:c.alloc(bits,name) for name,bits in graph.inputs.items()};flag=c.alloc(1,'positive interval flag')[0]
    snap=c.snapshot();oracle.old.emit_body(c,graph,tables,inputs);end=c.gate_count
    regs=c.node_registers;s=graph.semantic_nodes
    emit_interval(c,regs[s['lower_endpoint']],regs[s['threshold']],regs[s['upper_endpoint']],flag)
    c.uncompute(snap,end_gate_count=end)
    c.uniform_preparation_gates=[('H',wire) for reg in inputs.values() for wire in reg]
    c.reflection_core=tuple(w for reg in inputs.values() for w in reg)+(flag,)
    return c,inputs,flag

def counted_profile(g,tables,cache,mode):
    # These emitted body blocks are reused; the old objective comparator and
    # old q-bit uniform register are deliberately absent from this profile.
    old=oracle.counted_composition(g,tables,cache,mode)
    c=vo.ValueCircuit(store_gates=False,mode=mode)
    a,b,d=(c.alloc(g.width) for _ in range(3));flag=c.alloc(1)[0]
    emit_interval(c,a,b,d,flag);cm=c.resources()
    incoming=3*g.width+1;extra=cm['wires']-incoming
    forward=Counter(old['primitive_counts_forward_body'])
    counts=Counter({k:2*v for k,v in forward.items()});counts.update(cm['clifford_T_expanded_gate_counts'])
    inputs=sum(g.inputs.values());counts['H']+=inputs
    retained_wires=inputs+1;peak=retained_wires
    for block in old['blocks']:
        peak=max(peak,retained_wires+block['extra_qubits_upper']);retained_wires+=block['output_bits']
    peak=max(peak,retained_wires+extra)
    body_exp=sum(b['explicit_T_depth_serial_upper'] for b in old['blocks'])
    body_cons=sum(b['conservative_5layer_T_depth_serial_upper'] for b in old['blocks'])
    return {'schema':'qbp-shared-threshold-financial-profile-v1','architecture':mode,
      'configuration':g.configuration,'primitive_counts':dict(counts),
      'primitive_counts_forward_body':dict(forward),'logical_T_count':counts['T']+counts['Tdg'],
      'reflection_core_qubits':inputs+1,'uniform_H':inputs,'logical_qubits_upper':peak,
      'active_qubits_after_clean_prep':inputs+1,
      'explicit_T_depth_serial_upper':2*body_exp+cm['explicit_7T_sequence_T_depth'],
      'conservative_5layer_T_depth_serial_upper':2*body_cons+cm['conservative_5layer_CCX_T_depth_upper'],
      'explicit_all_gate_depth_serial_upper':old['explicit_all_gate_depth_serial_upper']+cm['explicit_Clifford_T_all_gate_dependency_depth'],
      'lookup_invocations_forward_and_cleanup':old['lookup_invocations_forward_and_cleanup'],
      'blocks':old['blocks'],'predicate_primitive_counts':cm['clifford_T_expanded_gate_counts'],
      'predicate_extra_qubits_upper':extra,'integer_fractional_bits':g.q,'integer_word_width':g.width,
      'scalar_table_entries':old['scalar_table_entries'],'scalar_table_program_bits':old['scalar_table_program_bits'],
      'source_sha256':sha(__file__),
      'scope':'Fully emitted compositional coherent body and exact clean interval predicate, retaining only actual uniform inputs and objective flag. Logical resource caps, not physical runtime.'}

def certificate(g,body,ordinates,prep,nu):
    if body['status']!='passed':return {'status':'rejected','reason':'Numeric graph domain/range guard failure','body_certificate':body}
    records={x['node']:x for x in body['nodes']};s=g.semantic_nodes;lam=F(prep['model']['lambda'])
    quantile=ordinates[s['normal_plus_L']];eN=F(quantile['maximum_node_error_upper'])
    massiv=quantile['mass_interval'];a_lo=F(massiv['lower_exact_fraction']);a_hi=F(massiv['upper_exact_fraction'])
    # tail/2+a/(2N) is decreasing in a; use its certified lower endpoint.
    D= (1-a_lo)/2+a_lo/F(2*(1<<10))+nu[0]*eN
    el=F(records[s['lower_endpoint']]['absolute_output_error_upper'])
    eu=F(records[s['upper_endpoint']]['absolute_output_error_upper'])
    ez=F(records[s['threshold']]['absolute_output_error_upper'])
    arithmetic=nu[0]*(el+eu+2*ez)/lam
    strike=F(1,100)**2*(nu[1]/lam**2+nu[0]/lam)/(12*(1<<8))
    eta=min(F(1),2*D+arithmetic+strike)
    # Actual integer ordering is established structurally. Numeric independent
    # interval guards establish overflow cannot break that integer identity.
    return {'schema':'qbp-shared-threshold-probability-certificate-v1','status':'passed',
      'configuration':g.configuration,'nodes':body['nodes'],'body_arithmetic_certificate':body,
      'probability_error_upper':str(eta),'monetary_value_error_upper':str(eta),
      'probability_actual_bounds':['0','1'],'inner_normal_CDF_error_upper':str(D),
      'inner_normal_node_rounding_error_upper':str(eN),'inner_normal_mass_interval':massiv,
      'endpoint_arithmetic_errors':{'lower':str(el),'upper':str(eu),'threshold':str(ez)},
      'arithmetic_probability_sandwich_upper':str(arithmetic),
      'strike_midpoint_bias_upper':str(strike),'inner_normal_interval_bias_upper':str(2*D),
      'implemented_upper_minus_lower':'The same represented lower register plus the represented nonnegative logT register, without overflow.',
      'decoded_observable':'H times native Boolean interval flag',
      'scope':'For every exact actual outer financial coordinate, conditional expectation of the binary implemented flag differs from R1/H by at most eta. This is not a pointwise error between a Boolean flag and a numerical payoff.',
      'normalisation':'Only original outer financial box mass is applied externally. Inner Gaussian truncation/quantisation is entirely charged through the CDF discrepancy.'}

def run():
    protocol=json.loads((RESULTS/'shared_threshold_protocol.json').read_text())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Frozen sources or inputs changed')
    data=json.loads(FACTOR.read_text());norms=json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())
    nu={r:F(next(x['norm_upper'] for x in norms['unit_norms'] if x['unit_derivative_order']==r and x['passed'])) for r in (0,1)}
    cache={};table_cache={};rows=[];archive=[];seen=set();start=time.perf_counter()
    for original,candidate in protocol['pairs']:
        prep=next(x for x in data['preparations'] if x['original_model_id']==original and x['candidate_id']==candidate)
        old=next(x for x in data['records'] if x['original_model_id']==original and x['factor_basis_candidate_id']==candidate and x['kind']=='positive_r1')
        g=build_graph(prep['model'],prep['domain'])
        tables,ordinates=oracle.make_tables(g,table_cache,True)
        body=bodycert.certify_graph(g,tables,ordinates)
        cert=certificate(g,body,ordinates,prep,nu)
        row={'original_model_id':original,'candidate_id':candidate,'model':prep['model'],
          'graph_semantic_nodes':g.semantic_nodes,'graph_configuration':g.configuration,
          'graph_nodes':[{ 'name':n,'kind':v.kind,'args':v.args,'meta':v.meta.as_dict() if isinstance(v.meta,oracle.ScalarSpec) else str(v.meta),
            'bits':v.bits,'output_bits':v.output_bits,'ideal_bounds':list(map(str,g.bounds[n]))} for n,v in g.nodes.items()],
          'probability_certificate':cert,'profiles':{},'strategies':[],
          'retained_direct_selected_strategy':old['selected_strategy']}
        for n,cc in ordinates.items():
            h=digest(cc)
            if h not in seen:
                seen.add(h);key=oracle.table_key(n,g.nodes[n]);archive.append({'certificate_sha256':h,'certificate':cc,
                    'packed_table':tables[key],'packed_table_sha256':digest(tables[key]),'graph_word_width':g.width})
        if cert['status']=='passed':
            eta=F(cert['probability_error_upper']);bud=copy.deepcopy(old['complete_deterministic_budget'])
            mass=F(bud['external_factor_upper']);mr=F(bud['external_factor_radius']);ec=F(bud['coefficient_radius_upper']);beta=F(bud['financial_transfer_upper'])
            det=ec+beta+mass*eta+mr
            bud.update({'arithmetic_before_mass_upper':cert['arithmetic_probability_sandwich_upper'],
                'auxiliary_components':{'inner_normal_interval_bias_upper':cert['inner_normal_interval_bias_upper'],
                    'strike_midpoint_bias_upper':cert['strike_midpoint_bias_upper'],
                    'bias_upper':str(F(cert['inner_normal_interval_bias_upper'])+F(cert['strike_midpoint_bias_upper']))},
                'total_deterministic_undiscounted_price_error_upper':str(det)})
            previous=old['coefficient_informed_prior'];lo=F(previous['finite_ideal_mean_lower']);hi=F(previous['finite_ideal_mean_upper'])
            prior={'probability_lower':str(max(F(0),lo-eta)),'probability_upper':str(min(F(1),hi+eta)),
                'native_positive':True,'finite_ideal_mean_lower':str(lo),'finite_ideal_mean_upper':str(hi),
                'information_scope':'Existing price-independent coefficient/curvature/U2 finite-mean interval, enlarged only by new conditional probability error.',
                'previous_mean_information_sha256':digest(previous)}
            row.update(complete_deterministic_budget=bud,coefficient_informed_prior=prior)
            if det>=F(1,100):row['rejection']={'reason':'Deterministic monetary allowance exhausted','deterministic_error_upper':str(det)}
            else:
                for mode in protocol['architectures']:
                    profile=counted_profile(g,tables,cache,mode);row['profiles'][mode]={'centred':profile}
                    strategy=informed.informed_complete_plan('positive',bud,cert,prior,row['profiles'][mode])
                    row['strategies'].append({'architecture':mode,'result':strategy})
                row['selected_strategy']=min(row['strategies'],key=lambda x:(x['result']['T_count'],x['architecture']))
        else:row['rejection']={'reason':cert['reason']}
        rows.append(row)
        selected=row.get('selected_strategy');print(original,candidate,cert['status'],
              'eta',cert.get('probability_error_upper'),'T',None if selected is None else selected['result']['T_count'],flush=True)
    archive_path=RESULTS/'shared_threshold_scalar_ordinates.json.gz'
    if archive_path.exists():raise FileExistsError('Preserve scalar archive')
    with archive_path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as stream:
            stream.write(json.dumps(archive,separators=(',',':'),default=str).encode())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Sources or inputs changed during execution')
    write(RESULTS/'shared_threshold_resources.json',{'schema':'qbp-shared-threshold-resources-v1',
       'status':'Completed, retaining every failure','protocol_sha256':sha(RESULTS/'shared_threshold_protocol.json'),
       **signatures(),'rows':rows,'scalar_archive_file':archive_path.name,'scalar_archive_sha256':sha(archive_path),
       'unique_scalar_certificates':len(archive),'emitted_block_profiles':len(cache),
       'generation_seconds_diagnostic':time.perf_counter()-start,
       'scope':'Full coherent financial preparations, conditional probability and full cent-price/QPE resource certificates. No quantum price outcomes or physical speed-up.'})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run'));p.add_argument('--repo',required=True)
    p.add_argument('--factor-resources',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    configure(a.repo,a.factor_resources,a.output)
    register() if a.mode=='register' else run()
