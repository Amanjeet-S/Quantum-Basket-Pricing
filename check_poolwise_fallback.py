"""Independent all-pool fallback comparison without a charging frontend."""
from pathlib import Path
from fractions import Fraction as F
import argparse, copy, gzip, hashlib, importlib.util, json, sys
sys.dont_write_bytecode=True
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v): return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p): return json.loads(gzip.decompress(Path(p).read_bytes()) if str(p).endswith('.gz') else Path(p).read_bytes())

def main(args):
    repo,results=Path(args.repo),Path(args.results)
    panelpath=results/'common_shift_poolwise_fallback.json'
    if not panelpath.exists(): panelpath=results/'common_shift_poolwise_fallback.json.gz'
    panel=read(panelpath);data=read(results/'common_shift_resources.json.gz')
    assert panel['common_shift_resources_sha256']==sha(results/'common_shift_resources.json.gz')
    assert panel['protocol_sha256']==sha(results/'common_shift_poolwise_fallback_protocol.json')
    spec=importlib.util.spec_from_file_location('retained_own_formula',repo/'check_expanded_precision_resources.py')
    own=importlib.util.module_from_spec(spec);spec.loader.exec_module(own);own.ALPHA=F(1,7200)
    from flint import arb,ctx
    ctx.prec=512
    rotations=read(repo/'results/qft_rotation_tapes.json')['rotations']
    bounds={x['key']:F(x['operator_norm_upper_rational']) for x in read(repo/'results/qft_rotation_certificates.json')['certificates']}
    seen=set();pairs=0
    for r in panel['records']:
        row=data['records'][r['record_index']]
        assert digest(row)==r['record_sha256']
        assert (r['record_index'],r['architecture']) not in seen
        seen.add((r['record_index'],r['architecture']))
        current=next(s['result'] for s in row['strategies'] if s['architecture']==r['architecture'])
        assert current['kind']=='complete_quantum_resource_plan'
        flag=current['flags'][0];plan=r['half_plan'];oldplan=flag['plan'];qft=flag['fourier_profile']
        assert F(plan['alpha'])==F(1,7200) and plan['repetitions']==27 and oldplan['repetitions']==25
        for k in oldplan:
            if k in ['alpha','repetitions','binomial_failure_upper','hoeffding_repetitions_upper']: continue
            a,b=oldplan[k],plan[k]
            try: assert F(a)==F(b),k
            except (ValueError,TypeError): assert a==b,k
        stat,failure=own.own_plan(plan,F(r['probability_tolerance']),F(r['prior_upper']),qft,rotations,bounds)
        gap=F(plan['good_probability_lower'])-F(1,2)
        gap_ball=arb(gap.numerator)/gap.denominator
        diagnostic=arb(7200).log()/(2*gap_ball**2)
        lower_ceiling=diagnostic.lower().ceil().unique_fmpz()
        upper_ceiling=diagnostic.upper().ceil().unique_fmpz()
        assert lower_ceiling==upper_ceiling
        odd=max(1,int(upper_ceiling));odd+=odd%2==0
        assert plan['hoeffding_repetitions_upper']==odd
        assert failure==F(r['half_failure_probability_upper'])
        assert digest(plan)==r['half_plan_sha256'] and digest(qft)==r['half_fourier_profile_sha256']
        original_charge=flag['complete_charge'];half=copy.deepcopy(original_charge)
        half['repetitions']=27
        half['metrics']={k:(v//25)*27 for k,v in original_charge['metrics'].items()}
        assert all(v%25==0 for v in original_charge['metrics'].values())
        for k in ['controlled_Q_calls','phase_register_measurements','complete_reinitialisations']:
            assert original_charge[k]%25==0;half[k]=original_charge[k]//25*27
        assert digest(half)==r['half_full_charge_sha256']
        assert r['T_full_alpha']==current['T_count']==original_charge['metrics']['T']
        assert r['T_half_alpha_fallback']==half['metrics']['T']==r['T_full_alpha']//25*27
        assert r['pilot_T_budget']==r['T_half_alpha_fallback']//10
        assert r['T_worst_history_upper']==r['T_half_alpha_fallback']+r['pilot_T_budget']
        budget=row['complete_deterministic_budget'];monetary=F(budget['total_deterministic_undiscounted_price_error_upper'])+F(budget['external_factor_upper'])*stat
        assert monetary==F(r['complete_monetary_error_upper'])<=F(1,100)
        assert r['M_full']==r['M_half']==plan['M'] and r['QFT_unchanged'] and r['per_repetition_charge_unchanged']
        pairs+=1
    assert pairs==808
    assert len(panel['rejections'])==12 and all(x['method']=='numerical_conditional_price' for x in panel['rejections'])
    for x in panel['rejections']:
        row=data['records'][x['record_index']];assert digest(row)==x['record_sha256']
        assert not row['strategies'] and row['rejection']==x['rejection']
    assert all(panel['invariance'].values())
    for c in panel['selected_comparisons']:
        for method,choice in c['methods'].items():
            group=[r for r in panel['records'] if r['model_id']==c['model_id'] and r['method']==method]
            best=min(r['T_worst_history_upper'] for r in group)
            assert choice['T_worst_history_upper']==best
            selected=next(r for r in group if digest(r)==choice['selected_record_sha256'])
            assert selected['T_worst_history_upper']==best
            assert min(r['T_full_alpha'] for r in group)==selected['T_full_alpha']
        rival=min(c['methods'][x]['T_worst_history_upper'] for x in ['raw_payoff','numerical_conditional_price'])
        event=c['methods']['shared_threshold']['T_worst_history_upper']
        assert F(c['saving_against_strongest_exact'])==F(rival-event,rival)
    out={'schema':'qbp-poolwise-fallback-root-audit-v1','status':'passed','auditor_sha256':sha(__file__),
         'independent_library_sha256':sha(repo/'check_expanded_precision_resources.py'),
         'panel_sha256':sha(panelpath),'amended_resources_sha256':sha(results/'common_shift_resources.json.gz'),
         'admitted_architecture_plans':pairs,'retained_rejections':12,'complete_charge_hash_reconstructions':pairs,
         'fallback_repetition_change':'25 to27; M,QFT andper-repetitionchargeunchanged',
         'protected_cap_formula':'F(T)=27T/25+floor(27T/250) for every admitted canonical capT',
         'all_pool_minima_verified':True,'source_or_confidence_frontend_imported':False,
         'scope':'Independently checked exact QFT/phase/binomial minima, full existing charge scaling hashes, monetary guarantees, rejections and protected-cap minimisation over all416settings/twoarchitectures. Caps do not certify adaptive expected costs or an optimal quantum algorithm.'}
    p=Path(args.output);assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print('Independent all-pool audit passed:808fallbacks andcompletechargehashes;12retainedrejects.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--results',required=True);p.add_argument('--output',required=True)
    main(p.parse_args())
