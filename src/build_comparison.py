"""Charge unchanged published envelopes and my separately proved controller.

The imported financial compiler and resource charger are used read-only.
I retain all source candidate prefixes rather than select a favourable case.
"""
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
import sys
from full_batch_miqae import plan,rational_record
from published_envelopes import theorem_envelopes

ROOT=Path(__file__).resolve().parent


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


@lru_cache(maxsize=2000)
def cached_envelopes(epsilon,alpha,prior):return theorem_envelopes(epsilon,alpha,prior)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',type=Path,default=ROOT)
    parser.add_argument('--config',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'results/published_comparison.json')
    args=parser.parse_args();repo=args.repo.resolve()
    config=args.config or repo/'results/architecture_resource_candidates.json'
    source_paths={'complete_estimation_study.py':repo/'complete_estimation_study.py',
                  'estimation_prior.py':repo/'estimation_prior.py',
                  'results/'+config.name:config}
    before={name:hashlib.sha256(p.read_bytes()).hexdigest() for name,p in source_paths.items()}
    sys.path.insert(0,str(repo))
    charge=load_module(repo/'complete_estimation_study.py','qbp_read_only_charge')
    prior_provider=load_module(repo/'estimation_prior.py','qbp_read_only_prior')
    data=json.loads(config.read_text());models={m['model_id']:m for m in data['models']}
    records=[]
    for index,row in enumerate(data['configurations']):
        prior=F(Decimal(prior_provider.prior_for_row(models[row['model_id']],row)['finite_mean_probability_upper_decimal']))
        output={'case_id':row['case_id'],'model_id':row['model_id'],'method':row['method'],
                'candidate_id':row.get('candidate_id'),'architecture':row.get('architecture',row['resources'].get('architecture')),
                'cap':row['cap'],'prefix':{k:row[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits')},
                'proved_finite_probability_upper':rational_record(prior),'confidence_profiles':{}}
        for scope,alpha in (('one_price_99_per_cent',F(1,100)),('selected_27_simultaneous_99_per_cent',F(1,2700))):
            profiles={}
            for label,bound in (('financial_prior',prior),('full_range',F(1))):
                setup=plan(F(1,200*row['cap']),alpha,bound)
                charge.plan=lambda e,a,p,retained=setup:retained
                charge.decimal_record=rational_record
                own=charge.architecture_budget(row,alpha,bound)
                own['controller_label']='my certified fixed-radius full-batch MIQAE-informed implementation; published62 coefficient is not applied'
                own['coupled_linear_T_cost_upper']=setup.cost_upper(row['resources']['logical_T_count'],row['resources']['grover_iterate_T_count'])
                envelopes=cached_envelopes(F(1,200*row['cap']),alpha,bound)
                charged={}
                for key,caps in envelopes.items():
                    if key=='financial_postprocessing':continue
                    setup_source=SimpleNamespace(**caps)
                    charge.plan=lambda e,a,p,retained=setup_source:retained
                    charge.decimal_record=lambda value:rational_record(vars(value))
                    charged[key]=charge.architecture_budget(row,alpha,bound)
                    charged[key]['published_scope']=caps['algorithm_scope'] if 'algorithm_scope' in caps else 'zero-shot proved prior'
                profiles[label]={'certified_full_batch_variant':own,'published_theorem_envelopes':charged,
                    'published_financial_postprocessing':rational_record(envelopes['financial_postprocessing'])}
            output['confidence_profiles'][scope]=profiles
        records.append(output)
        if (index+1)%200==0:print(json.dumps({'charged_candidates':index+1}),flush=True)
    for name,p in source_paths.items():assert hashlib.sha256(p.read_bytes()).hexdigest()==before[name],f'Source changed during retained run: {name}'
    selected=[]
    for architecture in sorted(set(r['architecture'] for r in records)):
        for case in sorted(set(r['case_id'] for r in records)):
            pool=[r for r in records if r['architecture']==architecture and r['case_id']==case]
            profiles={}
            for scope in ('one_price_99_per_cent','selected_27_simultaneous_99_per_cent'):
                for comparator in ('certified_full_batch_variant','iqae_published','miqae_full_batch_integer_corollary'):
                    def value(r):
                        profile=r['confidence_profiles'][scope]['financial_prior']
                        return profile['certified_full_batch_variant']['T_count_upper'] if comparator=='certified_full_batch_variant' else profile['published_theorem_envelopes'][comparator]['T_count_upper']
                    minimum=min(value(r) for r in pool)
                    profiles[scope+'__'+comparator]={'minimum_T_count_upper':minimum,
                        'tied_methods':sorted(set(r['method'] for r in pool if value(r)==minimum)),
                        'representative_prefix':next(r['prefix'] for r in pool if value(r)==minimum)}
            selected.append({'architecture':architecture,'case_id':case,'profiles':profiles})
    result={'status':'complete logical budgets over retained finite compiler pool; source theorem envelopes separated from my certified numerical variant',
        'source_hashes':before,'own_source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'full_batch_miqae.py',ROOT/'certified_scalar.py',ROOT/'published_envelopes.py',Path(__file__))},
        'confidence_scope':'one price99 per cent; selected27 profile jointly99 per cent for one preselected method/architecture percase; not all candidates simultaneously',
        'scope_limit':'ideal coherent finite oracle; source IQAE envelope assumes exact published postprocessing; final MIQAE chapter text inaccessible and authorv4 actually consulted; no hardware runtime advantage claim',
        'configurations':records,'selected_cases':selected}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'charged_candidates':len(records),'selected_case_architectures':len(selected),'output':str(args.output)},indent=2))


if __name__=='__main__':main()
