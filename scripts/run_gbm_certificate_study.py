"""Retain supplied-interval and rounded synthetic Gaussian diagnostics.

No empirical data are used; numerical references are not interval proofs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
from scipy.stats import chi2
import scipy



def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def fraction_bounds(record):
    return (Fraction(record['lower_exact_fraction']), Fraction(record['upper_exact_fraction']))


def exact_contains(record, value):
    lo, hi = fraction_bounds(record)
    return mp_number(lo) <= value <= mp_number(hi)


def serial(value):
    return None if value is None else mp.nstr(value, 175)


def prepare_output_directory(output, project_root):
    """Refuse repository writes and existing evidence before creating anything."""
    output, project_root = Path(output).resolve(), Path(project_root).resolve()
    if output == project_root or project_root in output.parents:
        raise ValueError('Output must be external to the resolved project root.')
    if output.exists():
        if not output.is_dir():
            raise ValueError('Output must be a new or empty directory, not an existing file.')
        if any(output.iterdir()):
            raise ValueError('Use a new or empty output directory; retained evidence is not overwritten.')
    output.mkdir(parents=True, exist_ok=True)
    return output


def numerical_verification_passes(manifest):
    """Statistical containment frequencies are deliberately not pass criteria."""
    if not all(manifest[key] for key in ('deterministic_reference_checks_pass',
                                        'quantile_cdf_checks_pass', 'sources_unchanged')):
        return False
    for record in manifest['synthetic_summary']:
        counts = record['counts']
        if (counts['unexpected_error'] != 0
                or counts['rounding_bound_checks'] != counts['attempts']
                or counts['reference_image_checks'] != counts['certified']
                or counts['certified'] + counts['abstained'] != counts['attempts']):
            return False
    return True


def run(args):
    output = prepare_output_directory(args.output, args.project_root)
    module_path = args.module.resolve()
    module = load_module(module_path, 'gbm_price_certificate_study_engine')
    here = Path(__file__).resolve()
    references = module_path.with_name('gbm_certificate_references.py')
    reference_module = load_module(references, 'gbm_certificate_independent_references')
    for name in ('capped_price','chi_square_cdf','digital_integral','exact_composition_example',
                 'exact_decimal_scatter','fixture_cases','full_image','mp_number',
                 'rounded_decimal_observations'):
        globals()[name] = getattr(reference_module,name)
    source_paths = [module_path, here, references]
    hashes_before = {str(path):sha(path) for path in source_paths}
    protocol = {
      'schema':'qbp-gbm-certificate-protocol-v1',
      'source_hashes_before':hashes_before,
      'study_type':'fixed synthetic numerical diagnostic; no empirical sample or exact RNG coverage claim',
      'seed':7362018, 'repetitions_per_case_size':args.repetitions,
      'sample_sizes':[25,60], 'observation_interval_exact':'1/252',
      'rounding':'exact nearest decimal6 with radius1/2000000 for each observation',
      'alpha_exact':'1/20', 'absolute_price_width_relative_cap':'1e-22',
      'standalone_quantile_width_exact':'1e-32', 'bridge_quantile_width_exact':'1e-24', 'reference_digits':[100,160],
      'dense_scan_points':1025,
      'comparison_design':'Full image vs endpoint-only and dense scan; SciPy plug-in endpoints uncertified.',
      'selection_record':'Adverse inputs chosen after inspecting existing binary64 K+H and extremum formulas; informed cases, not prospective random discovery.',
      'coverage_denominator':'Every attempted dataset, including abstentions and unexpected errors.',
      'input_semantics':'Decimal strings and fractions exact; binary64 inputs denote their exact stored value.',
      'environment':{'python':platform.python_version(),'numpy':np.__version__,
          'scipy':scipy.__version__,'mpmath':mp.__version__},
    }
    write(output/'protocol.json', protocol)
    print('Protocol retained.', flush=True)
    deterministic = []
    with mp.workdps(180):
        for case in fixture_cases():
            contract = module.ExactCappedGBM(**case['contract'])
            tolerance = Fraction(case['contract']['cap']) * Fraction('1e-22')
            try:
                cert = module.certified_price_image(contract,case['lower'],case['upper'],
                            absolute_tolerance=tolerance)
            except Exception as error:
                cert = {'status':'abstained', 'error':{'type':type(error).__name__,
                        'message':str(error)}}
            refs = {}
            for digits in (100,160):
                with mp.workdps(digits):
                    refs[digits] = full_image(case['contract'],case['lower'],case['upper'])
            ref = refs[160]
            with mp.workdps(100):
                integral_values = [digital_integral(case['contract'],case['lower']),
                                   digital_integral(case['contract'],case['upper'])]
                if ref['interior']:
                    integral_values.append(digital_integral(case['contract'],ref['peak']))
            with mp.workdps(80):
                lo, hi = mp_number(case['lower']), mp_number(case['upper'])
                grid = [capped_price(case['contract'],lo+(hi-lo)*j/1024) for j in range(1025)]
            with mp.workdps(160):
                reference_candidates = [capped_price(case['contract'],case['lower']),
                                        capped_price(case['contract'],case['upper'])]
                if ref['interior']:
                    reference_candidates.append(capped_price(case['contract'],ref['peak']))
                integral_relative_errors = [abs(value-target)/abs(target) if target else abs(value-target)
                    for value,target in zip(integral_values,reference_candidates)]
            image = cert.get('image_outer')
            lower_ok = image is not None and exact_contains(cert['minimum_price'], ref['lower'])
            upper_ok = image is not None and exact_contains(cert['maximum_price'], ref['upper'])
            diagnostic = {'case':case, 'certificate':cert,
                'reference160':{key:serial(value) if key!='interior' else value
                                for key,value in ref.items()},
                'reference100_difference':{'lower':serial(ref['lower']-refs[100]['lower']),
                                          'upper':serial(ref['upper']-refs[100]['upper'])},
                'digital_integral_values100':[serial(value) for value in integral_values],
                'digital_integral_relative_errors':[serial(value) for value in integral_relative_errors],
                'dense_scan':{'points':1025,'minimum':serial(min(grid)),'maximum':serial(max(grid)),
                              'maximum_gap_to_true_supremum':serial(ref['upper']-max(grid)),
                              'classification':'attained sample values, not a rigorous optimiser'},
                'endpoint_only_gap':serial(ref['upper']-ref['endpoint_only_upper']),
                'checks':{'minimum_encloses_reference':lower_ok,
                          'maximum_encloses_reference':upper_ok,
                          'sampled_prices_enclosed':image is not None and all(exact_contains(image,x) for x in grid),
                          'digital_integral_relative_error_below_1e_70':max(integral_relative_errors)<=mp.mpf('1e-70')}}
            deterministic.append(diagnostic)
            print('Deterministic:',case['name'],diagnostic['checks'],flush=True)
    write(output/'deterministic_images.json', deterministic)
    quantile_rows=[]
    with mp.workdps(160):
        for nu in (1,24,59):
            for p in ('1e-12','0.025','0.975','0.999999999999'):
                try:
                    cert=module.certified_chi_square_quantile(nu,p,absolute_tolerance='1e-32')
                except Exception as error:
                    quantile_rows.append({'degrees':nu,'probability_exact':p,
                        'status':'abstained','error':{'type':type(error).__name__,'message':str(error)},
                        'cdf_brackets_probability':False})
                    continue
                lo,hi=fraction_bounds(cert['quantile'])
                flo,fhi=chi_square_cdf(lo,nu),chi_square_cdf(hi,nu)
                quantile_rows.append({'degrees':nu,'probability_exact':p,'certificate':cert,
                    'mpmath_cdf_lower160':serial(flo),'mpmath_cdf_upper160':serial(fhi),
                    'cdf_brackets_probability':flo<=mp_number(p)<=fhi})
    write(output/'quantile_brackets.json',quantile_rows)
    print('Quantile brackets retained.',flush=True)
    base={'spot':'100','strike':'100','cap':'20','maturity':'1','rate':'0.03','dividend':'0.01'}
    with mp.workdps(160):
        peak=full_image(base,'0.2','0.6')['peak']
        near_peak=mp.nstr(peak,80)
    cases=[('ordinary',base,'0.2'),('near_stationary',base,near_peak),
           ('monotone_decreasing',{**base,'spot':'160'},'0.3')]
    summary=[]
    delta=1/252
    alpha=0.05
    for case_index,(name,contract,sigma_text) in enumerate(cases):
        with mp.workdps(160):
            true_price=capped_price(contract,sigma_text)
        sigma=float(sigma_text)
        for n in (25,60):
            counts=dict(attempts=args.repetitions,certified=0,abstained=0,unexpected_error=0,
                true_volatility_contained=0,true_price_contained=0,
                scipy_parameter_contained=0,scipy_full_image_contained=0,
                scipy_endpoint_only_contained=0,rounding_bound_checks=0,
                reference_image_checks=0)
            records=[]
            generator=np.random.default_rng(np.random.SeedSequence([7362018,case_index,n]))
            for repetition in range(args.repetitions):
                draws=generator.normal(0.04*delta,sigma*math.sqrt(delta),size=n)
                observations,error2=rounded_decimal_observations(draws,6)
                rounding_ok=error2<=n*Fraction(1,2_000_000)**2
                counts['rounding_bound_checks']+=int(rounding_ok)
                scatter=exact_decimal_scatter(observations)
                values=np.array([float(x) for x in observations])
                variance=float(np.var(values,ddof=1))/delta
                slo=math.sqrt((n-1)*variance/chi2.ppf(1-alpha/2,n-1))
                shi=math.sqrt((n-1)*variance/chi2.ppf(alpha/2,n-1))
                with mp.workdps(160):
                    simple=full_image(contract,slo,shi)
                    endpoints_lower=simple['lower']
                    counts['scipy_parameter_contained']+=int(slo<=sigma<=shi)
                    counts['scipy_full_image_contained']+=int(simple['lower']<=true_price<=simple['upper'])
                    counts['scipy_endpoint_only_contained']+=int(endpoints_lower<=true_price<=simple['endpoint_only_upper'])
                row={'repetition':repetition,'scatter_exact':str(scatter),
                     'rounding_error_squared_exact':str(error2),'rounding_bound_holds':rounding_ok,
                     'scipy_volatility_interval_uncertified':[slo,shi]}
                try:
                    cert=module.certified_gaussian_price_interval(module.ExactCappedGBM(**contract),
                        observations,'1/252',alpha='1/20',
                        rounding_radii=['1/2000000']*n,quantile_tolerance='1e-24',
                        absolute_tolerance=Fraction(contract['cap'])*Fraction('1e-22'))
                    row['certificate']=cert
                    if cert.get('status','certified')!='certified':
                        counts['abstained']+=1
                    else:
                        counts['certified']+=1
                        lo,hi=map(Fraction,cert['volatility_interval_exact'])
                        sigma_exact=Fraction(sigma_text)
                        counts['true_volatility_contained']+=int(lo<=sigma_exact<=hi)
                        with mp.workdps(160):
                            image=cert['price_image']
                            counts['true_price_contained']+=int(exact_contains(image['image_outer'],true_price))
                            reference=full_image(contract,lo,hi)
                            image_ok=(exact_contains(image['minimum_price'],reference['lower'])
                                      and exact_contains(image['maximum_price'],reference['upper']))
                        counts['reference_image_checks']+=int(image_ok)
                        row['independent_image_check']=image_ok
                except Exception as error:
                    counts['unexpected_error']+=1
                    row['error']={'type':type(error).__name__,'message':str(error)}
                records.append(row)
            summary.append({'case':name,'contract':contract,'true_sigma_exact':sigma_text,
                'observations':n,'counts':counts,'denominator':args.repetitions,
                'records':records,'interpretation':'Numerical coverage diagnostic under approximate RNG, not exact Gaussian sampling evidence.'})
            print('Synthetic:',name,n,counts,flush=True)
    write(output/'synthetic_bridge.json',summary)
    hashes_after={str(path):sha(path) for path in source_paths}
    manifest={'schema':'qbp-gbm-certificate-validation-v1',
      'deterministic_cases':len(deterministic),'quantile_brackets':len(quantile_rows),
      'synthetic_datasets':sum(x['counts']['attempts'] for x in summary),
      'deterministic_reference_checks_pass':all(all(x['checks'].values()) for x in deterministic),
      'quantile_cdf_checks_pass':all(x['cdf_brackets_probability'] for x in quantile_rows),
      'synthetic_summary':[{key:value for key,value in x.items() if key!='records'} for x in summary],
      'exact_composition_example':exact_composition_example(),
      'source_hashes_before':hashes_before,'source_hashes_after':hashes_after,
      'sources_unchanged':hashes_before==hashes_after,
      'source_consultation':'Local calibration module and companion image/composition proofs; no new publication consultation.',
      'certification_scope':'Exact declared arithmetic and GBM model, conditional on successful result. Supplied intervals do not certify own statistical coverage; rounded-data theorem assumes ideal iid Gaussian values within supplied radii. RNG draws approximate that law.',
    }
    manifest['overall_numerical_verification'] = numerical_verification_passes(manifest)
    write(output/'validation.json',manifest)
    print(json.dumps({key:value for key,value in manifest.items() if key not in ('synthetic_summary','source_hashes_before','source_hashes_after')},indent=2),flush=True)
    return manifest['overall_numerical_verification']


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--module',type=Path,help='Optional explicit engine source; defaults to project-root/src/gbm_price_certificate.py.')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repetitions',type=int,default=30)
    args=parser.parse_args(argv)
    if args.module is None:
        args.module=args.project_root/'src'/'gbm_price_certificate.py'
    if not 1<=args.repetitions<=1000:
        parser.error('repetitions must be in 1..1000')
    return 0 if run(args) else 1


if __name__=='__main__':
    raise SystemExit(main())
