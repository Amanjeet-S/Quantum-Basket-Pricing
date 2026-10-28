"""Independent mathematical/regression checks for exact price-image arithmetic."""
from fractions import Fraction as F
import importlib.util
import copy
import itertools
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import mpmath as mp

HERE = Path(__file__).resolve().parent
SOURCE = Path(os.environ.get('QBP_CERTIFICATE_MODULE', str(HERE.parents[0] / 'src' / 'gbm_price_certificate.py')))
spec = importlib.util.spec_from_file_location('gbm_price_certificate_tests_engine', SOURCE)
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)


def bounds(record):
    return F(record['lower_exact_fraction']), F(record['upper_exact_fraction'])


def mpf(f):
    value = F(f)
    return mp.mpf(value.numerator)/value.denominator


def price(c, sig):
    with mp.workdps(160):
        s,k,h,t,r,q=map(mpf,(c.spot,c.strike,c.cap,c.maturity,c.rate,c.dividend))
        v=mpf(sig)*mp.sqrt(t)
        phi=lambda z:mp.erfc(-z/mp.sqrt(2))/2
        def call(strike):
            d=(mp.log(s/strike)+(r-q)*t+v*v/2)/v
            return s*mp.exp(-q*t)*phi(d)-strike*mp.exp(-r*t)*phi(d-v)
        return call(k)-call(k+h)


def scatter(xs):
    values=list(map(F,xs))
    return sum((x*x for x in values),F())-sum(values,F())**2/len(values)


class CertificateTests(unittest.TestCase):
    def assertEncloses(self, interval, point):
        lo,hi=bounds(interval)
        with mp.workdps(180):
            self.assertLessEqual(mpf(lo),point)
            self.assertGreaterEqual(mpf(hi),point)

    def test_exact_decimal_and_binary_input_semantics(self):
        self.assertEqual(engine.rational('0.1'),F(1,10))
        self.assertEqual(engine.rational(0.1),F.from_float(0.1))
        self.assertNotEqual(engine.rational('0.1'),engine.rational(0.1))
        c=engine.ExactCappedGBM(cap='1e-15')
        self.assertEqual(c.strike+c.cap,F(100)+F(1,10**15))

    def test_exact_scatter_translation_and_scale(self):
        xs=['0.000007','-0.000002','0.000009','0.000003']
        expected=scatter(xs)
        result=engine.exact_scatter_enclosure(xs)
        self.assertEqual(F(result['observed_scatter_exact']),expected)
        self.assertEqual(bounds(result['latent_scatter']),(expected,expected))
        shifted=[F(x)+10**20 for x in xs]
        self.assertEqual(F(engine.exact_scatter_enclosure(shifted)['observed_scatter_exact']),expected)
        self.assertEqual(F(engine.exact_scatter_enclosure([7*F(x) for x in xs])['observed_scatter_exact']),49*expected)

    def test_rounded_scatter_encloses_all_vertices(self):
        xs=[F('-0.2'),F('0.1'),F('0.4')]
        radii=[F('0.01'),F('0.02'),F('0.03')]
        record=engine.exact_scatter_enclosure(xs,radii)
        lo,hi=bounds(record['latent_scatter'])
        self.assertEqual(F(record['rounding_l2_radius_squared_exact']),sum((x*x for x in radii),F()))
        for signs in itertools.product((-1,1),repeat=3):
            value=scatter([x+sign*r for x,sign,r in zip(xs,signs,radii)])
            self.assertLessEqual(lo,value)
            self.assertGreaterEqual(hi,value)
        self.assertLessEqual(lo,scatter(xs))
        self.assertGreaterEqual(hi,scatter(xs))

    def test_constant_sample_abstains_with_all_attempts_counted(self):
        result=engine.certified_gaussian_price_interval(engine.ExactCappedGBM(),['0.1']*5,'1/252')
        self.assertEqual(result['status'],'abstained')
        self.assertNotIn('price_image',result)
        self.assertIn('scatter',result)

    def test_large_rounding_radius_abstains(self):
        result=engine.certified_gaussian_price_interval(engine.ExactCappedGBM(),['0','0.001'],'1/252',rounding_radii=['1','1'])
        self.assertEqual(result['status'],'abstained')
        self.assertEqual(F(result['scatter']['latent_scatter']['lower_exact_fraction']),0)

    def test_binary_rounding_radius_adjacent_half_spacings(self):
        for value in (0.0,-0.0,1.0,-1.0,math.nextafter(0,1),sys.float_info.max,-sys.float_info.max):
            with self.subTest(value=value):
                radius=engine.binary64_rounding_radius(value)
                self.assertGreater(radius,0)
                for neighbour in (math.nextafter(value,-math.inf),math.nextafter(value,math.inf)):
                    if math.isfinite(neighbour):
                        self.assertGreaterEqual(radius,abs(F.from_float(neighbour)-F.from_float(value))/2)

    def test_quantile_matches_exact_exponential_case(self):
        with mp.workdps(160):
            for p in ('0.025','0.5','0.975'):
                record=engine.certified_chi_square_quantile(2,p,absolute_tolerance='1e-32')
                expected=-2*mp.log1p(-mpf(p))
                self.assertEncloses(record['quantile'],expected)
                lo,hi=bounds(record['quantile'])
                self.assertLessEqual(hi-lo,F('1e-32'))
                self.assertLess(1-mp.exp(-mpf(lo)/2),mpf(p))
                self.assertGreater(1-mp.exp(-mpf(hi)/2),mpf(p))
                self.assertLess(bounds(record['cdf_at_lower'])[1],F(p))
                self.assertGreater(bounds(record['cdf_at_upper'])[0],F(p))

    def test_nearly_equal_cdf_decision_retains_strict_rational_witness(self):
        p=('0.63212055882855767840447622983853913255418886896823216549216319830253850425510019664285272565408035625337267472315600479175302420720987099137334641050590121690781')
        record=engine.certified_chi_square_quantile(2,p)
        self.assertGreaterEqual(record['precision_bits'],512)
        self.assertLess(bounds(record['cdf_at_lower'])[1],F(p))
        self.assertGreater(bounds(record['cdf_at_upper'])[0],F(p))
        self.assertLessEqual(bounds(record['quantile'])[1]-bounds(record['quantile'])[0],F('1e-24'))

    def test_quantile_increases_in_probability(self):
        low=engine.certified_chi_square_quantile(24,'0.025')['quantile']
        high=engine.certified_chi_square_quantile(24,'0.975')['quantile']
        self.assertLess(bounds(low)[1],bounds(high)[0])

    def test_scalar_price_reference(self):
        c=engine.ExactCappedGBM()
        for sigma in ('0.01','0.2','0.6','2'):
            self.assertEncloses(engine.certified_price(c,sigma)['price'],price(c,sigma))

    def test_nonmonotone_image_requires_interior_maximum(self):
        c=engine.ExactCappedGBM()
        record=engine.certified_price_image(c,'0.2','0.6',absolute_tolerance='1e-28')
        with mp.workdps(160):
            critical=(mp.log(mpf(c.strike)*mpf(c.strike+c.cap)/mpf(c.spot)**2)/mpf(c.maturity)-2*mpf(c.rate-c.dividend))
            peak=mp.sqrt(critical)
            top=price(c,F(str(peak)))
            self.assertEncloses(record['maximum_price'],top)
            self.assertGreater(top,max(price(c,'0.2'),price(c,'0.6')))
        self.assertEqual(record['statistical_coverage_claim'],False)

    def test_tiny_cap_refinement_regression(self):
        c=engine.ExactCappedGBM(cap='1e-15',rate=0,dividend=0)
        record=engine.certified_price_image(c,'1e-9','6e-9',absolute_tolerance='1e-37')
        with mp.workdps(160):
            peak=mp.sqrt(mp.log1p(mp.mpf('1e-17')))
            self.assertEncloses(record['maximum_price'],price(c,F(str(peak))))
        self.assertLessEqual(bounds(record['maximum_price'])[1]-bounds(record['maximum_price'])[0],F('1e-37'))

    def test_singleton_and_monotone_image(self):
        c=engine.ExactCappedGBM()
        single=engine.certified_price_image(c,'0.2','0.2')
        self.assertEncloses(single['minimum_price'],price(c,'0.2'))
        self.assertEncloses(single['maximum_price'],price(c,'0.2'))
        c=engine.ExactCappedGBM(spot=160)
        decreasing=engine.certified_price_image(c,'0.2','0.6')
        self.assertEncloses(decreasing['minimum_price'],price(c,'0.6'))
        self.assertEncloses(decreasing['maximum_price'],price(c,'0.2'))

    def test_resource_exhaustion_is_explicit(self):
        with self.assertRaises(ArithmeticError):
            engine.certified_chi_square_quantile(24,'0.025',max_iterations=1)
        result=engine.certified_gaussian_price_interval(engine.ExactCappedGBM(),['-0.01','0.01','0.02'],'1/252',max_iterations=1)
        self.assertEqual(result['status'],'abstained')
        self.assertIn('budget',result['abstention_reason'])
        with self.assertRaises(ArithmeticError):
            engine.certified_price(engine.ExactCappedGBM(),'0.2',absolute_tolerance='1e-100',start_bits=64,max_bits=64)

    def test_bad_inputs_rejected(self):
        for value in (True,math.inf,math.nan,'NaN','Infinity'):
            with self.subTest(value=value),self.assertRaises(ValueError):
                engine.ExactCappedGBM(spot=value)
        for lo,hi in (('0','0.2'),('-0.1','0.2'),('0.3','0.2')):
            with self.assertRaises(ValueError):
                engine.certified_price_image(engine.ExactCappedGBM(),lo,hi)
        for observations,radii in ((['0'],None),(['0','1'],[-1,0]),(['0','1'],[0])):
            with self.assertRaises(ValueError):
                engine.exact_scatter_enclosure(observations,radii)

    def test_composition_does_not_require_independence(self):
        alpha,beta=F(1,10),F(1,20)
        # Failure probabilities depend on the dataset; each conditional budget
        # remains at most beta. This deliberately is a dependent finite law.
        rows=[(alpha,False,F(0)),(1-alpha,True,beta)]
        numerical_failure=sum((p*fail for p,_,fail in rows),F())
        simultaneous_success=sum((p*(1-fail) for p,parameter_ok,fail in rows if parameter_ok),F())
        self.assertLessEqual(numerical_failure,beta)
        self.assertGreaterEqual(simultaneous_success,1-alpha-beta)
        self.assertNotEqual(F(0),beta)  # data-dependent conditional probabilities

    def test_abstention_does_not_allow_conditional_coverage_claim(self):
        alpha=F(1,20)
        # A valid procedure may return only on the parameter-failure event.
        success_probability=alpha
        success_and_miss=alpha
        self.assertLessEqual(success_and_miss,alpha)
        self.assertEqual(success_and_miss/success_probability,1)
        self.assertGreater(success_and_miss/success_probability,alpha)


class StudyRunnerGuardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(os.environ.get('QBP_STUDY_RUNNER',str(HERE.parents[0]/'scripts'/'run_gbm_certificate_study.py')))
        spec=importlib.util.spec_from_file_location('gbm_certificate_runner_guard_tests',path)
        cls.runner=importlib.util.module_from_spec(spec)
        sys.modules[spec.name]=cls.runner
        spec.loader.exec_module(cls.runner)

    def test_project_outputs_rejected_before_creation(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name)/'project'
            root.mkdir()
            nested=root/'results'/'new'
            for output in (root,nested):
                with self.assertRaises(ValueError):
                    self.runner.prepare_output_directory(output,root)
            self.assertEqual(list(root.iterdir()),[])

    def test_symlink_into_project_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name)
            root=directory/'project'
            root.mkdir()
            alias=directory/'alias'
            alias.symlink_to(root,target_is_directory=True)
            with self.assertRaises(ValueError):
                self.runner.prepare_output_directory(alias/'new',root)
            self.assertEqual(list(root.iterdir()),[])

    def test_existing_evidence_and_file_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name)
            root=directory/'project'
            root.mkdir()
            evidence=directory/'evidence'
            evidence.mkdir()
            record=evidence/'retained.json'
            record.write_text('retained evidence')
            file=directory/'file'
            file.write_text('existing file')
            for output in (evidence,file):
                with self.assertRaises(ValueError):
                    self.runner.prepare_output_directory(output,root)
            self.assertEqual(record.read_text(),'retained evidence')
            self.assertEqual(file.read_text(),'existing file')
            self.assertEqual(list(evidence.iterdir()),[record])

    def test_external_new_or_empty_directory_is_allowed(self):
        with tempfile.TemporaryDirectory() as name:
            directory=Path(name)
            root=directory/'project'
            root.mkdir()
            new=directory/'new'
            self.assertEqual(self.runner.prepare_output_directory(new,root),new.resolve())
            self.assertTrue(new.is_dir())
            self.assertEqual(self.runner.prepare_output_directory(new,root),new.resolve())

    def test_numerical_status_excludes_statistical_misses(self):
        manifest={'deterministic_reference_checks_pass':True,'quantile_cdf_checks_pass':True,
          'sources_unchanged':True,'synthetic_summary':[{'counts':{'attempts':30,
          'certified':29,'abstained':1,'unexpected_error':0,'rounding_bound_checks':30,
          'reference_image_checks':29,'true_price_contained':0}}]}
        self.assertTrue(self.runner.numerical_verification_passes(manifest))
        for key in ('deterministic_reference_checks_pass','quantile_cdf_checks_pass','sources_unchanged'):
            changed=copy.deepcopy(manifest)
            changed[key]=False
            self.assertFalse(self.runner.numerical_verification_passes(changed))
        for key,value in (('unexpected_error',1),('rounding_bound_checks',29),
                          ('reference_image_checks',28),('abstained',0)):
            changed=copy.deepcopy(manifest)
            changed['synthetic_summary'][0]['counts'][key]=value
            self.assertFalse(self.runner.numerical_verification_passes(changed))

    def test_cli_status_fails_for_failed_verification(self):
        for result,status in ((False,1),(True,0)):
            with patch.object(self.runner,'run',return_value=result):
                self.assertEqual(self.runner.main(['--output','unused-external-output']),status)


if __name__=='__main__':
    unittest.main()
