"""Independent mathematical and numerical guards for drift-robust hedging."""
from fractions import Fraction as F
import importlib.util
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

import mpmath as mp

ROOT=Path(os.environ.get('QBP_DRIFT_PROJECT_ROOT',Path(__file__).resolve().parents[1]))
SOURCE=Path(os.environ.get('QBP_DRIFT_MODULE',ROOT/'src/gbm_drift_hedging.py'))
RUNNER=Path(os.environ.get('QBP_DRIFT_RUNNER',ROOT/'scripts/run_drift_hedge_study.py'))
for directory in (str(ROOT/'src'),str(SOURCE.parent)):
    if directory not in sys.path:sys.path.insert(0,directory)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module


e=load(SOURCE,'drift_test_engine')
r=load(SOURCE.with_name('gbm_drift_references.py'),'drift_test_references')
study=load(RUNNER,'drift_test_study')
import gbm_hedge_certificate as fixed
import gbm_hedge_references as static


def pair(record):return F(record['lower_exact_fraction']),F(record['upper_exact_fraction'])


def check_contains(test,record,value):
    lo,hi=pair(record);test.assertLessEqual(r.number(lo),value);test.assertGreaterEqual(r.number(hi),value)


class DriftRobustTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=e.ExactCappedGBM(spot='100',strike='100',cap='20',maturity='1',rate='.03',dividend='.01')
        cls.sigma=F('1/5');cls.lower=F('3/100');cls.upper=F('9/100')
        cls.cert=e.certify_drift_robust_hedge(cls.c,cls.sigma,cls.lower,cls.upper,F(1))

    def test_feasible_holdings_and_attained_bounds(self):
        cert=self.cert
        self.assertEqual(cert['status'],'certified');self.assertTrue(cert['robust_feasibility'])
        self.assertTrue(all(0<=F(x)<=1 for x in cert['holdings_exact']))
        self.assertLessEqual(F(cert['robust_spend_upper_exact']),F(cert['robust_budget_lower_exact']))
        lower=F(cert['optimal_loss_lower_bound_exact']);upper=F(cert['continuum_upper_exact'])
        self.assertGreater(lower,0);self.assertEqual(F(cert['optimality_gap_upper_exact']),upper-lower)
        self.assertLessEqual((upper-lower)/lower,F(cert['relative_gap_tolerance_exact']))

    def test_node_risks_contain_independent_piecewise_values(self):
        cert=self.cert;theta=tuple(F(x) for x in cert['holdings_exact'])
        with mp.workdps(160):
            for node,record in zip(cert['node_grid_exact'],cert['node_risks']):
                check_contains(self,record,static.fixed_risk(static.parameters(self.c,self.sigma,F(node)),theta))

    def test_continuum_upper_covers_independent_interior_diagnostics(self):
        theta=tuple(F(x) for x in self.cert['holdings_exact'])
        with mp.workdps(160):
            bound=r.number(F(self.cert['continuum_upper_exact']))
            for i in range(65):
                mu=self.lower+(self.upper-self.lower)*F(i,64)
                value=static.fixed_risk(static.parameters(self.c,self.sigma,mu),theta)
                self.assertLessEqual(value,bound)

    def test_exact_mixture_simplex_and_nonnegative_multipliers(self):
        weights=tuple(F(x) for x in self.cert['mixture_weights_exact'])
        self.assertEqual(sum(weights,F(0)),1);self.assertTrue(all(x>=0 for x in weights))
        self.assertTrue(all(F(x)>=0 for x in self.cert['financial_multipliers_exact']))
        self.assertEqual(len(weights),len(self.cert['node_grid_exact']))

    def test_mixture_dual_matches_independent_high_precision_formula(self):
        cert=self.cert;weights=tuple(F(x) for x in cert['mixture_weights_exact'])
        lam=tuple(F(x) for x in cert['financial_multipliers_exact'])
        with mp.workdps(160):
            laws=[static.parameters(self.c,self.sigma,F(mu)) for mu in cert['node_grid_exact']]
            coeff=[static.coefficients(model) for model in laws]
            matrix=mp.matrix(2);cross=mp.matrix(2,1);constant=mp.mpf(0)
            for w,(M,h,c) in zip(weights,coeff):
                matrix+=r.number(w)*M;cross+=r.number(w)*h;constant+=r.number(w)*c
            A=((-1,0),(0,-1),(1,0),(0,1),(laws[0]['s'],laws[0]['ck']))
            b=(0,0,1,1,r.number(F(cert['budget_multiplier_exact']))*laws[0]['premium'])
            z=mp.matrix([cross[j]-sum(r.number(lam[i])*A[i][j] for i in range(5))/2 for j in range(2)])
            dual=constant-sum(r.number(lam[i])*b[i] for i in range(5))-(z.T*(matrix**-1)*z)[0]
            check_contains(self,cert['raw_mixture_dual'],dual)
            for i in range(2):
                check_contains(self,cert['mixture_quadratic']['cross'][i],cross[i])
                for j in range(2):check_contains(self,cert['mixture_quadratic']['matrix'][i][j],matrix[i,j])
            check_contains(self,cert['mixture_quadratic']['constant'],constant)

    def test_cell_inflation_uses_maturity_sigma_and_mesh(self):
        cert=self.cert;width=F(cert['max_cell_width_exact'])
        exponent=self.c.maturity*width*width/(8*self.sigma*self.sigma)
        self.assertEqual(F(cert['cell_exponent_exact']),exponent)
        with mp.workdps(160):check_contains(self,cert['cell_inflation'],mp.exp(r.number(exponent)))
        self.assertEqual(F(cert['continuum_upper_exact']),max(pair(x)[1] for x in cert['node_risks'])*F(cert['max_cell_inflation_upper_exact']))

    def test_independent_finite_grid_optimum_is_bracketed(self):
        cert=self.cert
        with mp.workdps(160):
            qp=r.finite_grid_optimum(self.c,self.sigma,cert['node_grid_exact'],F(1))
            self.assertTrue(qp['success'])
            self.assertLessEqual(r.number(F(cert['optimal_loss_lower_bound_exact'])),qp['risk'])
            self.assertGreaterEqual(r.number(F(cert['continuum_upper_exact'])),qp['risk'])

    def test_singleton_drift_matches_static_model(self):
        cert=e.certify_drift_robust_hedge(self.c,self.sigma,F('.06'),F('.06'),1)
        self.assertEqual(cert['node_grid_exact'],['3/50']);self.assertEqual(F(cert['cell_exponent_exact']),0)
        self.assertEqual(F(cert['max_cell_inflation_upper_exact']),1)
        with mp.workdps(160):
            qp=static.optimum(static.parameters(self.c,self.sigma,F('.06')),F(1))
            self.assertLessEqual(r.number(F(cert['optimal_loss_lower_bound_exact'])),qp['risk'])
            self.assertGreaterEqual(r.number(F(cert['continuum_upper_exact'])),qp['risk'])

    def test_zero_budget_is_exact_no_hedge(self):
        cert=e.certify_drift_robust_hedge(self.c,self.sigma,self.lower,self.upper,0)
        self.assertTrue(cert['robust_feasibility']);self.assertEqual(tuple(F(x) for x in cert['holdings_exact']),(0,0))
        with mp.workdps(160):
            for mu,record in zip(cert['node_grid_exact'],cert['node_risks']):
                check_contains(self,record,static.fixed_risk(static.parameters(self.c,self.sigma,F(mu)),(0,0)))

    def test_endpoints_can_miss_a_large_interior_capped_claim_risk(self):
        with mp.workdps(160):
            values=[static.fixed_risk(static.parameters(self.c,F('.01'),F(mu)),(0,F('.5')))
                    for mu in ('.105','.19','.275')]
            self.assertGreater(values[1],4*max(values[0],values[2]))
            # This is an actual financially admissible holding at budget V.
            model=static.parameters(self.c,F('.01'),F('.19'))
            self.assertLess(model['ck']/2,model['premium'])
            self.assertLess(abs(values[1]-mp.mpf('68.08541970641572760856660083362491081449')),mp.mpf('1e-37'))

    def test_max_nodes_exhaustion_returns_unresolved_not_false_success(self):
        cert=e.certify_drift_robust_hedge(self.c,self.sigma,F('-.5'),F('.5'),1,
            initial_cells=1,max_nodes=2,relative_gap_tolerance='1e-30',cell_inflation_tolerance='1e-30')
        self.assertEqual(cert['status'],'unresolved');self.assertIn('unresolved_reason',cert)
        self.assertTrue(cert['robust_feasibility']);self.assertGreater(F(cert['continuum_upper_exact']),0)

    def test_invalid_hedge_inputs_are_rejected(self):
        for args in ((0,0,1,1),(-.1,0,1,1),(.2,1,0,1),(.2,0,1,-1),(.2,math.nan,1,1),(.2,0,math.inf,1)):
            with self.subTest(args=args),self.assertRaises((ValueError,TypeError)):
                e.certify_drift_robust_hedge(self.c,*args)
        for kwargs in ({'initial_cells':0},{'max_nodes':1},{'precision_bits':64},
                       {'relative_gap_tolerance':0},{'cell_inflation_tolerance':math.nan}):
            with self.subTest(kwargs=kwargs),self.assertRaises((ValueError,TypeError)):
                e.certify_drift_robust_hedge(self.c,self.sigma,0,1,1,**kwargs)

    def test_precision_context_restored_on_success_and_failure(self):
        previous=e.ctx.prec
        try:
            e.ctx.prec=111
            e.certify_drift_robust_hedge(self.c,self.sigma,self.lower,self.upper,1)
            self.assertEqual(e.ctx.prec,111)
            with self.assertRaises(ValueError):e.certify_drift_robust_hedge(self.c,self.sigma,1,0,1)
            self.assertEqual(e.ctx.prec,111)
        finally:e.ctx.prec=previous


class KnownSigmaInferenceTests(unittest.TestCase):
    def test_drift_centre_includes_dividend_ito_correction_and_spacing(self):
        cert=e.certified_known_sigma_drift_interval(('0.01','0.02'),'.2','.1',dividend='.03')
        self.assertEqual(F(cert['drift_centre_exact']),F('.03')+F('.2')**2/2+F('.015')/F('.1'))
        self.assertEqual(F(cert['observation_spacing_exact']),F('.1'))

    def test_quantile_cdf_witnesses_and_independent_gaussian_half_width(self):
        cert=e.certified_known_sigma_drift_interval(('0.01','0.02'),'.2','.1',alpha='.05')
        quantile=cert['chi_square_quantile'];p=F(quantile['probability_exact'])
        self.assertLess(pair(quantile['cdf_at_lower'])[1],p)
        self.assertGreater(pair(quantile['cdf_at_upper'])[0],p)
        with mp.workdps(160):
            z=mp.sqrt(2)*mp.erfinv(mp.mpf('.95'))
            expected=mp.mpf('.2')*z/mp.sqrt(mp.mpf('2')*mp.mpf('.1'))
            self.assertGreaterEqual(r.number(F(cert['gaussian_half_width_upper_exact'])),expected)
            self.assertLess(r.number(F(cert['gaussian_half_width_upper_exact']))-expected,mp.mpf('1e-35'))

    def test_rounding_cells_expand_mean_and_drift_interval_exactly(self):
        raw=e.certified_known_sigma_drift_interval(('0.01','0.02'),'.2','.1')
        rounded=e.certified_known_sigma_drift_interval(('0.01','0.02'),'.2','.1',rounding_radii=('.001','.003'))
        self.assertEqual(F(rounded['mean_rounding_radius_exact']),F('.002'))
        self.assertEqual(F(rounded['drift_centre_rounding_radius_exact']),F('.02'))
        lo,hi=pair(raw['drift_interval']);rlo,rhi=pair(rounded['drift_interval'])
        self.assertEqual(rlo,lo-F('.02'));self.assertEqual(rhi,hi+F('.02'))

    def test_one_observation_is_allowed_with_known_variance(self):
        cert=e.certified_known_sigma_drift_interval(('0',),'.2','1')
        self.assertEqual(cert['observations'],1);self.assertEqual(cert['status'],'certified')

    def test_invalid_confidence_inputs_are_rejected(self):
        for observations,sigma,spacing,kwargs in (([],'.2','1',{}),('00','.2','1',{}),
            (['0'],'0','1',{}),(['0'],'.2','0',{}),(['nan'],'.2','1',{}),
            (['0'],'.2','1',{'alpha':1}),(['0'],'.2','1',{'rounding_radii':[-1]}),
            (['0'],'.2','1',{'rounding_radii':[]})):
            with self.subTest(observations=observations,sigma=sigma,kwargs=kwargs),self.assertRaises((ValueError,TypeError)):
                e.certified_known_sigma_drift_interval(observations,sigma,spacing,**kwargs)


class EvidenceProtectionTests(unittest.TestCase):
    def test_internal_and_symlink_outputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();link=Path(tmp)/'link';link.symlink_to(root,target_is_directory=True)
            for out in (root/'new',link/'new'):
                with self.assertRaises(ValueError):study.output_destination(out,root)
            self.assertFalse((root/'new').exists())

    def test_existing_evidence_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();out=Path(tmp)/'out';out.mkdir();(out/'evidence').write_text('preserve')
            with self.assertRaises(FileExistsError):study.output_destination(out,root)
            self.assertEqual((out/'evidence').read_text(),'preserve')

    def test_missing_rows_failure_or_source_mutation_cannot_pass(self):
        good={'primary_rows':12,'confidence_rows':30,'dense_checks':12*5*257,
              'finite_grid_reference_rows':12,'bridge_rows':3,'bridge_dense_checks':3*4*257,'sources_unchanged':True,'failures':[]}
        self.assertTrue(study.overall_status(good))
        for key,value in (('primary_rows',11),('confidence_rows',29),('dense_checks',1),
                          ('finite_grid_reference_rows',11),('bridge_rows',2),('bridge_dense_checks',1),('sources_unchanged',False),('failures',['x'])):
            bad={**good,key:value};self.assertFalse(study.overall_status(bad))


if __name__=='__main__':unittest.main()
