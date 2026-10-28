"""Independent numerical and boundary checks of validated static hedging."""
from fractions import Fraction as F
import importlib.util
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import mpmath as mp

ROOT=Path(os.environ.get('QBP_HEDGE_PROJECT_ROOT',Path(__file__).resolve().parents[1]))
SOURCE=Path(os.environ.get('QBP_HEDGE_CERTIFICATE_MODULE',ROOT/'src/gbm_hedge_certificate.py'))
RUNNER=Path(os.environ.get('QBP_HEDGE_STUDY_RUNNER',ROOT/'scripts/run_hedge_certificate_study.py'))
for directory in (str(ROOT/'src'),str(SOURCE.parent)):
    if directory not in sys.path:sys.path.insert(0,directory)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module


e=load(SOURCE,'hedge_test_engine')
r=load(SOURCE.with_name('gbm_hedge_references.py'),'hedge_test_references')
study=load(RUNNER,'hedge_test_study')


def pair(record):return F(record['lower_exact_fraction']),F(record['upper_exact_fraction'])


def check_contains(test,record,value):
    lo,hi=pair(record)
    test.assertLessEqual(r.number(lo),value)
    test.assertGreaterEqual(r.number(hi),value)


class ValidatedHedgeTests(unittest.TestCase):
    def setUp(self):
        self.c=e.ExactCappedGBM()
        self.short=e.ExactCappedGBM(spot=100.,strike=100.,cap=20.,maturity=.05,rate=.03,dividend=.01)

    def test_fixed_risk_contains_independent_piecewise_moments(self):
        cert=e.certified_hedge_risk(self.c,.2,.06,(F(1,10),F(1,5)))
        with mp.workdps(160):
            model=r.parameters(self.c,F.from_float(.2),F.from_float(.06))
            check_contains(self,cert['risk'],r.fixed_risk(model,(F(1,10),F(1,5))))

    def test_unrestricted_signed_holdings_risk_is_supported(self):
        cert=e.certified_hedge_risk(self.c,.2,.06,(-.2,1.3))
        self.assertTrue(cert['risk']['strictly_positive'])
        with mp.workdps(160):
            model=r.parameters(self.c,F.from_float(.2),F.from_float(.06))
            check_contains(self,cert['risk'],r.fixed_risk(model,(F.from_float(-.2),F.from_float(1.3))))

    def test_zero_budget_is_exactly_cash_only(self):
        cert=e.certify_hedge(self.c,.2,.06,0,(.1,.2))
        self.assertEqual([F(x) for x in cert['holdings_exact']],[0,0])
        self.assertTrue(cert['robust_feasibility'])
        with mp.workdps(160):
            check_contains(self,cert['risk'],r.fixed_risk(r.parameters(self.c,F.from_float(.2),F.from_float(.06)),(0,0)))

    def test_original_outside_box_is_preserved_and_repaired(self):
        cert=e.certify_hedge(self.c,.2,.06,1,(-.1,1.1))
        self.assertEqual(tuple(F(x) for x in cert['original_holdings_exact']),(F.from_float(-.1),F.from_float(1.1)))
        self.assertGreater(F(cert['original_box_violation_exact']),0)
        self.assertTrue(all(0<=F(x)<=1 for x in cert['repaired_original_holdings_exact']))
        self.assertTrue(cert['robust_feasibility'])

    def test_unit_call_exact_fair_budget_violation_is_positive(self):
        cert=e.certify_hedge(self.short,.06,.06,1,(0,1))
        self.assertEqual(cert['original_budget_status'],'certified_infeasible')
        self.assertGreater(pair(cert['original_budget_violation'])[0],0)
        self.assertLess(F(cert['repaired_original_holdings_exact'][1]),1)

    def test_unit_call_is_funding_feasible_at_double_budget(self):
        cert=e.certify_hedge(self.short,.06,.06,2,(0,1))
        self.assertEqual(cert['original_feasibility_status'],'certified_feasible')
        self.assertEqual(tuple(F(x) for x in cert['repaired_original_holdings_exact']),(0,1))

    def test_short_dated_unit_call_retains_positive_tail_risk(self):
        cert=e.certified_hedge_risk(self.short,.06,.06,(0,1))
        self.assertGreater(pair(cert['risk'])[0],0)
        self.assertLess(pair(cert['risk'])[1],F('1e-40'))
        with mp.workdps(160):
            model=r.parameters(self.short,F.from_float(.06),F.from_float(.06))
            check_contains(self,cert['risk'],r.fixed_risk(model,(0,1)))

    def test_refined_primal_dual_bounds_contain_independent_optimum(self):
        cert=e.certify_hedge(self.c,.2,.06,1,(0,0))
        with mp.workdps(160):
            model=r.parameters(self.c,F.from_float(.2),F.from_float(.06));qp=r.optimum(model,F(1))
            self.assertLessEqual(r.number(F(cert['optimal_loss_lower_bound_exact'])),qp['risk'])
            self.assertGreaterEqual(r.number(F(cert['risk']['upper_exact_fraction'])),qp['risk'])
            self.assertLess(r.number(F(cert['risk']['upper_exact_fraction'])),r.fixed_risk(model,(0,0)))

    def test_eigenvalue_and_distance_bounds_are_verified(self):
        cert=e.certify_hedge(self.c,.2,.12,F(1,2),(.1,.2))
        self.assertGreater(F(cert['minimum_eigenvalue_lower_exact']),0)
        gap=F(cert['optimality_gap_upper_exact']);m=F(cert['minimum_eigenvalue_lower_exact'])
        self.assertEqual(F(cert['squared_distance_to_unique_optimum_upper_exact']),gap/m)
        with mp.workdps(160):
            model=r.parameters(self.c,F.from_float(.2),F.from_float(.12));qp=r.optimum(model,F(1,2))
            distance=sum((r.number(F(cert['holdings_exact'][i]))-qp['holdings'][i])**2 for i in range(2))
            self.assertLessEqual(distance,r.number(gap/m))

    def test_negative_proposed_multipliers_are_explicitly_preserved(self):
        cert=e.certify_hedge(self.c,.2,.06,1,(.1,.2),(-1,0,0,0,0),refine=False)
        self.assertEqual(F(cert['raw_multipliers_exact'][0]),-1)
        self.assertEqual(F(cert['safe_multipliers_exact'][0]),0)
        self.assertTrue(cert['robust_feasibility'])

    def test_directed_float_summaries_enclose_exact_endpoints(self):
        cert=e.certified_hedge_risk(self.short,.06,.06,(0,1))['risk']
        lo,hi=pair(cert)
        self.assertLessEqual(F.from_float(cert['lower_float']),lo)
        self.assertGreaterEqual(F.from_float(cert['upper_float']),hi)

    def test_malformed_inputs_are_rejected(self):
        for sigma,mu,holdings in ((0,.06,(0,0)),(-1,.06,(0,0)),(True,.06,(0,0)),(.2,math.nan,(0,0)),(.2,.06,(0,)),(.2,.06,'00')):
            with self.subTest(sigma=sigma,mu=mu,holdings=holdings),self.assertRaises((ValueError,TypeError)):
                e.certified_hedge_risk(self.c,sigma,mu,holdings)
        with self.assertRaises(ValueError):e.certify_hedge(self.c,.2,.06,-1,(0,0))

    def test_budget_exhaustion_is_not_a_false_certificate(self):
        with self.assertRaises(ArithmeticError):
            e.certify_hedge(self.short,.06,.06,1,(0,1),absolute_tolerance='1e-200',start_bits=64,max_bits=64)

    def test_precision_context_is_restored_on_success_and_failure(self):
        previous=e.ctx.prec
        try:
            e.ctx.prec=111
            e.certified_hedge_risk(self.c,.2,.06,(0,0))
            self.assertEqual(e.ctx.prec,111)
            with self.assertRaises(ArithmeticError):
                e.certified_hedge_risk(self.c,.2,.06,(0,0),absolute_tolerance='1e-200',start_bits=64,max_bits=64)
            self.assertEqual(e.ctx.prec,111)
        finally:e.ctx.prec=previous

    def test_invalid_precision_and_relative_tolerance_are_rejected(self):
        for kwargs in ({'start_bits':True},{'max_bits':32},{'relative_tolerance':0},{'relative_tolerance':math.inf}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                e.certified_hedge_risk(self.c,.2,.06,(0,0),**kwargs)


class EvidenceProtectionTests(unittest.TestCase):
    def test_internal_output_is_rejected_before_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();output=root/'new'
            with self.assertRaises(ValueError):study.output_destination(output,root)
            self.assertFalse(output.exists())

    def test_symlink_output_into_project_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();link=Path(tmp)/'link';link.symlink_to(root,target_is_directory=True)
            with self.assertRaises(ValueError):study.output_destination(link/'new',root)
            self.assertFalse((root/'new').exists())

    def test_existing_evidence_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();output=Path(tmp)/'out';output.mkdir();sentinel=output/'saved';sentinel.write_text('unchanged')
            with self.assertRaises(FileExistsError):study.output_destination(output,root)
            self.assertEqual(sentinel.read_text(),'unchanged')

    def test_existing_output_file_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();output=Path(tmp)/'out';output.write_text('evidence')
            with self.assertRaises(FileExistsError):study.output_destination(output,root)
            self.assertEqual(output.read_text(),'evidence')

    def test_empty_external_destination_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'project';root.mkdir();output=Path(tmp)/'out';output.mkdir()
            self.assertEqual(study.output_destination(output,root),output.resolve())

    def test_failed_numerics_produce_failure_status(self):
        base={'failures':[],'sources_unchanged':True,'fixed_risk_rows':168,'refined_rows':126,'slsqp_rows':126,'positive_integral_rows':60}
        self.assertTrue(study.overall_status(base))
        self.assertFalse(study.overall_status({**base,'failures':['a:failed']}))
        self.assertFalse(study.overall_status({**base,'sources_unchanged':False}))
        with patch.object(study,'run',return_value={'overall_numerical_verification':False}):
            self.assertEqual(study.main(['--output','unused']),1)
        with patch.object(study,'run',return_value={'overall_numerical_verification':True}):
            self.assertEqual(study.main(['--output','unused']),0)


if __name__=='__main__':unittest.main()
