import unittest
import math
from calibrate_mixed_mode_loading_v3_5 import (
    phase_ratio_support_gate, phase_spread_diagnostic, calibration_accept,
)


def record(**kw):
    d=dict(
        KI_MPa_sqrt_m=1.0,
        KII_MPa_sqrt_m=1.0,
        achieved_psi_deg=45.002359170380934,
        projection_n=50,
        projection_fit_count=4,
        projection_condition=100.0,
        projection_psi_spread_deg=24.515744698936196,
    )
    d.update(kw)
    return d


class CalibrationTests(unittest.TestCase):
    def test_reported_45_degree_case_is_accepted_with_warning(self):
        r=record()
        ok,err,reasons=calibration_accept(r,45.0,0.75)
        spread,warning=phase_spread_diagnostic(r,20.0)
        self.assertTrue(ok)
        self.assertLess(abs(err),0.01)
        self.assertEqual(reasons,[])
        self.assertTrue(warning)
        self.assertAlmostEqual(spread,24.515744698936196)

    def test_large_phase_error_still_rejected(self):
        r=record(achieved_psi_deg=49.0,projection_psi_spread_deg=1.0)
        ok,err,reasons=calibration_accept(r,45.0,0.75)
        self.assertFalse(ok)
        self.assertAlmostEqual(err,4.0)
        self.assertIn('phase_error',reasons)

    def test_insufficient_support_rejected(self):
        r=record(projection_fit_count=1)
        support,reasons=phase_ratio_support_gate(r)
        self.assertFalse(support)
        self.assertIn('fits<2',reasons)

    def test_nonfinite_mode_rejected(self):
        r=record(KI_MPa_sqrt_m=math.nan)
        support,reasons=phase_ratio_support_gate(r)
        self.assertFalse(support)
        self.assertIn('nonfinite_mode_ratio',reasons)

if __name__=='__main__':unittest.main()
