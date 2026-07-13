from __future__ import annotations
import unittest
from run_mixed_mode_fem_czm_v2_campaign import angle_error, next_alpha

class ControlTests(unittest.TestCase):
    def test_angle_error_wrap(self):
        self.assertAlmostEqual(angle_error(-179, 179), 2.0)
        self.assertAlmostEqual(angle_error(179, -179), -2.0)

    def test_newton_first_update(self):
        h=[{'alpha_deg':0.0,'error_deg':20.0}]
        self.assertAlmostEqual(next_alpha(h,2.0,-70,70,15),-10.0)

    def test_bracketed_update(self):
        h=[{'alpha_deg':-5.0,'error_deg':-2.0},{'alpha_deg':5.0,'error_deg':2.0}]
        self.assertAlmostEqual(next_alpha(h,1.0,-70,70,15),0.0)

    def test_step_is_safeguarded(self):
        h=[{'alpha_deg':0.0,'error_deg':100.0}]
        self.assertAlmostEqual(next_alpha(h,0.1,-70,70,15),-15.0)

if __name__=='__main__':unittest.main()
