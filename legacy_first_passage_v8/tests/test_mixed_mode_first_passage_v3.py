from __future__ import annotations
import math, unittest
from arrhenius_fracture.mixed_mode_first_passage_v3 import MixedModeContext, MixedModeFrontEngineMixin, maximum_hoop_drive

class Dummy(MixedModeFrontEngineMixin):
    pass

class Tests(unittest.TestCase):
    def drives(self,psi,K=10.0):
        d=Dummy();d._mm=MixedModeContext(0.0,target_mode_phase_deg=psi,shear_emission_weight=1.0)
        return d._mm_drives(K)
    def test_mode_I(self):
        ko,ke,m=self.drives(0.0)
        self.assertAlmostEqual(m['KI'],10.0)
        self.assertAlmostEqual(m['KII'],0.0)
        self.assertAlmostEqual(ko,10.0,places=7)
    def test_KJ_identity(self):
        for psi in (-60,-30,0,30,60):
            _,_,m=self.drives(psi,13.0)
            self.assertAlmostEqual(math.hypot(m['KI'],m['KII']),13.0,places=12)
    def test_phase_sign(self):
        _,_,a=self.drives(35.0)
        _,_,b=self.drives(-35.0)
        self.assertAlmostEqual(a['KI'],b['KI'],places=12)
        self.assertAlmostEqual(a['KII'],-b['KII'],places=12)
    def test_kink_sign_symmetry(self):
        _,pa=maximum_hoop_drive(8.0,3.0)
        _,pb=maximum_hoop_drive(8.0,-3.0)
        self.assertAlmostEqual(pa,-pb,places=8)
    def test_deterministic_threshold(self):
        c=MixedModeContext(0.0,target_mode_phase_deg=0.0,stochastic_first_passage=False)
        self.assertEqual(c.draw_threshold(),1.0)

if __name__=='__main__':unittest.main()
