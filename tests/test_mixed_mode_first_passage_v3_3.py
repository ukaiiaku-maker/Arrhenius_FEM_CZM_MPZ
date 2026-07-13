from __future__ import annotations
import math, unittest
import numpy as np
from arrhenius_fracture.mixed_mode_first_passage_v3_3 import (
    MixedModeContext, MixedModeFrontEngineMixin, maximum_hoop_drive,
    loading_angle_from_mode_basis, phase_projection_gate, angle_error_deg,
    mode_signs_from_basis, apply_mode_signs, circular_robust_spread_deg)

class Dummy(MixedModeFrontEngineMixin):
    pass

class Tests(unittest.TestCase):
    def drives(self,psi,K=10.0):
        d=Dummy();d._mm=MixedModeContext(0.0,target_mode_phase_deg=psi,shear_emission_weight=1.0)
        return d._mm_drives(K)
    def test_mode_I(self):
        ko,ke,m=self.drives(0.0)
        self.assertAlmostEqual(m['KI'],10.0);self.assertAlmostEqual(m['KII'],0.0);self.assertAlmostEqual(ko,10.0,places=7)
    def test_KJ_identity(self):
        for psi in (-60,-30,0,30,60):
            _,_,m=self.drives(psi,13.0);self.assertAlmostEqual(math.hypot(m['KI'],m['KII']),13.0,places=12)
    def test_kink_sign_symmetry(self):
        _,a=maximum_hoop_drive(8.0,3.0);_,b=maximum_hoop_drive(8.0,-3.0);self.assertAlmostEqual(a,-b,places=8)
    def test_deterministic_threshold(self):
        c=MixedModeContext(0.0,target_mode_phase_deg=0.0,stochastic_first_passage=False);self.assertEqual(c.draw_threshold(),1.0)
    def test_basis_diagonal(self):
        M=np.array([[2.0,0.0],[0.0,1.0]])
        self.assertAlmostEqual(loading_angle_from_mode_basis(M,0.0),0.0,places=12)
        self.assertAlmostEqual(loading_angle_from_mode_basis(M,45.0),math.degrees(math.atan(2.0)),places=12)
    def test_basis_cross_coupling(self):
        M=np.array([[2.0,0.2],[-0.1,1.3]])
        alpha=loading_angle_from_mode_basis(M,30.0)
        q=np.array([math.cos(math.radians(alpha)),math.sin(math.radians(alpha))])
        k=M@q
        self.assertAlmostEqual(math.degrees(math.atan2(k[1],k[0])),30.0,places=10)
    def test_phase_gate_ignores_amplitude_residual(self):
        r={'KI_MPa_sqrt_m':2.0,'KII_MPa_sqrt_m':1.0,'achieved_psi_deg':26.565,
           'projection_n':40,'projection_fit_count':4,'projection_condition':1e5,
           'projection_psi_spread_deg':2.0,'projection_rel_rmse':0.95,
           'amplitude_projection_reliable':False}
        ok,reasons=phase_projection_gate(r);self.assertTrue(ok);self.assertEqual(reasons,[])
    def test_phase_gate_rejects_spread(self):
        r={'KI_MPa_sqrt_m':2.0,'KII_MPa_sqrt_m':1.0,'achieved_psi_deg':26.565,
           'projection_n':40,'projection_fit_count':4,'projection_condition':1e5,
           'projection_psi_spread_deg':40.0}
        ok,reasons=phase_projection_gate(r,max_phase_spread_deg=20.0);self.assertFalse(ok);self.assertIn('phase_spread',reasons)
    def test_angle_wrap(self):
        self.assertAlmostEqual(angle_error_deg(-179,179),2.0)
    def test_reported_basis_sign_normalization(self):
        raw=np.array([[-3.13374403e-2,2.28086067e-4],[3.55930199e-6,1.01476363e-3]])
        signs=mode_signs_from_basis(raw)
        np.testing.assert_allclose(signs,[-1.0,1.0])
        M=np.diag(signs)@raw
        for target in (-60,-30,0,30,60):
            alpha=loading_angle_from_mode_basis(M,target,max_abs_alpha_deg=89.9)
            q=np.array([math.cos(math.radians(alpha)),math.sin(math.radians(alpha))])
            k=M@q
            psi=math.degrees(math.atan2(k[1],k[0]))
            self.assertAlmostEqual(psi,target,places=9)
    def test_circular_spread_wraps_branch_cut(self):
        vals=np.array([179.2,-179.4,178.8,-178.9])
        self.assertLess(circular_robust_spread_deg(vals),3.0)
    def test_circular_spread_rejects_opposite_modes(self):
        vals=np.array([0.0,0.0,180.0,180.0])
        self.assertGreater(circular_robust_spread_deg(vals),100.0)
    def test_apply_mode_signs(self):
        KI,KII,psi=apply_mode_signs(-2.0,1.0,np.array([-1.0,1.0]))
        self.assertEqual(KI,2.0);self.assertEqual(KII,1.0)
        self.assertAlmostEqual(psi,math.degrees(math.atan2(1,2)))

if __name__=='__main__':unittest.main()
