import math,unittest
import numpy as np
from arrhenius_fracture.mixed_mode_first_passage_v4 import (traction_phase_deg,shear_sign_from_basis,loading_angle_from_response_basis,energy_matrix_from_basis,energy_phase_from_matrix,process_zone_traction_probe)
class Mesh:
 def __init__(self):
  self.nodes=np.array([[0,0],[1,0],[0,1],[1,1],[2,0],[2,1]],float)*1e-5
  self.elems=np.array([[0,1,2],[1,3,2],[1,4,3],[4,5,3]])
  self.area_e=np.ones(4)*.5e-10
class Tests(unittest.TestCase):
 def test_phase(self):self.assertAlmostEqual(traction_phase_deg(1,1),45)
 def test_shear_sign(self):self.assertEqual(shear_sign_from_basis([[2,.1],[.1,-3]]),-1)
 def test_basis_diagonal(self):
  a=loading_angle_from_response_basis([[2,0],[0,1]],45);self.assertAlmostEqual(a,math.degrees(math.atan2(1,.5)))
 def test_basis_cross(self):
  M=np.array([[2,.2],[.1,1.5]]);a=loading_angle_from_response_basis(M,30);q=np.array([math.cos(math.radians(a)),math.sin(math.radians(a))]);r=M@q;self.assertAlmostEqual(math.degrees(math.atan2(r[1],r[0])),30,places=8)
 def test_energy_matrix(self):
  G=energy_matrix_from_basis(2,3,3.5,1);self.assertTrue(np.allclose(G,[[2,1],[1,3]]))
 def test_energy_phase_finite(self):self.assertTrue(np.isfinite(energy_phase_from_matrix([[2,.2],[.2,1]],20)))
 def test_probe(self):
  m=Mesh();sig=np.array([[1,1,1,1],[3,3,3,3],[.5,.5,.5,.5]],float)*1e9;d=np.zeros(6)
  r=process_zone_traction_probe(m,sig,d,[0,0],[1,0],1.5e-5,annulus_half_width=.9,sector_half_angle_deg=85,min_elements=2)
  self.assertTrue(r['reliable']);self.assertAlmostEqual(r['sigma_nn_Pa'],3e9)
if __name__=='__main__':unittest.main()
