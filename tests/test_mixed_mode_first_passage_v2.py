from __future__ import annotations
import unittest
import numpy as np
from arrhenius_fracture.mixed_mode_first_passage_v2 import (
    maximum_hoop_drive, project_near_tip_modes, MixedModeContext, _williams_shapes)

class DummyMesh: pass

def synthetic_mesh(KI=8e6,KII=-3e6,noise=0.0,outlier=False):
    mesh=DummyMesh();pts=[]
    for r in np.linspace(1.2e-5,5.5e-5,7):
        for th in np.radians([-90,-60,-30,0,30,60,90]):pts.append([r*np.cos(th),r*np.sin(th)])
    cent=np.asarray(pts);nodes=[];elems=[]
    for c in cent:
        j=len(nodes);nodes += [c+[-1e-7,-1e-7],c+[1e-7,-1e-7],c+[0,2e-7]];elems.append([j,j+1,j+2])
    mesh.nodes=np.asarray(nodes);mesh.elems=np.asarray(elems);mesh.area_e=np.ones(len(elems))*1e-14;mesh.hbar_tip=8e-6;mesh.hbar=8e-6
    r=np.linalg.norm(cent,axis=1);th=np.arctan2(cent[:,1],cent[:,0]);fI,fII=_williams_shapes(th);b=1/np.sqrt(2*np.pi*r)
    sl=b[:,None]*(KI*fI+KII*fII)+np.array([1.2e8,-0.4e8,0.2e8])
    rng=np.random.default_rng(2);sl += noise*rng.normal(size=sl.shape)
    if outlier:sl[4] += np.array([2e9,-1e9,1.5e9])
    sig=sl.T.copy();d=np.zeros(len(nodes));return mesh,sig,d

class Tests(unittest.TestCase):
    def test_mode_I_hoop(self):
        k,th=maximum_hoop_drive(10,0);self.assertAlmostEqual(k,10);self.assertAlmostEqual(th,0)
    def test_deterministic_default(self):
        c=MixedModeContext(0);self.assertFalse(c.stochastic_first_passage);self.assertEqual(c.draw_threshold(),1.0)
    def test_full_field_exact(self):
        m,s,d=synthetic_mesh();g=project_near_tip_modes(m,s,d,[0,0],[1,0],1e-5,6e-5,105,.85)
        self.assertAlmostEqual(g['KI_Pa_sqrt_m']/8e6,1,places=5);self.assertAlmostEqual(g['KII_Pa_sqrt_m']/(-3e6),1,places=5);self.assertTrue(g['mode_projection_reliable'])
    def test_robust_to_outlier(self):
        m,s,d=synthetic_mesh(noise=2e6,outlier=True);g=project_near_tip_modes(m,s,d,[0,0],[1,0],1e-5,6e-5,105,.85)
        self.assertLess(abs(g['KI_Pa_sqrt_m']/8e6-1),.08);self.assertLess(abs(g['KII_Pa_sqrt_m']/(-3e6)-1),.12)
    def test_sign_symmetry(self):
        m,s,d=synthetic_mesh(KII=3e6);a=project_near_tip_modes(m,s,d,[0,0],[1,0],1e-5,6e-5,105,.85)
        m,s,d=synthetic_mesh(KII=-3e6);b=project_near_tip_modes(m,s,d,[0,0],[1,0],1e-5,6e-5,105,.85)
        self.assertAlmostEqual(a['mode_phase_deg'],-b['mode_phase_deg'],places=5)
if __name__=='__main__':unittest.main()
