from __future__ import annotations
import unittest
import numpy as np
from arrhenius_fracture.mixed_mode_first_passage_v1 import (
    maximum_hoop_drive, project_near_tip_modes, MixedModeContext,
    _mixed_solve_factory)

class DummyMesh:
    pass

class MixedModeTests(unittest.TestCase):
    def test_maximum_hoop_mode_I(self):
        k,th=maximum_hoop_drive(10.0,0.0)
        self.assertAlmostEqual(k,10.0,places=8)
        self.assertAlmostEqual(th,0.0,places=6)

    def test_maximum_hoop_changes_sign_with_mode_II(self):
        _,a=maximum_hoop_drive(0.0,5.0)
        _,b=maximum_hoop_drive(0.0,-5.0)
        self.assertAlmostEqual(a,-b,places=6)
        self.assertGreater(abs(a),60.0)

    def test_seeded_threshold_reproducible(self):
        a=MixedModeContext(0.0,solver_seed=123)
        b=MixedModeContext(0.0,solver_seed=123)
        self.assertAlmostEqual(a.draw_threshold(),b.draw_threshold(),places=14)

    def test_projection_recovers_synthetic_modes(self):
        # Small triangular centroids in a forward wedge. Set all element stresses
        # to the exact leading Williams theta=0 values.
        mesh=DummyMesh()
        cent=np.array([[1.0,0.0],[2.0,.05],[3.0,-.05],[4.0,.08],[5.0,-.08]])*1e-5
        # triangles whose centroids equal cent approximately
        nodes=[]; elems=[]
        for i,c in enumerate(cent):
            j=len(nodes); nodes += [c+[-1e-7,-1e-7],c+[1e-7,-1e-7],c+[0,2e-7]]; elems.append([j,j+1,j+2])
        mesh.nodes=np.array(nodes);mesh.elems=np.array(elems);mesh.area_e=np.ones(len(elems))*1e-14
        KI=8e6;KII=-3e6;r=np.linalg.norm(cent,axis=1);b=1/np.sqrt(2*np.pi*r)
        sig=np.zeros((3,len(elems)));sig[1]=KI*b;sig[2]=KII*b
        d=np.zeros(len(nodes))
        got=project_near_tip_modes(mesh,sig,d,np.zeros(2),np.array([1.,0.]),5e-6,7e-5,20)
        self.assertAlmostEqual(got['KI_Pa_sqrt_m']/KI,1.0,places=5)
        self.assertAlmostEqual(got['KII_Pa_sqrt_m']/KII,1.0,places=5)

if __name__=='__main__':unittest.main()
