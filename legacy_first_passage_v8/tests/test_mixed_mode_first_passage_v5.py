import math
import unittest

import numpy as np

from arrhenius_fracture.mixed_mode_first_passage_v5 import (
    directional_drive_factors,
    directional_shape_metrics,
    energy_matrix_from_basis,
    loading_angle_from_response_basis,
    phase_derivative_deg_per_deg,
    safeguarded_alpha_update,
    shear_sign_from_basis,
    traction_phase_deg,
)


class MechanicsTests(unittest.TestCase):
    def test_mode_I_normalization(self):
        f = directional_drive_factors(.87, .04, .87, .04, 1.0, 5.0)
        self.assertAlmostEqual(f["cleavage_factor"], 1.0)
        self.assertAlmostEqual(f["emission_factor"], 1.0)

    def test_shear_augments_emission_not_cleavage(self):
        f = directional_drive_factors(.87, .60, .87, .04, 1.0, 5.0)
        self.assertAlmostEqual(f["cleavage_factor"], 1.0)
        self.assertGreater(f["emission_factor"], 1.0)

    def test_factor_cap_is_audited(self):
        f = directional_drive_factors(20.0, 20.0, 1.0, 0.0, 1.0, 3.0)
        self.assertEqual(f["cleavage_factor"], 3.0)
        self.assertEqual(f["emission_factor"], 3.0)
        self.assertTrue(f["directional_factor_cap_active"])

    def test_directional_metric_finite(self):
        S = np.array([[0.2e9, 0.3e9], [0.3e9, 1.0e9]])
        m = directional_shape_metrics(S, 45.0, .3, np.array([1.0, 0.0]))
        self.assertTrue(m["reliable"])
        self.assertGreater(m["cleavage_shape"], 0.0)
        self.assertGreaterEqual(m["slip_shape"], 0.0)

    def test_basis_cross(self):
        M = np.array([[2.0, .2], [.1, 1.5]])
        alpha = loading_angle_from_response_basis(M, 30.0)
        q = np.array([math.cos(math.radians(alpha)), math.sin(math.radians(alpha))])
        r = M @ q
        self.assertAlmostEqual(math.degrees(math.atan2(r[1], r[0])), 30.0, places=8)

    def test_phase_derivative(self):
        M = np.eye(2)
        self.assertAlmostEqual(phase_derivative_deg_per_deg(M, 12.0), 1.0)

    def test_safeguarded_update(self):
        M = np.eye(2)
        self.assertAlmostEqual(safeguarded_alpha_update(10.0, 10.0, 0.0, M), 0.0)

    def test_shear_sign(self):
        self.assertEqual(shear_sign_from_basis([[2.0, .1], [.1, -3.0]]), -1.0)

    def test_phase(self):
        self.assertAlmostEqual(traction_phase_deg(1.0, 1.0), 45.0)

    def test_energy_matrix(self):
        G = energy_matrix_from_basis(2.0, 3.0, 3.5, 1.0)
        self.assertTrue(np.allclose(G, [[2.0, 1.0], [1.0, 3.0]]))


if __name__ == "__main__":
    unittest.main()
