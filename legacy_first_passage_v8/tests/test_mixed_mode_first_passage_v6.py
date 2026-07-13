import math
import unittest

import numpy as np

from arrhenius_fracture.mixed_mode_first_passage_v6 import (
    directional_drive_factors,
    loading_coefficients_from_response_basis,
    loading_coefficients_from_z,
    loading_z_from_coefficients,
    normalize_loading_coefficients,
    phase_from_response_z,
    safeguarded_event_z_update,
)


class MechanicsTests(unittest.TestCase):
    def test_loading_coefficients_normalized(self):
        q = normalize_loading_coefficients(3.0, 4.0)
        self.assertAlmostEqual(math.hypot(*q), 1.0)

    def test_response_basis_recovers_target(self):
        M = np.array([[2.0, 0.2], [0.1, 1.5]])
        q = loading_coefficients_from_response_basis(M, 30.0)
        r = M @ np.array(q)
        self.assertAlmostEqual(math.degrees(math.atan2(r[1], r[0])), 30.0, places=8)

    def test_z_round_trip_large_ratio(self):
        q = normalize_loading_coefficients(1.0, 5000.0)
        z = loading_z_from_coefficients(*q)
        q2 = loading_coefficients_from_z(z)
        self.assertTrue(np.allclose(q, q2, rtol=1e-10, atol=1e-12))

    def test_z_has_no_89p9_saturation(self):
        q = loading_coefficients_from_z(12.0)
        angle = math.degrees(math.atan2(q[1], q[0]))
        self.assertGreater(angle, 89.9)
        self.assertLess(angle, 90.0)

    def test_event_update_bracketed(self):
        M = np.eye(2)
        samples = [
            {"loading_z": -1.0, "achieved_psi_deg": -40.0},
            {"loading_z": 1.0, "achieved_psi_deg": 40.0},
        ]
        z = safeguarded_event_z_update(samples, 0.0, M)
        self.assertAlmostEqual(z, 0.0)

    def test_event_update_first_sample_uses_backend_slope(self):
        M = np.eye(2)
        samples = [{"loading_z": 0.0, "achieved_psi_deg": 0.0}]
        z = safeguarded_event_z_update(samples, 20.0, M, max_step=2.0)
        self.assertGreater(z, 0.0)

    def test_phase_from_response_z(self):
        M = np.eye(2)
        z = np.arcsinh(1.0)
        self.assertAlmostEqual(phase_from_response_z(M, z), 45.0)

    def test_mode_I_factor_normalization(self):
        f = directional_drive_factors(.8, .1, .8, .1, 1.0, 5.0)
        self.assertAlmostEqual(f["cleavage_factor"], 1.0)
        self.assertAlmostEqual(f["emission_factor"], 1.0)


if __name__ == "__main__":
    unittest.main()
