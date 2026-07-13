import math
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from arrhenius_fracture.mixed_mode_first_passage_v8 import (
    directional_drive_factors,
    loading_alpha_deg_from_coefficients,
    loading_coefficients_from_alpha_deg,
    loading_coefficients_from_response_basis,
    normalize_loading_coefficients,
    phase_from_response_alpha,
    safeguarded_event_alpha_update,
)
import calibrate_mixed_mode_loading_v8 as cal


class MechanicsTests(unittest.TestCase):
    def test_loading_coefficients_normalized(self):
        q = normalize_loading_coefficients(3.0, 4.0)
        self.assertAlmostEqual(math.hypot(*q), 1.0)

    def test_response_basis_recovers_target(self):
        M = np.array([[2.0, 0.2], [0.1, 1.5]])
        q = loading_coefficients_from_response_basis(M, 30.0)
        r = M @ np.array(q)
        self.assertAlmostEqual(math.degrees(math.atan2(r[1], r[0])), 30.0, places=8)

    def test_reported_basis_retains_negative_opening(self):
        M = np.array([[3.81221929e7, 5.20675615e5],
                      [1.91804414e4, 2.01374690e5]])
        q = loading_coefficients_from_response_basis(M, 30.0)
        self.assertLess(q[0], 0.0)
        r = M @ np.array(q)
        self.assertAlmostEqual(math.degrees(math.atan2(r[1], r[0])), 30.0, places=6)

    def test_full_circle_round_trip_negative_opening(self):
        q = normalize_loading_coefficients(-0.0045, 0.99999)
        a = loading_alpha_deg_from_coefficients(*q)
        q2 = loading_coefficients_from_alpha_deg(a)
        self.assertTrue(np.allclose(q, q2, atol=1e-12))
        self.assertGreater(a, 90.0)

    def test_phase_from_response_alpha(self):
        M = np.eye(2)
        self.assertAlmostEqual(phase_from_response_alpha(M, 135.0), 135.0)

    def test_event_update_bracketed(self):
        M = np.eye(2)
        samples = [
            {"loading_alpha_unwrapped_deg": -20.0, "achieved_psi_deg": -20.0},
            {"loading_alpha_unwrapped_deg": 20.0, "achieved_psi_deg": 20.0},
        ]
        a = safeguarded_event_alpha_update(samples, 0.0, M)
        self.assertAlmostEqual(a, 0.0)

    def test_mode_I_factor_normalization(self):
        f = directional_drive_factors(.8, .1, .8, .1, 1.0, 5.0)
        self.assertAlmostEqual(f["cleavage_factor"], 1.0)
        self.assertAlmostEqual(f["emission_factor"], 1.0)

    def test_exact_backend_root_recovers_negative_branch(self):
        # Mimic the reported v6 failure: the basis guess near -88.5 deg gives
        # only about -4.3 deg, while the exact root is near -91.2 deg.
        root = -91.2

        def fake_run_probe(py, a, qo, qs, target, out, ref_c=1.0, ref_s=0.0, shear_sign=1.0):
            alpha = math.degrees(math.atan2(qs, qo))
            # unwrap near the desired negative branch
            if alpha > 0:
                alpha -= 360.0
            phase = target + 9.5 * (alpha - root)
            phase = max(min(phase, 170.0), -170.0)
            mag = 1.0e6
            pr = math.radians(phase)
            return {
                "traction_phase_probe_reliable": True,
                "traction_probe_reliable": False,  # directional gate must not reject phase
                "reference_sigma_nn_Pa": mag * math.cos(pr),
                "reference_tau_tn_Pa": mag * math.sin(pr),
                "cleavage_shape": 0.8,
                "slip_shape": 0.1,
            }

        a = types.SimpleNamespace(
            psi_tol_deg=0.75,
            max_root_iters=20,
            max_alpha_step_deg=20.0,
        )
        with tempfile.TemporaryDirectory() as td, patch.object(cal, "run_probe", fake_run_probe):
            selected, hist = cal.solve_exact_target(
                "python", a, -30.0, -88.5, Path(td), 0.8, 0.1, 1.0)
        self.assertTrue(selected["phase_sample_reliable"])
        self.assertLess(abs(selected["traction_phase_error_deg"]), 0.75)
        self.assertAlmostEqual(selected["loading_alpha_unwrapped_deg"], root, delta=0.2)
        self.assertGreaterEqual(len(hist), 2)


if __name__ == "__main__":
    unittest.main()
