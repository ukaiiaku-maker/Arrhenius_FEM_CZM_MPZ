import unittest
from pathlib import Path


class CampaignTests(unittest.TestCase):
    def test_calibration_invokes_exact_backend(self):
        s = Path("calibrate_mixed_mode_loading_v6.py").read_text()
        self.assertIn('"--crack-backend", "adaptive_czm"', s)
        self.assertIn('"--steps", "1"', s)
        self.assertIn("first_production_step_verified", s)
        self.assertNotIn("d[bnd.notch_nodes]", s)

    def test_direct_coefficients_are_passed(self):
        s = Path("run_mixed_mode_fem_czm_v6_campaign.py").read_text()
        self.assertIn("--mixity-open-coeff", s)
        self.assertIn("--mixity-shear-coeff", s)
        self.assertNotIn("--mixity-loading-angle-deg", s)

    def test_event_controller_uses_empirical_z(self):
        s = Path("run_mixed_mode_fem_czm_v6_campaign.py").read_text()
        self.assertIn("safeguarded_event_z_update", s)
        self.assertIn("mixed_mode_control_history_v6.csv", s)
        self.assertIn("selected_loading_z", s)

    def test_barriers_and_censoring_retained(self):
        s = Path("run_mixed_mode_fem_czm_v6_campaign.py").read_text()
        self.assertIn("barrier_fingerprint_sha256", s)
        self.assertIn("right_censored_phase_mismatch", s)
        self.assertIn("event_phase_mismatch", s)


if __name__ == "__main__":
    unittest.main()
