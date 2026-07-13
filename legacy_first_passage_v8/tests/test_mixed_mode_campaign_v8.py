import unittest
from pathlib import Path


class CampaignTests(unittest.TestCase):
    def test_calibration_invokes_exact_backend(self):
        s = Path("calibrate_mixed_mode_loading_v8.py").read_text()
        self.assertIn('"--crack-backend", "adaptive_czm"', s)
        self.assertIn('"--steps", "1"', s)
        self.assertIn("solve_exact_target", s)
        self.assertIn("first_production_step_verified", s)

    def test_phase_reliability_is_separate(self):
        s = Path("arrhenius_fracture/mixed_mode_first_passage_v8.py").read_text()
        self.assertIn("traction_phase_probe_reliable", s)
        self.assertIn("directional_metrics_reliable", s)

    def test_full_circle_coefficients_are_passed(self):
        s = Path("run_mixed_mode_fem_czm_v8_campaign.py").read_text()
        self.assertIn("loading_coefficients_from_alpha_deg", s)
        self.assertIn("selected_loading_alpha_deg", s)
        self.assertNotIn("loading_coefficients_from_z", s)
        self.assertNotIn("selected_loading_z", s)

    def test_event_controller_uses_empirical_alpha(self):
        s = Path("run_mixed_mode_fem_czm_v8_campaign.py").read_text()
        self.assertIn("safeguarded_event_alpha_update", s)
        self.assertIn("mixed_mode_control_history_v8.csv", s)

    def test_barriers_and_censoring_retained(self):
        s = Path("run_mixed_mode_fem_czm_v8_campaign.py").read_text()
        self.assertIn("barrier_fingerprint_sha256", s)
        self.assertIn("right_censored_phase_mismatch", s)
        self.assertIn("event_phase_mismatch", s)


if __name__ == "__main__":
    unittest.main()
