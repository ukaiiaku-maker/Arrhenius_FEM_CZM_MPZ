import unittest
from pathlib import Path


class CampaignTests(unittest.TestCase):
    def test_calibrated_tip_not_raw_probe_drives_barrier(self):
        s = Path("arrhenius_fracture/mixed_mode_first_passage_v5.py").read_text()
        self.assertIn("existing_sharp_front_sigma_tip_from_directionally_scaled_KJ", s)
        self.assertIn("finite_radius_traction_role", s)
        self.assertNotIn("class DirectTractionEngineMixin", s)

    def test_command_passes_reference_normalization(self):
        s = Path("run_mixed_mode_fem_czm_v5_campaign.py").read_text()
        for flag in ("--reference-cleavage-shape", "--reference-slip-shape",
                     "--target-traction-phase-deg", "--crystal-aniso",
                     "--crystal-compete"):
            self.assertIn(flag, s)

    def test_event_state_controller_present(self):
        s = Path("run_mixed_mode_fem_czm_v5_campaign.py").read_text()
        self.assertIn("max-control-iters", s)
        self.assertIn("mixed_mode_control_history_v5.csv", s)
        self.assertIn("event_phase_control_converged", s)

    def test_censor_aware_statuses(self):
        s = Path("run_mixed_mode_fem_czm_v5_campaign.py").read_text()
        self.assertIn("right_censored_phase_mismatch", s)
        self.assertIn("event_phase_mismatch", s)

    def test_barrier_audit_fingerprint(self):
        s = Path("run_mixed_mode_fem_czm_v5_campaign.py").read_text()
        self.assertIn("barrier_fingerprint_sha256", s)
        self.assertIn("barrier_audit.json", s)


if __name__ == "__main__":
    unittest.main()
