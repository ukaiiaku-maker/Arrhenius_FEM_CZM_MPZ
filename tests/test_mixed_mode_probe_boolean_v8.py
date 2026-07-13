import tempfile
import unittest
from pathlib import Path

import pandas as pd

from calibrate_mixed_mode_loading_v8 import bool_col, read_probe, phase_sample_reliable


class ProbeBooleanTests(unittest.TestCase):
    def test_numeric_and_text_true_values(self):
        for value in (True, 1, 1.0, "1", "1.0", "true", "YES"):
            with self.subTest(value=value):
                self.assertTrue(bool_col(value))

    def test_false_values(self):
        for value in (False, 0, 0.0, "0", "0.0", "false", "no", "nan", None):
            with self.subTest(value=value):
                self.assertFalse(bool_col(value))

    def test_read_probe_numeric_boolean_remains_reliable(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            pd.DataFrame([
                {
                    "traction_phase_probe_reliable": 1.0,
                    "traction_probe_reliable": 0.0,
                    "reference_sigma_nn_Pa": 2.0,
                    "reference_tau_tn_Pa": -1.0,
                }
            ]).to_csv(out / "anisotropic_calibrated_tip_calls.csv", index=False)
            row = read_probe(out)
            self.assertTrue(phase_sample_reliable(row))

    def test_unreliable_zero_flag_is_rejected(self):
        row = {
            "traction_phase_probe_reliable": 0.0,
            "reference_sigma_nn_Pa": 2.0,
            "reference_tau_tn_Pa": -1.0,
        }
        self.assertFalse(phase_sample_reliable(row))


if __name__ == "__main__":
    unittest.main()
