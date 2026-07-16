from __future__ import annotations

from pathlib import Path

import run_mpz_v9_18_5_5_mode_i_rcurve as mode_runner
import run_mpz_v9_18_5_5_persistent_plastic_wake as campaign_runner


def test_mode_i_runner_routes_to_v91855():
    text = Path(mode_runner.__file__).read_text()
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_5" in text


def test_campaign_runner_routes_to_v91855():
    text = Path(campaign_runner.__file__).read_text()
    assert "run_mpz_v9_18_5_5_mode_i_rcurve.py" in text
