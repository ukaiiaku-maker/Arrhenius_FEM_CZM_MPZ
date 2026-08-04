"""Real-engine integration tests for the PF KJ-target controller.

These tests do not require the (currently unavailable, see
CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER" and PF_REFERENCE_REGENERATION_CONTRACT.md)
real PF reference artifacts. They exercise the actual `sharp_front.run_2d`
solver -- real mesh, real FEM assembly/solve, real plasticity update -- via
a small, fast, synthetic target trajectory, proving the controller wiring
(FEM_CZM_HANDOFF.md section 6, Gate 1 of FEM_PF_PARITY_SCORECARD.md) is
transactional and numerically sane against the real engine, not just
against the synthetic callbacks in test_pf_theta0_j_controlled_loading_v10051840.py.

This is Gate 1 evidence only (numerical controller correctness). It says
nothing about physical FEM/PF parity, which requires the real PF reference
data (Gate 2 onward).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from arrhenius_fracture import sharp_front

MESH_ARGS = ["--nx", "6", "--ny", "10"]


def _write_synthetic_target(path: Path, *, n: int = 10, dt_s: float = 8.4) -> None:
    step = np.arange(1, n + 1, dtype=float)
    KJ = np.linspace(1.0e5, 5.0e6, n)
    J = np.linspace(10.0, 200.0, n)
    frame = pd.DataFrame(
        {
            "step": step,
            "dt_cur_s": np.full(n, dt_s),
            "J_effective_direct_J_per_m2": J,
            "J_signed_direct_J_per_m2": J,
            "KJ_Pa_sqrtm": KJ,
            "B": np.zeros(n),
            "N_em": np.zeros(n),
            "crack_extension_m": np.zeros(n),
        }
    )
    frame.to_csv(path, index=False)


def _run(tmp_path: Path, *, steps: int, extra_args: list[str] | None = None) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    csv_path = tmp_path / "steps_1000K.csv"
    _write_synthetic_target(csv_path)
    out_dir = tmp_path / "out"
    argv = [
        "--mode", "2d",
        "--steps", str(steps),
        *MESH_ARGS,
        "--temperatures", "1000",
        "--out", str(out_dir),
        "--print-every", "1",
        "--no-plots",
        "--dt", "8.4",
    ]
    if extra_args is None:
        extra_args = ["--pf-kj-target-csv", str(csv_path)]
    sharp_front.main(argv + extra_args)
    return out_dir


def _load_audit(out_dir: Path) -> dict:
    audit_path = out_dir / "pf_kj_target_controller_audit_1000K.json"
    assert audit_path.is_file(), f"expected controller audit file at {audit_path}"
    return json.loads(audit_path.read_text())


def test_real_engine_controller_tracks_synthetic_target_within_tolerance(tmp_path: Path):
    out_dir = _run(tmp_path, steps=3)
    audit = _load_audit(out_dir)

    assert audit["schema"] == "pf_kj_target_controller_audit_v10_0_5_18_4_0"
    records = audit["records"]
    assert len(records) == 3

    for i, record in enumerate(records, start=1):
        assert record["step"] == i
        assert record["converged"] is True
        assert 0 <= record["iterations"] <= 25  # default max_iterations
        target = record["KJ_target_Pa_sqrtm"]
        achieved = record["KJ_achieved_Pa_sqrtm"]
        assert achieved == pytest.approx(target, rel=2.0e-3)
        # physical time is inherited from the PF interval (8.4 s/row), not
        # from solver step count -- both fields must agree exactly with the
        # cumulative dt, matching FEM_CZM_HANDOFF.md section 6's requirement
        # that physical time intervals are preserved.
        assert record["physical_time_target_s"] == pytest.approx(8.4 * i, rel=1.0e-12)
        assert record["physical_time_accepted_s"] == pytest.approx(8.4 * i, rel=1.0e-12)


def test_real_engine_rejected_trials_are_recorded_and_precede_the_accepted_value(tmp_path: Path):
    out_dir = _run(tmp_path, steps=1)
    audit = _load_audit(out_dir)
    record = audit["records"][0]

    # A synthetic target starting from Uapp=0 requires more than one trial
    # (the elastic predictor alone under/overshoots off a cold start), so
    # this configuration is expected to exercise the safeguarded-secant
    # rejection path for real, not just converge in one shot.
    assert record["iterations"] >= 2
    assert len(record["rejected_trials"]) == record["iterations"] - 1
    for trial in record["rejected_trials"]:
        assert trial["KJ_achieved_Pa_sqrtm"] != pytest.approx(
            record["KJ_target_Pa_sqrtm"], rel=2.0e-3
        )


def test_real_engine_run_is_deterministic_across_repeats(tmp_path: Path):
    """No cross-trial or cross-run state leakage: an identical invocation
    must reproduce the identical accepted controller trajectory bit-for-bit
    (this front-engine/mesh configuration has no RNG in its path, so this
    is a legitimate exact-equality check, not merely a tolerance check)."""
    out_a = _run(tmp_path / "a", steps=3)
    out_b = _run(tmp_path / "b", steps=3)
    audit_a = _load_audit(out_a)
    audit_b = _load_audit(out_b)

    for record_a, record_b in zip(audit_a["records"], audit_b["records"]):
        assert record_a["Uapp_accepted_m"] == record_b["Uapp_accepted_m"]
        assert record_a["KJ_achieved_Pa_sqrtm"] == record_b["KJ_achieved_Pa_sqrtm"]
        assert record_a["iterations"] == record_b["iterations"]
        assert record_a["rejected_trials"] == record_b["rejected_trials"]


def test_flag_off_path_writes_no_controller_audit_file(tmp_path: Path):
    out_dir = _run(tmp_path, steps=3, extra_args=[])
    audit_path = out_dir / "pf_kj_target_controller_audit_1000K.json"
    assert not audit_path.exists()
