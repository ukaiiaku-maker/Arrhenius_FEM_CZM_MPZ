from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.mode_i_first_passage_v9_14 import _forward_mode_i_plane
from arrhenius_fracture.remesh_audit_v914 import audit_case
import run_mpz_v9_14_event_driven_remesh as runner


def test_forward_plane_is_exact_mode_i():
    plane = _forward_mode_i_plane(45.0)[0]
    assert list(plane["t"]) == [1.0, 0.0]
    assert list(plane["n"]) == [0.0, 1.0]


def test_v914_command_selects_event_remesh_driver():
    class Args:
        parameter_root = Path("mpz_v9_11_parameters")
        T_K = 700.0
        target_extension_um = 25.0
        steps = 100
        nx = 18; ny = 36
        tip_h_fine = 1e-6; tip_ratio = 1.2
        dU = 2e-7; dt = 8.4
        n_stagger = 2; print_every = 25
        adaptive_event_target = 0.15; da_phys_um = 5.0
        mpz_length_um = 100.0; mpz_n_bins = 200
        crystal_theta_deg = 45.0
        save_snapshots = 2; snapshot_cols = 2; snapshot_by_extension_um = 10.0
    cmd = runner._build_command_v914(Args(), "DBTT", Path("runs/x"), True)
    assert cmd[1] == "run_mpz_v9_14_mode_i_rcurve.py"


def test_audit_requires_same_load_post_event_equilibrium(tmp_path):
    pd.DataFrame({
        "step": [1, 2, 3],
        "Uapp_m": [1.0e-4, 1.1e-4, 1.2e-4],
        "n_fire": [0, 1, 0],
        "adaptive_frac": [1.0, 1e-8, 1.0],
    }).to_csv(tmp_path / "steps_0700K.csv", index=False)
    (tmp_path / "czm_advance_log.json").write_text(json.dumps([
        {"length_m": 5e-6, "x0": 0.0, "x1": 5e-6}
    ]))
    (tmp_path / "field_snapshot_manifest_700K.json").write_text("{}")
    (tmp_path / "cohesive_elements.csv").write_text("index\n")
    out = audit_case(tmp_path, 700.0)
    assert out["requirements_1_to_4_passed"]
    assert not out["requirement_5_passed"]
    assert not out["all_five_requirements_passed"]
