import csv
import importlib.util
import math
from pathlib import Path

import numpy as np

MODULE = Path(__file__).parents[1] / "scripts" / "analyze_dbtt1000_k_decomposition.py"
spec = importlib.util.spec_from_file_location("dbtt_kdiag", MODULE)
diag = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diag)


def _row(step, u, f, k, ext=0.0, fire=0.0):
    return {
        "step": float(step),
        "Uapp_m": float(u),
        "Ftop_N": float(f),
        "KJ_Pa_sqrtm": float(k),
        "sigma_tip_Pa": 2.0e9,
        "sigma_back_Pa": 1.0e9,
        "crack_extension_m": float(ext),
        "n_fire": float(fire),
        "N_em": 5.0,
        "mpz_K_shield_Pa_sqrt_m": 0.0,
    }


def test_shared_fresh_reference_and_remote_apparent_k():
    rows = [_row(1, 1e-6, 2e5, 4e6), _row(2, 2e-6, 4e5, 4.1e6, 1e-5, 1)]
    match = diag.first_row_match(rows[0], dict(rows[0]))
    assert max(match.values()) == 0.0
    out = diag.process_model(
        "x", rows, a0_m=5e-4, width_m=2e-3,
        k_over_f0=20.0, reference_calibration=None,
    )
    assert math.isclose(out[0]["K_app_initial_geometry_Pa_sqrtm"], 4e6)
    assert math.isclose(out[1]["K_app_initial_geometry_Pa_sqrtm"], 8e6)
    assert out[1]["K_app_SENT_Pa_sqrtm"] > out[1]["K_app_initial_geometry_Pa_sqrtm"]


def test_event_mask_uses_fire_or_extension_change():
    rows = [
        _row(1, 1e-6, 2e5, 4e6),
        _row(2, 1e-6, 2e5, 4e6, 1e-5, 0),
        _row(3, 1e-6, 2e5, 4e6, 1e-5, 1),
    ]
    mask = diag.event_mask(rows)
    assert mask.tolist() == [False, True, True]


def test_tip_radius_proxy_matches_analytic_relation():
    radius = 3e-6
    sigma = 5e9
    k = sigma*np.sqrt(2*np.pi*radius)
    r = _row(1, 1e-6, 2e5, k)
    r["sigma_tip_Pa"] = sigma
    assert math.isclose(diag.infer_tip_radius(r), radius, rel_tol=1e-12)


def test_reference_calibration_interpolation(tmp_path):
    path = tmp_path / "cal.csv"
    with path.open("w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["projected_crack_length_m", "K_over_F_Pa_sqrtm_per_N"])
        w.writerow([5e-4, 20.0])
        w.writerow([1.5e-3, 40.0])
    cal = diag.load_reference_calibration(path)
    values = diag.interp_reference(np.asarray([5e-4, 1e-3, 1.5e-3]), cal)
    assert np.allclose(values, [20.0, 30.0, 40.0])
