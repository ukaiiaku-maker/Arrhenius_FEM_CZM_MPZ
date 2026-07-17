from __future__ import annotations

import audit_v91858_state_coupled_differentiation as audit


def _case(name: str, K_values, U_values):
    curve = []
    K0 = float(K_values[0])
    for i, (K, U) in enumerate(zip(K_values, U_values), start=1):
        curve.append({
            "x_um": 5.0 * i,
            "K": float(K),
            "U": float(U),
            "Kshield": 0.0,
            "K_normalized": float(K) / K0,
            "geometry_factor": float(K) / float(U),
            "shield_fraction": 0.0,
        })
    return {"class": name, "T_K": 700.0, "curve": curve}


def test_identical_scaled_curves_are_rejected():
    a = _case("weakT", [10.0, 11.0, 12.0], [1.0, 1.0, 1.0])
    b = _case("DBTT", [20.0, 22.0, 24.0], [2.0, 2.0, 2.0])
    row = audit._pair_metrics(a, b, 0.02, 0.01)
    assert row["normalized_shape_collapse"] is True
    assert row["maximum_normalized_K_separation"] == 0.0
    assert row["maximum_relative_geometry_factor_separation"] == 0.0


def test_state_dependent_shape_is_accepted():
    a = _case("weakT", [10.0, 11.0, 12.0], [1.0, 1.0, 1.0])
    b = _case("DBTT", [20.0, 24.0, 30.0], [2.0, 2.0, 2.0])
    row = audit._pair_metrics(a, b, 0.02, 0.01)
    assert row["normalized_shape_collapse"] is False
    assert row["maximum_normalized_K_separation"] > 0.02


def test_geometry_factor_difference_is_accepted_even_when_normK_is_close():
    a = _case("weakT", [10.0, 11.0, 12.0], [1.0, 1.0, 1.0])
    b = _case("DBTT", [20.0, 22.1, 24.1], [2.0, 2.2, 2.4])
    row = audit._pair_metrics(a, b, 0.02, 0.01)
    assert row["maximum_relative_geometry_factor_separation"] > 0.01
    assert row["normalized_shape_collapse"] is False
