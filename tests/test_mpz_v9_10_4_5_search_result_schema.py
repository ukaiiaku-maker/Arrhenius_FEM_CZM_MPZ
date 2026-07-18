from __future__ import annotations

import numpy as np
from scipy.optimize import OptimizeResult

import optimize_mpz_v9_10_4_5_narrow_dbtt as patched


def test_incomplete_detailed_result_supplies_empty_event_detail() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    result = patched.stabilize_detailed_result(
        {
            "objective": 1.1e6,
            "completion_loss": 1.0e5,
            "parameters": {"source_sites_per_system": 1.0},
            "temperature_detail": [{"T_K": 300.0, "completed": False}],
        },
        x,
        details=True,
    )
    assert result["event_detail"] == []
    assert result["temperature_detail"][0]["T_K"] == 300.0
    assert result["evaluation_status"] == "INCOMPLETE_CANDIDATE"


def test_early_rejection_has_stable_detailed_schema() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    result = patched.stabilize_detailed_result(
        {
            "objective": 1.0e8,
            "parameters": {},
        },
        x,
        details=True,
    )
    assert result["parameters"] == {}
    assert result["temperature_detail"] == []
    assert result["event_detail"] == []
    assert result["evaluation_status"] == "EARLY_REJECTED_CANDIDATE"


def test_non_detailed_result_is_not_expanded() -> None:
    x = np.zeros(len(patched._base.PARAMETER_NAMES), dtype=float)
    original = {"objective": 42.0}
    result = patched.stabilize_detailed_result(original, x, details=False)
    assert result == original
    assert "event_detail" not in result


def test_de_checkpoint_roundtrip(tmp_path) -> None:
    state = tmp_path / "restart_000_de_state.npz"
    original = OptimizeResult(
        x=np.asarray([1.0, 2.0]),
        fun=3.5,
        population=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
        population_energies=np.asarray([3.5, 7.5]),
        nit=2,
        nfev=16,
    )
    patched._save_de_result(state, original)
    loaded = patched._load_de_result(state)
    assert np.allclose(loaded.x, original.x)
    assert np.allclose(loaded.population, original.population)
    assert np.allclose(loaded.population_energies, original.population_energies)
    assert loaded.fun == original.fun
    assert loaded.nit == original.nit
    assert loaded.nfev == original.nfev
