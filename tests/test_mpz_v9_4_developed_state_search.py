import numpy as np
import pandas as pd

from search_mpz_v9_4_developed_state import (
    SHARED_COLUMNS,
    apply_shared,
    convergence_length,
    parse_temperatures,
    sample_shared_parameters,
)


def _base_row():
    return pd.Series({
        "mpz_source_sites_per_system": 40.0,
        "mpz_trap_barrier_eV": 0.65,
        "mpz_retained_recovery_barrier_eV": 1.50,
        "c_blunt": 0.50,
        "mpz_source_recovery_rate_s": 1.0e-5,
        "mpz_source_refresh_length_m": 2.0e-6,
        "mpz_length_m": 5.0e-5,
        "mpz_pair_annihilation_rate_per_count_s": 0.0,
    })


def test_temperature_parser_requires_all_three_classes():
    out = parse_temperatures(
        "ceramic:300,1200;weakT:300,700,1200;DBTT:300,900,1200"
    )
    assert out["ceramic"] == [300.0, 1200.0]
    assert out["weakT"][-1] == 1200.0
    assert out["DBTT"][1] == 900.0


def test_shared_samples_only_change_active_production_state_parameters():
    base = _base_row()
    rows = sample_shared_parameters(8, 94131, base)
    assert len(rows) == 13
    assert set(rows.columns) == set(SHARED_COLUMNS)
    assert np.all(rows["mpz_source_sites_per_system"] > 0)
    assert np.all(rows["mpz_length_m"] > 0)

    updated = apply_shared(base, rows.iloc[-1].to_dict())
    for name in SHARED_COLUMNS:
        assert float(updated[name]) == float(rows.iloc[-1][name])
    assert updated["mpz_pair_annihilation_rate_per_count_s"] == 0.0


def test_convergence_length_detects_plateau():
    events = []
    for i in range(30):
        value = 1.0 - np.exp(-i / 4.0)
        events.append({"a_um": 5.0 * i, "shield_fraction": value})
    length = convergence_length(
        events, "shield_fraction", relative_tolerance=0.12,
        absolute_tolerance=0.01,
    )
    assert np.isfinite(length)
    assert 20.0 <= length <= 145.0
