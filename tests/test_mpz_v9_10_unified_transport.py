import numpy as np
import pandas as pd

from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.moving_process_zone_v910 import MovingProcessZoneState
from optimize_mpz_v9_10_unified_global import (
    BOUNDS,
    PARAMETER_NAMES,
    UnifiedObjective,
    ZeroDSettings,
    bounds_list,
    decode,
)


def midpoint_vector():
    return np.array([0.5 * (lo + hi) for lo, hi in bounds_list()], dtype=float)


def test_ordered_absolute_barriers_and_full_dimension():
    x = midpoint_vector()
    p = decode(x)
    assert len(PARAMETER_NAMES) == len(BOUNDS) == len(x) == 25
    assert p["peierls_H0_eV"] < p["taylor_H0_eV"]
    assert p["peierls_nu0_s"] == 1.0e12
    assert p["taylor_nu0_s"] == 1.0e11


def test_encounter_requires_peierls_transport():
    rho = np.array([1.0e12, 1.0e16])
    zero = MovingProcessZoneState.encounter_rate_s(0.0, 1.0e-8, rho, 10.0)
    assert np.all(zero == 0.0)
    driven = MovingProcessZoneState.encounter_rate_s(2.0, 1.0e-8, rho, 1.0)
    assert np.all(driven > 0.0)
    assert driven[1] > driven[0]


def test_exact_mobile_retained_exchange_conserves_population():
    mobile = np.array([[10.0, 4.0]])
    retained = np.array([[0.0, 6.0]])
    m2, r2, trapped, released = MovingProcessZoneState._exchange_mobile_retained(
        mobile, retained,
        np.array([2.0, 0.5]),
        np.array([0.5, 2.0]),
        0.25,
    )
    assert np.allclose(m2 + r2, mobile + retained)
    assert trapped >= 0.0
    assert released >= 0.0


def make_cfg(trap_barrier):
    cfg = MovingProcessZoneConfig(length_m=5.0e-6, n_bins=20, n_systems=2)
    cfg.pt_emit_G00_eV = 1.0
    cfg.pt_emit_gT_eV_per_K = 0.0
    cfg.pt_emit_sigc0_Pa = 2.0e9
    cfg.pt_emit_sT_Pa_per_K = 0.0
    cfg.pt_peierls_energy_ratio = 0.2
    cfg.pt_peierls_entropy_ratio = 0.0
    cfg.pt_peierls_nu0_s = 1.0e12
    cfg.pt_taylor_energy_ratio = 2.0
    cfg.pt_taylor_entropy_ratio = 0.0
    cfg.pt_taylor_nu0_s = 1.0e11
    cfg.pt_taylor_m_cap = float("inf")
    cfg.pt_mobile_saturation_density_m2 = float("inf")
    cfg.pt_mobile_density_floor_m2 = 0.0
    cfg.pt_jump_length_min_m = 0.0
    cfg.pt_taylor_phi_max = float("inf")
    cfg.pt_encounter_efficiency = 1.0
    cfg.trap_barrier_eV = trap_barrier
    cfg.retained_recovery_nu0_s = 0.0
    cfg.mobile_recovery_rate_s = 0.0
    return cfg


def test_legacy_trap_barrier_is_inactive():
    a = MovingProcessZoneState(make_cfg(0.01))
    b = MovingProcessZoneState(make_cfg(20.0))
    a.mobile[:, :2] = 2.0
    b.mobile[:, :2] = 2.0
    out_a = a.evolve(1.0e-4, 900.0, 2.0e9, 2.74e-10)
    out_b = b.evolve(1.0e-4, 900.0, 2.0e9, 2.74e-10)
    assert np.allclose(a.mobile, b.mobile)
    assert np.allclose(a.retained, b.retained)
    assert out_a["legacy_trap_barrier_active"] == 0.0
    assert out_b["legacy_trap_barrier_active"] == 0.0


def test_unified_objective_returns_finite_penalty_or_score():
    targets = pd.DataFrame({
        "target_class": ["weakT", "weakT"],
        "T_K": [300.0, 700.0],
        "K_init_target": [15.0, 15.0],
        "K_init_scale": [2.0, 2.0],
        "K_plateau_target": [20.0, 20.0],
        "K_plateau_scale": [3.0, 3.0],
        "early_rise_per_100um_target": [1.0, 1.0],
        "early_rise_scale": [0.5, 0.5],
        "plateau_rise_per_100um_target": [0.0, 0.0],
        "plateau_rise_scale": [0.3, 0.3],
        "delta_KR_min": [2.0, 2.0],
        "delta_KR_max": [9.0, 9.0],
        "weight": [1.0, 1.0],
    })
    objective = UnifiedObjective(ZeroDSettings(
        target_class="weakT",
        temperatures=np.array([300.0, 700.0]),
        targets=targets,
        dK=2.0,
        Kmax=20.0,
        target_extension_um=20.0,
        da_um=5.0,
    ))
    value = objective(midpoint_vector())
    assert np.isfinite(value)
    assert value >= 0.0
