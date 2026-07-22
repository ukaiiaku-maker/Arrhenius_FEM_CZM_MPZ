"""Regression tests for v9.12 reduced work/energy diagnostics."""
from __future__ import annotations

import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CandidateParameters,
    CommonPhysics,
    EmergentGNDState,
    ExpFloorSurface,
    PTMechanism,
    ProtocolSegment,
    run_temperature_protocol,
)


def fast_candidate() -> CandidateParameters:
    cleavage = ExpFloorSurface(2.0, 0.0, 4.0e9, 0.0, 1.0, 1.0, 0.05)
    emission = ExpFloorSurface(0.05, 0.0, 1.0e9, 0.0, 2.0, 1.0, 0.01)
    return CandidateParameters(
        candidate_id="energy_test",
        cleavage=cleavage,
        emission=emission,
        peierls=PTMechanism(0.05, 0.0, 2.0, 1.0, 1.0e12),
        taylor=PTMechanism(0.10, 0.0, 2.0, 1.0, 1.0e11),
        rho_source0_m2=1.0e14,
        source_refresh_length_m=5.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
        recovery_nu0_s=1.0e12,
        recovery_H0_eV=0.05,
    )


def test_power_snapshot_obeys_stress_decomposition():
    physics = CommonPhysics(
        n_bins=3,
        n_systems=2,
        mpz_length_m=3.0e-6,
        source_zone_length_m=1.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    state.mobile_m2[:] = 2.0e13

    velocity = np.asarray(
        [[2.0e-3, 1.0e-3, -1.0e-3], [1.5e-3, 0.5e-3, -0.5e-3]]
    )
    tau_external = np.asarray(
        [[100.0e6, 80.0e6, -60.0e6], [90.0e6, 70.0e6, -50.0e6]]
    )
    tau_shielding = np.asarray(
        [[-12.0e6, -8.0e6, 6.0e6], [-10.0e6, -7.0e6, 5.0e6]]
    )
    tau_gnd = np.asarray(
        [[10.0e6, -5.0e6, 4.0e6], [-8.0e6, 3.0e6, -2.0e6]]
    )
    snapshot = state._plastic_power_snapshot(
        {
            "velocity_m_s": velocity,
            "tau_external_Pa": tau_external,
            "tau_nonlocal_shielding_Pa": tau_shielding,
            "tau_gnd_Pa": tau_gnd,
            "tau_eff_Pa": tau_external + tau_shielding + tau_gnd,
        }
    )

    assert np.isclose(
        snapshot["effective_plastic_power_J_per_m_s"],
        snapshot["external_plastic_power_J_per_m_s"]
        + snapshot["nonlocal_shielding_power_J_per_m_s"]
        + snapshot["internal_stress_power_J_per_m_s"],
    )
    assert snapshot["effective_plastic_dissipation_J_per_m_s"] >= 0.0


def test_spatial_advance_accumulates_finite_diagnostic_work():
    physics = CommonPhysics(
        n_bins=8,
        n_systems=2,
        mpz_length_m=8.0e-6,
        source_zone_length_m=2.0e-6,
    )
    state = EmergentGNDState(fast_candidate(), physics)
    state.max_feedback_substep_s = 0.025
    state.coupled_operator_substeps = 2

    state.advance_time(0.10, 20.0, 900.0)
    state.translate_tip(1.0e-6)
    diag = state.diagnostics(0.8, 20.0, 900.0)

    for key in (
        "external_plastic_work_J_per_m",
        "nonlocal_shielding_work_J_per_m",
        "internal_stress_work_J_per_m",
        "effective_plastic_work_J_per_m",
        "effective_plastic_dissipation_J_per_m",
        "external_plastic_work_per_crack_area_J_m2",
        "effective_plastic_dissipation_per_crack_area_J_m2",
        "mobile_line_energy_J_per_m",
        "retained_line_energy_J_per_m",
        "total_line_energy_J_per_m",
    ):
        assert np.isfinite(diag[key])
    assert diag["effective_plastic_dissipation_J_per_m"] >= 0.0
    assert diag["mobile_line_energy_J_per_m"] >= 0.0
    assert diag["retained_line_energy_J_per_m"] >= 0.0
    assert np.isclose(
        diag["total_line_energy_J_per_m"],
        diag["mobile_line_energy_J_per_m"]
        + diag["retained_line_energy_J_per_m"],
    )
    assert np.isclose(
        diag["effective_plastic_work_J_per_m"],
        diag["external_plastic_work_J_per_m"]
        + diag["nonlocal_shielding_work_J_per_m"]
        + diag["internal_stress_work_J_per_m"],
        rtol=1.0e-10,
        atol=1.0e-20,
    )


def test_temperature_result_serializes_energy_histories():
    physics = CommonPhysics(
        n_bins=8,
        n_systems=2,
        mpz_length_m=8.0e-6,
        source_zone_length_m=2.0e-6,
    )
    protocol = [
        ProtocolSegment(0.0, 1.0e-6, 10.0, 12.0, 0.05),
        ProtocolSegment(1.0e-6, 2.0e-6, 12.0, 14.0, 0.05),
    ]
    result = run_temperature_protocol(
        fast_candidate(),
        physics,
        protocol,
        900.0,
    )
    payload = result.as_dict()

    for key in (
        "external_plastic_work_J_per_m",
        "nonlocal_shielding_work_J_per_m",
        "internal_stress_work_J_per_m",
        "effective_plastic_work_J_per_m",
        "effective_plastic_dissipation_J_per_m",
        "external_plastic_work_per_crack_area_J_m2",
        "effective_plastic_dissipation_per_crack_area_J_m2",
        "mobile_line_energy_J_per_m",
        "retained_line_energy_J_per_m",
        "total_line_energy_J_per_m",
    ):
        assert key in payload
        assert len(payload[key]) == len(protocol)
        assert np.all(np.isfinite(payload[key]))

    metadata = payload["numerical_integration"]
    assert metadata["energy_bookkeeping"] == (
        "reduced_1d_orowan_power_and_log_line_energy_v1"
    )
    assert metadata["energy_bookkeeping_feedback_active"] is False
