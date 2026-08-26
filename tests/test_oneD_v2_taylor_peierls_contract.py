from __future__ import annotations

import math

import numpy as np
import pytest

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from arrhenius_fracture.emergent_gnd_types_v912 import EV_TO_J, KB_J_PER_K
from arrhenius_fracture.zero_d_persistent_v913 import (
    ZeroDState,
    reduction_geometry,
    source_kinetic_diagnostics,
)
from reduced_fracture_v2.taylor_peierls_contract import (
    CLEAVAGE_FIELDS,
    EMISSION_FIELDS,
    FIXED_ATTEMPT_FREQUENCIES,
    FROZEN_ACTIVE_FIELDS,
    SEARCH_WHITELIST,
    SOURCE_FIXED_PROJECTION_FACTORS,
    barrier_hashes,
    full_material_hash,
    validate_candidate,
)
from tests.test_zero_d_persistent_v913 import candidate, physics


def _row() -> dict[str, float]:
    return {
        field: float(index + 1)
        for index, field in enumerate(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
    }


def test_only_taylor_peierls_whitelist_may_change() -> None:
    control = _row()
    varied = dict(control)
    for index, field in enumerate(SEARCH_WHITELIST):
        varied[field] = float(varied[field]) + 0.125 * (index + 1)
    validate_candidate(varied, control)
    assert set(FROZEN_ACTIVE_FIELDS) == (
        set(ACTIVE_CANDIDATE_PARAMETER_FIELDS) - set(SEARCH_WHITELIST)
    )
    assert set(FIXED_ATTEMPT_FREQUENCIES) <= set(FROZEN_ACTIVE_FIELDS)

    varied["rho_source0_m2"] = np.nextafter(control["rho_source0_m2"], np.inf)
    with pytest.raises(ValueError, match="non-whitelisted"):
        validate_candidate(varied, control)


def test_barrier_hashes_and_material_identity_are_bit_stable() -> None:
    control = _row()
    same = dict(control)
    assert barrier_hashes(same) == barrier_hashes(control)
    assert full_material_hash(same) == full_material_hash(control)
    same[CLEAVAGE_FIELDS[0]] = np.nextafter(control[CLEAVAGE_FIELDS[0]], np.inf)
    assert barrier_hashes(same)["cleavage_barrier_sha256"] != barrier_hashes(control)["cleavage_barrier_sha256"]
    same = dict(control)
    same[EMISSION_FIELDS[0]] = np.nextafter(control[EMISSION_FIELDS[0]], np.inf)
    assert barrier_hashes(same)["emission_barrier_sha256"] != barrier_hashes(control)["emission_barrier_sha256"]


def test_attempt_frequencies_and_projection_factors_are_source_fixed() -> None:
    model = candidate()
    assert model.peierls.nu0_s == 1.0e12
    assert model.taylor.nu0_s == 1.0e11
    assert model.peierls.stress_fraction == pytest.approx(1.0 / math.sqrt(3.0))
    assert model.taylor.stress_fraction == pytest.approx(1.0 / math.sqrt(3.0))
    assert SOURCE_FIXED_PROJECTION_FACTORS == {
        "peierls_stress_fraction": 1.0 / math.sqrt(3.0),
        "taylor_stress_fraction": 1.0 / math.sqrt(3.0),
    }


def test_timescale_diagnostics_reproduce_source_owned_rates() -> None:
    model = candidate()
    common = physics()
    reduced = reduction_geometry(common)
    state = ZeroDState(
        mobile_m2=np.asarray([2.0e12, 3.0e12]),
        retained_m2=np.asarray([5.0e12, 7.0e12]),
        local_slip_count_by_system=np.asarray([2.0, 3.0]),
        cumulative_activations=np.zeros(2),
    )
    temperature = 1000.0
    diagnostic = source_kinetic_diagnostics(
        model,
        common,
        reduced,
        state,
        K_MPa_sqrt_m=20.0,
        temperature_K=temperature,
    )
    external = np.asarray(diagnostic["resolved_emission_drive_Pa_by_system"])
    forest = common.rho_forest_floor_m2 + np.sum(state.mobile_m2 + state.retained_m2)
    spacing = 1.0 / (2.0 * math.sqrt(forest))
    jump = common.jump_fraction_of_forest_spacing * spacing
    p_stress = model.peierls.stress_fraction * external
    p_barrier = model.peierls.surface(model.emission).barrier_eV(p_stress, temperature)
    p_rate = model.peierls.nu0_s * np.exp(
        np.clip(-np.asarray(p_barrier) * EV_TO_J / (KB_J_PER_K * temperature), -700.0, 0.0)
    )
    velocity = jump * p_rate
    mfp = common.mean_free_path_coefficient / math.sqrt(forest)
    encounter = common.encounter_efficiency * np.abs(velocity) / mfp
    assert np.asarray(diagnostic["peierls_rate_s_by_system"]) == pytest.approx(p_rate)
    assert np.asarray(diagnostic["peierls_velocity_m_s_by_system"]) == pytest.approx(velocity)
    assert np.asarray(diagnostic["encounter_rate_s_by_system"]) == pytest.approx(encounter)
    assert np.asarray(diagnostic["chi_ret_by_system"]) == pytest.approx(
        np.asarray(diagnostic["transport_time_s_by_system"]) * encounter
    )
    assert np.asarray(diagnostic["chi_taylor_completion_by_system"]) == pytest.approx(
        np.asarray(diagnostic["transport_time_s_by_system"])
        * np.asarray(diagnostic["taylor_completion_rate_s_by_system"])
    )
