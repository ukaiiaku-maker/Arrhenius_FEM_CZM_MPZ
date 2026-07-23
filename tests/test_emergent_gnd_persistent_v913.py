import math

import numpy as np
import pytest

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_state_v913 import (
    EmergentGNDState,
    effective_front_width_m,
    persistent_site_multiplicity,
    solve_backstress_limited_activations,
)
from arrhenius_fracture.emergent_gnd_types_v912 import ExpFloorSurface, PTMechanism
from arrhenius_fracture.emergent_gnd_types_v913 import CandidateParameters, CommonPhysics


def candidate() -> CandidateParameters:
    surface = ExpFloorSurface(
        G00_eV=0.4,
        gT_eV_per_K=0.0,
        sigc0_Pa=2.0e9,
        sT_Pa_per_K=0.0,
        exp_a=1.0,
        exp_n=1.0,
        floor_fraction=0.01,
    )
    return CandidateParameters(
        candidate_id="test",
        cleavage=surface,
        emission=surface,
        peierls=PTMechanism(0.05, 0.0, 1.0, 1.0, 1.0e8),
        taylor=PTMechanism(0.10, 0.0, 1.0, 1.0, 1.0e8),
        rho_source0_m2=1.4115242646890916e16,
        source_refresh_length_m=1.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
        recovery_nu0_s=0.0,
        c_blunt=1.411283192139077,
    )


def physics(**updates) -> CommonPhysics:
    values = dict(
        n_bins=8,
        mpz_length_m=8.0e-6,
        source_zone_length_m=2.0e-6,
        blunting_length_m=1.0e-6,
        activation_to_line_content_per_system=(1.0, 1.0),
    )
    values.update(updates)
    result = CommonPhysics(**values)
    result.validate()
    return result


def test_reference_width_and_multiplicity_match_v10221():
    width = effective_front_width_m(
        5.0e12,
        reference_width_m=10.0e-6,
        reference_density_m2=5.0e12,
        minimum_width_m=0.625e-6,
        maximum_width_m=50.0e-6,
    )
    assert width == pytest.approx(10.0e-6)
    arc = 25.0e-12 / (1.0e-6 * 10.0e-6)
    multiplicity = persistent_site_multiplicity(
        1.4115242646890916e16,
        1.0e-6,
        width,
        arc,
    )
    assert multiplicity == pytest.approx(3.528810661722729e5)


def test_implicit_root_blocks_without_inventory():
    prefactor = 50.0
    value = solve_backstress_limited_activations(
        multiplicity=3.5e5,
        dt_s=1.0,
        drive_stress_Pa=1.0e9,
        rho_initial_m2=1.0e14,
        rho_increment_per_activation_m2=1.0e11,
        backstress_prefactor_Pa_sqrt_m2=prefactor,
        rate_function=lambda sigma: 1.0e4 * math.exp(sigma / 1.0e9),
    )
    block = (((1.0e9 / prefactor) ** 2) - 1.0e14) / 1.0e11
    assert 0.0 < value < block


def test_implicit_root_accepts_exact_mechanical_blocking_boundary():
    prefactor = 100.0
    rho_increment = 1.0e12
    drive = 1.0e9
    block = (drive / prefactor) ** 2 / rho_increment
    value = solve_backstress_limited_activations(
        multiplicity=1.0e6,
        dt_s=1.0,
        drive_stress_Pa=drive,
        rho_initial_m2=0.0,
        rho_increment_per_activation_m2=rho_increment,
        backstress_prefactor_Pa_sqrt_m2=prefactor,
        # A finite zero-stress thermal rate used to destroy the bracket through
        # floating-point leakage just below the exact blocking state.
        rate_function=lambda sigma: 1.0,
    )
    assert value == pytest.approx(block, rel=1.0e-8)


def test_front_width_floor_and_multiplicity_are_independent_of_mpz_dx():
    coarse = EmergentGNDState(
        candidate(),
        physics(
            n_bins=8,
            mpz_length_m=8.0e-6,
            minimum_front_width_m=50.0e-9,
        ),
    )
    fine = EmergentGNDState(
        candidate(),
        physics(
            n_bins=64,
            mpz_length_m=8.0e-6,
            minimum_front_width_m=50.0e-9,
        ),
    )
    for state in (coarse, fine):
        state.mobile_m2[...] = 1.0e20
    coarse_geometry = coarse.source_geometry()
    fine_geometry = fine.source_geometry()
    assert coarse.dx != fine.dx
    assert coarse_geometry["front_width_m"] == pytest.approx(50.0e-9)
    assert fine_geometry["front_width_m"] == pytest.approx(50.0e-9)
    assert coarse_geometry["multiplicity_per_system"] == pytest.approx(
        fine_geometry["multiplicity_per_system"]
    )
    assert coarse_geometry["front_width_at_minimum"] == 1.0
    assert fine_geometry["front_width_at_minimum"] == 1.0


def test_persistent_emission_conserves_line_content_and_does_not_deplete():
    state = EmergentGNDState(candidate(), physics())
    rates = state.local_rates(40.0, 900.0)
    before = state.source_available_fraction()
    state._emit_exact(rates, 1.0e-6)
    after = state.source_available_fraction()
    deposited = float(np.sum(state.mobile_m2) * state.cell_area_m2)
    assert before == 1.0
    assert after == 1.0
    assert np.sum(state.last_source_activations) > 0.0
    assert deposited == pytest.approx(float(np.sum(state.last_line_content)))
    assert np.all(state.source_available_m2 == candidate().rho_source0_m2)


def test_crack_advance_resharpens_accumulated_slip_without_source_refresh():
    state = EmergentGNDState(candidate(), physics())
    state.accumulated_slip_m2[0, 1, 0] = 1.0e22
    before = state.tip_radius_m()
    source_before = state.source_available_m2.copy()
    state.translate_tip(state.dx)
    after = state.tip_radius_m()
    assert after < before
    assert np.array_equal(state.source_available_m2, source_before)
    assert state.source_available_fraction() == 1.0


def test_coupled_moving_tip_is_invariant_to_reporting_interval_split():
    common = physics(
        coupled_moving_tip_enabled=True,
        moving_tip_cfl=0.5,
    )
    whole = EmergentGNDState(candidate(), common)
    split = EmergentGNDState(candidate(), common)
    whole.advance_coupled_segment(
        duration_s=2.0e-9,
        da_m=1.0e-6,
        K_start_MPa_sqrt_m=20.0,
        K_end_MPa_sqrt_m=40.0,
        T_K=900.0,
    )
    split.advance_coupled_segment(
        duration_s=1.0e-9,
        da_m=0.5e-6,
        K_start_MPa_sqrt_m=20.0,
        K_end_MPa_sqrt_m=30.0,
        T_K=900.0,
    )
    split.advance_coupled_segment(
        duration_s=1.0e-9,
        da_m=0.5e-6,
        K_start_MPa_sqrt_m=30.0,
        K_end_MPa_sqrt_m=40.0,
        T_K=900.0,
    )
    assert whole.extension_m == pytest.approx(split.extension_m)
    assert whole.time_s == pytest.approx(split.time_s)
    assert np.allclose(whole.mobile_m2, split.mobile_m2)
    assert np.allclose(whole.retained_m2, split.retained_m2)
    assert np.allclose(whole.accumulated_slip_m2, split.accumulated_slip_m2)
    assert np.allclose(
        whole.cumulative_source_activations,
        split.cumulative_source_activations,
    )


def test_encounter_efficiency_scales_geometric_storage_rate():
    unit = EmergentGNDState(candidate(), physics(encounter_efficiency=1.0))
    scaled = EmergentGNDState(candidate(), physics(encounter_efficiency=9.0))
    unit_rate = unit.local_rates(20.0, 900.0)["encounter_s"]
    scaled_rate = scaled.local_rates(20.0, 900.0)["encounter_s"]
    assert np.allclose(scaled_rate, 9.0 * unit_rate)


def test_taylor_local_stress_uses_shared_2d_phi_limit():
    uncapped = EmergentGNDState(
        candidate(),
        physics(taylor_phi_max=float("inf")),
    )
    capped = EmergentGNDState(candidate(), physics(taylor_phi_max=20.0))
    uncapped_rate = uncapped.local_rates(
        1.0,
        900.0,
    )["taylor_completion_s"]
    capped_rate = capped.local_rates(
        1.0,
        900.0,
    )["taylor_completion_s"]
    assert np.all(capped_rate <= uncapped_rate)
    assert np.any(capped_rate < uncapped_rate)


def test_zero_spatial_transport_preserves_peierls_encounter_storage():
    state = EmergentGNDState(
        candidate(),
        physics(
            encounter_efficiency=9.0,
            mobile_transport_velocity_scale=0.0,
        ),
    )
    rates = state.local_rates(20.0, 900.0)
    assert np.all(rates["velocity_m_s"] == 0.0)
    assert np.any(np.abs(rates["peierls_velocity_m_s"]) > 0.0)
    assert np.any(rates["encounter_s"] > 0.0)


def test_extension_resolved_emission_geometry_replaces_constant_projection():
    state = EmergentGNDState(
        candidate(),
        physics(
            emission_schmid_factors=(1.0, 1.0),
            emission_geometry_extension_m=(0.0, 2.0e-6),
            emission_geometry_factors=((0.1, 0.01), (0.5, 0.02)),
        ),
    )
    assert state.emission_drive_factors() == pytest.approx((0.1, 0.01))
    state.translate_tip(2.0e-6)
    assert state.emission_drive_factors() == pytest.approx((0.5, 0.02))
    rates = state.local_rates(20.0, 900.0)
    expected = 0.5 / 0.02
    ratio = (
        rates["emission_drive_Pa"][0]
        / rates["emission_drive_Pa"][1]
    )
    assert ratio == pytest.approx(expected)


def test_empty_emission_geometry_preserves_constant_projection():
    state = EmergentGNDState(
        candidate(),
        physics(emission_schmid_factors=(0.4, 0.2)),
    )
    assert state.emission_drive_factors() == pytest.approx((0.4, 0.2))


def test_recovery_is_forced_off_and_metadata_declares_persistent_contract():
    row = {
        "candidate_id": "row",
        "cleave_G00_eV": "1", "cleave_gT_eV_per_K": "0",
        "cleave_sigc0_GPa": "2", "cleave_sT_GPa_per_K": "0",
        "cleave_exp_a": "1", "cleave_exp_n": "1", "cleave_floor_frac": "0.01",
        "emit_G00_eV": "1", "emit_gT_eV_per_K": "0",
        "emit_sigc0_GPa": "2", "emit_sT_GPa_per_K": "0",
        "emit_exp_a": "1", "emit_exp_n": "1", "emit_floor_frac": "0.01",
        "peierls_H0_eV": "0.1", "peierls_activation_entropy_kB": "0",
        "peierls_exp_a": "1", "peierls_exp_n": "1", "peierls_nu0_s": "1e8",
        "taylor_H0_eV": "0.2", "taylor_activation_entropy_kB": "0",
        "taylor_exp_a": "1", "taylor_exp_n": "1", "taylor_nu0_s": "1e8",
        "rho_source0_m2": "1e16", "source_refresh_length_um": "5",
        "taylor_corr_rho_c_m2": "1e14", "taylor_corr_scale": "1",
        "recovery_nu0_s": "1e20", "recovery_H0_eV": "1",
        "recovery_activation_entropy_kB": "2", "c_blunt": "1",
    }
    parsed = candidate_from_registry_row(row)
    assert parsed.recovery_nu0_s == 0.0
    state = EmergentGNDState(parsed, physics())
    assert float(state.local_rates(20.0, 900.0)["recovery_rate_s"]) == 0.0
    metadata = state.integration_metadata()
    assert metadata["finite_source_inventory"] is False
    assert metadata["source_depletion_on_emission"] is False
    assert metadata["source_refresh_on_crack_advance"] is False
    assert metadata["explicit_recovery_active"] is False
    assert metadata["front_width_grid_spacing_coupling_active"] is False


def test_front_width_requires_an_explicit_positive_physical_floor():
    with pytest.raises(ValueError, match="positive physical length"):
        physics(minimum_front_width_m=0.0)
