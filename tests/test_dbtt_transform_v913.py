from __future__ import annotations

import math

import numpy as np
import pytest

from arrhenius_fracture.dbtt_transform_v913 import (
    anchored_cleavage_pivot_row,
    scale_cleavage_stress_axis_row,
    surface_linear_values,
    temperature_scale_candidate_row,
    validate_positive_barrier_domain,
)
from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from arrhenius_fracture.emergent_gnd_types_v912 import KB_EV_PER_K


@pytest.fixture
def candidate_row() -> dict[str, float | str]:
    return {
        "candidate_id": "v912_targeted_local_peak_013476_0083",
        "Tref_K": 481.33,
        "cleave_G00_eV": 3.412309742899366,
        "cleave_gT_eV_per_K": 0.0085569133975694,
        "cleave_sigc0_GPa": 4.178054527802145,
        "cleave_sT_GPa_per_K": 0.004399175167258,
        "cleave_exp_a": 0.4846030,
        "cleave_exp_n": 2.115325,
        "cleave_floor_frac": 0.05438703,
        "emit_G00_eV": 2.216504,
        "emit_gT_eV_per_K": 0.000762811,
        "emit_sigc0_GPa": 5.224678,
        "emit_sT_GPa_per_K": -0.004530284,
        "emit_exp_a": 0.2784131,
        "emit_exp_n": 0.8357097,
        "emit_floor_frac": 0.06870673,
        "peierls_H0_eV": 3.767854,
        "peierls_activation_entropy_kB": -8.446862,
        "peierls_exp_a": 1.143917,
        "peierls_exp_n": 1.764089,
        "peierls_nu0_s": 1.0e12,
        "taylor_H0_eV": 0.2980938,
        "taylor_activation_entropy_kB": 27.12420,
        "taylor_exp_a": 0.1601883,
        "taylor_exp_n": 1.529218,
        "taylor_nu0_s": 1.0e11,
        "rho_source0_m2": 6.348286e15,
        "taylor_corr_rho_c_m2": 9.491109e16,
        "taylor_corr_scale": 1.327004,
        "c_blunt": 1.411283,
    }


def test_temperature_scale_preserves_arrhenius_factors(candidate_row):
    scale = 0.7
    original = candidate_from_registry_row(candidate_row)
    transformed_row = temperature_scale_candidate_row(candidate_row, scale)
    transformed = candidate_from_registry_row(transformed_row)

    assert transformed.cleavage.Tref_K == pytest.approx(
        scale * original.cleavage.Tref_K
    )
    assert transformed.peierls.nu0_s == original.peierls.nu0_s
    assert transformed.taylor.nu0_s == original.taylor.nu0_s

    for temperature in (700.0, 900.0, 1000.0, 1200.0):
        for stress_GPa in (0.0, 1.0, 3.0, 6.0, 10.0):
            stress = stress_GPa * 1.0e9
            for old_surface, new_surface in (
                (original.cleavage, transformed.cleavage),
                (original.emission, transformed.emission),
                (
                    original.peierls.surface(original.emission),
                    transformed.peierls.surface(transformed.emission),
                ),
                (
                    original.taylor.surface(original.emission),
                    transformed.taylor.surface(transformed.emission),
                ),
            ):
                old_exponent = float(
                    old_surface.barrier_eV(stress, temperature)
                ) / (KB_EV_PER_K * temperature)
                new_exponent = float(
                    new_surface.barrier_eV(stress, scale * temperature)
                ) / (KB_EV_PER_K * scale * temperature)
                assert math.exp(-new_exponent) == pytest.approx(
                    math.exp(-old_exponent),
                    rel=2.0e-11,
                    abs=2.0e-11,
                )


def test_anchored_pivot_preserves_complete_anchor_surface(candidate_row):
    anchor = 1000.0
    shelf = 700.0
    pivoted_row = anchored_cleavage_pivot_row(
        candidate_row,
        shelf_temperature_K=shelf,
        anchor_temperature_K=anchor,
        shelf_energy_factor=0.9,
        shelf_stress_factor=0.8,
    )
    original = candidate_from_registry_row(candidate_row)
    pivoted = candidate_from_registry_row(pivoted_row)

    for stress_GPa in np.linspace(0.0, 12.0, 25):
        stress = float(stress_GPa) * 1.0e9
        assert float(pivoted.cleavage.barrier_eV(stress, anchor)) == pytest.approx(
            float(original.cleavage.barrier_eV(stress, anchor)),
            rel=2.0e-13,
            abs=2.0e-13,
        )

    old_energy, old_stress = surface_linear_values(candidate_row, "cleave", shelf)
    new_energy, new_stress = surface_linear_values(pivoted_row, "cleave", shelf)
    assert new_energy == pytest.approx(0.9 * old_energy)
    assert new_stress == pytest.approx(0.8 * old_stress)

    changed = {
        "cleave_G00_eV",
        "cleave_gT_eV_per_K",
        "cleave_sigc0_GPa",
        "cleave_sT_GPa_per_K",
    }
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        if field not in changed:
            assert float(pivoted_row[field]) == pytest.approx(
                float(candidate_row[field])
            )


def test_global_cleavage_stress_scale_changes_only_cleavage_axis(candidate_row):
    scaled = scale_cleavage_stress_axis_row(candidate_row, 0.6)
    assert scaled["cleave_sigc0_GPa"] == pytest.approx(
        0.6 * float(candidate_row["cleave_sigc0_GPa"])
    )
    assert scaled["cleave_sT_GPa_per_K"] == pytest.approx(
        0.6 * float(candidate_row["cleave_sT_GPa_per_K"])
    )
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        if field not in ("cleave_sigc0_GPa", "cleave_sT_GPa_per_K"):
            assert float(scaled[field]) == pytest.approx(
                float(candidate_row[field])
            )


def test_positive_domain_rejects_overaggressive_anchor_pivot(candidate_row):
    pivoted = anchored_cleavage_pivot_row(
        candidate_row,
        shelf_temperature_K=700.0,
        anchor_temperature_K=1000.0,
        shelf_energy_factor=0.2,
        shelf_stress_factor=0.2,
    )
    with pytest.raises(ValueError, match="positive domain"):
        validate_positive_barrier_domain(
            pivoted,
            (300.0, 700.0, 1000.0),
            minimum_zero_stress_energy_eV=0.05,
            minimum_characteristic_stress_GPa=0.1,
        )
