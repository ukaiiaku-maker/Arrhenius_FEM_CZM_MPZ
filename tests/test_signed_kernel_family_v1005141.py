from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
)
from arrhenius_fracture.signed_kernel_family_v1005141 import (
    FAMILY_SCHEMA,
    SignedShieldingKernelFamilyV1005141,
    load_signed_shielding_artifact_v1005141,
)


EXTENSIONS = [0.0, 200e-6, 500e-6, 800e-6]
CONVERSION = [0.9124087591240877, 0.9124087591240877]


def family_payload() -> dict:
    x = ((np.arange(40, dtype=float) + 0.5) * 2.5e-6).tolist()
    states = []
    for index, extension in enumerate(EXTENSIONS):
        base = 10.0 * (index + 1)
        active_I = np.vstack(
            (
                base + np.arange(40, dtype=float),
                -0.5 * base - np.arange(40, dtype=float),
            )
        )
        active_II = 0.25 * active_I
        zeros = np.zeros_like(active_I)
        states.append(
            {
                "state_id": f"E{int(round(extension * 1e6)):03d}",
                "crack_extension_m": extension,
                "opening_strength_fraction": 0.0,
                "r_eff_over_r0": 1.0,
                "active_kernel_I_Pa_sqrt_m_per_signed_line": active_I.tolist(),
                "active_kernel_II_Pa_sqrt_m_per_signed_line": active_II.tolist(),
                "wake_kernel_I_Pa_sqrt_m_per_signed_line": zeros.tolist(),
                "wake_kernel_II_Pa_sqrt_m_per_signed_line": zeros.tolist(),
            }
        )
    return {
        "schema": FAMILY_SCHEMA,
        "candidate_independent": True,
        "counts_are_signed_burgers_lines": True,
        "normalization_is_mechanically_derived": True,
        "active_kernel_mechanically_measured": True,
        "kernel_from_signed_interaction_integral": True,
        "signed_burgers_population_required": True,
        "production_parameterization_allowed": True,
        "constitutive_K_shield_cap": False,
        "constitutive_K_shield_cap_present": False,
        "wake_kernel_forced_zero": True,
        "wake_shielding_supported": False,
        "crack_extension_m_semantics": "cumulative_crack_path_extension_m",
        "v10_2_13_state_semantics": {
            "cumulative_crack_path_extension_used_for_interpolation": True,
            "crack_extension_extrapolation_allowed": False,
            "analytical_r_eff_used_for_interpolation": False,
            "opening_strength_fraction_used_for_interpolation": False,
        },
        "active_x_m": x,
        "wake_x_m": x,
        "activation_to_line_content_by_system": CONVERSION,
        "interpolation": {
            "method": "inverse_distance",
            "neighbors": 4,
            "power": 2.0,
            "envelope_relative_tolerance": 1e-10,
            "extrapolation_allowed": False,
        },
        "states": states,
    }


def write_family(tmp_path: Path) -> Path:
    path = tmp_path / "v10_2_14_active_only_campaign_family.json"
    path.write_text(json.dumps(family_payload()))
    return path


def runtime_grid() -> np.ndarray:
    return (np.arange(80, dtype=float) + 0.5) * (50e-6 / 80.0)


def make_state(family: SignedShieldingKernelFamilyV1005141):
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    return PersistentSiteSignedMPZStateV100514(
        candidate,
        family,
        G_Pa=160e9,
        nu=0.28,
        b_m=2.74e-10,
        r0_m=1e-6,
        blunting_length_m=0.5e-6,
    )


def test_loads_actual_pf_family_schema_and_conversion(tmp_path):
    path = write_family(tmp_path)
    artifact = load_signed_shielding_artifact_v1005141(path)
    assert isinstance(artifact, SignedShieldingKernelFamilyV1005141)
    assert artifact.schema == FAMILY_SCHEMA
    assert [state.crack_extension_m for state in artifact.states] == EXTENSIONS
    assert artifact.activation_to_line_content_by_system == pytest.approx(CONVERSION)
    assert artifact.audit_payload()["wake_kernel_forced_zero"] is True


def test_exact_state_and_spatial_projection(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    x = runtime_grid()
    snapshot = family.snapshot(0.0, x, x)
    assert snapshot.active_kernel_Pa_sqrt_m_per_signed_line.shape == (2, 80)
    assert snapshot.metadata["state_weights"] == {"E000": 1.0}
    # Runtime centers below the first measured station use exact endpoint hold.
    assert snapshot.active_kernel_Pa_sqrt_m_per_signed_line[0, 0] == pytest.approx(10.0)
    assert np.allclose(snapshot.wake_kernel_Pa_sqrt_m_per_signed_line, 0.0)


def test_inverse_distance_extension_interpolation(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    x = runtime_grid()
    snapshot = family.snapshot(100e-6, x, x)
    weights = snapshot.metadata["state_weights"]
    distances = np.abs(np.asarray(EXTENSIONS) - 100e-6)
    raw = 1.0 / (distances / (800e-6)) ** 2
    expected = raw / raw.sum()
    for state_id, value in zip(("E000", "E200", "E500", "E800"), expected):
        assert weights[state_id] == pytest.approx(value)
    expected_first = sum(value * base for value, base in zip(expected, (10, 20, 30, 40)))
    assert snapshot.active_kernel_Pa_sqrt_m_per_signed_line[0, 0] == pytest.approx(
        expected_first
    )


def test_family_forbids_crack_extension_extrapolation(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    x = runtime_grid()
    with pytest.raises(ValueError, match="exceeds"):
        family.snapshot(801e-6, x, x)


def test_state_uses_family_conversion_and_dynamic_mode_I_kernel(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    state = make_state(family)
    assert state.activation_to_line_content_by_system == pytest.approx(CONVERSION)
    state.retained_positive[0, 0] = 1.0
    K0 = state.shielding_K()
    state.advance_total_m = 100e-6
    K100 = state.shielding_K()
    assert K0 != pytest.approx(K100)
    audit = state.kernel_artifact_audit()
    assert audit["artifact_kind"] == "crack_extension_kernel_family"
    assert audit["current_cumulative_crack_path_extension_m"] == pytest.approx(100e-6)


def test_advance_preflights_family_envelope_without_mutating_state(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    state = make_state(family)
    state.mobile_positive[0, 0] = 3.0
    before = state.mobile_positive.copy()
    with pytest.raises(ValueError, match="exceeds"):
        state.advance(801e-6)
    assert state.advance_total_m == 0.0
    assert np.array_equal(state.mobile_positive, before)


def test_split_inherits_parent_kernel_coordinate(tmp_path):
    family = SignedShieldingKernelFamilyV1005141.from_json(write_family(tmp_path))
    state = make_state(family)
    state.advance_total_m = 100e-6
    child = state.split(0.4)
    assert child.advance_total_m == pytest.approx(100e-6)
    assert child.kernel_artifact_audit()[
        "current_cumulative_crack_path_extension_m"
    ] == pytest.approx(100e-6)


def test_family_restart_round_trip(tmp_path):
    path = write_family(tmp_path)
    family = SignedShieldingKernelFamilyV1005141.from_json(path)
    original = make_state(family)
    original.mobile_positive[0, 0] = 2.0
    original.retained_negative[1, 3] = 5.0
    original.advance_total_m = 100e-6
    payload = original.state_dict()
    restored = PersistentSiteSignedMPZStateV100514.from_state_dict(payload)
    assert restored.kernel_family is not None
    assert restored.advance_total_m == pytest.approx(100e-6)
    assert np.array_equal(restored.mobile_positive, original.mobile_positive)
    assert np.array_equal(restored.retained_negative, original.retained_negative)
    assert restored.kernel_artifact_audit()["current_state_weights"] == pytest.approx(
        original.kernel_artifact_audit()["current_state_weights"]
    )


def test_rejects_constitutive_shielding_cap(tmp_path):
    payload = family_payload()
    payload["constitutive_K_shield_cap"] = True
    path = tmp_path / "capped.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="no constitutive shielding cap"):
        SignedShieldingKernelFamilyV1005141.from_json(path)
