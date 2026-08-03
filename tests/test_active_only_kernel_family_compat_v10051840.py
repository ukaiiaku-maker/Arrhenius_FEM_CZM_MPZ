from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.active_only_kernel_family_compat_v10051840 import (
    AUDIT_FILE,
    NORMALIZED_FILE,
    load_active_only_kernel_family_compat,
)
from arrhenius_fracture.signed_kernel_family_v1005141 import (
    FAMILY_SCHEMA,
    SignedShieldingKernelFamilyV1005141,
)


def _payload() -> dict:
    common = {
        "r_eff_over_r0": 1.0,
        "opening_strength_fraction": 0.0,
    }
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
        "active_x_m": [1.0e-6, 2.0e-6],
        "wake_x_m": [],
        "activation_to_line_content_by_system": [2.0, 3.0],
        "states": [
            {
                "state_id": "s0",
                "crack_extension_m": 0.0,
                "active_kernel_I_Pa_sqrt_m_per_signed_line": [
                    [1.0, 2.0],
                    [3.0, 4.0],
                ],
                "active_kernel_II_Pa_sqrt_m_per_signed_line": [
                    [0.1, 0.2],
                    [0.3, 0.4],
                ],
                "wake_kernel_I_Pa_sqrt_m_per_signed_line": [[], []],
                "wake_kernel_II_Pa_sqrt_m_per_signed_line": [[], []],
                **common,
            },
            {
                "state_id": "s1",
                "crack_extension_m": 1.0e-3,
                "active_kernel_I_Pa_sqrt_m_per_signed_line": [
                    [2.0, 3.0],
                    [4.0, 5.0],
                ],
                "active_kernel_II_Pa_sqrt_m_per_signed_line": [
                    [0.2, 0.3],
                    [0.4, 0.5],
                ],
                **common,
            },
        ],
        "interpolation": {
            "method": "inverse_distance",
            "neighbors": 2,
            "power": 2.0,
            "envelope_relative_tolerance": 1.0e-10,
            "extrapolation_allowed": False,
        },
    }


def test_empty_wake_family_loads_without_changing_active_kernel(tmp_path: Path):
    source = tmp_path / "family.json"
    source.write_text(json.dumps(_payload()))
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    out = tmp_path / "audit"

    artifact = load_active_only_kernel_family_compat(source, out)

    assert isinstance(artifact, SignedShieldingKernelFamilyV1005141)
    assert artifact.source_path == str(source.resolve())
    assert artifact.wake_x_m.tolist() == [0.0]
    assert np.array_equal(
        artifact.states[0].active_I,
        np.asarray([[1.0, 2.0], [3.0, 4.0]]),
    )
    assert np.array_equal(artifact.states[0].wake_I, np.zeros((2, 1)))
    assert np.array_equal(artifact.states[1].wake_II, np.zeros((2, 1)))

    snapshot = artifact.snapshot(
        0.0,
        np.asarray([1.0e-6, 2.0e-6]),
        np.asarray([1.0e-6, 2.0e-6]),
    )
    assert np.array_equal(
        snapshot.active_kernel_Pa_sqrt_m_per_signed_line,
        np.asarray([[1.0, 2.0], [3.0, 4.0]]),
    )
    assert np.array_equal(
        snapshot.wake_kernel_Pa_sqrt_m_per_signed_line,
        np.zeros((2, 2)),
    )

    audit = json.loads((out / AUDIT_FILE).read_text())
    assert audit["source_sha256"] == source_sha
    assert audit["compatibility_applied"] is True
    assert audit["active_kernel_modified"] is False
    assert audit["runtime_wake_shielding_enabled"] is False
    assert audit["physics_modified"] is False
    assert (out / NORMALIZED_FILE).is_file()


def test_nonzero_wake_data_fails_closed_when_source_grid_is_empty(tmp_path: Path):
    payload = _payload()
    payload["states"][0][
        "wake_kernel_I_Pa_sqrt_m_per_signed_line"
    ] = [[1.0], [0.0]]
    source = tmp_path / "family.json"
    source.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="nonzero or nonfinite wake data"):
        load_active_only_kernel_family_compat(source, tmp_path / "audit")


def test_empty_wake_requires_explicit_active_only_flags(tmp_path: Path):
    payload = _payload()
    payload["wake_shielding_supported"] = True
    source = tmp_path / "family.json"
    source.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="wake_shielding_supported=false"):
        load_active_only_kernel_family_compat(source, tmp_path / "audit")
