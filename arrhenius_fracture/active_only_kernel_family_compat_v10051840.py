"""Audited schema compatibility for PF active-only kernel families.

Some direct prescribed-geometry PF families encode the disabled wake channel as
an empty grid and omit (or store empty) wake matrices.  The inherited FEM/CZM
v10.0.5.14.1 loader predates that representation and requires a nonempty source
wake grid even though it subsequently forces the runtime wake kernel to exactly
zero.

This adapter is permitted only when the source artifact explicitly declares
``wake_kernel_forced_zero=true`` and ``wake_shielding_supported=false``.  It
creates a one-point, exactly-zero *representation* for the disabled source wake
channel, leaves every active-kernel coefficient and interpolation state
unchanged, and retains the original source path and SHA-256 in the returned
artifact and audit.  No runtime wake shielding is enabled.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from . import signed_kernel_family_v1005141 as _kernel

MODEL_ID = "active_only_empty_wake_schema_compatibility_v10_0_5_18_4_0"
AUDIT_FILE = "active_only_kernel_family_compatibility_v10_0_5_18_4_0.json"
NORMALIZED_FILE = "active_only_kernel_family_compatibility_input_v10_0_5_18_4_0.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _array_size(value: Any) -> int:
    try:
        return int(np.asarray(value, dtype=float).size)
    except (TypeError, ValueError):
        return -1


def _empty_or_missing(value: Any) -> bool:
    return value is None or _array_size(value) == 0


def _require_exact_active_only_contract(payload: dict[str, Any]) -> None:
    if payload.get("schema") != _kernel.FAMILY_SCHEMA:
        raise ValueError(
            "empty-wake compatibility supports only schema "
            f"{_kernel.FAMILY_SCHEMA!r}"
        )
    if payload.get("wake_kernel_forced_zero") is not True:
        raise ValueError(
            "empty-wake compatibility requires wake_kernel_forced_zero=true"
        )
    if payload.get("wake_shielding_supported") is not False:
        raise ValueError(
            "empty-wake compatibility requires wake_shielding_supported=false"
        )
    if payload.get("constitutive_K_shield_cap") is not False:
        raise ValueError(
            "empty-wake compatibility requires constitutive_K_shield_cap=false"
        )


def _normalized_payload(payload: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    """Return a loader-compatible copy and whether normalization was needed."""
    wake_x = payload.get("wake_x_m")
    if not _empty_or_missing(wake_x):
        return copy.deepcopy(payload), False

    _require_exact_active_only_contract(payload)
    conversion = np.asarray(
        payload.get("activation_to_line_content_by_system"), dtype=float
    ).reshape(-1)
    if conversion.size == 0 or not np.all(np.isfinite(conversion)):
        raise ValueError(
            "active-only kernel family has no valid activation-to-line conversion"
        )
    n_systems = int(conversion.size)

    raw_states = payload.get("states")
    if not isinstance(raw_states, list) or len(raw_states) < 2:
        raise ValueError("active-only kernel family requires at least two states")

    normalized = copy.deepcopy(payload)
    normalized["wake_x_m"] = [0.0]
    normalized["fem_empty_wake_schema_compatibility"] = {
        "model_id": MODEL_ID,
        "source_wake_grid_was_empty": True,
        "normalized_source_wake_grid_points": 1,
        "normalized_source_wake_value": 0.0,
        "runtime_wake_shielding_enabled": False,
        "active_kernel_modified": False,
        "physics_modified": False,
    }

    wake_keys = (
        "wake_kernel_I_Pa_sqrt_m_per_signed_line",
        "wake_kernel_II_Pa_sqrt_m_per_signed_line",
    )
    for index, state in enumerate(normalized["states"]):
        if not isinstance(state, dict):
            raise ValueError(f"states[{index}] must be a mapping")
        for key in wake_keys:
            raw = state.get(key)
            if not _empty_or_missing(raw):
                array = np.asarray(raw, dtype=float)
                if not np.all(np.isfinite(array)) or not np.allclose(array, 0.0):
                    raise ValueError(
                        f"states[{index}].{key} contains nonzero or nonfinite wake data"
                    )
            state[key] = np.zeros((n_systems, 1), dtype=float).tolist()
    return normalized, True


def load_active_only_kernel_family_compat(
    path: str | Path,
    audit_root: str | Path,
):
    source = Path(path).expanduser().resolve()
    out = Path(audit_root).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    payload = json.loads(source.read_text())
    normalized, applied = _normalized_payload(payload)

    if not applied:
        artifact = _kernel.load_signed_shielding_artifact_v1005141(source)
        audit = {
            "schema": MODEL_ID,
            "source_path": str(source),
            "source_sha256": _sha256(source),
            "compatibility_applied": False,
            "reason": "source already contains a nonempty wake representation",
            "active_kernel_modified": False,
            "runtime_wake_shielding_enabled": False,
            "physics_modified": False,
        }
        (out / AUDIT_FILE).write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n"
        )
        return artifact

    normalized_path = out / NORMALIZED_FILE
    normalized_path.write_text(
        json.dumps(normalized, indent=2, sort_keys=True) + "\n"
    )
    artifact = _kernel.load_signed_shielding_artifact_v1005141(normalized_path)
    if not isinstance(artifact, _kernel.SignedShieldingKernelFamilyV1005141):
        raise TypeError("normalized active-only artifact did not load as a family")

    metadata = copy.deepcopy(artifact.metadata)
    metadata.update(
        {
            "fem_empty_wake_schema_compatibility_model": MODEL_ID,
            "fem_empty_wake_schema_compatibility_applied": True,
            "original_source_path": str(source),
            "original_source_sha256": _sha256(source),
            "normalized_representation_path": str(normalized_path),
            "normalized_representation_sha256": _sha256(normalized_path),
            "active_kernel_modified": False,
            "runtime_wake_shielding_enabled": False,
            "physics_modified": False,
        }
    )
    artifact = replace(
        artifact,
        metadata=metadata,
        source_path=str(source),
    )

    audit = {
        "schema": MODEL_ID,
        "source_path": str(source),
        "source_sha256": _sha256(source),
        "normalized_representation_path": str(normalized_path),
        "normalized_representation_sha256": _sha256(normalized_path),
        "compatibility_applied": True,
        "source_wake_grid_was_empty": True,
        "normalized_source_wake_grid_points": 1,
        "normalized_source_wake_value": 0.0,
        "active_kernel_modified": False,
        "interpolation_states_modified": False,
        "activation_normalization_modified": False,
        "runtime_wake_shielding_enabled": False,
        "physics_modified": False,
    }
    (out / AUDIT_FILE).write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n"
    )
    return artifact


@contextmanager
def installed_active_only_kernel_family_compat(
    audit_root: str | Path,
) -> Iterator[None]:
    """Patch only the v10.0.5.14.1 local loader symbol for one run."""
    from . import mode_i_first_passage_v10_0_5_14_1_persistent_site_family as entry

    original = entry.load_signed_shielding_artifact_v1005141

    def compatible_loader(path: str | Path):
        return load_active_only_kernel_family_compat(path, audit_root)

    entry.load_signed_shielding_artifact_v1005141 = compatible_loader
    try:
        yield
    finally:
        entry.load_signed_shielding_artifact_v1005141 = original


__all__ = [
    "AUDIT_FILE",
    "MODEL_ID",
    "NORMALIZED_FILE",
    "installed_active_only_kernel_family_compat",
    "load_active_only_kernel_family_compat",
]
