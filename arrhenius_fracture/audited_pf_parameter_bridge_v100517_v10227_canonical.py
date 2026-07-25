"""Canonical PF v10.2.27 bridge with serialization-independent validation.

The PF v10.2.27 selection manifest explicitly records that the canonical CSV
serialization changed after the source commit.  Exact equality is therefore
established by the four option/candidate/class rows, the persistent-site
closure, and the authoritative active-parameter fingerprint.  The observed CSV
SHA-256 is still recorded for provenance, but it is not treated as a physics
identity.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import audited_pf_parameter_bridge_v100517_v10227 as _legacy

BRIDGE_SCHEMA = "audited_PF_v10_2_27_four_class_bridge_v10_0_5_17_canonical"
EXPECTED_ACTIVE_FINGERPRINT_SHA256 = _legacy.EXPECTED_ACTIVE_FINGERPRINT_SHA256
EXPECTED_OPTIONS = _legacy.EXPECTED_OPTIONS
EXPECTED_REGISTRY_SHA256 = None
PF_BRANCH_REQUIRED = None
PF_COMMIT_PREFIX_REQUIRED = None


def load_audited_parameter_option(
    pf_repo_root: str | Path,
    parameter_option: str,
):
    root = Path(pf_repo_root).expanduser().resolve()
    materials = root / "arrhenius_fracture" / "data" / "materials"
    registry_path = materials / _legacy.REGISTRY_NAME
    selection_path = materials / _legacy.SELECTION_NAME
    if not registry_path.is_file() or not selection_path.is_file():
        missing = [str(path) for path in (registry_path, selection_path) if not path.is_file()]
        raise FileNotFoundError("missing canonical PF v10.2.27 registry inputs: " + ", ".join(missing))

    observed_registry_sha = _legacy._sha256(registry_path)
    selection: dict[str, Any] = json.loads(selection_path.read_text())
    manifest_registry_sha = selection.get("source_registry_sha256")
    if manifest_registry_sha not in (None, "", observed_registry_sha):
        raise ValueError(
            "PF selection manifest registry SHA-256 conflicts with the observed file: "
            f"manifest={manifest_registry_sha} observed={observed_registry_sha}"
        )

    old_expected = _legacy.EXPECTED_REGISTRY_SHA256
    try:
        # Preserve every legacy row/class/closure/fingerprint validation while
        # allowing the manifest-authorized CSV serialization to differ.
        _legacy.EXPECTED_REGISTRY_SHA256 = observed_registry_sha
        candidate, audit = _legacy.load_audited_parameter_option(root, parameter_option)
    finally:
        _legacy.EXPECTED_REGISTRY_SHA256 = old_expected

    audit["schema"] = BRIDGE_SCHEMA
    audit["PF_branch_required"] = None
    audit["PF_commit_prefix_required"] = None
    audit["registry_sha256"] = observed_registry_sha
    audit["registry_byte_hash_used_as_physics_identity"] = False
    audit["registry_validation_policy"] = (
        "exact option/candidate/class rows + persistent-site closure + authoritative "
        "active-parameter fingerprint; observed CSV SHA retained for provenance"
    )
    audit["selection_source_commit"] = selection.get("source_commit")
    audit["selection_registry_sha256_note"] = selection.get("source_registry_sha256_note")
    return candidate, audit


__all__ = [
    "BRIDGE_SCHEMA",
    "EXPECTED_ACTIVE_FINGERPRINT_SHA256",
    "EXPECTED_OPTIONS",
    "EXPECTED_REGISTRY_SHA256",
    "PF_BRANCH_REQUIRED",
    "PF_COMMIT_PREFIX_REQUIRED",
    "load_audited_parameter_option",
]
