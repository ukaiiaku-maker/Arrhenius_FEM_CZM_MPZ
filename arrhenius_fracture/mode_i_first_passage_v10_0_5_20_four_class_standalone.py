"""Standalone four-class FEM/CZM with persistent geometry-failure diagnostics."""
from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
import sys

from . import mode_i_first_passage_v10_0_5_18_four_class_standalone as _base18
from . import mode_i_first_passage_v10_0_5_19_four_class_standalone as _base19
from .adaptive_czm_quality_partition_v100519 import (
    installed_adaptive_quality_partition_v100519,
)
from .adaptive_czm_partition_failure_audit_v100520 import (
    MODEL_ID as FAILURE_AUDIT_MODEL,
    installed_partition_failure_audit_v100520,
)

POINT_RELEASE = "10.0.5.20-frozen-four-class-partition-failure-audit"
MODEL_ID = "FEM_CZM_full_2D_frozen_four_class_partition_failure_audit_v10_0_5_20"
PRODUCTION_MANIFEST = _base19.PRODUCTION_MANIFEST
SELECTION_MANIFEST = _base19.SELECTION_MANIFEST
TRANSFER_MANIFEST = _base19.TRANSFER_MANIFEST


def _out_path(argv: list[str]) -> Path | None:
    if "--out" not in argv:
        return None
    index = argv.index("--out")
    if index + 1 >= len(argv):
        return None
    return Path(argv[index + 1]).expanduser().resolve()


@contextmanager
def _installed_v100520_geometry_path():
    with installed_adaptive_quality_partition_v100519():
        with installed_partition_failure_audit_v100520():
            yield


def _rewrite_metadata(out: Path | None) -> None:
    _base19._rewrite_metadata(out)
    if out is None:
        return

    manifest = out / PRODUCTION_MANIFEST
    if manifest.is_file():
        payload = json.loads(manifest.read_text())
        payload["schema"] = "persistent_site_production_manifest_v10_0_5_20_frozen_four_class"
        payload["model"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        physics = dict(payload.get("physics_contract", {}) or {})
        physics.update(
            {
                "adaptive_czm_partition_failure_audit_model": FAILURE_AUDIT_MODEL,
                "raw_backend_rejection_reasons_persisted_after_rollback": True,
                "geometry_behavior_changed_by_failure_audit": False,
                "quality_floors_changed_by_failure_audit": False,
                "constitutive_physics_changed_by_failure_audit": False,
            }
        )
        payload["physics_contract"] = physics
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")

    selection = out / SELECTION_MANIFEST
    if selection.is_file():
        payload = json.loads(selection.read_text())
        payload["schema"] = MODEL_ID
        payload["point_release"] = POINT_RELEASE
        policy = dict(payload.get("policy", {}) or {})
        policy.update(
            {
                "adaptive_czm_partition_failure_audit_model": FAILURE_AUDIT_MODEL,
                "raw_backend_rejection_reasons_persisted_after_rollback": True,
            }
        )
        payload["policy"] = policy
        selection.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out = _out_path(user_args)
    try:
        with _installed_v100520_geometry_path():
            return _base18.main(user_args)
    finally:
        _rewrite_metadata(out)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
