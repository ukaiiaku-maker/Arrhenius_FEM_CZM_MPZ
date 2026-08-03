"""Diagnostics-only physical-refinement metadata lifecycle repair.

The inherited v10.0.5.13 wrapper verifies the final refinement radius after the
2-D solver returns. The protected solver can expose either no runtime mesh or a
later runtime mesh whose refinement annotations were not propagated, even though
an earlier construction/rebuild mesh was correctly annotated. This context
captures the latest verified annotated mesh and uses it only for the post-run
metadata query when the runtime object cannot verify the requested radius.

The v10.0.5.13.5 long-corridor wrapper substitutes its robust physical mesh
constructor after this lifecycle context is entered. That constructor is also
captured here so the post-run audit cannot lose valid refinement metadata merely
because of wrapper-installation order.

Mesh construction, refinement, mechanics, constitutive state, and quality gates
are unchanged.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

from . import mode_i_first_passage_v10_0_5_13_barrier_only as _entry
from . import mode_i_first_passage_v10_0_5_13_5_barrier_only as _long_corridor

AUDIT_MODEL = "physical_refinement_metadata_lifecycle_v10_0_5_16_1"
_LAST_REFINED_MESH: Any = None


def clear_refinement_audit_v100516() -> None:
    global _LAST_REFINED_MESH
    _LAST_REFINED_MESH = None


def _mesh_has_refinement_metadata(mesh: Any) -> bool:
    if mesh is None:
        return False
    radius = getattr(mesh, "production_refinement_radius_m", None)
    return radius is not None


def _capture_if_verified(mesh: Any) -> Any:
    """Retain an annotated mesh for diagnostics and return it unchanged."""
    global _LAST_REFINED_MESH
    if _mesh_has_refinement_metadata(mesh):
        _LAST_REFINED_MESH = mesh
    return mesh


def refinement_audit_payload_v100516() -> dict[str, Any]:
    mesh = _LAST_REFINED_MESH
    return {
        "model_id": AUDIT_MODEL,
        "captured_mesh_available": mesh is not None,
        "captured_mesh_has_refinement_metadata": _mesh_has_refinement_metadata(mesh),
        "production_refinement_radius_m": (
            None
            if mesh is None
            else getattr(mesh, "production_refinement_radius_m", None)
        ),
        "production_refinement_policy": (
            None
            if mesh is None
            else getattr(mesh, "production_refinement_policy", None)
        ),
        "physics_or_mesh_modified": False,
        "post_run_metadata_lifecycle_only": True,
        "long_corridor_constructor_capture_active": True,
    }


@contextmanager
def installed_refinement_audit_v100516() -> Iterator[None]:
    clear_refinement_audit_v100516()
    old_make = _entry.make_physical_refinement_mesh_v100510
    old_make_long_corridor = (
        _long_corridor.make_physical_refinement_mesh_v1005135
    )
    old_annotate = _entry._annotate_mesh
    old_payload = _entry._mesh_payload

    def make_and_capture(*args, **kwargs):
        return _capture_if_verified(old_make(*args, **kwargs))

    def make_long_corridor_and_capture(*args, **kwargs):
        return _capture_if_verified(
            old_make_long_corridor(*args, **kwargs)
        )

    def annotate_and_capture(*args, **kwargs):
        return _capture_if_verified(old_annotate(*args, **kwargs))

    def payload_with_capture(mesh, requested_radius_m):
        global _LAST_REFINED_MESH
        runtime_payload = old_payload(mesh, requested_radius_m)
        runtime_verified = bool(runtime_payload.get("actual_radius_verified", False))

        captured_payload = None
        captured_verified = False
        if _LAST_REFINED_MESH is not None:
            captured_payload = old_payload(_LAST_REFINED_MESH, requested_radius_m)
            captured_verified = bool(
                captured_payload.get("actual_radius_verified", False)
            )

        if runtime_verified:
            selected = runtime_payload
            if _mesh_has_refinement_metadata(mesh):
                _LAST_REFINED_MESH = mesh
            captured_used = False
            fallback_reason = "runtime_mesh_verified"
        elif captured_verified:
            selected = captured_payload
            captured_used = True
            fallback_reason = (
                "runtime_mesh_missing"
                if mesh is None
                else "runtime_mesh_unannotated_or_unverified"
            )
        else:
            selected = runtime_payload if mesh is not None else (
                captured_payload if captured_payload is not None else runtime_payload
            )
            captured_used = mesh is None and captured_payload is not None
            fallback_reason = "no_verified_refinement_metadata_available"

        selected = dict(selected)
        selected["metadata_lifecycle_model"] = AUDIT_MODEL
        selected["runtime_mesh_pointer_missing"] = mesh is None
        selected["runtime_mesh_metadata_verified"] = runtime_verified
        selected["captured_mesh_metadata_verified"] = captured_verified
        selected["captured_mesh_used"] = captured_used
        selected["metadata_fallback_reason"] = fallback_reason
        selected["long_corridor_constructor_capture_active"] = True
        return selected

    _entry.make_physical_refinement_mesh_v100510 = make_and_capture
    _long_corridor.make_physical_refinement_mesh_v1005135 = (
        make_long_corridor_and_capture
    )
    _entry._annotate_mesh = annotate_and_capture
    _entry._mesh_payload = payload_with_capture
    try:
        yield
    finally:
        _entry.make_physical_refinement_mesh_v100510 = old_make
        _long_corridor.make_physical_refinement_mesh_v1005135 = (
            old_make_long_corridor
        )
        _entry._annotate_mesh = old_annotate
        _entry._mesh_payload = old_payload


__all__ = [
    "AUDIT_MODEL",
    "clear_refinement_audit_v100516",
    "installed_refinement_audit_v100516",
    "refinement_audit_payload_v100516",
]
