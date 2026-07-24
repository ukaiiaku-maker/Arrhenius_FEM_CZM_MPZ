"""Diagnostics-only physical-refinement metadata lifecycle repair.

The inherited v10.0.5.13 wrapper verifies the final refinement radius after the
2-D solver returns.  The protected solver may clear its runtime mesh reference
before that post-run query, causing a completed calculation to be marked failed.
This context captures the already annotated construction/rebuild mesh and uses it
only when the post-run runtime pointer is absent.  Mesh construction, refinement,
mechanics, and quality gates are unchanged.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

from . import mode_i_first_passage_v10_0_5_13_barrier_only as _entry

AUDIT_MODEL = "physical_refinement_metadata_lifecycle_v10_0_5_16"
_LAST_REFINED_MESH: Any = None


def clear_refinement_audit_v100516() -> None:
    global _LAST_REFINED_MESH
    _LAST_REFINED_MESH = None


def refinement_audit_payload_v100516() -> dict[str, Any]:
    mesh = _LAST_REFINED_MESH
    return {
        "model_id": AUDIT_MODEL,
        "captured_mesh_available": mesh is not None,
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
    }


@contextmanager
def installed_refinement_audit_v100516() -> Iterator[None]:
    global _LAST_REFINED_MESH
    clear_refinement_audit_v100516()
    old_make = _entry.make_physical_refinement_mesh_v100510
    old_annotate = _entry._annotate_mesh
    old_payload = _entry._mesh_payload

    def make_and_capture(*args, **kwargs):
        global _LAST_REFINED_MESH
        mesh = old_make(*args, **kwargs)
        _LAST_REFINED_MESH = mesh
        return mesh

    def annotate_and_capture(*args, **kwargs):
        global _LAST_REFINED_MESH
        mesh = old_annotate(*args, **kwargs)
        _LAST_REFINED_MESH = mesh
        return mesh

    def payload_with_capture(mesh, requested_radius_m):
        global _LAST_REFINED_MESH
        selected = mesh if mesh is not None else _LAST_REFINED_MESH
        if selected is not None:
            _LAST_REFINED_MESH = selected
        payload = old_payload(selected, requested_radius_m)
        payload["metadata_lifecycle_model"] = AUDIT_MODEL
        payload["runtime_mesh_pointer_missing"] = mesh is None
        payload["captured_mesh_used"] = mesh is None and selected is not None
        return payload

    _entry.make_physical_refinement_mesh_v100510 = make_and_capture
    _entry._annotate_mesh = annotate_and_capture
    _entry._mesh_payload = payload_with_capture
    try:
        yield
    finally:
        _entry.make_physical_refinement_mesh_v100510 = old_make
        _entry._annotate_mesh = old_annotate
        _entry._mesh_payload = old_payload


__all__ = [
    "AUDIT_MODEL",
    "clear_refinement_audit_v100516",
    "installed_refinement_audit_v100516",
    "refinement_audit_payload_v100516",
]
