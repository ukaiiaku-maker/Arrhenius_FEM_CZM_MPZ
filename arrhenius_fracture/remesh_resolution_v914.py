"""Resolution acceptance for the v9.14 event-centered remesher.

The remesher must not bisect an existing cohesive edge.  Therefore a strict
maximum-edge test over every element in the tip patch can reject a valid mesh
solely because the crack topology contains a protected long edge.  Acceptance
requires both:

1. the actual tip-local mean edge length ``hbar_tip`` is at or below the target;
2. no remaining oversized candidate has a non-cohesive, non-duplicate longest
   edge that the current conforming refinement operator could still split.

Protected and duplicate-blocked candidates remain explicit diagnostics.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def classify_remaining_candidates(backend, mesh, tip, direction) -> dict[str, Any]:
    candidates = list(backend._candidate_elements(mesh, tip, direction))
    cohesive = backend._cohesive_edges()
    refinable = []
    protected = []
    duplicate = []
    for elem_id in candidates:
        i, j, edge_length = backend._longest_edge(mesh, int(elem_id))
        edge = (min(int(i), int(j)), max(int(i), int(j)))
        row = {
            "element": int(elem_id),
            "edge_i": int(i),
            "edge_j": int(j),
            "edge_length_m": float(edge_length),
        }
        if edge in cohesive:
            protected.append(row)
            continue
        q = 0.5 * (mesh.nodes[int(i)] + mesh.nodes[int(j)])
        tol = max(1.0e-14, 1.0e-8 * float(backend.target_h_m))
        if float(np.min(np.linalg.norm(mesh.nodes - q[None, :], axis=1))) <= tol:
            duplicate.append(row)
            continue
        refinable.append(row)
    return {
        "remaining_oversized_candidate_count": int(len(candidates)),
        "remaining_refinable_noncohesive_count": int(len(refinable)),
        "remaining_protected_cohesive_count": int(len(protected)),
        "remaining_duplicate_midpoint_count": int(len(duplicate)),
        "remaining_refinable_noncohesive": refinable,
        "remaining_protected_cohesive": protected,
        "remaining_duplicate_midpoint": duplicate,
    }


def install_resolution_acceptance(backend_class):
    """Patch ``_refine_forward_patch`` and return the original method."""
    original = backend_class._refine_forward_patch

    def wrapped(self, mesh, displacement, damage, tip, direction):
        result = original(self, mesh, displacement, damage, tip, direction)
        new_mesh, new_boundary, new_damage, new_u, parent, audit = result
        remaining = classify_remaining_candidates(
            self, new_mesh, np.asarray(tip, float), np.asarray(direction, float)
        )
        strict_max_edge = bool(audit.get("patch_target_satisfied", False))
        hbar_tip = float(new_mesh.hbar_tip)
        hbar_target = bool(
            np.isfinite(hbar_tip)
            and hbar_tip <= float(self.target_h_m) * (1.0 + 1.0e-10)
        )
        no_refinable = remaining["remaining_refinable_noncohesive_count"] == 0
        audit.update({
            "strict_all_edges_target_satisfied": strict_max_edge,
            "tip_mean_edge_target_satisfied": hbar_target,
            **remaining,
            "patch_target_satisfied": bool(
                strict_max_edge or (hbar_target and no_refinable)
            ),
            "patch_target_interpretation": (
                "strict_all_edges"
                if strict_max_edge
                else (
                    "tip_mean_resolved_remaining_edges_are_topology_protected"
                    if hbar_target and no_refinable
                    else "unresolved_refinable_noncohesive_edges"
                )
            ),
        })
        return new_mesh, new_boundary, new_damage, new_u, parent, audit

    backend_class._refine_forward_patch = wrapped
    return original


def restore_resolution_acceptance(backend_class, original) -> None:
    backend_class._refine_forward_patch = original


__all__ = [
    "classify_remaining_candidates",
    "install_resolution_acceptance",
    "restore_resolution_acceptance",
]
