"""Local crack-tip mesh patch planning utilities.

The initial production backend uses topology-only edge splitting so element
history remains exactly attached to unchanged bulk elements.  This module
provides the patch-selection boundary for the next backend: local retriangulation
can be added here without touching FrontEngine, fatigue cycle jumping, or branch
bookkeeping.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class PatchPlan:
    element_ids: np.ndarray
    node_ids: np.ndarray
    center: np.ndarray
    radius: float


def select_tip_patch(mesh, tip_xy, radius):
    tip = np.asarray(tip_xy, float)
    cent = mesh.nodes[mesh.elems].mean(axis=1)
    elem_ids = np.where(np.linalg.norm(cent - tip[None, :], axis=1) <= float(radius))[0]
    node_ids = np.unique(mesh.elems[elem_ids].ravel()) if len(elem_ids) else np.array([], dtype=int)
    return PatchPlan(elem_ids, node_ids, tip.copy(), float(radius))
