"""Arrhenius cohesive-interface support for sharp-front and kinetic-CZM modes.

Cohesive interfaces remain a mechanical representation of Arrhenius cleavage
progress.  No critical traction, opening, fracture energy, or second failure
criterion is evaluated here.  The v10 kinetic mode may create a trial interface
with ``damage=B``; the legacy modes retain abrupt ``damage=1`` insertion.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from scipy import sparse


@dataclass
class CohesiveElement:
    """Zero-thickness two-node line interface controlled by an Arrhenius clock."""

    plus_nodes: tuple[int, int]
    minus_nodes: tuple[int, int]
    normal: np.ndarray
    tangent: np.ndarray
    length: float
    damage: float = 1.0
    clock: float = 0.0
    status: str = "committed"
    front_id: int = -1
    event_index: int = 0
    barrier_kind: str = "exp_floor"
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.normal = np.asarray(self.normal, dtype=float)
        self.tangent = np.asarray(self.tangent, dtype=float)
        nrm = float(np.linalg.norm(self.normal))
        trm = float(np.linalg.norm(self.tangent))
        if nrm <= 0.0 or trm <= 0.0:
            raise ValueError("cohesive normal/tangent must be non-zero")
        self.normal /= nrm
        self.tangent /= trm
        self.damage = float(np.clip(self.damage, 0.0, 1.0))
        self.clock = float(np.clip(self.clock, 0.0, 1.0))
        self.length = float(max(self.length, 0.0))
        self.status = str(self.status)

    @property
    def nodes4(self) -> tuple[int, int, int, int]:
        return (*self.plus_nodes, *self.minus_nodes)

    def set_clock_damage(self, B: float, coupling: str = "clock_linear") -> None:
        q = float(np.clip(B, 0.0, 1.0))
        if q + 1.0e-14 < self.clock:
            raise ValueError("trial cohesive clock/damage must be monotonic")
        mode = str(coupling).strip().lower()
        if mode == "clock_linear":
            damage = q
        elif mode == "abrupt":
            damage = 1.0 if q >= 1.0 - 1.0e-12 else 0.0
        else:
            raise ValueError("cohesive opening coupling must be abrupt or clock_linear")
        self.clock = q
        self.damage = float(np.clip(damage, self.damage, 1.0))
        self.status = "committed" if self.damage >= 1.0 - 1.0e-12 else "trial"


@dataclass
class CohesiveNetwork:
    """Collection of Arrhenius-controlled cohesive interfaces."""

    elements: List[CohesiveElement] = field(default_factory=list)
    penalty_normal_Pa_per_m: float = 1.0e18
    penalty_tangent_Pa_per_m: float = 1.0e18
    compression_penalty_factor: float = 1.0

    def add(self, elem: CohesiveElement) -> None:
        self.elements.append(elem)

    def active_count(self) -> int:
        return sum(1 for e in self.elements if e.damage < 1.0 - 1e-14)

    def failed_count(self) -> int:
        return sum(1 for e in self.elements if e.damage >= 1.0 - 1e-14)

    def to_rows(self) -> np.ndarray:
        rows = []
        for i, e in enumerate(self.elements):
            rows.append([
                i, e.front_id, e.event_index,
                e.plus_nodes[0], e.plus_nodes[1],
                e.minus_nodes[0], e.minus_nodes[1],
                e.length, e.damage, e.clock,
                e.tangent[0], e.tangent[1],
                e.normal[0], e.normal[1],
            ])
        return np.asarray(rows, dtype=float) if rows else np.zeros((0, 14), dtype=float)


def _jump_operator() -> np.ndarray:
    """Map 8 interface dofs to midpoint displacement jump [ux, uy]."""
    A = np.zeros((2, 8), dtype=float)
    A[0, 0] = 0.5; A[1, 1] = 0.5
    A[0, 2] = 0.5; A[1, 3] = 0.5
    A[0, 4] = -0.5; A[1, 5] = -0.5
    A[0, 6] = -0.5; A[1, 7] = -0.5
    return A


def cohesive_element_diagnostics(
    elem: CohesiveElement,
    u: np.ndarray,
    network: CohesiveNetwork,
) -> dict[str, float | str]:
    A = _jump_operator()
    p0, p1 = elem.plus_nodes
    m0, m1 = elem.minus_nodes
    nodes = np.array([p0, p1, m0, m1], dtype=int)
    edofs = np.empty(8, dtype=int)
    edofs[0::2] = 2 * nodes
    edofs[1::2] = 2 * nodes + 1
    jump = A @ np.asarray(u, dtype=float)[edofs]
    dn = float(elem.normal @ jump)
    dt = float(elem.tangent @ jump)
    intact = max(1.0 - float(elem.damage), 0.0)
    kn = float(network.penalty_normal_Pa_per_m)
    kt = float(network.penalty_tangent_Pa_per_m)
    kn_eff = kn * (
        max(float(network.compression_penalty_factor), 0.0) if dn < 0.0 else intact
    )
    kt_eff = kt * intact
    tn = kn_eff * dn
    tt = kt_eff * dt
    return {
        "normal_jump_m": dn,
        "tangential_jump_m": dt,
        "normal_traction_Pa": tn,
        "tangential_traction_Pa": tt,
        "damage": float(elem.damage),
        "clock": float(elem.clock),
        "status": str(elem.status),
    }


def cohesive_contribution(
    network: Optional[CohesiveNetwork],
    u: np.ndarray,
    ndof: int,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Assemble cohesive tangent and internal force with unilateral contact."""
    if network is None or not network.elements:
        return sparse.csr_matrix((ndof, ndof)), np.zeros(ndof)

    A = _jump_operator()
    rows = []
    cols = []
    vals = []
    R = np.zeros(ndof, dtype=float)

    kn = float(network.penalty_normal_Pa_per_m)
    kt = float(network.penalty_tangent_Pa_per_m)
    kc_factor = float(max(network.compression_penalty_factor, 0.0))

    for elem in network.elements:
        p0, p1 = elem.plus_nodes
        m0, m1 = elem.minus_nodes
        nodes = np.array([p0, p1, m0, m1], dtype=int)
        edofs = np.empty(8, dtype=int)
        edofs[0::2] = 2 * nodes
        edofs[1::2] = 2 * nodes + 1
        ue = u[edofs]
        jump = A @ ue

        n = elem.normal
        t = elem.tangent
        dn = float(n @ jump)
        dtt = float(t @ jump)
        intact = max(1.0 - float(elem.damage), 0.0)
        kn_eff = kn * (kc_factor if dn < 0.0 else intact)
        kt_eff = kt * intact
        Kglob = kn_eff * np.outer(n, n) + kt_eff * np.outer(t, t)

        L = max(float(elem.length), 0.0)
        Ke = A.T @ Kglob @ A * L
        Re = A.T @ (Kglob @ jump) * L
        np.add.at(R, edofs, Re)

        elem.metadata.update({
            "normal_jump_m": dn,
            "tangential_jump_m": dtt,
            "normal_traction_Pa": kn_eff * dn,
            "tangential_traction_Pa": kt_eff * dtt,
            "damage": float(elem.damage),
            "clock": float(elem.clock),
            "status": str(elem.status),
        })

        ii = np.repeat(edofs[:, None], 8, axis=1)
        jj = np.repeat(edofs[None, :], 8, axis=0)
        rows.append(ii.ravel()); cols.append(jj.ravel()); vals.append(Ke.ravel())

    if not rows:
        return sparse.csr_matrix((ndof, ndof)), R
    rr = np.concatenate(rows); cc = np.concatenate(cols); vv = np.concatenate(vals)
    K = sparse.csr_matrix((vv, (rr, cc)), shape=(ndof, ndof))
    return K, R


__all__ = [
    "CohesiveElement",
    "CohesiveNetwork",
    "cohesive_element_diagnostics",
    "cohesive_contribution",
]
