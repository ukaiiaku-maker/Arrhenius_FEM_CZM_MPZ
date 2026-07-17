"""Transactional trial-cohesive lifecycle for ``kinetic_campaign_czm``.

The adaptive backend remains responsible for exact-angle topology construction.
This module adds an active trial interface whose damage is controlled only by
the Arrhenius cleavage clock.  It also supplies a full rollback object spanning
backend bookkeeping, mesh-facing arrays, cohesive state, and front-local
kinetics.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from .cohesive import CohesiveElement
from .crack_backend import AdaptiveCZMBackend, CrackAdvanceResult

VALID_OPENING_COUPLINGS = ("abrupt", "clock_linear")


def normalize_opening_coupling(value: str) -> str:
    mode = str(value).strip().lower().replace("-", "_")
    aliases = {"linear": "clock_linear", "clock": "clock_linear"}
    mode = aliases.get(mode, mode)
    if mode not in VALID_OPENING_COUPLINGS:
        raise ValueError(
            f"unknown CZM opening coupling {value!r}; expected {VALID_OPENING_COUPLINGS}"
        )
    return mode


@dataclass
class KineticCZMTransactionSnapshot:
    """Complete state required to undo one active trial event."""

    mesh: Any
    boundary: Any
    displacement: np.ndarray
    damage: np.ndarray
    bulk_history: dict[str, Any]
    backend_state: dict[str, Any]
    cohesive_elements: list[CohesiveElement]
    advance_log: list[dict[str, Any]]
    front_state: Any
    front_position: np.ndarray | None = None
    front_path: Any = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def capture(
        cls,
        *,
        backend: AdaptiveCZMBackend,
        mesh: Any,
        boundary: Any,
        displacement: np.ndarray,
        damage: np.ndarray,
        front_engine: Any | None,
        bulk_history: Mapping[str, Any] | None = None,
        front_position: np.ndarray | None = None,
        front_path: Any = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> "KineticCZMTransactionSnapshot":
        if front_engine is None:
            front_state = None
        elif hasattr(front_engine, "snapshot_kinetic_state"):
            front_state = front_engine.snapshot_kinetic_state()
        else:
            front_state = copy.deepcopy(front_engine)
        return cls(
            mesh=copy.deepcopy(mesh),
            boundary=copy.deepcopy(boundary),
            displacement=np.asarray(displacement, dtype=float).copy(),
            damage=np.asarray(damage, dtype=float).copy(),
            bulk_history=copy.deepcopy(dict(bulk_history or {})),
            backend_state=copy.deepcopy(backend._transaction_snapshot()),
            cohesive_elements=copy.deepcopy(backend.cohesive_network.elements),
            advance_log=copy.deepcopy(list(backend.advance_log)),
            front_state=copy.deepcopy(front_state),
            front_position=(
                None
                if front_position is None
                else np.asarray(front_position, dtype=float).copy()
            ),
            front_path=copy.deepcopy(front_path),
            metadata=copy.deepcopy(dict(metadata or {})),
        )

    def restore_backend_and_front(
        self,
        backend: AdaptiveCZMBackend,
        front_engine: Any | None,
    ) -> None:
        backend.cohesive_network.elements = copy.deepcopy(self.cohesive_elements)
        backend.tip_nodes = {
            int(fid): (int(a), int(b), np.asarray(xy, dtype=float).copy())
            for fid, (a, b, xy) in self.backend_state["tip_nodes"].items()
        }
        backend.event_counter = int(self.backend_state["event_counter"])
        backend.advance_log = copy.deepcopy(self.advance_log)
        if front_engine is None or self.front_state is None:
            return
        if hasattr(front_engine, "restore_kinetic_state") and isinstance(
            self.front_state, Mapping
        ):
            front_engine.restore_kinetic_state(copy.deepcopy(self.front_state))
        else:
            restored = copy.deepcopy(self.front_state)
            front_engine.__dict__.clear()
            front_engine.__dict__.update(restored.__dict__)

    def restored_payload(self) -> dict[str, Any]:
        return {
            "mesh": copy.deepcopy(self.mesh),
            "boundary": copy.deepcopy(self.boundary),
            "displacement": self.displacement.copy(),
            "damage": self.damage.copy(),
            "bulk_history": copy.deepcopy(self.bulk_history),
            "front_position": (
                None if self.front_position is None else self.front_position.copy()
            ),
            "front_path": copy.deepcopy(self.front_path),
            "metadata": copy.deepcopy(self.metadata),
        }


@dataclass
class ActiveTrialSegment:
    front_id: int
    event_index: int
    element_indices: tuple[int, ...]
    log_indices: tuple[int, ...]
    direction: np.ndarray
    requested_length_m: float
    coupling: str
    transaction: KineticCZMTransactionSnapshot
    committed: bool = False

    @property
    def progress(self) -> float:
        return float(self.transaction.metadata.get("last_progress", 0.0))


class KineticTrialAdaptiveCZMBackend(AdaptiveCZMBackend):
    """Adaptive CZM backend with one active trial interface per front."""

    name = "kinetic_trial_adaptive_czm"

    def __init__(self, *args: Any, opening_coupling: str = "clock_linear", **kwargs: Any):
        kwargs["event_damage"] = 0.0
        super().__init__(*args, **kwargs)
        self.opening_coupling = normalize_opening_coupling(opening_coupling)
        self.active_trials: dict[int, ActiveTrialSegment] = {}

    def _trial_elements(self, trial: ActiveTrialSegment) -> list[CohesiveElement]:
        return [self.cohesive_network.elements[i] for i in trial.element_indices]

    def begin_trial_segment(
        self,
        *,
        mesh: Any,
        boundary: Any,
        damage: np.ndarray,
        displacement: np.ndarray,
        p0: np.ndarray,
        p1: np.ndarray,
        direction: np.ndarray,
        front_id: int,
        front_engine: Any | None = None,
        bulk_history: Mapping[str, Any] | None = None,
        front_position: np.ndarray | None = None,
        front_path: Any = None,
        **kwargs: Any,
    ) -> CrackAdvanceResult:
        fid = int(front_id)
        if fid in self.active_trials:
            raise RuntimeError(f"front {fid} already owns an active trial interface")
        transaction = KineticCZMTransactionSnapshot.capture(
            backend=self,
            mesh=mesh,
            boundary=boundary,
            displacement=displacement,
            damage=damage,
            front_engine=front_engine,
            bulk_history=bulk_history,
            front_position=front_position,
            front_path=front_path,
            metadata={"last_progress": 0.0},
        )
        n_elem = len(self.cohesive_network.elements)
        n_log = len(self.advance_log)
        result = super().advance(
            mesh=mesh,
            boundary=boundary,
            damage=damage,
            displacement=displacement,
            p0=p0,
            p1=p1,
            direction=direction,
            front_id=fid,
            **kwargs,
        )
        if not result.inserted:
            transaction.restore_backend_and_front(self, front_engine)
            return result

        element_indices = tuple(range(n_elem, len(self.cohesive_network.elements)))
        log_indices = tuple(range(n_log, len(self.advance_log)))
        if not element_indices:
            transaction.restore_backend_and_front(self, front_engine)
            return CrackAdvanceResult(
                mesh, boundary, damage, displacement, 0.0, False,
                reason="trial_topology_inserted_without_cohesive_elements",
            )
        event_index = int(self.cohesive_network.elements[element_indices[0]].event_index)
        trial = ActiveTrialSegment(
            front_id=fid,
            event_index=event_index,
            element_indices=element_indices,
            log_indices=log_indices,
            direction=np.asarray(direction, dtype=float).copy(),
            requested_length_m=float(np.linalg.norm(np.asarray(p1) - np.asarray(p0))),
            coupling=self.opening_coupling,
            transaction=transaction,
        )
        self.active_trials[fid] = trial
        for elem in self._trial_elements(trial):
            elem.damage = 0.0
            elem.clock = 0.0
            elem.status = "trial"
            elem.metadata.update({
                "kinetic_campaign_czm_trial": True,
                "trial_direction_frozen": True,
                "trial_event_index": event_index,
                "opening_coupling": self.opening_coupling,
            })
        for i in log_indices:
            self.advance_log[i].update({
                "damage": 0.0,
                "clock": 0.0,
                "status": "trial",
                "kinetic_campaign_czm_trial": True,
                "trial_direction_frozen": True,
                "physical_event_index": event_index,
            })
        return result

    def update_trial_segment(
        self,
        front_id: int,
        B: float,
        *,
        coupling: str | None = None,
    ) -> list[CohesiveElement]:
        fid = int(front_id)
        if fid not in self.active_trials:
            raise RuntimeError(f"front {fid} has no active trial interface")
        trial = self.active_trials[fid]
        mode = normalize_opening_coupling(coupling or trial.coupling)
        q = float(np.clip(B, 0.0, 1.0))
        previous = float(trial.transaction.metadata.get("last_progress", 0.0))
        if q + 1.0e-14 < previous:
            raise ValueError("trial cohesive progress must be monotonic")
        for elem in self._trial_elements(trial):
            elem.set_clock_damage(q, mode)
            elem.metadata["kinetic_campaign_czm_trial"] = elem.status == "trial"
        for i in trial.log_indices:
            self.advance_log[i].update({
                "damage": float(self.cohesive_network.elements[
                    trial.element_indices[0]
                ].damage),
                "clock": q,
                "status": "committed" if q >= 1.0 - 1.0e-12 else "trial",
            })
        trial.transaction.metadata["last_progress"] = q
        return self._trial_elements(trial)

    def commit_trial_segment(self, front_id: int) -> list[CohesiveElement]:
        fid = int(front_id)
        elements = self.update_trial_segment(fid, 1.0)
        trial = self.active_trials[fid]
        trial.committed = True
        for elem in elements:
            elem.damage = 1.0
            elem.clock = 1.0
            elem.status = "committed"
            elem.metadata.update({
                "kinetic_campaign_czm_trial": False,
                "kinetic_campaign_czm_committed": True,
                "mpz_advance_on_commit_m": 0.0,
            })
        for i in trial.log_indices:
            self.advance_log[i].update({
                "damage": 1.0,
                "clock": 1.0,
                "status": "committed",
                "kinetic_campaign_czm_trial": False,
                "kinetic_campaign_czm_committed": True,
                "mpz_advance_on_commit_m": 0.0,
            })
        del self.active_trials[fid]
        return elements

    def rollback_trial_segment(
        self,
        front_id: int,
        *,
        front_engine: Any | None = None,
    ) -> dict[str, Any]:
        fid = int(front_id)
        trial = self.active_trials.pop(fid, None)
        if trial is None:
            raise RuntimeError(f"front {fid} has no active trial interface")
        trial.transaction.restore_backend_and_front(self, front_engine)
        return trial.transaction.restored_payload()

    def active_trial(self, front_id: int) -> ActiveTrialSegment | None:
        return self.active_trials.get(int(front_id))

    def active_trial_diagnostics(self, front_id: int) -> dict[str, Any]:
        trial = self.active_trial(front_id)
        if trial is None:
            return {}
        elements = self._trial_elements(trial)
        return {
            "front_id": int(trial.front_id),
            "trial_event_id": int(trial.event_index),
            "trial_element_count": len(elements),
            "trial_cohesive_damage": max(float(e.damage) for e in elements),
            "cleavage_clock_B": max(float(e.clock) for e in elements),
            "trial_status": "committed" if trial.committed else "trial",
            "trial_direction": trial.direction.tolist(),
            "trial_direction_frozen": True,
            "requested_length_m": float(trial.requested_length_m),
        }


__all__ = [
    "VALID_OPENING_COUPLINGS",
    "normalize_opening_coupling",
    "KineticCZMTransactionSnapshot",
    "ActiveTrialSegment",
    "KineticTrialAdaptiveCZMBackend",
]
