"""Event-centered conservative remeshing for the v9.14 FEM/CZM solver.

The existing :class:`AdaptiveCZMBackend` realizes one Arrhenius renewal as one
physical crack increment. This backend retains that event law, refines a forward
patch around the *new* crack tip, transfers every piecewise-constant Gauss-point
field through an exact parent map, and immediately re-equilibrates the FEM at the
unchanged event displacement and physical time.

The remesh is refinement-only. Every new triangle has one old parent, every new
node is introduced at an existing non-cohesive edge midpoint, and pre-existing
cohesive node numbers/history objects are untouched. Children inherit their
parent state and their areas sum to the parent area, giving an exact conservative
operator for the solver's element-centered state representation.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .crack_backend import AdaptiveCZMBackend, CrackAdvanceResult
from .event_equilibrium_v914 import ACTIVE_CONTEXT
from .mesh import make_boundary_data, rebuild_tri_mesh


class EventRemeshCZMBackend(AdaptiveCZMBackend):
    """One physical CZM event, conservative tip remesh, same-load equilibrium."""

    name = "event_remesh_czm"

    def __init__(
        self,
        *,
        geom,
        target_h_m: float,
        patch_radius_m: float,
        max_edge_splits_per_event: int = 256,
        target_edge_factor: float = 1.25,
        forward_back_margin_m: float = 0.0,
        min_remesh_triangle_quality: float = 0.02,
        require_post_event_equilibrium: bool = True,
        fail_fast_on_event_error: bool = True,
        **kwargs: Any,
    ) -> None:
        super().__init__(geom=geom, **kwargs)
        self.target_h_m = max(float(target_h_m), 1.0e-12)
        self.patch_radius_m = max(float(patch_radius_m), 2.0 * self.target_h_m)
        self.max_edge_splits_per_event = max(int(max_edge_splits_per_event), 0)
        self.target_edge_factor = max(float(target_edge_factor), 1.0)
        self.forward_back_margin_m = max(float(forward_back_margin_m), 0.0)
        self.min_remesh_triangle_quality = max(
            float(min_remesh_triangle_quality), 1.0e-8
        )
        self.require_post_event_equilibrium = bool(require_post_event_equilibrium)
        self.fail_fast_on_event_error = bool(fail_fast_on_event_error)
        self.remesh_audit: list[dict[str, Any]] = []
        self.failed_event_attempts: list[dict[str, Any]] = []
        self.physical_event_counter = 0

    @staticmethod
    def _max_edge_lengths(mesh) -> np.ndarray:
        x = mesh.nodes[mesh.elems]
        return np.maximum.reduce([
            np.linalg.norm(x[:, 1] - x[:, 0], axis=1),
            np.linalg.norm(x[:, 2] - x[:, 1], axis=1),
            np.linalg.norm(x[:, 0] - x[:, 2], axis=1),
        ])

    @staticmethod
    def _longest_edge(mesh, elem_id: int) -> tuple[int, int, float]:
        conn = [int(v) for v in mesh.elems[int(elem_id)]]
        candidates = []
        for a, b in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
            length = float(np.linalg.norm(mesh.nodes[b] - mesh.nodes[a]))
            candidates.append((length, min(a, b), max(a, b)))
        length, i, j = max(candidates, key=lambda z: z[0])
        return int(i), int(j), float(length)

    def _cohesive_edges(self) -> set[tuple[int, int]]:
        edges: set[tuple[int, int]] = set()
        if self.cohesive_network is None:
            return edges
        for ce in self.cohesive_network.elements:
            for a, b in (ce.plus_nodes, ce.minus_nodes):
                edges.add((min(int(a), int(b)), max(int(a), int(b))))
        return edges

    @staticmethod
    def _cohesive_signature(network) -> list[tuple[Any, ...]]:
        if network is None:
            return []
        out = []
        for ce in network.elements:
            out.append((
                tuple(int(x) for x in ce.plus_nodes),
                tuple(int(x) for x in ce.minus_nodes),
                float(ce.damage),
                float(ce.clock),
                int(ce.front_id),
                int(ce.event_index),
                str(ce.barrier_kind),
                json.dumps(ce.metadata, sort_keys=True, default=str),
            ))
        return out

    def _patch_mask(self, mesh, tip: np.ndarray, direction: np.ndarray) -> np.ndarray:
        cent = mesh.nodes[mesh.elems].mean(axis=1)
        rel = cent - tip[None, :]
        tangent = direction / max(float(np.linalg.norm(direction)), 1.0e-300)
        normal = np.array([-tangent[1], tangent[0]])
        xi = rel @ tangent
        eta = rel @ normal
        radius = self.patch_radius_m
        return (
            (xi >= -self.forward_back_margin_m)
            & (xi <= radius)
            & (np.abs(eta) <= radius)
            & (xi * xi + eta * eta <= radius * radius)
        )

    def _candidate_elements(self, mesh, tip: np.ndarray, direction: np.ndarray) -> list[int]:
        mask = self._patch_mask(mesh, tip, direction)
        lengths = self._max_edge_lengths(mesh)
        ids = np.where(mask & (lengths > self.target_edge_factor * self.target_h_m))[0]
        if ids.size == 0:
            return []
        cent = mesh.nodes[mesh.elems].mean(axis=1)
        distance = np.linalg.norm(cent[ids] - tip[None, :], axis=1)
        order = np.lexsort((distance, -lengths[ids]))
        return ids[order].astype(int).tolist()

    def _refine_forward_patch(
        self,
        mesh,
        displacement: np.ndarray,
        damage: np.ndarray,
        tip: np.ndarray,
        direction: np.ndarray,
    ):
        current_mesh = mesh
        current_u = np.asarray(displacement, float).copy()
        current_d = np.asarray(damage, float).copy()
        cumulative = np.arange(mesh.ne, dtype=int)
        nsplit = 0
        rejected_quality = 0
        rejected_cohesive = 0
        rejected_duplicate = 0
        split_rows: list[dict[str, Any]] = []

        while nsplit < self.max_edge_splits_per_event:
            candidates = self._candidate_elements(current_mesh, tip, direction)
            if not candidates:
                break
            cohesive_edges = self._cohesive_edges()
            accepted = False
            for e in candidates:
                i, j, edge_length = self._longest_edge(current_mesh, e)
                edge_key = (min(i, j), max(i, j))
                if edge_key in cohesive_edges:
                    rejected_cohesive += 1
                    continue
                q = 0.5 * (current_mesh.nodes[i] + current_mesh.nodes[j])
                tol = max(1.0e-14, 1.0e-8 * self.target_h_m)
                if float(np.min(np.linalg.norm(current_mesh.nodes - q[None, :], axis=1))) <= tol:
                    rejected_duplicate += 1
                    continue

                trial_mesh, trial_u, _reason, meta, parent_map = self._insert_point_on_edge(
                    current_mesh, current_u, q, i, j
                )
                if trial_mesh is None or parent_map is None:
                    continue
                qmin = float(meta.get("min_triangle_quality", np.nan))
                if np.isfinite(qmin) and qmin < self.min_remesh_triangle_quality:
                    rejected_quality += 1
                    continue
                if trial_mesh.nn != current_mesh.nn + 1:
                    raise RuntimeError("event remesh expected one midpoint node per edge split")

                dmid = max(float(current_d[i]), float(current_d[j]))
                trial_d = np.concatenate([current_d, np.array([dmid], dtype=float)])
                parent_map = np.asarray(parent_map, dtype=int)
                cumulative = cumulative[parent_map]
                split_rows.append({
                    "split_index": int(nsplit),
                    "parent_element": int(e),
                    "edge_i": int(i),
                    "edge_j": int(j),
                    "edge_length_before_m": edge_length,
                    "midpoint_x_m": float(q[0]),
                    "midpoint_y_m": float(q[1]),
                    "min_triangle_quality": qmin,
                    "n_elements_after": int(trial_mesh.ne),
                })
                current_mesh = trial_mesh
                current_u = np.asarray(trial_u, float)
                current_d = trial_d
                nsplit += 1
                accepted = True
                break
            if not accepted:
                break

        current_mesh = rebuild_tri_mesh(
            current_mesh.nodes, current_mesh.elems, tip_centers=[tip]
        )
        current_bnd = make_boundary_data(current_mesh, self.geom)
        post_area = np.asarray(mesh.area_e, float)
        final_area = np.asarray(current_mesh.area_e, float)
        inherited_area = np.bincount(
            cumulative, weights=final_area, minlength=mesh.ne
        )[: mesh.ne]
        area_error = inherited_area - post_area
        relative_area_error = np.abs(area_error) / np.maximum(post_area, 1.0e-300)
        patch_ids = np.where(self._patch_mask(current_mesh, tip, direction))[0]
        max_patch_edge = (
            float(np.max(self._max_edge_lengths(current_mesh)[patch_ids]))
            if patch_ids.size else 0.0
        )
        audit = {
            "n_edge_splits": int(nsplit),
            "n_elements_before_patch": int(mesh.ne),
            "n_elements_after_patch": int(current_mesh.ne),
            "n_nodes_before_patch": int(mesh.nn),
            "n_nodes_after_patch": int(current_mesh.nn),
            "target_h_m": float(self.target_h_m),
            "target_edge_factor": float(self.target_edge_factor),
            "patch_radius_m": float(self.patch_radius_m),
            "hbar_tip_before_m": float(mesh.hbar_tip),
            "hbar_tip_after_m": float(current_mesh.hbar_tip),
            "max_edge_in_patch_after_m": max_patch_edge,
            "patch_target_satisfied": bool(
                max_patch_edge <= self.target_edge_factor * self.target_h_m * (1.0 + 1.0e-10)
            ),
            "split_budget_exhausted": bool(
                nsplit >= self.max_edge_splits_per_event
                and bool(self._candidate_elements(current_mesh, tip, direction))
            ),
            "total_area_before_m2": float(np.sum(post_area)),
            "total_area_after_m2": float(np.sum(final_area)),
            "relative_total_area_error": float(
                abs(np.sum(final_area) - np.sum(post_area))
                / max(abs(np.sum(post_area)), 1.0e-300)
            ),
            "max_parent_area_conservation_error": (
                float(np.max(np.abs(area_error))) if area_error.size else 0.0
            ),
            "max_parent_relative_area_conservation_error": (
                float(np.max(relative_area_error)) if relative_area_error.size else 0.0
            ),
            "rejected_quality": int(rejected_quality),
            "rejected_cohesive_edge": int(rejected_cohesive),
            "rejected_duplicate_midpoint": int(rejected_duplicate),
            "split_rows": split_rows,
        }
        return current_mesh, current_bnd, current_d, current_u, cumulative, audit

    def _rollback_failed_event(
        self,
        transaction,
        kwargs,
        reason: str,
        details: dict[str, Any],
    ) -> CrackAdvanceResult:
        self._transaction_rollback(transaction)
        record = {"reason": str(reason), **details}
        self.failed_event_attempts.append(record)
        if self.fail_fast_on_event_error:
            raise RuntimeError(
                "v9.14 event remesh transaction rolled back and the case was "
                f"aborted: {reason}"
            )
        return CrackAdvanceResult(
            kwargs["mesh"], kwargs["boundary"], kwargs["damage"],
            kwargs["displacement"], 0.0, False,
            angle_error_deg=180.0, reason=str(reason), elem_parent_map=None,
        )

    def advance(self, **kwargs) -> CrackAdvanceResult:
        pre_mesh = kwargs["mesh"]
        transaction = self._transaction_snapshot()
        cohesive_before = self._cohesive_signature(self.cohesive_network)
        try:
            result = super().advance(**kwargs)
            if not result.inserted or result.moved <= 0.0:
                return result

            new_log_rows = self.advance_log[int(transaction["n_log"]):]
            if not new_log_rows:
                raise RuntimeError("physical event inserted no advance-log segment")
            tip = np.array([
                float(new_log_rows[-1]["x1"]), float(new_log_rows[-1]["y1"])
            ], dtype=float)
            start = np.array([
                float(new_log_rows[0]["x0"]), float(new_log_rows[0]["y0"])
            ], dtype=float)
            direction = tip - start
            direction /= max(float(np.linalg.norm(direction)), 1.0e-300)

            patch_mesh, patch_bnd, patch_d, patch_u, patch_map, patch_audit = \
                self._refine_forward_patch(
                    result.mesh, result.displacement, result.damage, tip, direction
                )
            event_map = (
                np.arange(result.mesh.ne, dtype=int)
                if result.elem_parent_map is None
                else np.asarray(result.elem_parent_map, dtype=int)
            )
            composed = event_map[np.asarray(patch_map, dtype=int)]
            if composed.shape != (patch_mesh.ne,):
                raise RuntimeError("composed remesh parent map has invalid shape")

            existing_count = len(cohesive_before)
            physical_event_index = int(self.physical_event_counter)
            new_cohesive = self.cohesive_network.elements[existing_count:]
            if not new_cohesive:
                raise RuntimeError("physical crack event created no cohesive topology")
            for ce in new_cohesive:
                ce.metadata = dict(ce.metadata)
                ce.metadata["physical_event_index_v914"] = physical_event_index
                ce.metadata["physical_event_subsegment_count_v914"] = len(new_cohesive)
            cohesive_after = self._cohesive_signature(self.cohesive_network)
            preexisting_unchanged = cohesive_before == cohesive_after[:existing_count]

            ACTIVE_CONTEXT.set_solver(__import__(
                "arrhenius_fracture.fem", fromlist=["solve_dirichlet"]
            ).solve_dirichlet)
            segments = [
                (np.array([r["x0"], r["y0"]], float),
                 np.array([r["x1"], r["y1"]], float))
                for r in self.advance_log
            ]
            equilibrium_record: dict[str, Any] = {}
            if self.require_post_event_equilibrium:
                patch_u, equilibrium_record = ACTIVE_CONTEXT.equilibrate(
                    pre_mesh=pre_mesh,
                    pre_boundary=kwargs["boundary"],
                    pre_displacement=kwargs["displacement"],
                    new_mesh=patch_mesh,
                    new_boundary=patch_bnd,
                    new_damage=patch_d,
                    new_displacement=patch_u,
                    parent_map=composed,
                    cohesive_network=self.cohesive_network,
                    new_tip=tip,
                    direction=direction,
                    crack_segments=segments,
                    event_index=physical_event_index,
                    front_id=int(kwargs.get("front_id", 0)),
                )

            event_audit = {
                "physical_event_index": physical_event_index,
                "front_id": int(kwargs.get("front_id", 0)),
                "x0_m": float(start[0]),
                "y0_m": float(start[1]),
                "x1_m": float(tip[0]),
                "y1_m": float(tip[1]),
                "physical_advance_m": float(result.moved),
                "n_topological_subsegments": int(len(new_cohesive)),
                "one_physical_cohesive_event": True,
                "n_elements_before_event": int(pre_mesh.ne),
                "n_elements_after_event_before_patch": int(result.mesh.ne),
                "n_elements_after_patch": int(patch_mesh.ne),
                "preexisting_cohesive_state_unchanged": bool(preexisting_unchanged),
                "cohesive_count_before_event": int(existing_count),
                "cohesive_count_after_event_and_patch": int(len(cohesive_after)),
                "parent_map_min": int(np.min(composed)) if composed.size else -1,
                "parent_map_max": int(np.max(composed)) if composed.size else -1,
                "parent_map_valid": bool(
                    composed.size == patch_mesh.ne
                    and np.min(composed, initial=0) >= 0
                    and np.max(composed, initial=-1) < pre_mesh.ne
                ),
                "post_event_equilibrium_required": self.require_post_event_equilibrium,
                "post_event_equilibrium_completed": bool(
                    (not self.require_post_event_equilibrium) or equilibrium_record
                ),
                **patch_audit,
                **{f"equilibrium_{k}": v for k, v in equilibrium_record.items()},
            }
            self.remesh_audit.append(event_audit)
            self.physical_event_counter += 1
            for row in new_log_rows:
                row.update({
                    "event_remesh_v914": True,
                    "physical_event_index_v914": physical_event_index,
                    "physical_event_subsegments_v914": len(new_cohesive),
                    "event_remesh_edge_splits": int(patch_audit["n_edge_splits"]),
                    "event_remesh_hbar_tip_after_m": float(patch_mesh.hbar_tip),
                    "event_remesh_parent_area_error": float(
                        patch_audit["max_parent_relative_area_conservation_error"]
                    ),
                    "post_event_same_load_equilibrium": bool(equilibrium_record),
                })
            return CrackAdvanceResult(
                patch_mesh, patch_bnd, patch_d, patch_u,
                float(result.moved), True,
                angle_error_deg=float(result.angle_error_deg),
                selected_edge_length=float(result.selected_edge_length),
                reason="ok_event_remeshed_and_equilibrated",
                elem_parent_map=np.ascontiguousarray(composed, dtype=int),
            )
        except Exception as exc:
            return self._rollback_failed_event(
                transaction,
                kwargs,
                f"event_remesh_or_equilibrium_error:{type(exc).__name__}:{exc}",
                {"physical_event_index": int(self.physical_event_counter)},
            )

    def write_diagnostics(self, out_dir: str) -> None:
        super().write_diagnostics(out_dir)
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": "event_remesh_czm_v914",
            "backend": self.name,
            "target_h_m": self.target_h_m,
            "target_edge_factor": self.target_edge_factor,
            "patch_radius_m": self.patch_radius_m,
            "max_edge_splits_per_event": self.max_edge_splits_per_event,
            "require_post_event_equilibrium": self.require_post_event_equilibrium,
            "fail_fast_on_event_error": self.fail_fast_on_event_error,
            "events": self.remesh_audit,
            "failed_event_attempts": self.failed_event_attempts,
            "n_events": len(self.remesh_audit),
            "n_failed_event_attempts": len(self.failed_event_attempts),
            "all_parent_maps_valid": bool(
                self.remesh_audit and all(x["parent_map_valid"] for x in self.remesh_audit)
            ),
            "all_patch_targets_satisfied": bool(
                self.remesh_audit and all(x["patch_target_satisfied"] for x in self.remesh_audit)
            ),
            "all_preexisting_cohesive_states_unchanged": bool(
                self.remesh_audit
                and all(x["preexisting_cohesive_state_unchanged"] for x in self.remesh_audit)
            ),
            "all_events_one_physical_cohesive_event": bool(
                self.remesh_audit
                and all(x["one_physical_cohesive_event"] for x in self.remesh_audit)
            ),
            "all_post_event_equilibria_completed": bool(
                self.remesh_audit
                and all(x["post_event_equilibrium_completed"] for x in self.remesh_audit)
            ),
            "max_parent_relative_area_conservation_error": max(
                (float(x["max_parent_relative_area_conservation_error"])
                 for x in self.remesh_audit), default=float("nan")
            ),
            "max_relative_total_area_error": max(
                (float(x["relative_total_area_error"])
                 for x in self.remesh_audit), default=float("nan")
            ),
            "max_same_load_boundary_drift": max(
                (float(x.get("equilibrium_max_relative_boundary_displacement_drift", np.nan))
                 for x in self.remesh_audit), default=float("nan")
            ),
            "max_post_equilibrium_free_residual_N": max(
                (float(x.get("equilibrium_free_residual_norm_after_N", np.nan))
                 for x in self.remesh_audit), default=float("nan")
            ),
        }
        (out / "event_remesh_audit_v914.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
        rows = [{k: v for k, v in event.items() if k != "split_rows"}
                for event in self.remesh_audit]
        if rows:
            with (out / "event_remesh_audit_v914.csv").open("w", newline="") as fp:
                writer = csv.DictWriter(fp, fieldnames=sorted({k for row in rows for k in row}))
                writer.writeheader()
                writer.writerows(rows)


def build_event_remesh_backend(args, geom) -> EventRemeshCZMBackend:
    """Construct the v9.14 backend from the existing sharp-front namespace."""
    target_h = float(getattr(args, "event_remesh_target_h_m", 0.0) or 0.0)
    if target_h <= 0.0:
        target_h = float(getattr(args, "tip_h_fine", 0.0) or 1.0e-6)
    da = float(getattr(args, "da_phys", 0.0) or 5.0e-6)
    patch_radius = float(getattr(args, "event_remesh_patch_radius_m", 0.0) or 0.0)
    if patch_radius <= 0.0:
        patch_radius = max(20.0 * target_h, 4.0 * da, 20.0e-6)
    ACTIVE_CONTEXT.set_solver(__import__(
        "arrhenius_fracture.fem", fromlist=["solve_dirichlet"]
    ).solve_dirichlet)
    return EventRemeshCZMBackend(
        geom=geom,
        penalty_normal_Pa_per_m=float(getattr(args, "czm_penalty_normal", 1.0e18)),
        penalty_tangent_Pa_per_m=float(getattr(args, "czm_penalty_tangent", 1.0e18)),
        max_angle_error_deg=float(getattr(args, "czm_max_angle_error_deg", 35.0)),
        event_damage=float(getattr(args, "czm_event_damage", 1.0)),
        min_area_ratio=float(getattr(args, "czm_min_area_ratio", 0.08)),
        min_triangle_quality=float(getattr(args, "czm_min_triangle_quality", 0.035)),
        max_node_move_factor=float(getattr(args, "czm_max_node_move_factor", 1.75)),
        max_hrefine_subsegments=int(getattr(args, "czm_max_hrefine_subsegments", 512)),
        target_h_m=target_h,
        patch_radius_m=patch_radius,
        max_edge_splits_per_event=int(
            getattr(args, "event_remesh_max_edge_splits", 256) or 256
        ),
        target_edge_factor=float(
            getattr(args, "event_remesh_target_edge_factor", 1.25) or 1.25
        ),
        forward_back_margin_m=float(
            getattr(args, "event_remesh_back_margin_m", target_h) or target_h
        ),
        min_remesh_triangle_quality=float(
            getattr(args, "event_remesh_min_quality", 0.02) or 0.02
        ),
        require_post_event_equilibrium=bool(
            getattr(args, "event_remesh_require_equilibrium", True)
        ),
        fail_fast_on_event_error=bool(
            getattr(args, "event_remesh_fail_fast", True)
        ),
    )


__all__ = ["EventRemeshCZMBackend", "build_event_remesh_backend"]
