"""Same-time, same-load post-event equilibrium support for v9.14.

The sharp-front loop commits an Arrhenius renewal after an equilibrium/J solve.
The v9.14 crack backend changes the mesh for that one event and then calls this
context before returning.  The context transfers the piecewise-constant bulk
state through the exact parent map, solves the remeshed FEM at the *existing*
Dirichlet displacement, recomputes stress/energy, refreshes the 2-D MPZ sampling
profile and evaluates J at the new tip.  No physical time, remote displacement
or hazard action is advanced here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np


@dataclass
class MechanicsSnapshot:
    mesh: Any
    displacement: np.ndarray
    ep_gp: np.ndarray
    rho_gp: np.ndarray
    damage: np.ndarray
    D: np.ndarray
    material: Any
    kappa: float
    cohesive_network: Any


@dataclass
class JCallSnapshot:
    callable: Callable[..., Any]
    material: Any
    ell: float
    cfg: Any
    exclude_radius: float
    J_before: float
    K_before: float
    coupling_context: Any = None


@dataclass
class EventEquilibriumContext:
    original_assemble: Callable[..., Any] | None = None
    solve_callback: Callable[..., Any] | None = None
    latest_mechanics: MechanicsSnapshot | None = None
    latest_j: JCallSnapshot | None = None
    records: list[dict[str, Any]] = field(default_factory=list)

    def clear(self) -> None:
        self.latest_mechanics = None
        self.latest_j = None
        self.records.clear()

    def set_solver(self, callback: Callable[..., Any]) -> None:
        self.solve_callback = callback

    def record_mechanics(
        self,
        mesh,
        displacement,
        ep_gp,
        rho_gp,
        damage,
        D,
        material,
        kappa,
        cohesive_network,
    ) -> None:
        # State is copied because the outer staggered loop may update its arrays
        # after assembly and before the crack event is committed.  v9.14 is gated
        # to tip_only production use, for which these fields are unchanged by the
        # intervening plasticity call; the copies make the contract explicit.
        self.latest_mechanics = MechanicsSnapshot(
            mesh=mesh,
            displacement=np.asarray(displacement, float).copy(),
            ep_gp=np.asarray(ep_gp, float).copy(),
            rho_gp=np.asarray(rho_gp, float).copy(),
            damage=np.asarray(damage, float).copy(),
            D=np.asarray(D, float).copy(),
            material=material,
            kappa=float(kappa),
            cohesive_network=cohesive_network,
        )

    def record_j_call(
        self,
        callback: Callable[..., Any],
        material,
        ell: float,
        cfg,
        exclude_radius: float,
        J: float,
        K: float,
        coupling_context=None,
    ) -> None:
        self.latest_j = JCallSnapshot(
            callable=callback,
            material=material,
            ell=float(ell),
            cfg=cfg,
            exclude_radius=float(exclude_radius),
            J_before=float(J),
            K_before=float(K),
            coupling_context=coupling_context,
        )

    @staticmethod
    def _prescribed_mask(mesh, boundary) -> np.ndarray:
        prescribed = np.zeros(mesh.ndof, dtype=bool)
        prescribed[2 * np.asarray(boundary.top_nodes, dtype=int) + 1] = True
        prescribed[2 * np.asarray(boundary.bot_nodes, dtype=int) + 1] = True
        prescribed[2 * int(boundary.left_bot)] = True
        prescribed[2 * int(boundary.left_bot) + 1] = True
        prescribed[2 * int(boundary.right_bot)] = True
        return prescribed

    @staticmethod
    def _boundary_values(displacement: np.ndarray, boundary) -> tuple[float, float]:
        u = np.asarray(displacement, float)
        top = np.asarray(boundary.top_nodes, dtype=int)
        bot = np.asarray(boundary.bot_nodes, dtype=int)
        if top.size == 0 or bot.size == 0:
            raise RuntimeError("post-event equilibrium requires non-empty top/bottom boundary sets")
        return (
            float(np.mean(u[2 * top + 1])),
            float(np.mean(u[2 * bot + 1])),
        )

    @staticmethod
    def _refresh_coupling_context(js: JCallSnapshot, rho_new: np.ndarray) -> None:
        """Expose transferred density to the v9.11 J/profile wrapper.

        The J wrapper also samples the finite-radius FEM stress/density profile
        used by the front MPZ.  Without this update it would see the old mesh-sized
        rho array during the same-load post-event J evaluation.
        """
        context = js.coupling_context
        if context is None:
            return
        context.bulk_rho_gp = np.asarray(rho_new, float).copy()
        if hasattr(context, "bulk_retained_rho_gp"):
            previous = getattr(context, "bulk_retained_rho_gp", None)
            if previous is not None:
                context.bulk_retained_rho_gp = np.asarray(rho_new, float).copy()
        if hasattr(context, "bulk_mobile_rho_gp"):
            previous = getattr(context, "bulk_mobile_rho_gp", None)
            if previous is not None:
                context.bulk_mobile_rho_gp = np.zeros_like(rho_new, dtype=float)

    def equilibrate(
        self,
        *,
        pre_mesh,
        pre_boundary,
        pre_displacement: np.ndarray,
        new_mesh,
        new_boundary,
        new_damage: np.ndarray,
        new_displacement: np.ndarray,
        parent_map: np.ndarray,
        cohesive_network,
        new_tip: np.ndarray,
        direction: np.ndarray,
        crack_segments,
        event_index: int,
        front_id: int,
    ) -> tuple[np.ndarray, dict[str, Any]]:
        if self.original_assemble is None or self.solve_callback is None:
            raise RuntimeError("v9.14 post-event equilibrium hooks are not installed")
        snap = self.latest_mechanics
        if snap is None:
            raise RuntimeError("no pre-event mechanics snapshot is available")
        if snap.mesh is not pre_mesh:
            raise RuntimeError("pre-event mechanics snapshot does not match crack-event mesh")

        pmap = np.asarray(parent_map, dtype=int)
        if pmap.shape != (new_mesh.ne,):
            raise RuntimeError(
                f"parent map has shape {pmap.shape}; expected ({new_mesh.ne},)"
            )
        if np.min(pmap, initial=0) < 0 or np.max(pmap, initial=-1) >= pre_mesh.ne:
            raise RuntimeError("parent map references an invalid pre-event element")

        ep_new = np.ascontiguousarray(snap.ep_gp[:, pmap])
        rho_new = np.ascontiguousarray(snap.rho_gp[pmap])
        u0 = np.asarray(new_displacement, float).copy()
        dnew = np.asarray(new_damage, float)
        Uy_top, Uy_bot = self._boundary_values(pre_displacement, pre_boundary)

        K0, R0, sigma0, seq0, s10, psi0 = self.original_assemble(
            new_mesh,
            u0,
            ep_new,
            rho_new,
            dnew,
            snap.D,
            snap.material,
            kappa=snap.kappa,
            cohesive_network=cohesive_network,
        )
        prescribed = self._prescribed_mask(new_mesh, new_boundary)
        free = ~prescribed
        residual_before = float(np.linalg.norm(R0[free])) if np.any(free) else 0.0

        ueq, Ftop = self.solve_callback(
            K0, R0, u0, new_boundary, Uy_top, Uy_bot
        )
        K1, R1, sigma1, seq1, s11, psi1 = self.original_assemble(
            new_mesh,
            ueq,
            ep_new,
            rho_new,
            dnew,
            snap.D,
            snap.material,
            kappa=snap.kappa,
            cohesive_network=cohesive_network,
        )
        residual_after = float(np.linalg.norm(R1[free])) if np.any(free) else 0.0
        top_after, bot_after = self._boundary_values(ueq, new_boundary)

        J_after = float("nan")
        K_after = float("nan")
        j_active_elements = 0
        j_reason = "no_recorded_pre_event_J_call"
        mpz_profile_recomputed = False
        mpz_profile_reliable = False
        if self.latest_j is not None:
            js = self.latest_j
            self._refresh_coupling_context(js, rho_new)
            try:
                J_after, K_after, jinfo = js.callable(
                    new_mesh,
                    ueq,
                    sigma1,
                    psi1,
                    dnew,
                    np.asarray(new_tip, float),
                    np.asarray(direction, float),
                    js.material,
                    js.ell,
                    cfg=js.cfg,
                    crack_segments=crack_segments,
                    exclude_radius=js.exclude_radius,
                )
                J_after = float(J_after)
                K_after = float(K_after)
                j_active_elements = int(jinfo.get("n_active_elements", 0))
                j_reason = "ok"
                if js.coupling_context is not None:
                    mpz_profile_recomputed = getattr(
                        js.coupling_context, "mpz_profile_2d", None
                    ) is not None
                    mpz_profile_reliable = bool(
                        jinfo.get("mpz_2d_profile_reliable", False)
                    )
            except Exception as exc:
                j_reason = f"{type(exc).__name__}:{exc}"

        total_area_old = float(np.sum(pre_mesh.area_e))
        total_area_new = float(np.sum(new_mesh.area_e))
        rho_integral_old = float(np.sum(snap.rho_gp * pre_mesh.area_e))
        rho_integral_new = float(np.sum(rho_new * new_mesh.area_e))
        ep_integral_old = np.sum(snap.ep_gp * pre_mesh.area_e[None, :], axis=1)
        ep_integral_new = np.sum(ep_new * new_mesh.area_e[None, :], axis=1)
        elastic_energy_before_solve = float(np.sum(psi0 * new_mesh.area_e))
        elastic_energy_after_solve = float(np.sum(psi1 * new_mesh.area_e))

        scale_u = max(abs(Uy_top), abs(Uy_bot), 1.0e-30)
        record = {
            "event_index": int(event_index),
            "front_id": int(front_id),
            "same_time_equilibrium": True,
            "physical_time_increment_s": 0.0,
            "hazard_action_increment": 0.0,
            "Uy_top_event_m": Uy_top,
            "Uy_bottom_event_m": Uy_bot,
            "Uy_top_after_equilibrium_m": top_after,
            "Uy_bottom_after_equilibrium_m": bot_after,
            "max_relative_boundary_displacement_drift": float(
                max(abs(top_after - Uy_top), abs(bot_after - Uy_bot)) / scale_u
            ),
            "free_residual_norm_before_N": residual_before,
            "free_residual_norm_after_N": residual_after,
            "free_residual_reduction": float(
                residual_after / max(residual_before, 1.0e-300)
            ),
            "reaction_top_after_N": float(Ftop),
            "J_before_event_N_per_m": (
                float(self.latest_j.J_before) if self.latest_j is not None else float("nan")
            ),
            "KJ_before_event_Pa_sqrt_m": (
                float(self.latest_j.K_before) if self.latest_j is not None else float("nan")
            ),
            "J_after_event_equilibrium_N_per_m": J_after,
            "KJ_after_event_equilibrium_Pa_sqrt_m": K_after,
            "J_after_event_active_elements": int(j_active_elements),
            "J_after_event_status": j_reason,
            "mpz_profile_recomputed_after_event": bool(mpz_profile_recomputed),
            "mpz_profile_reliable_after_event": bool(mpz_profile_reliable),
            "elastic_energy_before_equilibrium_J_per_m": elastic_energy_before_solve,
            "elastic_energy_after_equilibrium_J_per_m": elastic_energy_after_solve,
            "total_mesh_area_before_m2": total_area_old,
            "total_mesh_area_after_m2": total_area_new,
            "relative_total_mesh_area_error": float(
                abs(total_area_new - total_area_old) / max(abs(total_area_old), 1.0e-300)
            ),
            "rho_area_integral_before": rho_integral_old,
            "rho_area_integral_after": rho_integral_new,
            "relative_rho_area_integral_error": float(
                abs(rho_integral_new - rho_integral_old)
                / max(abs(rho_integral_old), 1.0e-300)
            ),
            "ep_area_integral_before": ep_integral_old.tolist(),
            "ep_area_integral_after": ep_integral_new.tolist(),
            "max_relative_ep_area_integral_error": float(
                np.max(
                    np.abs(ep_integral_new - ep_integral_old)
                    / np.maximum(np.abs(ep_integral_old), 1.0e-30)
                )
            ),
            "n_elements_before": int(pre_mesh.ne),
            "n_elements_after": int(new_mesh.ne),
            "n_nodes_before": int(pre_mesh.nn),
            "n_nodes_after": int(new_mesh.nn),
        }
        self.records.append(record)
        return np.asarray(ueq, float), record


ACTIVE_CONTEXT = EventEquilibriumContext()


def install_mechanics_recorder(fem_module) -> Callable[..., Any]:
    """Patch FEM assembly so the event backend sees the accepted pre-event state."""
    original = fem_module.assemble_mechanics
    ACTIVE_CONTEXT.original_assemble = original

    def wrapped(
        mesh,
        u,
        ep_gp,
        rho_gp,
        d,
        D,
        mat,
        kappa=1.0e-6,
        cohesive_network=None,
    ):
        out = original(
            mesh,
            u,
            ep_gp,
            rho_gp,
            d,
            D,
            mat,
            kappa=kappa,
            cohesive_network=cohesive_network,
        )
        ACTIVE_CONTEXT.record_mechanics(
            mesh,
            u,
            ep_gp,
            rho_gp,
            d,
            D,
            mat,
            kappa,
            cohesive_network,
        )
        return out

    fem_module.assemble_mechanics = wrapped
    return original


def restore_mechanics_recorder(fem_module, original: Callable[..., Any]) -> None:
    fem_module.assemble_mechanics = original


__all__ = [
    "ACTIVE_CONTEXT",
    "EventEquilibriumContext",
    "install_mechanics_recorder",
    "restore_mechanics_recorder",
]
