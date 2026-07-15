"""Explicit bulk mobile/retained state for the v9.11 full 2-D transfer gate.

Two physically distinct modes are supported without changing any calibrated
barrier parameter:

``tip_only``
    The continuum FEM is elastic. The moving crack-tip MPZ remains the sole
    source and owner of mobile/retained dislocations.

``bulk_same_pt_km``
    Every integration point carries explicit mobile and retained/forest density.
    The same class manifest supplies the independent Arrhenius Peierls and Taylor
    surfaces used by the crack-tip MPZ. Peierls controls mobile glide, geometric
    encounters retain mobile content, Taylor completion releases retained
    content, and the existing Kocks--Mecking storage/recovery coefficients evolve
    total bulk density. No fixed mobile fraction is used.

The initial validation uses the static-mesh ``adaptive_czm`` backend. A change in
integration-point count is rejected in ``bulk_same_pt_km`` rather than silently
reinitializing the explicit carrier inventory. Conservative remesh transfer is a
separate gate before this mode is used with tip-following remeshing.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .bulk_plasticity_v9102 import independent_config_from_dislocation_config
from .emission_derived_plasticity_v9102 import EmissionDerivedPeierlsTaylorModel
from .mpz_parameterization_v911 import apply_pt_dislocation_config

VALID_BULK_MODES = ("tip_only", "bulk_same_pt_km")


def normalize_bulk_mode(value: str) -> str:
    key = str(value).strip().lower().replace("-", "_")
    aliases = {
        "tip": "tip_only",
        "tip_source_only": "tip_only",
        "tip_sources_only": "tip_only",
        "elastic_bulk": "tip_only",
        "bulk": "bulk_same_pt_km",
        "bulk_uniform": "bulk_same_pt_km",
        "bulk_uniform_hardening": "bulk_same_pt_km",
        "bulk_pt_km": "bulk_same_pt_km",
    }
    key = aliases.get(key, key)
    if key not in VALID_BULK_MODES:
        raise ValueError(
            f"unknown bulk plasticity mode {value!r}; expected one of "
            f"{', '.join(VALID_BULK_MODES)}"
        )
    return key


def _von_mises_plane_strain(sigma_gp: np.ndarray, nu: float):
    sigma = np.asarray(sigma_gp, dtype=float)
    sx, sy, txy = sigma[0], sigma[1], sigma[2]
    szz = float(nu) * (sx + sy)
    p = (sx + sy + szz) / 3.0
    sd_xx, sd_yy, sd_zz = sx - p, sy - p, szz - p
    norm_s = np.sqrt(sd_xx**2 + sd_yy**2 + sd_zz**2 + 2.0 * txy**2)
    seq = np.sqrt(1.5) * norm_s
    return seq, norm_s, sd_xx, sd_yy, txy


def _zero_info(ne: int, mode: str) -> dict[str, Any]:
    z = np.zeros(int(ne), dtype=float)
    return {
        "dWp_requested_gp": z.copy(),
        "dWp_accepted_gp": z.copy(),
        "dep_eq_requested_gp": z.copy(),
        "dep_eq_accepted_gp": z.copy(),
        "thermo_scale_gp": np.ones(int(ne), dtype=float),
        "thermo_admissible_gp": z.copy(),
        "thermo_hazard_gp": z.copy(),
        "thermo_mode": "bulk_disabled" if mode == "tip_only" else "time_cone",
        "bulk_plasticity_mode": mode,
        "bulk_fixed_mobile_fraction_active": False,
        "bulk_explicit_mobile_retained_state": mode == "bulk_same_pt_km",
    }


class BulkPlasticityControllerV911:
    """Mode-aware wrapper around the full solver's plasticity update."""

    def __init__(self, mode: str, row: dict[str, Any], context: Any | None = None):
        self.mode = normalize_bulk_mode(mode)
        self.row = dict(row)
        self.context = context
        self.mobile_rho_m2: np.ndarray | None = None
        self.retained_rho_m2: np.ndarray | None = None
        self.last_temperature_K: float | None = None
        self.calls = 0
        self.reset_count = 0
        self.storage_added_mean_m2 = 0.0
        self.dynamic_recovery_removed_mean_m2 = 0.0
        self.static_recovery_removed_mean_m2 = 0.0
        self.exchange_trapped_mean_m2 = 0.0
        self.exchange_released_mean_m2 = 0.0
        self.accepted_dep_mean_acc = 0.0
        self.max_hazard = 0.0
        self.mesh_change_rejected = False

    @property
    def explicit_state_active(self) -> bool:
        return self.mode == "bulk_same_pt_km"

    def reset(self) -> None:
        self.mobile_rho_m2 = None
        self.retained_rho_m2 = None
        self.reset_count += 1

    def _ensure_state(self, rho_fallback: np.ndarray, T_K: float) -> None:
        rho = np.maximum(np.asarray(rho_fallback, dtype=float).reshape(-1), 0.0)
        if self.last_temperature_K is None or not math.isclose(
            float(T_K), float(self.last_temperature_K), rel_tol=0.0, abs_tol=1.0e-12
        ):
            self.mobile_rho_m2 = None
            self.retained_rho_m2 = None
            self.last_temperature_K = float(T_K)
        if self.mobile_rho_m2 is None or self.retained_rho_m2 is None:
            self.mobile_rho_m2 = np.zeros_like(rho)
            self.retained_rho_m2 = rho.copy()
            return
        if self.mobile_rho_m2.size != rho.size:
            self.mesh_change_rejected = True
            raise RuntimeError(
                "bulk_same_pt_km explicit mobile/retained state encountered a change "
                "in integration-point count. The current validation supports the "
                "static-mesh adaptive_czm backend; conservative remesh transfer must "
                "be implemented before using this mode with tip remeshing."
            )

    @staticmethod
    def _exchange(
        mobile: np.ndarray,
        retained: np.ndarray,
        encounter_rate_s: np.ndarray,
        taylor_release_rate_s: np.ndarray,
        dt_s: float,
    ) -> tuple[np.ndarray, np.ndarray, float, float]:
        """Exact two-state mobile/retained exchange over ``dt_s``."""
        if dt_s <= 0.0:
            return mobile.copy(), retained.copy(), 0.0, 0.0
        m0 = np.maximum(np.asarray(mobile, dtype=float), 0.0)
        r0 = np.maximum(np.asarray(retained, dtype=float), 0.0)
        ke = np.maximum(np.asarray(encounter_rate_s, dtype=float), 0.0)
        kt = np.maximum(np.asarray(taylor_release_rate_s, dtype=float), 0.0)
        total = m0 + r0
        rate = ke + kt
        frac_r_eq = np.divide(ke, rate, out=np.zeros_like(rate), where=rate > 0.0)
        r_eq = frac_r_eq * total
        decay = np.exp(-np.minimum(rate * float(dt_s), 700.0))
        r1 = r_eq + (r0 - r_eq) * decay
        r1 = np.clip(r1, 0.0, total)
        m1 = total - r1
        trapped = float(np.mean(np.maximum(r1 - r0, 0.0)))
        released = float(np.mean(np.maximum(r0 - r1, 0.0)))
        return m1, r1, trapped, released

    @staticmethod
    def _static_recovery_rate(
        rho_total: np.ndarray,
        mat: Any,
        disl_cfg: Any,
        T_K: float,
    ) -> np.ndarray:
        out = np.zeros_like(rho_total)
        if not bool(getattr(disl_cfg, "use_static_recovery", False)):
            return out
        if float(T_K) / max(float(getattr(mat, "Tm", 1.0)), 1.0) < float(
            getattr(disl_cfg, "Tfrac_on", 1.0)
        ):
            return out
        kev = 8.617333262e-5
        b = max(abs(float(mat.b)), 1.0e-30)
        Dl = (
            float(getattr(disl_cfg, "Dl0a", 0.0))
            * np.exp(-float(getattr(disl_cfg, "Ea_eV", 0.0)) / (kev * max(float(T_K), 1.0e-9)))
            + float(getattr(disl_cfg, "Dl0b", 0.0))
            * np.exp(-float(getattr(disl_cfg, "Eb_eV", 0.0)) / (kev * max(float(T_K), 1.0e-9)))
        )
        kB = 1.380649e-23
        x = (
            float(getattr(disl_cfg, "kpp", 0.0))
            * float(mat.E)
            * b**4
            * np.sqrt(np.maximum(rho_total, 1.0e-300))
            / max(kB * float(T_K), 1.0e-30)
        )
        x = np.clip(x, 0.0, 700.0)
        exm1 = np.maximum(np.exp(x) - 1.0, 1.0e-30)
        rho32 = rho_total * np.sqrt(np.maximum(rho_total, 1.0e-300))
        out = float(getattr(disl_cfg, "kprime", 0.0)) * (Dl / b) * rho32 / exm1
        return np.clip(out, 0.0, float(getattr(disl_cfg, "gamma_cap", np.inf)))

    def _apply_storage_recovery(
        self,
        dot_ep_gp: np.ndarray,
        mat: Any,
        disl_cfg: Any,
        T_K: float,
        dt_s: float,
    ) -> tuple[float, float, float]:
        """Use the repository's existing Kocks--Mecking coefficients."""
        assert self.mobile_rho_m2 is not None
        assert self.retained_rho_m2 is not None
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0:
            return 0.0, 0.0, 0.0
        mobile = np.maximum(self.mobile_rho_m2, 0.0)
        retained = np.maximum(self.retained_rho_m2, 0.0)
        total = mobile + retained
        forest = np.maximum(retained, float(getattr(disl_cfg, "pt_forest_density_floor_m2", 0.0)))
        dot_ep = np.maximum(np.asarray(dot_ep_gp, dtype=float), 0.0)
        b = max(abs(float(mat.b)), 1.0e-30)
        storage_rate = (
            float(getattr(disl_cfg, "k_store", 0.0))
            * np.sqrt(np.maximum(forest, 1.0e-300))
            / b
            * dot_ep
        )
        dynamic_rate = float(getattr(disl_cfg, "k_dyn", 0.0)) * total * dot_ep
        static_rate = self._static_recovery_rate(total, mat, disl_cfg, T_K)

        storage = dt * np.maximum(storage_rate, 0.0)
        mobile = mobile + storage
        total_after_storage = mobile + retained
        removal = dt * np.maximum(dynamic_rate + static_rate, 0.0)
        removal = np.minimum(removal, total_after_storage)
        frac = np.divide(
            removal,
            total_after_storage,
            out=np.zeros_like(removal),
            where=total_after_storage > 0.0,
        )
        mobile *= 1.0 - frac
        retained *= 1.0 - frac
        self.mobile_rho_m2 = np.maximum(mobile, 0.0)
        self.retained_rho_m2 = np.maximum(retained, 0.0)
        return (
            float(np.mean(storage)),
            float(np.mean(dt * np.maximum(dynamic_rate, 0.0))),
            float(np.mean(dt * np.maximum(static_rate, 0.0))),
        )

    def _state_rates(
        self,
        seq: np.ndarray,
        mat: Any,
        T_K: float,
        disl_cfg: Any,
    ) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
        assert self.mobile_rho_m2 is not None
        assert self.retained_rho_m2 is not None
        cfg = independent_config_from_dislocation_config(disl_cfg)
        model = EmissionDerivedPeierlsTaylorModel(cfg)
        forest = np.maximum(
            self.retained_rho_m2,
            float(getattr(disl_cfg, "pt_forest_density_floor_m2", 0.0)),
        )
        rates = model.rates(
            seq,
            forest,
            T_K,
            mat.b,
            rho_mobile_m2=self.mobile_rho_m2,
        )
        peierls = np.maximum(np.asarray(rates["peierls_rate_s"], dtype=float), 0.0)
        jump = np.maximum(np.asarray(rates["jump_length_m"], dtype=float), 0.0)
        encounter = (
            max(float(getattr(disl_cfg, "pt_encounter_efficiency", 0.0)), 0.0)
            * jump
            * peierls
            * np.sqrt(np.maximum(forest, 0.0))
        )
        return rates, forest, encounter

    def update(
        self,
        ep_gp: np.ndarray,
        rho_gp: np.ndarray,
        sigma_gp: np.ndarray,
        mat: Any,
        T_K: float,
        dt_s: float,
        plast_model: Any,
        disl_cfg: Any,
        *,
        return_info: bool = False,
    ):
        """Apply one bulk update while preserving the sharp-front call contract."""
        del plast_model
        apply_pt_dislocation_config(disl_cfg, self.row)
        disl_cfg.pt_mobile_fraction = 0.0
        self.calls += 1
        ne = int(np.asarray(rho_gp).size)

        if self.mode == "tip_only":
            dot_ep = np.zeros(ne, dtype=float)
            info = _zero_info(ne, self.mode)
            self._publish(np.asarray(rho_gp, dtype=float), None, info)
            if return_info:
                return ep_gp, rho_gp, dot_ep, info
            return ep_gp, rho_gp, dot_ep

        self._ensure_state(rho_gp, T_K)
        assert self.mobile_rho_m2 is not None
        assert self.retained_rho_m2 is not None

        dt = max(float(dt_s), 0.0)
        seq, norm_s, sd_xx, sd_yy, txy = _von_mises_plane_strain(sigma_gp, mat.nu)

        rates0, _, encounter0 = self._state_rates(seq, mat, T_K, disl_cfg)
        taylor0 = np.maximum(
            np.asarray(rates0["taylor_completion_rate_s"], dtype=float), 0.0
        )
        self.mobile_rho_m2, self.retained_rho_m2, trapped0, released0 = self._exchange(
            self.mobile_rho_m2,
            self.retained_rho_m2,
            encounter0,
            taylor0,
            0.5 * dt,
        )

        rates, _, _ = self._state_rates(seq, mat, T_K, disl_cfg)
        peierls = np.maximum(np.asarray(rates["peierls_rate_s"], dtype=float), 0.0)
        jump = np.maximum(np.asarray(rates["jump_length_m"], dtype=float), 0.0)
        eq_factor = max(float(getattr(disl_cfg, "pt_equivalent_strain_factor", 1.0)), 0.0)
        dot_ep_kin = (
            eq_factor
            * np.maximum(self.mobile_rho_m2, 0.0)
            * abs(float(mat.b))
            * jump
            * peierls
        )

        sqrt23 = math.sqrt(2.0 / 3.0)
        G = max(float(mat.G), 1.0e-30)
        dgamma_eq = 0.999 * np.maximum(seq, 0.0) / (3.0 * G)
        dep_event = max(float(getattr(disl_cfg, "thermo_event_strain", 1.0e-4)), 1.0e-16)
        hazard = np.clip(dot_ep_kin * dt / dep_event, 0.0, 80.0)
        p_event = 1.0 - np.exp(-hazard)
        dep_eq_max = dgamma_eq / sqrt23
        dep_eq = p_event * dep_eq_max
        dgamma = sqrt23 * np.maximum(dep_eq, 0.0)
        dot_ep_gp = np.where(dt > 0.0, dep_eq / dt, 0.0)

        safe_norm = np.maximum(norm_s, 1.0e-30)
        ep_gp[0, :] += dgamma * 1.5 * (sd_xx / safe_norm)
        ep_gp[1, :] += dgamma * 1.5 * (sd_yy / safe_norm)
        ep_gp[2, :] += dgamma * 1.5 * (txy / safe_norm)

        storage, dyn_removed, static_removed = self._apply_storage_recovery(
            dot_ep_gp, mat, disl_cfg, T_K, dt
        )

        rates1, _, encounter1 = self._state_rates(seq, mat, T_K, disl_cfg)
        taylor1 = np.maximum(
            np.asarray(rates1["taylor_completion_rate_s"], dtype=float), 0.0
        )
        self.mobile_rho_m2, self.retained_rho_m2, trapped1, released1 = self._exchange(
            self.mobile_rho_m2,
            self.retained_rho_m2,
            encounter1,
            taylor1,
            0.5 * dt,
        )

        seq_after = np.maximum(seq - 3.0 * G * dgamma, 0.0)
        dWp_requested = seq * np.maximum(dot_ep_kin * dt, 0.0)
        dWp_accepted = 0.5 * (seq + seq_after) * np.maximum(dep_eq, 0.0)
        scale = np.divide(
            dWp_accepted,
            np.maximum(dWp_requested, 1.0e-300),
            out=np.ones_like(dWp_accepted),
            where=dWp_requested > 0.0,
        )

        self.storage_added_mean_m2 += storage
        self.dynamic_recovery_removed_mean_m2 += dyn_removed
        self.static_recovery_removed_mean_m2 += static_removed
        self.exchange_trapped_mean_m2 += trapped0 + trapped1
        self.exchange_released_mean_m2 += released0 + released1
        self.accepted_dep_mean_acc += float(np.mean(dep_eq))
        self.max_hazard = max(self.max_hazard, float(np.max(hazard)) if hazard.size else 0.0)

        info = {
            "dWp_requested_gp": dWp_requested,
            "dWp_accepted_gp": dWp_accepted,
            "dep_eq_requested_gp": dot_ep_kin * dt,
            "dep_eq_accepted_gp": dep_eq,
            "thermo_scale_gp": scale,
            "thermo_admissible_gp": (seq > 0.0).astype(float),
            "thermo_hazard_gp": hazard,
            "thermo_mode": "time_cone_explicit_mobile_retained",
            "bulk_plasticity_mode": self.mode,
            "bulk_fixed_mobile_fraction_active": False,
            "bulk_explicit_mobile_retained_state": True,
            "bulk_mobile_rho_mean_m2": float(np.mean(self.mobile_rho_m2)),
            "bulk_mobile_rho_max_m2": float(np.max(self.mobile_rho_m2)),
            "bulk_retained_rho_mean_m2": float(np.mean(self.retained_rho_m2)),
            "bulk_retained_rho_max_m2": float(np.max(self.retained_rho_m2)),
            "bulk_storage_added_mean_m2": storage,
            "bulk_dynamic_recovery_removed_mean_m2": dyn_removed,
            "bulk_static_recovery_removed_mean_m2": static_removed,
            "bulk_exchange_trapped_mean_m2": trapped0 + trapped1,
            "bulk_exchange_released_mean_m2": released0 + released1,
            "bulk_peierls_rate_max_s": float(np.max(peierls)) if peierls.size else 0.0,
            "bulk_taylor_release_rate_max_s": float(np.max(taylor1)) if taylor1.size else 0.0,
        }

        rho_out = self.retained_rho_m2.copy()
        self._publish(rho_out, self.mobile_rho_m2, info)
        if return_info:
            return ep_gp, rho_out, dot_ep_gp, info
        return ep_gp, rho_out, dot_ep_gp

    def _publish(
        self,
        retained: np.ndarray,
        mobile: np.ndarray | None,
        info: dict[str, Any],
    ) -> None:
        if self.context is None:
            return
        self.context.bulk_rho_gp = np.asarray(retained, dtype=float).copy()
        self.context.bulk_retained_rho_gp = np.asarray(retained, dtype=float).copy()
        self.context.bulk_mobile_rho_gp = (
            None if mobile is None else np.asarray(mobile, dtype=float).copy()
        )
        compact = {
            k: v
            for k, v in info.items()
            if np.isscalar(v) and not isinstance(v, np.ndarray)
        }
        self.context.latest.update(compact)
        if self.context.records:
            self.context.records[-1].update(compact)

    def summary(self) -> dict[str, Any]:
        mobile = self.mobile_rho_m2
        retained = self.retained_rho_m2
        return {
            "bulk_plasticity_mode": self.mode,
            "bulk_same_manifest_peierls_taylor": True,
            "bulk_fixed_mobile_fraction_active": False,
            "bulk_explicit_mobile_retained_state": self.explicit_state_active,
            "bulk_initial_mobile_fraction_assumption": None,
            "bulk_state_update_calls": int(self.calls),
            "bulk_state_reset_count": int(self.reset_count),
            "bulk_mesh_change_rejected": bool(self.mesh_change_rejected),
            "bulk_mobile_rho_mean_m2": None if mobile is None else float(np.mean(mobile)),
            "bulk_mobile_rho_max_m2": None if mobile is None else float(np.max(mobile)),
            "bulk_retained_rho_mean_m2": None if retained is None else float(np.mean(retained)),
            "bulk_retained_rho_max_m2": None if retained is None else float(np.max(retained)),
            "bulk_storage_added_mean_acc_m2": float(self.storage_added_mean_m2),
            "bulk_dynamic_recovery_removed_mean_acc_m2": float(
                self.dynamic_recovery_removed_mean_m2
            ),
            "bulk_static_recovery_removed_mean_acc_m2": float(
                self.static_recovery_removed_mean_m2
            ),
            "bulk_exchange_trapped_mean_acc_m2": float(self.exchange_trapped_mean_m2),
            "bulk_exchange_released_mean_acc_m2": float(self.exchange_released_mean_m2),
            "bulk_accepted_dep_mean_acc": float(self.accepted_dep_mean_acc),
            "bulk_max_time_cone_hazard": float(self.max_hazard),
            "bulk_hardening_law": (
                "none_bulk_elastic" if self.mode == "tip_only" else
                "existing_Kocks_Mecking_k_store_k_dyn_plus_existing_static_recovery"
            ),
            "bulk_remesh_transfer_status": (
                "not_applicable_bulk_elastic" if self.mode == "tip_only" else
                "static_mesh_validation_only"
            ),
        }

    def write_outputs(self, out: str | Path) -> dict[str, Any]:
        root = Path(out)
        root.mkdir(parents=True, exist_ok=True)
        payload = self.summary()
        (root / "bulk_state_v9_11_summary.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
        np.savez_compressed(
            root / "bulk_state_v9_11_final.npz",
            mode=np.asarray(self.mode),
            mobile_rho_m2=(
                np.asarray([], dtype=float)
                if self.mobile_rho_m2 is None
                else self.mobile_rho_m2
            ),
            retained_rho_m2=(
                np.asarray([], dtype=float)
                if self.retained_rho_m2 is None
                else self.retained_rho_m2
            ),
        )
        return payload


def plasticity_wrapper_factory(
    original: Any,
    controller: BulkPlasticityControllerV911,
):
    """Return a function matching ``plasticity.update_plasticity``."""
    del original

    def wrapped(
        ep_gp,
        rho_gp,
        sigma_gp,
        mat,
        T,
        dt,
        plast_model,
        disl_cfg,
        *args,
        **kwargs,
    ):
        return_info = bool(kwargs.pop("return_info", False))
        if args:
            return_info = bool(args[0])
            if len(args) > 1:
                raise TypeError("unexpected positional arguments to update_plasticity")
        if kwargs:
            raise TypeError(f"unexpected keyword arguments: {sorted(kwargs)}")
        return controller.update(
            ep_gp,
            rho_gp,
            sigma_gp,
            mat,
            T,
            dt,
            plast_model,
            disl_cfg,
            return_info=return_info,
        )

    return wrapped


__all__ = [
    "VALID_BULK_MODES",
    "BulkPlasticityControllerV911",
    "normalize_bulk_mode",
    "plasticity_wrapper_factory",
]
