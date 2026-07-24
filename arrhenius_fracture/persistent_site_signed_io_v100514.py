"""Moving-frame, diagnostics, split, and restart mixin."""
from __future__ import annotations

import copy
from dataclasses import asdict
from typing import Any

import numpy as np

from .persistent_site_registry_v100514 import PersistentSiteRowV100514
from .persistent_site_signed_support_v100514 import (
    SignedShieldingKernelV100514,
    _shift_wake_forward,
    _translate_toward_tip,
)
from .signed_kernel_family_v1005141 import (
    SignedShieldingKernelFamilyV1005141,
    load_signed_shielding_artifact_v1005141,
)


class PersistentSiteSignedIOMixin:
    def advance(self, distance_m: float) -> dict[str, float]:
        distance = max(float(distance_m), 0.0)
        radius_before = self.blunted_radius()
        crossed_totals = {"mobile": 0.0, "retained": 0.0, "slip": 0.0}
        discarded = 0.0
        for active_name, wake_name, group in (
            ("mobile_positive", "wake_mobile_positive", "mobile"),
            ("mobile_negative", "wake_mobile_negative", "mobile"),
            ("retained_positive", "wake_retained_positive", "retained"),
            ("retained_negative", "wake_retained_negative", "retained"),
            ("accumulated_slip_positive", "wake_slip_positive", "slip"),
            ("accumulated_slip_negative", "wake_slip_negative", "slip"),
        ):
            old_wake, lost_wake = _shift_wake_forward(
                getattr(self, wake_name), distance, self.wake_dx
            )
            active, crossed, lost_active = _translate_toward_tip(
                getattr(self, active_name),
                distance,
                self.dx,
                self.wake_n_bins,
                self.wake_dx,
            )
            setattr(self, active_name, active)
            setattr(self, wake_name, old_wake + crossed)
            crossed_totals[group] += float(np.sum(crossed))
            discarded += lost_wake + lost_active
        self.available_sites = self.site_capacity.copy()
        self.tip_source_activity = np.ones(self.n_systems)
        self.advance_total_m += distance
        self.wake_discarded_total += discarded
        # Evaluate the new family state immediately.  This fails closed if an
        # accepted crack advance would require forbidden atlas extrapolation.
        kernel_audit = self.kernel_artifact_audit()
        radius_after = self.blunted_radius()
        geometry = self.source_geometry()
        out = {
            "wake_mobile": crossed_totals["mobile"],
            "wake_retained": crossed_totals["retained"],
            "wake_slip": crossed_totals["slip"],
            "source_sites_refreshed": 0.0,
            "available_site_fraction": 1.0,
            "persistent_source_inventory_active": 0.0,
            "tip_radius_before_advance_m": radius_before,
            "tip_radius_after_advance_m": radius_after,
            "tip_resharpening_by_advance_m": max(
                radius_before - radius_after, 0.0
            ),
            "persistent_site_front_width_m": float(geometry["front_width_m"]),
            "persistent_site_multiplicity_per_system": float(
                geometry["multiplicity_per_system"]
            ),
            "fractional_moving_frame": 1.0,
            "cumulative_crack_path_extension_m": self.advance_total_m,
            "signed_kernel_artifact": kernel_audit,
        }
        self.last_advance = copy.deepcopy(out)
        return out

    def split(self, daughter_fraction: float) -> "PersistentSiteSignedMPZStateV100514":
        fraction = float(np.clip(daughter_fraction, 0.0, 1.0))
        child = self.copy()
        for name in (
            "mobile_positive",
            "mobile_negative",
            "retained_positive",
            "retained_negative",
            "accumulated_slip_positive",
            "accumulated_slip_negative",
            "wake_mobile_positive",
            "wake_mobile_negative",
            "wake_retained_positive",
            "wake_retained_negative",
            "wake_slip_positive",
            "wake_slip_negative",
        ):
            original = getattr(self, name)
            setattr(child, name, original * fraction)
            setattr(self, name, original * (1.0 - fraction))
        child.available_sites = child.site_capacity.copy()
        self.available_sites = self.site_capacity.copy()
        child.advance_total_m = 0.0
        child.current_kernel_snapshot()
        return child

    def diagnostics(
        self, G=None, nu=None, b=None, r0=None, c_blunt=None
    ) -> dict[str, Any]:
        geometry = self.source_geometry()
        rho, sigma_back = self.backstress()
        signed_retained = self.retained_positive - self.retained_negative
        kernel_audit = self.kernel_artifact_audit()
        return {
            "state_model": self.state_model,
            "persistent_site_density_m2": self.candidate.rho_source0_m2,
            "persistent_site_multiplicity_per_system": geometry[
                "multiplicity_per_system"
            ],
            "persistent_site_source_area_m2": geometry["source_area_m2"],
            "persistent_site_front_width_m": geometry["front_width_m"],
            "persistent_site_width_density_m2": geometry["rho_width_m2"],
            "persistent_tip_radius_m": geometry["tip_radius_m"],
            "persistent_active_arc_factor": geometry["active_arc_factor"],
            "front_width_grid_independent": True,
            "ahead_of_tip_dx_used_as_front_width_floor": False,
            "available_site_fraction": 1.0,
            "persistent_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "source_sites_refreshed": 0.0,
            "legacy_source_sites_per_system_active": False,
            "legacy_source_refresh_length_active": False,
            "explicit_recovery_active": False,
            "mobile_positive_total": float(np.sum(self.mobile_positive)),
            "mobile_negative_total": float(np.sum(self.mobile_negative)),
            "retained_positive_total": float(np.sum(self.retained_positive)),
            "retained_negative_total": float(np.sum(self.retained_negative)),
            "accumulated_slip_positive_total": float(
                np.sum(self.accumulated_slip_positive)
            ),
            "accumulated_slip_negative_total": float(
                np.sum(self.accumulated_slip_negative)
            ),
            "signed_retained_total": float(np.sum(signed_retained)),
            "unsigned_backstress_density_mean_m2": float(np.mean(rho)),
            "sigma_back_mean_Pa": float(np.mean(sigma_back)),
            "mpz_K_shield_Pa_sqrt_m": self.shielding_K(),
            "wake_shielding_active": self.wake_shielding,
            "mobile_shield_fraction": self.candidate.mobile_shield_fraction,
            "time_s": self.time_s,
            "advance_total_m": self.advance_total_m,
            "cumulative_crack_path_extension_m": self.advance_total_m,
            "kernel_source": self.kernel_source_path,
            "signed_kernel_artifact": kernel_audit,
        }

    def _validate_state_arrays(self) -> None:
        active_shape = (self.n_systems, self.n_bins)
        wake_shape = (self.n_systems, self.wake_n_bins)
        for name in (
            "mobile_positive",
            "mobile_negative",
            "retained_positive",
            "retained_negative",
            "accumulated_slip_positive",
            "accumulated_slip_negative",
        ):
            array = np.asarray(getattr(self, name), dtype=float)
            if array.shape != active_shape or not np.all(np.isfinite(array)):
                raise ValueError(f"{name} shape/content is invalid for restart")
            setattr(self, name, np.maximum(array, 0.0))
        for name in (
            "wake_mobile_positive",
            "wake_mobile_negative",
            "wake_retained_positive",
            "wake_retained_negative",
            "wake_slip_positive",
            "wake_slip_negative",
        ):
            array = np.asarray(getattr(self, name), dtype=float)
            if array.shape != wake_shape or not np.all(np.isfinite(array)):
                raise ValueError(f"{name} shape/content is invalid for restart")
            setattr(self, name, np.maximum(array, 0.0))
        self.available_sites = self.site_capacity.copy()
        self.tip_source_activity = np.ones(self.n_systems, dtype=float)
        self.current_kernel_snapshot()

    def state_dict(self) -> dict[str, Any]:
        arrays = {
            name: getattr(self, name).tolist()
            for name in (
                "mobile_positive",
                "mobile_negative",
                "retained_positive",
                "retained_negative",
                "accumulated_slip_positive",
                "accumulated_slip_negative",
                "wake_mobile_positive",
                "wake_mobile_negative",
                "wake_retained_positive",
                "wake_retained_negative",
                "wake_slip_positive",
                "wake_slip_negative",
            )
        }
        return {
            "schema": "persistent_site_signed_mpz_v10_0_5_14_1",
            "candidate": asdict(self.candidate),
            "kernel_source": self.kernel_source_path,
            "kernel_artifact": self.kernel_artifact_audit(),
            "state_config": {
                "G_Pa": self.G_Pa,
                "nu": self.nu,
                "b_m": self.b_m,
                "r0_m": self.r0_m,
                "blunting_length_m": self.blunting_length_m,
                "wake_shielding": self.wake_shielding,
                "max_transport_cfl": self.max_transport_cfl,
                "max_transport_substeps": self.max_transport_substeps,
            },
            **arrays,
            "time_s": self.time_s,
            "advance_total_m": self.advance_total_m,
            "emitted_total": self.emitted_total,
            "escaped_total": self.escaped_total,
            "recovered_total": self.recovered_total,
            "wake_discarded_total": self.wake_discarded_total,
        }

    @classmethod
    def from_state_dict(
        cls,
        payload: dict[str, Any],
        kernel: (
            SignedShieldingKernelV100514
            | SignedShieldingKernelFamilyV1005141
            | None
        ) = None,
    ) -> "PersistentSiteSignedMPZStateV100514":
        if payload.get("schema") not in {
            "persistent_site_signed_mpz_v10_0_5_14",
            "persistent_site_signed_mpz_v10_0_5_14_1",
        }:
            raise ValueError("unsupported persistent-site restart schema")
        candidate = PersistentSiteRowV100514(
            **dict(payload["candidate"])
        ).validate()
        if kernel is None:
            source = str(payload.get("kernel_source", ""))
            if not source:
                raise ValueError("restart requires a signed shielding artifact")
            kernel = load_signed_shielding_artifact_v1005141(source)
        config = dict(payload.get("state_config", {}))
        required = ("G_Pa", "nu", "b_m", "r0_m")
        missing = [name for name in required if name not in config]
        if missing:
            raise ValueError(f"restart state_config lacks {missing}")
        obj = cls(
            candidate,
            kernel,
            G_Pa=float(config["G_Pa"]),
            nu=float(config["nu"]),
            b_m=float(config["b_m"]),
            r0_m=float(config["r0_m"]),
            blunting_length_m=float(
                config.get("blunting_length_m", 0.5e-6)
            ),
            wake_shielding=bool(config.get("wake_shielding", False)),
            max_transport_cfl=float(config.get("max_transport_cfl", 0.35)),
            max_transport_substeps=int(
                config.get("max_transport_substeps", 2000)
            ),
        )
        for name in (
            "mobile_positive",
            "mobile_negative",
            "retained_positive",
            "retained_negative",
            "accumulated_slip_positive",
            "accumulated_slip_negative",
            "wake_mobile_positive",
            "wake_mobile_negative",
            "wake_retained_positive",
            "wake_retained_negative",
            "wake_slip_positive",
            "wake_slip_negative",
        ):
            if name not in payload:
                raise ValueError(f"restart payload lacks {name}")
            setattr(obj, name, np.asarray(payload[name], dtype=float))
        for name in (
            "time_s",
            "advance_total_m",
            "emitted_total",
            "escaped_total",
            "recovered_total",
            "wake_discarded_total",
        ):
            setattr(obj, name, float(payload.get(name, 0.0)))
        obj._validate_state_arrays()
        return obj


__all__ = ["PersistentSiteSignedIOMixin"]
