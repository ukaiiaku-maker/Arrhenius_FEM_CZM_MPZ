"""PF v10.2.22 persistent-site signed moving-process-zone state for FEM/CZM.

The implementation is split into auditable geometry/emission, transport, and
restart/moving-frame mixins; this module exposes the production state class.
"""
from __future__ import annotations

import copy
from typing import Any

import numpy as np

from .persistent_site_registry_v100514 import PersistentSiteRowV100514
from . import persistent_site_signed_support_v100514 as _support
from .persistent_site_complementarity_v100514 import (
    solve_backstress_limited_activations,
)

# Install the audited endpoint treatment before the core mixin imports the
# support symbol. This mirrors the PF v10.2.22 audited wrapper while keeping the
# correction explicit in the versioned FEM/CZM state entry.
_support.solve_backstress_limited_activations = solve_backstress_limited_activations

from .persistent_site_signed_support_v100514 import (  # noqa: E402
    KERNEL_SCHEMA,
    MODEL_ID,
    SignedShieldingKernelV100514,
    effective_front_width_m,
    persistent_site_multiplicity,
)
from .signed_kernel_family_v1005141 import (  # noqa: E402
    SignedShieldingKernelFamilyV1005141,
)
from .persistent_site_signed_core_v100514 import (  # noqa: E402
    PersistentSiteSignedCoreMixin,
)
from .persistent_site_signed_io_v100514 import (  # noqa: E402
    PersistentSiteSignedIOMixin,
)
from .persistent_site_signed_transport_v100514 import (  # noqa: E402
    PersistentSiteSignedTransportMixin,
)


class PersistentSiteSignedMPZStateV100514(
    PersistentSiteSignedIOMixin,
    PersistentSiteSignedTransportMixin,
    PersistentSiteSignedCoreMixin,
):
    """Signed active/wake MPZ state with persistent, nondepleting sources."""

    state_model = MODEL_ID

    def __init__(
        self,
        candidate: PersistentSiteRowV100514,
        kernel: SignedShieldingKernelV100514 | SignedShieldingKernelFamilyV1005141,
        *,
        G_Pa: float,
        nu: float,
        b_m: float,
        r0_m: float,
        blunting_length_m: float = 0.5e-6,
        wake_shielding: bool = False,
        max_transport_cfl: float = 0.35,
        max_transport_substeps: int = 2000,
    ) -> None:
        self.candidate = copy.deepcopy(candidate.validate())
        self.n_systems = int(candidate.n_slip_channels)
        self.n_bins = int(candidate.n_bins_recommended)
        self.length_m = float(candidate.L_pz_um_recommended) * 1.0e-6
        self.dx = self.length_m / self.n_bins
        self.x = (np.arange(self.n_bins, dtype=float) + 0.5) * self.dx
        self.wake_n_bins = self.n_bins
        self.wake_dx = self.dx
        self.wake_x = (
            np.arange(self.wake_n_bins, dtype=float) + 0.5
        ) * self.wake_dx
        self.G_Pa = float(G_Pa)
        self.nu = float(nu)
        self.b_m = abs(float(b_m))
        self.r0_m = max(float(r0_m), self.b_m, 1.0e-12)
        self.blunting_length_m = max(float(blunting_length_m), self.b_m)
        self.wake_shielding = bool(wake_shielding)
        self.max_transport_cfl = max(float(max_transport_cfl), 1.0e-6)
        self.max_transport_substeps = max(int(max_transport_substeps), 1)

        self.kernel_family: SignedShieldingKernelFamilyV1005141 | None = None
        if isinstance(kernel, SignedShieldingKernelFamilyV1005141):
            self.kernel_family = copy.deepcopy(kernel)
            self.kernel_family.validate(
                self.n_systems,
                self.n_bins,
                active_x_m=self.x,
                wake_x_m=self.wake_x,
            )
            self.kernel = self.kernel_family.snapshot(0.0, self.x, self.wake_x)
            self.activation_to_line_content_by_system = (
                self.kernel_family.activation_to_line_content_by_system.copy()
            )
        else:
            self.kernel = copy.deepcopy(kernel)
            self.kernel.validate(self.n_systems, self.n_bins)
            if self.kernel.active_x_m is not None and not np.allclose(
                np.asarray(self.kernel.active_x_m, dtype=float),
                self.x,
                rtol=1.0e-12,
                atol=1.0e-18,
            ):
                raise ValueError(
                    "signed shielding kernel coordinates do not match the MPZ grid"
                )
            self.activation_to_line_content_by_system = (
                self.kernel.activation_to_line_content_by_system.copy()
            )
        if self.activation_to_line_content_by_system.shape != (self.n_systems,):
            raise ValueError("activation-to-line conversion system count mismatch")
        if np.any(self.activation_to_line_content_by_system <= 0.0):
            raise ValueError("activation-to-line conversion must be positive")

        shape = (self.n_systems, self.n_bins)
        wake_shape = (self.n_systems, self.wake_n_bins)
        for name in (
            "mobile_positive",
            "mobile_negative",
            "retained_positive",
            "retained_negative",
            "accumulated_slip_positive",
            "accumulated_slip_negative",
        ):
            setattr(self, name, np.zeros(shape, dtype=float))
        for name in (
            "wake_mobile_positive",
            "wake_mobile_negative",
            "wake_retained_positive",
            "wake_retained_negative",
            "wake_slip_positive",
            "wake_slip_negative",
        ):
            setattr(self, name, np.zeros(wake_shape, dtype=float))
        self.site_capacity = np.ones(self.n_systems, dtype=float)
        self.available_sites = self.site_capacity.copy()
        self.tip_source_activity = np.ones(self.n_systems, dtype=float)
        self.reference_area_m2 = (
            float(candidate.reference_source_area_um2) * 1.0e-12
        )
        self.reference_width_m = (
            float(candidate.reference_front_width_um) * 1.0e-6
        )
        self.reference_density_m2 = float(candidate.rho_forest_floor_m2)
        self.source_zone_length_m = (
            float(candidate.source_zone_length_um) * 1.0e-6
        )
        self.minimum_front_width_m = (
            float(candidate.minimum_front_width_um) * 1.0e-6
        )
        self.maximum_front_width_m = self.length_m
        self.active_arc_factor = self.reference_area_m2 / (
            self.r0_m * self.reference_width_m
        )
        self.time_s = 0.0
        self.advance_total_m = 0.0
        self.emitted_total = 0.0
        self.escaped_total = 0.0
        self.recovered_total = 0.0
        self.wake_discarded_total = 0.0
        self.last_emission: dict[str, Any] = {}
        self.last_transport: dict[str, Any] = {}
        self.last_advance: dict[str, Any] = {}

    def current_kernel_snapshot(self) -> SignedShieldingKernelV100514:
        """Return the mode-I signed kernel at the committed crack-path extension."""
        if self.kernel_family is None:
            return self.kernel
        self.kernel = self.kernel_family.snapshot(
            self.advance_total_m,
            self.x,
            self.wake_x,
        )
        return self.kernel

    def advance(self, distance_m: float) -> dict[str, float]:
        """Preflight the family envelope before mutating the moving-frame state."""
        distance = max(float(distance_m), 0.0)
        if self.kernel_family is not None:
            self.kernel_family.snapshot(
                self.advance_total_m + distance,
                self.x,
                self.wake_x,
            )
        return super().advance(distance)

    def split(self, daughter_fraction: float) -> "PersistentSiteSignedMPZStateV100514":
        """A daughter inherits the parent's cumulative crack-path coordinate."""
        parent_extension = self.advance_total_m
        child = super().split(daughter_fraction)
        child.advance_total_m = parent_extension
        child.current_kernel_snapshot()
        return child

    @property
    def kernel_source_path(self) -> str:
        if self.kernel_family is not None:
            return self.kernel_family.source_path
        return self.kernel.source_path

    def kernel_artifact_audit(self) -> dict[str, Any]:
        if self.kernel_family is not None:
            payload = self.kernel_family.audit_payload()
            payload["current_cumulative_crack_path_extension_m"] = (
                self.advance_total_m
            )
            payload["current_state_weights"] = self.current_kernel_snapshot().metadata.get(
                "state_weights", {}
            )
            return payload
        return {
            "schema": self.kernel.metadata.get("schema", KERNEL_SCHEMA),
            "artifact_kind": "fixed_kernel",
            "source_path": self.kernel.source_path,
            "active_shape": list(
                self.kernel.active_kernel_Pa_sqrt_m_per_signed_line.shape
            ),
            "wake_shape": list(
                self.kernel.wake_kernel_Pa_sqrt_m_per_signed_line.shape
            ),
            "activation_to_line_content_by_system": (
                self.activation_to_line_content_by_system.tolist()
            ),
        }


__all__ = [
    "MODEL_ID",
    "KERNEL_SCHEMA",
    "SignedShieldingKernelV100514",
    "SignedShieldingKernelFamilyV1005141",
    "PersistentSiteSignedMPZStateV100514",
    "effective_front_width_m",
    "persistent_site_multiplicity",
    "solve_backstress_limited_activations",
]
