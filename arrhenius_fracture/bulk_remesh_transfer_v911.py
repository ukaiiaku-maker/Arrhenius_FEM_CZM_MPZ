"""Runtime remesh-transfer completion for the explicit v9.11 bulk state.

The full solver already transfers its scalar integration-point density field when
adaptive CZM topology changes the element inventory. The v9.11 controller keeps
an additional mobile field. This patch treats the transferred scalar field as
the retained/forest state and transfers the mobile state with the pre-remesh
global mobile-to-retained ratio. It preserves the independently calibrated PT
kinetics and avoids silently reinitializing carriers.

This is a first scalar-state transfer gate. A later slip-system-resolved model
should transfer each system with element-overlap weights.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def install_bulk_remesh_transfer_patch() -> None:
    from .bulk_state_v911 import BulkPlasticityControllerV911

    cls = BulkPlasticityControllerV911
    if bool(getattr(cls, "_v911_remesh_transfer_installed", False)):
        return

    original_ensure = cls._ensure_state
    original_summary = cls.summary

    def _ensure_state_with_transfer(self, rho_fallback, T_K):
        rho = np.maximum(np.asarray(rho_fallback, dtype=float).reshape(-1), 0.0)
        old_mobile = None if self.mobile_rho_m2 is None else self.mobile_rho_m2.copy()
        old_retained = None if self.retained_rho_m2 is None else self.retained_rho_m2.copy()
        try:
            return original_ensure(self, rho, T_K)
        except RuntimeError as exc:
            text = str(exc)
            if "change in integration-point count" not in text:
                raise
            if old_mobile is None or old_retained is None:
                raise

            retained_sum = float(np.sum(np.maximum(old_retained, 0.0)))
            mobile_sum = float(np.sum(np.maximum(old_mobile, 0.0)))
            if retained_sum > 0.0:
                mobile_to_retained = mobile_sum / retained_sum
            else:
                total = mobile_sum + retained_sum
                mobile_to_retained = 0.0 if total <= 0.0 else mobile_sum / max(total, 1.0e-300)

            # rho is the full solver's already-transferred scalar forest field.
            # Reconstruct the unresolved mobile field without inventing a fixed
            # material fraction: use the actual pre-remesh state ratio.
            self.retained_rho_m2 = rho.copy()
            self.mobile_rho_m2 = np.maximum(mobile_to_retained * rho, 0.0)
            self.mesh_change_rejected = False
            self.mesh_transfer_count = int(getattr(self, "mesh_transfer_count", 0)) + 1
            self.last_mesh_transfer_old_size = int(old_retained.size)
            self.last_mesh_transfer_new_size = int(rho.size)
            self.last_mesh_transfer_mobile_to_retained = float(mobile_to_retained)
            self.bulk_remesh_transfer_status = (
                "solver_transferred_retained_field_plus_pre_remesh_mobile_ratio"
            )
            return None

    def _summary_with_transfer(self) -> dict[str, Any]:
        out = original_summary(self)
        out.update({
            "bulk_mesh_change_rejected": bool(self.mesh_change_rejected),
            "bulk_remesh_transfer_count": int(getattr(self, "mesh_transfer_count", 0)),
            "bulk_remesh_transfer_old_size": getattr(self, "last_mesh_transfer_old_size", None),
            "bulk_remesh_transfer_new_size": getattr(self, "last_mesh_transfer_new_size", None),
            "bulk_remesh_mobile_to_retained_ratio": getattr(
                self, "last_mesh_transfer_mobile_to_retained", None
            ),
            "bulk_remesh_transfer_status": (
                "not_applicable_bulk_elastic"
                if self.mode == "tip_only"
                else getattr(
                    self,
                    "bulk_remesh_transfer_status",
                    "not_exercised_solver_scalar_transfer_available",
                )
            ),
        })
        return out

    cls._ensure_state = _ensure_state_with_transfer
    cls.summary = _summary_with_transfer
    cls._v911_remesh_transfer_installed = True


__all__ = ["install_bulk_remesh_transfer_patch"]
