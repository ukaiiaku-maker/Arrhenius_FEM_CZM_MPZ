"""First-class material history container for adaptive crack geometry.

The current production driver still stores arrays directly for backward
compatibility.  This object provides the migration boundary so remeshing and
future persistent background-state carriers can be introduced without changing
hazard, plasticity, or fatigue kernels.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class MaterialState:
    ep_gp: np.ndarray
    rho_gp: np.ndarray
    pz_store_gp: Optional[np.ndarray] = None
    pz_mobile_gp: Optional[np.ndarray] = None
    pz_escape_gp: Optional[np.ndarray] = None
    pz_emit_gp: Optional[np.ndarray] = None

    def validate(self, ne: int) -> None:
        if self.ep_gp.shape != (3, ne):
            raise ValueError(f"ep_gp shape {self.ep_gp.shape} != (3,{ne})")
        if self.rho_gp.shape != (ne,):
            raise ValueError(f"rho_gp shape {self.rho_gp.shape} != ({ne},)")
        for name in ("pz_store_gp", "pz_mobile_gp", "pz_escape_gp", "pz_emit_gp"):
            arr = getattr(self, name)
            if arr is not None and arr.shape != (ne,):
                raise ValueError(f"{name} shape {arr.shape} != ({ne},)")
