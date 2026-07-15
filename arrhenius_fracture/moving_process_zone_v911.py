"""v9.11 MPZ state with 2-D profiles and optional stochastic source events."""
from __future__ import annotations

import copy
import math
import os

import numpy as np

from .moving_process_zone_v9102 import MovingProcessZoneState as _V9102State
from .process_zone_2d_v911 import ProcessZone2DProfile
from .stochastic_kinetics_v911 import (
    make_rng,
    normalize_event_statistics,
    sample_effective_binomial,
)


class MovingProcessZoneState(_V9102State):
    """Independent-shape MPZ with optional 2-D forest/stress-shape inputs.

    The 2-D scalar density augments the local forest density used by Peierls/Taylor
    transport. It is not interpreted as signed shielding. The direct K shield
    remains the unresolved retained-line integral inherited from v9.10.2.

    In stochastic mode, finite source sites emit through Bernoulli/binomial
    realizations of the integrated Arrhenius hazard. Deterministic expected-value
    source depletion remains available for regression.
    """

    state_model = "moving_pz_v911_independent_shapes_2d_profile"

    def __init__(self, cfg):
        super().__init__(cfg)
        self._profile_2d: ProcessZone2DProfile | None = None
        self.event_statistics = normalize_event_statistics(
            getattr(
                cfg,
                "event_statistics",
                os.environ.get("ARRHENIUS_EVENT_STATISTICS", "deterministic"),
            )
        )
        self.stochastic_emission = bool(
            getattr(
                cfg,
                "stochastic_emission",
                os.environ.get("ARRHENIUS_STOCHASTIC_EMISSION", "1") != "0"
                and self.event_statistics == "stochastic",
            )
        )
        self.stochastic_seed = int(
            getattr(
                cfg,
                "stochastic_seed",
                os.environ.get("ARRHENIUS_STOCHASTIC_SEED", "1"),
            )
        )
        self.stochastic_stream = int(getattr(cfg, "stochastic_emission_stream", 17011))
        self._emission_rng = make_rng(self.stochastic_seed, self.stochastic_stream)
        self.stochastic_emission_events = 0

    def set_2d_profile(self, profile: ProcessZone2DProfile | None) -> None:
        self._profile_2d = profile

    @staticmethod
    def _resample(values, n):
        q = np.asarray(values, float).reshape(-1)
        if q.size == n:
            return q.copy()
        old = (np.arange(q.size, dtype=float) + 0.5) / max(q.size, 1)
        new = (np.arange(n, dtype=float) + 0.5) / n
        return np.interp(new, old, q)

    def local_forest_density_m2(self) -> np.ndarray:
        local = np.asarray(super().local_forest_density_m2(), float)
        if self._profile_2d is None:
            return local
        floor = float(self.cfg.pt_forest_density_floor_m2)
        bulk = self._resample(self._profile_2d.forest_density_m2, self.n_bins)
        retained_excess = np.maximum(local - floor, 0.0)
        return np.maximum(bulk, floor) + retained_excess

    def local_stress_profile_Pa(self, tip_stress_Pa: float) -> np.ndarray:
        if self._profile_2d is None or not self._profile_2d.reliable:
            return super().local_stress_profile_Pa(tip_stress_Pa)
        shape = self._resample(self._profile_2d.stress_shape, self.n_bins)
        return max(float(tip_stress_Pa), 0.0) * np.maximum(shape, 0.0)

    def _source_commit_from_hazard(
        self,
        hazard_integral: float,
        system_weights: np.ndarray | None = None,
    ) -> np.ndarray:
        if self.event_statistics != "stochastic" or not self.stochastic_emission:
            return super()._source_commit_from_hazard(hazard_integral, system_weights)

        H = max(float(hazard_integral), 0.0)
        p = 1.0 - math.exp(-min(H, 700.0))
        if system_weights is None:
            system_weights = np.ones(self.n_systems, dtype=float)
        w = np.maximum(np.asarray(system_weights, dtype=float).reshape(-1), 0.0)
        if w.size < self.n_systems:
            w = np.pad(w, (0, self.n_systems - w.size), mode="edge")
        w = w[: self.n_systems]
        if float(np.sum(w)) <= 0.0:
            w[:] = 1.0
        active_fraction = w / max(float(np.max(w)), 1.0e-300)

        emitted = np.zeros(self.n_systems, dtype=float)
        for i in range(self.n_systems):
            # Preserve the deterministic model's mean exactly.
            pi = float(np.clip(active_fraction[i] * p, 0.0, 1.0))
            emitted[i] = sample_effective_binomial(
                self._emission_rng,
                float(self.available_sites[i]),
                pi,
            )
        self.available_sites = np.maximum(self.available_sites - emitted, 0.0)
        nsrc = max(min(int(self.cfg.source_bin_count), self.n_bins), 1)
        self.mobile[:, :nsrc] += emitted[:, None] / nsrc
        self.accumulated_slip[:, :nsrc] += emitted[:, None] / nsrc
        self.emitted_total += float(np.sum(emitted))
        self.stochastic_emission_events += int(np.count_nonzero(emitted > 0.0))
        return emitted

    def split(self, daughter_fraction: float) -> "MovingProcessZoneState":
        child = super().split(daughter_fraction)
        if self.event_statistics == "stochastic":
            seeds = self._emission_rng.integers(
                0, np.iinfo(np.uint64).max, size=2, dtype=np.uint64
            )
            self._emission_rng = make_rng(int(seeds[0]), self.stochastic_stream + 2)
            child._emission_rng = make_rng(int(seeds[1]), self.stochastic_stream + 3)
            self.stochastic_stream += 2
            child.stochastic_stream = self.stochastic_stream + 1
        return child

    def state_dict(self):
        out = super().state_dict()
        out["stochastic_v911"] = {
            "event_statistics": self.event_statistics,
            "stochastic_emission": self.stochastic_emission,
            "stochastic_seed": self.stochastic_seed,
            "stochastic_stream": self.stochastic_stream,
            "stochastic_emission_events": self.stochastic_emission_events,
            "emission_rng_state": copy.deepcopy(self._emission_rng.bit_generator.state),
        }
        return out

    @classmethod
    def from_state_dict(cls, payload):
        obj = super().from_state_dict(payload)
        stochastic = dict(payload.get("stochastic_v911", {}))
        obj.event_statistics = normalize_event_statistics(
            stochastic.get("event_statistics", "deterministic")
        )
        obj.stochastic_emission = bool(stochastic.get("stochastic_emission", False))
        obj.stochastic_seed = int(stochastic.get("stochastic_seed", 1))
        obj.stochastic_stream = int(stochastic.get("stochastic_stream", 17011))
        obj.stochastic_emission_events = int(
            stochastic.get("stochastic_emission_events", 0)
        )
        obj._emission_rng = make_rng(obj.stochastic_seed, obj.stochastic_stream)
        state = stochastic.get("emission_rng_state")
        if state is not None:
            obj._emission_rng.bit_generator.state = copy.deepcopy(state)
        return obj

    def diagnostics(self, G_shear, nu, b, r0, c_blunt):
        out = super().diagnostics(G_shear, nu, b, r0, c_blunt)
        if self._profile_2d is None:
            out.update({
                "mpz_2d_profile_active": 0.0,
                "bulk_scalar_rho_used_for_signed_shielding": 0.0,
            })
        else:
            out.update(self._profile_2d.diagnostics())
            out["mpz_2d_profile_active"] = 1.0
            out["bulk_scalar_rho_used_for_signed_shielding"] = 0.0
        out.update({
            "mpz_shielding_source": "unresolved_retained_line_integral",
            "bulk_plastic_shielding_source": "FEM_stress_redistribution_already_in_J",
            "bulk_K_shield_subtracted_again": 0.0,
            "explicit_GND_backstress_active": 0.0,
            "event_statistics_stochastic": float(self.event_statistics == "stochastic"),
            "stochastic_emission_active": float(
                self.event_statistics == "stochastic" and self.stochastic_emission
            ),
            "stochastic_emission_events": float(self.stochastic_emission_events),
        })
        return out


__all__ = ["MovingProcessZoneState"]
