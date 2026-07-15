"""Reproducible stochastic event statistics for the v9.11 Arrhenius transfer.

Arrhenius rates define integrated hazards. A physical first-passage event occurs
when the integrated hazard crosses an exponentially distributed unit-mean action
threshold. A deterministic threshold of one is retained as an explicit
regression/mean-field ablation.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

import numpy as np

VALID_EVENT_STATISTICS = ("deterministic", "stochastic")


def normalize_event_statistics(value: str) -> str:
    key = str(value).strip().lower().replace("-", "_")
    aliases = {
        "mean": "deterministic",
        "mean_field": "deterministic",
        "fixed": "deterministic",
        "random": "stochastic",
        "poisson": "stochastic",
        "gumbel": "stochastic",
        "first_passage": "stochastic",
    }
    key = aliases.get(key, key)
    if key not in VALID_EVENT_STATISTICS:
        raise ValueError(
            f"unknown event statistics {value!r}; expected one of "
            f"{', '.join(VALID_EVENT_STATISTICS)}"
        )
    return key


def make_rng(seed: int, stream: int) -> np.random.Generator:
    """Build a stable independent PCG64 stream from integer identifiers."""
    ss = np.random.SeedSequence([int(seed), int(stream)])
    return np.random.default_rng(ss)


@dataclass
class HazardThresholdSnapshot:
    target: float
    event_index: int
    rng_state: dict[str, Any]


class HazardThresholdStream:
    """Sequence of renewal thresholds in integrated-hazard coordinates.

    For ``stochastic`` statistics each threshold is Exp(1). This is the exact
    time-change construction of a non-homogeneous Poisson first-passage process.
    Under a monotonic stress/load ramp, transforming these action thresholds
    through the Arrhenius rate naturally produces extreme-value/Gumbel-like
    failure-load statistics; no ad-hoc scatter in barrier energy is added.
    """

    def __init__(self, mode: str = "deterministic", seed: int = 1, stream: int = 0):
        self.mode = normalize_event_statistics(mode)
        self.seed = int(seed)
        self.stream = int(stream)
        self.rng = make_rng(self.seed, self.stream)
        self.event_index = 0
        self.target = self._draw()

    def _draw(self) -> float:
        if self.mode == "deterministic":
            return 1.0
        return max(float(self.rng.exponential(1.0)), 1.0e-14)

    def consume(self, accumulated_action: float, max_events: int = 1) -> tuple[float, int, list[float]]:
        """Consume crossed thresholds and return residual action and event count."""
        action = max(float(accumulated_action), 0.0)
        limit = max(int(max_events), 0)
        crossed: list[float] = []
        while limit > 0 and action + 1.0e-15 >= self.target:
            used = float(self.target)
            action = max(action - used, 0.0)
            crossed.append(used)
            self.event_index += 1
            self.target = self._draw()
            limit -= 1
        return action, len(crossed), crossed

    def snapshot(self) -> HazardThresholdSnapshot:
        return HazardThresholdSnapshot(
            target=float(self.target),
            event_index=int(self.event_index),
            rng_state=copy.deepcopy(self.rng.bit_generator.state),
        )

    def restore(self, snapshot: HazardThresholdSnapshot) -> None:
        self.target = float(snapshot.target)
        self.event_index = int(snapshot.event_index)
        self.rng.bit_generator.state = copy.deepcopy(snapshot.rng_state)

    def fork(self, child_stream_offset: int = 1) -> "HazardThresholdStream":
        """Create an independent child renewal stream for a new branch."""
        seeds = self.rng.integers(0, np.iinfo(np.uint64).max, size=2, dtype=np.uint64)
        self.seed = int(seeds[0])
        self.stream += 2 * int(child_stream_offset)
        self.rng = make_rng(self.seed, self.stream)
        self.event_index = 0
        self.target = self._draw()
        return HazardThresholdStream(
            self.mode,
            seed=int(seeds[1]),
            stream=self.stream + 1,
        )

    def state_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "seed": self.seed,
            "stream": self.stream,
            "event_index": self.event_index,
            "target": self.target,
            "rng_state": copy.deepcopy(self.rng.bit_generator.state),
        }

    @classmethod
    def from_state_dict(cls, payload: dict[str, Any]) -> "HazardThresholdStream":
        obj = cls(
            payload.get("mode", "deterministic"),
            int(payload.get("seed", 1)),
            int(payload.get("stream", 0)),
        )
        obj.event_index = int(payload.get("event_index", 0))
        obj.target = float(payload.get("target", obj.target))
        state = payload.get("rng_state")
        if state is not None:
            obj.rng.bit_generator.state = copy.deepcopy(state)
        return obj


def sample_effective_binomial(
    rng: np.random.Generator,
    capacity: float,
    probability: float,
) -> float:
    """Binomial sample for a possibly non-integer effective site capacity."""
    cap = max(float(capacity), 0.0)
    p = float(np.clip(probability, 0.0, 1.0))
    n = int(np.floor(cap))
    frac = cap - n
    out = float(rng.binomial(n, p)) if n > 0 else 0.0
    if frac > 0.0 and float(rng.random()) < p:
        out += frac
    return min(out, cap)


__all__ = [
    "VALID_EVENT_STATISTICS",
    "HazardThresholdSnapshot",
    "HazardThresholdStream",
    "make_rng",
    "normalize_event_statistics",
    "sample_effective_binomial",
]
