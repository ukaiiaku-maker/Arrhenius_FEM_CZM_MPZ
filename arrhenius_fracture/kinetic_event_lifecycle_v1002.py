"""Transactional physical-time lifecycle for progressive kinetic CZM events.

The controller is geometry agnostic.  The production driver supplies callbacks
that create the next trial interface, advance its kinetic/cohesive state, and
update the accepted crack geometry after a commit.  Rejected damage increments
consume no physical time.  Unused time after a commit remains at the same
applied load and is offered to a newly equilibrated trial event.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Callable


@dataclass
class EventLifecycleConfig:
    min_retry_dt_s: float = 1.0e-18
    max_retries_per_substep: int = 64
    max_accepted_substeps_per_interval: int = 10000
    time_tolerance_s: float = 0.0

    def validate(self) -> "EventLifecycleConfig":
        if self.min_retry_dt_s <= 0.0:
            raise ValueError("min_retry_dt_s must be positive")
        if self.max_retries_per_substep < 1:
            raise ValueError("max_retries_per_substep must be positive")
        if self.max_accepted_substeps_per_interval < 1:
            raise ValueError("max_accepted_substeps_per_interval must be positive")
        if self.time_tolerance_s < 0.0:
            raise ValueError("time_tolerance_s must be nonnegative")
        return self


@dataclass
class AcceptedLifecycleStep:
    context: Any
    result: Any
    requested_dt_s: float
    retry_count: int


@dataclass
class EventLifecycleResult:
    requested_dt_s: float
    consumed_dt_s: float
    unused_dt_s: float
    committed_events: int
    rejected_attempts: int
    accepted_steps: list[AcceptedLifecycleStep] = field(default_factory=list)
    stopped_at_target: bool = False


class KineticEventLifecycleController:
    """Consume one outer physical-time interval without losing event time."""

    def __init__(self, config: EventLifecycleConfig | None = None):
        self.config = (config or EventLifecycleConfig()).validate()

    @staticmethod
    def _accepted(result: Any) -> bool:
        return bool(getattr(result, "accepted"))

    @staticmethod
    def _committed(result: Any) -> bool:
        return bool(getattr(result, "committed"))

    @staticmethod
    def _consumed(result: Any) -> float:
        return float(getattr(result, "dt_consumed_s"))

    @staticmethod
    def _recommended(result: Any) -> float | None:
        value = getattr(result, "recommended_dt_s", None)
        if value is None:
            return None
        value = float(value)
        return value if math.isfinite(value) and value > 0.0 else None

    def consume_interval(
        self,
        *,
        total_dt_s: float,
        ensure_trial: Callable[[], Any],
        advance_trial: Callable[[float], Any],
        on_accepted: Callable[[Any, Any, float, int], None] | None = None,
        on_commit: Callable[[Any, Any], None] | None = None,
        target_reached: Callable[[], bool] | None = None,
    ) -> EventLifecycleResult:
        requested = max(float(total_dt_s), 0.0)
        remaining = requested
        rejected = 0
        committed = 0
        accepted_steps: list[AcceptedLifecycleStep] = []
        stopped_at_target = False
        ulp_scale = requested if requested > 0.0 else 1.0
        tol = max(
            float(self.config.time_tolerance_s),
            32.0 * math.ulp(ulp_scale),
        )

        while remaining > tol:
            if target_reached is not None and bool(target_reached()):
                stopped_at_target = True
                break
            if len(accepted_steps) >= int(
                self.config.max_accepted_substeps_per_interval
            ):
                raise RuntimeError(
                    "kinetic event lifecycle exceeded the accepted-substep "
                    "safety limit without consuming the outer interval"
                )

            context = ensure_trial()
            dt_try = remaining
            retry_count = 0

            while True:
                result = advance_trial(dt_try)
                if self._accepted(result):
                    break

                rejected += 1
                retry_count += 1
                if retry_count > int(self.config.max_retries_per_substep):
                    raise RuntimeError(
                        "kinetic event lifecycle exceeded the retry limit for "
                        "one transactional cohesive increment"
                    )

                recommended = self._recommended(result)
                if recommended is None or recommended >= dt_try * (
                    1.0 - 1.0e-12
                ):
                    recommended = 0.5 * dt_try
                dt_try = min(recommended, remaining)
                if dt_try < float(self.config.min_retry_dt_s):
                    raise RuntimeError(
                        "kinetic event lifecycle retry interval fell below "
                        f"min_retry_dt_s={self.config.min_retry_dt_s:.6e}"
                    )

            consumed = self._consumed(result)
            if not math.isfinite(consumed) or consumed <= 0.0:
                raise RuntimeError(
                    "accepted kinetic cohesive step consumed no finite physical time"
                )
            if consumed > dt_try + max(tol, 1.0e-12 * dt_try):
                raise RuntimeError(
                    "accepted kinetic cohesive step consumed more time than requested"
                )
            if (
                not self._committed(result)
                and consumed < dt_try - max(tol, 1.0e-12 * dt_try)
            ):
                raise RuntimeError(
                    "noncommitted kinetic cohesive step returned unexplained unused time"
                )

            remaining = max(remaining - consumed, 0.0)
            accepted = AcceptedLifecycleStep(
                context=context,
                result=result,
                requested_dt_s=dt_try,
                retry_count=retry_count,
            )
            accepted_steps.append(accepted)
            if on_accepted is not None:
                on_accepted(context, result, dt_try, retry_count)

            if self._committed(result):
                committed += 1
                if on_commit is not None:
                    on_commit(context, result)

        consumed_total = max(requested - remaining, 0.0)
        return EventLifecycleResult(
            requested_dt_s=requested,
            consumed_dt_s=consumed_total,
            unused_dt_s=remaining,
            committed_events=committed,
            rejected_attempts=rejected,
            accepted_steps=accepted_steps,
            stopped_at_target=stopped_at_target,
        )


__all__ = [
    "EventLifecycleConfig",
    "AcceptedLifecycleStep",
    "EventLifecycleResult",
    "KineticEventLifecycleController",
]
