"""Compact stochastic-identity fingerprint for Gate 1 transactionality audits.

The direct-step mechanical/plastic state-equality proof
(tests/test_pf_theta0_j_controlled_loading_real_engine_v10051840.py) shows
the FEM mechanics are transactional across rejected KJ-target controller
trials. It says nothing about the *stochastic* Arrhenius first-passage
state -- the hazard RNG streams, the currently-drawn cleavage threshold,
the materialized event-length factor, B, N_em, and the moving-tip MPZ
population -- which live on the front-engine object, not in the mesh
arrays, and are committed only by `front_engine.step(...)` after the
controller accepts a trial (see the "sharp_front.py architecture" notes
in CLAUDE_PROGRESS.md: front-engine state is untouched during trial
evaluation for the non-deflect path).

This module captures a SHORT, comparison-ready fingerprint of that state
without ever dumping a full RNG bit-generator state or a full MPZ array,
and without ever triggering a draw merely to audit it: any field that
does not exist on the given engine (e.g. the minimal `legacy_scalar`
front engine has no RNG, threshold, event-length factor, or MPZ state)
is reported as `NOT_MATERIALIZED`, never fabricated.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any

MODEL_ID = "pf_theta0_stochastic_fingerprint_v10_0_5_18_4_0"

NOT_MATERIALIZED = "not_materialized"


def _short_digest(payload: str) -> str:
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _rng_digest(rng: Any) -> str:
    if rng is None:
        return NOT_MATERIALIZED
    state = rng.bit_generator.state
    return _short_digest(json.dumps(state, sort_keys=True, default=str))


def _emission_rng_digest(rngs: Any) -> str:
    if not rngs:
        return NOT_MATERIALIZED
    payload = json.dumps(
        [r.bit_generator.state for r in rngs], sort_keys=True, default=str
    )
    return _short_digest(payload)


def _scalar_or_not_materialized(value: Any) -> float | str:
    if value is None:
        return NOT_MATERIALIZED
    return float(value)


@dataclass(frozen=True)
class MPZFingerprint:
    mobile_total: float | str
    retained_total: float | str
    emitted_total: float | str
    escaped_total: float | str
    recovered_total: float | str
    wake_mobile_total: float | str
    wake_retained_total: float | str
    digest: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "mobile_total": self.mobile_total,
            "retained_total": self.retained_total,
            "emitted_total": self.emitted_total,
            "escaped_total": self.escaped_total,
            "recovered_total": self.recovered_total,
            "wake_mobile_total": self.wake_mobile_total,
            "wake_retained_total": self.wake_retained_total,
            "digest": self.digest,
        }


def _capture_mpz(mpz_state: Any) -> MPZFingerprint:
    if mpz_state is None:
        return MPZFingerprint(
            mobile_total=NOT_MATERIALIZED,
            retained_total=NOT_MATERIALIZED,
            emitted_total=NOT_MATERIALIZED,
            escaped_total=NOT_MATERIALIZED,
            recovered_total=NOT_MATERIALIZED,
            wake_mobile_total=NOT_MATERIALIZED,
            wake_retained_total=NOT_MATERIALIZED,
            digest=NOT_MATERIALIZED,
        )
    import numpy as np

    mobile_total = float(np.sum(mpz_state.mobile))
    retained_total = float(np.sum(mpz_state.retained))
    emitted_total = float(mpz_state.emitted_total)
    escaped_total = float(mpz_state.escaped_total)
    recovered_total = float(mpz_state.recovered_total)
    wake_mobile_total = float(mpz_state.wake_mobile_total)
    wake_retained_total = float(mpz_state.wake_retained_total)
    payload = "|".join(
        f"{v:.17g}"
        for v in (
            mobile_total,
            retained_total,
            emitted_total,
            escaped_total,
            recovered_total,
            wake_mobile_total,
            wake_retained_total,
        )
    )
    return MPZFingerprint(
        mobile_total=mobile_total,
        retained_total=retained_total,
        emitted_total=emitted_total,
        escaped_total=escaped_total,
        recovered_total=recovered_total,
        wake_mobile_total=wake_mobile_total,
        wake_retained_total=wake_retained_total,
        digest=_short_digest(payload),
    )


@dataclass(frozen=True)
class StochasticFingerprint:
    hazard_rng_digest: str
    emission_rng_digest: str
    cleave_threshold_action: float | str
    stochastic_event_length_factor: float | str
    B: float | str
    N_em: float | str
    mpz: MPZFingerprint

    def matches(self, other: "StochasticFingerprint") -> bool:
        """Exact equality on every field: the fail-closed check that a
        rejected (or in-progress) controller trial left no residue."""
        return (
            self.hazard_rng_digest == other.hazard_rng_digest
            and self.emission_rng_digest == other.emission_rng_digest
            and self.cleave_threshold_action == other.cleave_threshold_action
            and self.stochastic_event_length_factor
            == other.stochastic_event_length_factor
            and self.B == other.B
            and self.N_em == other.N_em
            and self.mpz == other.mpz
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "hazard_rng_digest": self.hazard_rng_digest,
            "emission_rng_digest": self.emission_rng_digest,
            "cleave_threshold_action": self.cleave_threshold_action,
            "stochastic_event_length_factor": self.stochastic_event_length_factor,
            "B": self.B,
            "N_em": self.N_em,
            "mpz": self.mpz.as_dict(),
        }


def capture_stochastic_fingerprint(eng: Any) -> StochasticFingerprint:
    """Read-only capture of `eng`'s current stochastic/hazard/MPZ identity.

    Never triggers a draw or advances any state; only reads attributes
    that already exist on the engine. Attributes absent on this engine
    class (e.g. the minimal legacy_scalar front engine has none of the
    persistent-site stochastic attributes) are reported as
    NOT_MATERIALIZED rather than fabricated or defaulted to zero.
    """
    hazard_rng = getattr(eng, "_hazard_rng", None)
    emission_rngs = getattr(eng, "_emission_rngs", None)
    cleave_threshold = getattr(eng, "hazard_threshold_action", None)
    event_length_factor = getattr(eng, "stochastic_event_length_factor", None)
    B = getattr(eng, "B", None)
    N_em = getattr(eng, "N_em", None)
    mpz_state = getattr(eng, "mpz_state", None)

    return StochasticFingerprint(
        hazard_rng_digest=_rng_digest(hazard_rng),
        emission_rng_digest=_emission_rng_digest(emission_rngs),
        cleave_threshold_action=_scalar_or_not_materialized(cleave_threshold),
        stochastic_event_length_factor=_scalar_or_not_materialized(event_length_factor),
        B=_scalar_or_not_materialized(B),
        N_em=_scalar_or_not_materialized(N_em),
        mpz=_capture_mpz(mpz_state),
    )


__all__ = [
    "MODEL_ID",
    "NOT_MATERIALIZED",
    "MPZFingerprint",
    "StochasticFingerprint",
    "capture_stochastic_fingerprint",
]
