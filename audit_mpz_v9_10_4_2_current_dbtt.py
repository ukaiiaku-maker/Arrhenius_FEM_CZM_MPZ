#!/usr/bin/env python3
"""Crash-safe wrapper for the v9.10.4 historical DBTT audit.

The transition feature extractor intentionally returns a compact invalid result
when one or more temperatures have no finite toughness.  The v9.10.4/v9.10.4.1
audit assumed a valid transition unconditionally and raised ``KeyError`` after
all solver calls had finished.  This wrapper preserves that invalid status,
adds NaN placeholders for optional feature fields, prints a clear diagnostic,
and then delegates to the existing progress/checkpointing audit.

No solver, constitutive, timestep, or transition-loss equation is changed.
"""
from __future__ import annotations

from typing import Any

import numpy as np

import audit_mpz_v9_10_4_current_dbtt as audit
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    best_adjacent_transition as _best_adjacent_transition,
)


_OPTIONAL_TRANSITION_FIELDS = (
    "split_index",
    "transition_low_K",
    "transition_high_K",
    "low_shelf",
    "high_shelf",
    "shelf_ratio",
    "robust_shelf_ratio",
    "main_jump",
    "jump_concentration",
    "low_span_fraction",
    "high_span_fraction",
    "secondary_jump_ratio",
    "plasticity_off_ratio",
)


def guarded_best_adjacent_transition(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Return a stable transition schema even when toughness is non-finite."""
    result = dict(_best_adjacent_transition(*args, **kwargs))
    if bool(result.get("valid", False)):
        return result

    mode = "new_full_vs_off" if kwargs.get("plasticity_off_toughness") is not None else "legacy"
    reason = str(result.get("reason", "invalid_transition"))
    print(
        f"[transition-invalid] mode={mode} reason={reason} "
        "feature metrics will be written as NaN",
        flush=True,
    )
    for key in _OPTIONAL_TRANSITION_FIELDS:
        result.setdefault(key, np.nan)
    result.setdefault("penalties", {})
    return result


def transition_summary_fields(prefix: str, result: dict[str, Any]) -> dict[str, Any]:
    """Expose the stable fields used by tests and downstream inspection."""
    valid = bool(result.get("valid", False))
    return {
        f"{prefix}_transition_valid": valid,
        f"{prefix}_transition_reason": "" if valid else str(result.get("reason", "invalid_transition")),
        f"{prefix}_transition_loss": float(result.get("loss", np.nan)),
        f"{prefix}_shelf_ratio": float(result.get("shelf_ratio", np.nan)),
        f"{prefix}_jump_concentration": float(result.get("jump_concentration", np.nan)),
        f"{prefix}_transition_low_K": float(result.get("transition_low_K", np.nan)),
        f"{prefix}_transition_high_K": float(result.get("transition_high_K", np.nan)),
        f"{prefix}_plasticity_off_ratio": float(result.get("plasticity_off_ratio", np.nan)),
    }


def main() -> None:
    # The existing audit imported the transition function directly.  Replace
    # that module-local reference before delegating to its checkpointed main.
    audit.best_adjacent_transition = guarded_best_adjacent_transition
    audit.main()


if __name__ == "__main__":
    main()
