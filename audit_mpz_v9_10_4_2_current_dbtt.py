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

import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

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
_RECORDED_TRANSITIONS: list[tuple[str, dict[str, Any]]] = []


def guarded_best_adjacent_transition(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Return a stable transition schema even when toughness is non-finite."""
    result = dict(_best_adjacent_transition(*args, **kwargs))
    mode = "new" if kwargs.get("plasticity_off_toughness") is not None else "old"
    if not bool(result.get("valid", False)):
        reason = str(result.get("reason", "invalid_transition"))
        print(
            f"[transition-invalid] mode={mode} reason={reason} "
            "feature metrics will be written as NaN",
            flush=True,
        )
        for key in _OPTIONAL_TRANSITION_FIELDS:
            result.setdefault(key, np.nan)
        result.setdefault("penalties", {})
    _RECORDED_TRANSITIONS.append((mode, dict(result)))
    return result


def transition_summary_fields(prefix: str, result: dict[str, Any]) -> dict[str, Any]:
    """Return stable validity, reason, and feature fields for one transition."""
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


def _argument_path(flag: str) -> Path | None:
    try:
        index = sys.argv.index(flag)
    except ValueError:
        return None
    if index + 1 >= len(sys.argv):
        return None
    return Path(sys.argv[index + 1]).resolve()


def _augment_final_summary(out: Path) -> None:
    summary_path = out / "current_dbtt_old_vs_new_summary.csv"
    if not summary_path.exists():
        return
    summary = pd.read_csv(summary_path)
    expected = 2 * len(summary)
    if len(_RECORDED_TRANSITIONS) != expected:
        print(
            f"[transition-summary-warning] recorded={len(_RECORDED_TRANSITIONS)} "
            f"expected={expected}; validity columns were not added",
            flush=True,
        )
        return

    augmented_rows: list[dict[str, Any]] = []
    for candidate_index, (_, row) in enumerate(summary.iterrows()):
        old_mode, old_result = _RECORDED_TRANSITIONS[2 * candidate_index]
        new_mode, new_result = _RECORDED_TRANSITIONS[2 * candidate_index + 1]
        if old_mode != "old" or new_mode != "new":
            raise RuntimeError(
                f"unexpected transition recording order: {old_mode!r}, {new_mode!r}"
            )
        data = row.to_dict()
        data.update(transition_summary_fields("old", old_result))
        data.update(transition_summary_fields("new", new_result))
        augmented_rows.append(data)

    pd.DataFrame(augmented_rows).to_csv(summary_path, index=False)
    invalid_new = sum(not bool(row["new_transition_valid"]) for row in augmented_rows)

    report_path = out / "current_dbtt_audit.json"
    if report_path.exists():
        report = json.loads(report_path.read_text())
        report.update(
            {
                "transition_guard_version": "v9.10.4.2",
                "n_invalid_old_transitions": sum(
                    not bool(row["old_transition_valid"]) for row in augmented_rows
                ),
                "n_invalid_new_transitions": invalid_new,
            }
        )
        report_path.write_text(json.dumps(report, indent=2))

    print(
        f"[transition-summary] wrote validity/reason fields; "
        f"invalid_new_transitions={invalid_new}",
        flush=True,
    )


def main() -> None:
    # The existing audit imported the transition function directly. Replace
    # that module-local reference before delegating to its checkpointed main.
    _RECORDED_TRANSITIONS.clear()
    audit.best_adjacent_transition = guarded_best_adjacent_transition
    audit.main()
    out = _argument_path("--out")
    if out is not None:
        _augment_final_summary(out)


if __name__ == "__main__":
    main()
