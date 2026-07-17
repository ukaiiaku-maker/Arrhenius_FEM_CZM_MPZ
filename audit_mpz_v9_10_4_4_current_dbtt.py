#!/usr/bin/env python3
"""DBTT audit wrapper with terminal state and timestep-limiter diagnostics."""
from __future__ import annotations

import json
from typing import Any

import audit_mpz_v9_10_4_current_dbtt as audit
import audit_mpz_v9_10_4_2_current_dbtt as transition_guard
from arrhenius_fracture.reduced_campaign_front_v9104 import (
    simulate_reduced_response as _simulate_reduced_response,
)

_RESULTS: list[dict[str, Any]] = []
_ORIGINAL_WRITE_PARTIAL = audit.write_partial_outputs


def _recording_simulate(parameters, T_K, settings, *, mode="full"):
    result = dict(
        _simulate_reduced_response(parameters, T_K, settings, mode=mode)
    )
    _RESULTS.append({"T_K": float(T_K), "mode": str(mode), "result": result})
    counts = result.get("dt_limiter_counts", {})
    print(
        "[terminal] "
        f"T={float(T_K):g}K mode={mode} "
        f"reason={result.get('termination_reason')} "
        f"K={result.get('terminal_K_MPa_sqrt_m')} "
        f"B={result.get('terminal_B')} a_um={result.get('terminal_a_um')} "
        f"available={result.get('terminal_available_sites')} "
        f"mobile={result.get('terminal_mobile_count')} "
        f"retained={result.get('terminal_retained_count')} "
        f"limiter={result.get('dominant_dt_limiter')} "
        f"min_dt={result.get('minimum_selected_dt_s')} "
        f"limiter_counts={json.dumps(counts, sort_keys=True)}",
        flush=True,
    )
    return result


def _diagnostic_fields(prefix: str, result: dict[str, Any]) -> dict[str, Any]:
    names = (
        "termination_reason",
        "terminal_K_MPa_sqrt_m",
        "terminal_time_s",
        "terminal_B",
        "terminal_a_um",
        "terminal_available_sites",
        "terminal_mobile_count",
        "terminal_retained_count",
        "terminal_slip_count",
        "terminal_K_shield_MPa_sqrt_m",
        "terminal_sigma_back_max_Pa",
        "terminal_lambda_c_s",
        "terminal_lambda_e_max_s",
        "dominant_dt_limiter",
        "minimum_selected_dt_s",
        "dt_selected_s",
        "dt_limiter",
    )
    fields = {f"{prefix}_{name}": result.get(name) for name in names}
    fields[f"{prefix}_dt_limiter_counts_json"] = json.dumps(
        result.get("dt_limiter_counts", {}), sort_keys=True
    )
    for key, value in result.items():
        if key.startswith("dt_limit_"):
            fields[f"{prefix}_{key}"] = value
    return fields


def _write_partial_with_diagnostics(out, rows, event_rows, summaries):
    # Two updated reduced-front calls (full and plasticity_off) are recorded for
    # each completed temperature row. Legacy zero-D calls do not enter _RESULTS.
    for index, row in enumerate(rows):
        full_index = 2 * index
        off_index = full_index + 1
        if off_index >= len(_RESULTS):
            break
        full = _RESULTS[full_index]
        off = _RESULTS[off_index]
        if full["mode"] != "full" or off["mode"] != "plasticity_off":
            raise RuntimeError(
                f"unexpected reduced-run order at row {index}: "
                f"{full['mode']!r}, {off['mode']!r}"
            )
        row.update(_diagnostic_fields("new_full", full["result"]))
        row.update(_diagnostic_fields("new_off", off["result"]))
    _ORIGINAL_WRITE_PARTIAL(out, rows, event_rows, summaries)


def main() -> None:
    _RESULTS.clear()
    audit.simulate_reduced_response = _recording_simulate
    audit.write_partial_outputs = _write_partial_with_diagnostics
    transition_guard.main()


if __name__ == "__main__":
    main()
