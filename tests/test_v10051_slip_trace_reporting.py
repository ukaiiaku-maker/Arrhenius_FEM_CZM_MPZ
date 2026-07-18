from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from arrhenius_fracture.slip_trace_reporting_v10051 import (
    CHANNEL_TABLE,
    NORMALIZED_AUDIT,
    NORMALIZED_RESULTS,
    normalize_output,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_source_output(root: Path, *, include_channels: bool = True) -> None:
    audit = {
        "implementation_certified": True,
        "capture_count": 4,
        "nonzero_emission_drive_capture_count": 4,
        "drive_factor_normalization_or_clipping_active": False,
        "directional_multiplier_applied_after_hazard": False,
        "fit_derived_shielding_cap_active": False,
        "latest": {"slip_system_names": ["(110)", "(1-10)"]},
    }
    (root / "parallel_opening_emission_v10_0_5_audit.json").write_text(
        json.dumps(audit, indent=2)
    )
    (root / "mode_i_v10_0_5_results.json").write_text(
        json.dumps([{"T_K": 700.0, "N_em_final": 0.0}], indent=2)
    )
    fields = ["step", "T_K", "anisotropic_KJ_Pa_sqrt_m", "B", "N_em"]
    if include_channels:
        for index in range(2):
            fields.extend([
                f"slip_drive_factor_{index}",
                f"slip_tau_signed_Pa_{index}",
                f"sigma_emission_effective_Pa_{index}",
                f"sigma_emission_backstress_Pa_{index}",
                f"lambda_emit_s-1_{index}",
                f"dN_emit_{index}",
            ])
    rows = []
    for step in range(2):
        row = {
            "step": step,
            "T_K": 700.0,
            "anisotropic_KJ_Pa_sqrt_m": 1.0e6 * (step + 1),
            "B": 0.1 * step,
            "N_em": 0.0,
        }
        if include_channels:
            for index in range(2):
                row.update({
                    f"slip_drive_factor_{index}": 0.05 * (index + 1),
                    f"slip_tau_signed_Pa_{index}": (-1.0) ** index * 1.0e7,
                    f"sigma_emission_effective_Pa_{index}": 1.0e8,
                    f"sigma_emission_backstress_Pa_{index}": 0.0,
                    f"lambda_emit_s-1_{index}": 0.0,
                    f"dN_emit_{index}": 0.0,
                })
        rows.append(row)
    with (root / "steps_0700K.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_zero_emission_is_valid_and_source_outputs_are_immutable(tmp_path):
    _write_source_output(tmp_path)
    source_names = [
        "parallel_opening_emission_v10_0_5_audit.json",
        "mode_i_v10_0_5_results.json",
        "steps_0700K.csv",
    ]
    before = {name: _sha256(tmp_path / name) for name in source_names}

    payload = normalize_output(tmp_path)

    after = {name: _sha256(tmp_path / name) for name in source_names}
    assert before == after
    assert payload["implementation_certified"] is True
    assert payload["plastic_channel_representation"] == "reduced_2d_slip_trace_channels"
    assert payload["full_3d_bcc_slip_system_model_active"] is False
    assert payload["emission_observation_required_for_implementation_certification"] is False
    assert payload["emission_observed_in_this_run"] is False
    assert payload["zero_emission_is_valid_implementation_outcome"] is True
    assert payload["physics_recomputed"] is False
    assert payload["source_outputs_modified"] is False
    assert payload["channel_rows_written"] == 4
    assert (tmp_path / NORMALIZED_AUDIT).is_file()
    assert (tmp_path / NORMALIZED_RESULTS).is_file()
    assert (tmp_path / CHANNEL_TABLE).is_file()

    rows = list(csv.DictReader((tmp_path / CHANNEL_TABLE).open(newline="")))
    assert len(rows) == 4
    assert rows[0]["slip_trace_channel_name"].startswith("2D_slip_trace_channel_0")
    normalized = json.loads((tmp_path / NORMALIZED_RESULTS).read_text())
    assert normalized[0]["emission_observed_in_this_run"] is False
    assert normalized[0]["slip_trace_channel_count"] == 2


def test_missing_per_channel_diagnostics_fails_closed(tmp_path):
    _write_source_output(tmp_path, include_channels=False)
    with pytest.raises(RuntimeError, match="per-channel drive columns"):
        normalize_output(tmp_path)
