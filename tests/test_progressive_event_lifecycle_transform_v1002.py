from __future__ import annotations

from types import SimpleNamespace

import pytest

from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1002 import (
    _merge_progressive_infos,
    build_progressive_run_2d_v1002,
    progressive_runtime_payload_v1002,
    reset_progressive_runtime_v1002,
)


def test_v1002_transform_compiles_against_actual_production_run_2d():
    reset_progressive_runtime_v1002()
    transformed = build_progressive_run_2d_v1002(sharp_front.run_2d)
    assert transformed is not sharp_front.run_2d
    assert transformed._v10_progressive_source_transform is True
    assert transformed._v1002_event_lifecycle is True
    payload = progressive_runtime_payload_v1002()
    assert payload["anchor_counts"] == {
        "backend_construction": 1,
        "adaptive_target": 1,
        "single_front_step": 1,
    }
    assert payload["rejected_step_retry_active"] is True
    assert payload["unused_time_carry_active"] is True
    assert payload["same_load_re_equilibration_after_commit"] is True


def test_interval_info_merges_increment_fields_and_final_state():
    infos = [
        {
            "dB_block": 0.2,
            "dN_emit_block": 3.0,
            "micro_advance_step_m": 1.0e-6,
            "B": 0.2,
            "trial_status": "trial",
        },
        {
            "dB_block": 0.8,
            "dN_emit_block": 4.0,
            "micro_advance_step_m": 4.0e-6,
            "B": 0.0,
            "trial_status": "committed",
        },
    ]
    lifecycle = SimpleNamespace(
        committed_events=1,
        consumed_dt_s=0.5,
        unused_dt_s=0.0,
        rejected_attempts=2,
        stopped_at_target=False,
    )
    out = _merge_progressive_infos(infos, lifecycle)
    assert out["dB_block"] == pytest.approx(1.0)
    assert out["dN_emit_block"] == pytest.approx(7.0)
    assert out["micro_advance_step_m"] == pytest.approx(5.0e-6)
    assert out["B"] == pytest.approx(0.0)
    assert out["trial_status"] == "committed"
    assert out["n_fire"] == 1
    assert out["event_lifecycle_rejected_attempts"] == 2
    assert out["v_crack"] == pytest.approx(1.0e-5)
