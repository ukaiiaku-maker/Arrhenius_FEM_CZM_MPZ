from __future__ import annotations

import json

import arrhenius_fracture.mode_i_first_passage_v9_14 as mode_i_v914
from arrhenius_fracture.mode_i_first_passage_v9_14 import (
    _LoadingProxy,
    _PostEventEquilibriumController,
)


def test_loading_proxy_controls_existing_values():
    class Loading:
        dt = 8.4
        dU_top = 2.0e-7

    ctl = _PostEventEquilibriumController()
    loading = _LoadingProxy(Loading(), ctl)
    assert loading.dt * 0.5 == 4.2
    assert loading.dU_top * 0.5 == 1.0e-7
    ctl.schedule_physical_event()
    assert loading.dt * 1.0 == 0.0
    assert loading.dU_top * 1.0 == 0.0
    assert ctl.corrected_event_ids == [1]


def test_v914_hooks_v911_runtime_step_drives(monkeypatch, tmp_path):
    def fake_step_drives(self, *args, **kwargs):
        return {"n_fire": 1, "fired": True}

    monkeypatch.setattr(
        mode_i_v914.MovingProcessZone2DFrontEngine,
        "step_drives",
        fake_step_drives,
    )

    def fake_main(argv):
        dummy = object.__new__(mode_i_v914.MovingProcessZone2DFrontEngine)
        out = mode_i_v914.MovingProcessZone2DFrontEngine.step_drives(
            dummy, 1.0, 1.0, 700.0, 1.0
        )
        assert out["v914_physical_event_id"] == 1
        return []

    monkeypatch.setattr(mode_i_v914._base, "main", fake_main)
    mode_i_v914.main(["--out", str(tmp_path)])
    payload = json.loads(
        (tmp_path / "post_event_equilibrium_audit_v914.json").read_text()
    )
    assert payload["runtime_engine_hook"] == (
        "MovingProcessZone2DFrontEngine.step_drives"
    )
    assert payload["events_scheduled"] == 1
    assert payload["scheduled_event_ids"] == [1]
