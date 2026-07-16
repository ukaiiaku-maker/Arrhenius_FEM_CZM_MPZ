from __future__ import annotations

import audit_matched_stress_classes_v917 as base
import audit_matched_stress_classes_v9171 as routed
from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState


def test_v9171_audit_routes_to_same_state_as_fem(monkeypatch):
    seen = {}

    def fake_main():
        seen["state"] = base.MovingProcessZoneState
        return "ok"

    monkeypatch.setattr(base, "main", fake_main)
    original = base.MovingProcessZoneState
    result = routed.main()

    assert result == "ok"
    assert seen["state"] is MovingProcessZoneState
    assert base.MovingProcessZoneState is original
