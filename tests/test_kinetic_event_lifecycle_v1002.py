from __future__ import annotations

from types import SimpleNamespace

import pytest

from arrhenius_fracture.kinetic_event_lifecycle_v1002 import (
    EventLifecycleConfig,
    KineticEventLifecycleController,
)


def result(*, accepted, committed=False, consumed=0.0, recommended=None):
    return SimpleNamespace(
        accepted=accepted,
        committed=committed,
        dt_consumed_s=consumed,
        recommended_dt_s=recommended,
    )


def test_rejected_increment_retries_without_consuming_time():
    controller = KineticEventLifecycleController()
    calls = []

    def ensure_trial():
        return "trial-0"

    def advance(dt):
        calls.append(dt)
        if len(calls) == 1:
            return result(accepted=False, recommended=0.25)
        return result(accepted=True, consumed=dt)

    out = controller.consume_interval(
        total_dt_s=1.0,
        ensure_trial=ensure_trial,
        advance_trial=advance,
    )
    assert out.rejected_attempts == 1
    assert out.committed_events == 0
    assert out.consumed_dt_s == pytest.approx(1.0)
    assert out.unused_dt_s == pytest.approx(0.0)
    assert calls == pytest.approx([1.0, 0.25, 0.75])


def test_unused_time_after_commit_is_offered_to_next_trial():
    controller = KineticEventLifecycleController()
    current_event = {"value": 0}
    contexts = []

    def ensure_trial():
        context = f"trial-{current_event['value']}"
        contexts.append(context)
        return context

    def advance(dt):
        if current_event["value"] == 0:
            return result(accepted=True, committed=True, consumed=0.4)
        return result(accepted=True, committed=False, consumed=dt)

    def on_commit(context, step):
        assert context == "trial-0"
        current_event["value"] += 1

    out = controller.consume_interval(
        total_dt_s=1.0,
        ensure_trial=ensure_trial,
        advance_trial=advance,
        on_commit=on_commit,
    )
    assert out.committed_events == 1
    assert out.consumed_dt_s == pytest.approx(1.0)
    assert out.unused_dt_s == pytest.approx(0.0)
    assert contexts == ["trial-0", "trial-1"]
    assert [x.requested_dt_s for x in out.accepted_steps] == pytest.approx([1.0, 0.6])


def test_target_stop_preserves_remaining_interval_after_commit():
    controller = KineticEventLifecycleController()
    target = {"reached": False}

    def ensure_trial():
        return "trial-0"

    def advance(dt):
        return result(accepted=True, committed=True, consumed=0.4)

    def on_commit(context, step):
        target["reached"] = True

    out = controller.consume_interval(
        total_dt_s=1.0,
        ensure_trial=ensure_trial,
        advance_trial=advance,
        on_commit=on_commit,
        target_reached=lambda: target["reached"],
    )
    assert out.stopped_at_target
    assert out.committed_events == 1
    assert out.consumed_dt_s == pytest.approx(0.4)
    assert out.unused_dt_s == pytest.approx(0.6)


def test_invalid_accepted_zero_time_is_rejected():
    controller = KineticEventLifecycleController(
        EventLifecycleConfig(max_retries_per_substep=2)
    )
    with pytest.raises(RuntimeError, match="consumed no finite physical time"):
        controller.consume_interval(
            total_dt_s=1.0,
            ensure_trial=lambda: "trial",
            advance_trial=lambda dt: result(accepted=True, consumed=0.0),
        )
