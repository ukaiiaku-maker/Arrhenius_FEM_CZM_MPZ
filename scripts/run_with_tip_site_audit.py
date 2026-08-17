#!/usr/bin/env python3
"""Run the normal sharp-front driver while serializing accepted tip-site history.

The physics path is unchanged.  The wrapper only registers engine instances
through their existing reset method and, after the normal driver returns,
serializes the already-existing ``emission_event_history`` from the instance
with the longest accepted history.

The production stochastic-emission engine already records for each accepted
emission event:

* multiplicity_at_event
* rate_per_site_at_event_s
* aggregate_hazard_at_event_s
* tip_radius_before_event_m / tip_radius_after_event_m
* front_width_at_event_m
* line_content_added and threshold information

Rollback restores this history transactionally, so the final retained history
is an accepted-state diagnostic rather than a new constitutive state variable.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import weakref
from pathlib import Path


def _install_registry():
    from arrhenius_fracture.persistent_site_joint_K_ramp_emission_v10051832 import (
        PersistentSiteJointKRampEventDrivenFrontEngineV10051832 as Engine,
    )

    registry: list[weakref.ReferenceType] = []
    original_reset = Engine.reset

    def audited_reset(self, *args, **kwargs):
        result = original_reset(self, *args, **kwargs)
        registry.append(weakref.ref(self))
        return result

    Engine.reset = audited_reset
    return Engine, original_reset, registry


def _live_instances(registry):
    seen = set()
    out = []
    for ref in registry:
        obj = ref()
        if obj is None or id(obj) in seen:
            continue
        seen.add(id(obj)); out.append(obj)
    return out


def select_longest_history(instances):
    if not instances:
        return None
    return max(instances, key=lambda x: len(getattr(x, "emission_event_history", [])))


def _safe_geometry(engine):
    try:
        return copy.deepcopy(engine.mpz_state.source_geometry())
    except Exception as exc:
        return {"unavailable": True, "error": repr(exc)}


def _dump(path: Path, registry, *, argv, status, error=None):
    instances = _live_instances(registry)
    selected = select_longest_history(instances)
    histories = []
    for idx, engine in enumerate(instances):
        history = copy.deepcopy(getattr(engine, "emission_event_history", []))
        histories.append({
            "engine_index": idx,
            "python_object_id": id(engine),
            "history_length": len(history),
            "final_source_geometry": _safe_geometry(engine),
            "history": history,
        })
    payload = {
        "schema": "dbtt_tip_site_accepted_history_v1",
        "diagnostic_only": True,
        "physics_modified": False,
        "driver_argv": argv,
        "status": status,
        "error": error,
        "registered_engine_instances": len(instances),
        "selected_history_length": (
            len(getattr(selected, "emission_event_history", [])) if selected is not None else 0
        ),
        "selected_history": (
            copy.deepcopy(getattr(selected, "emission_event_history", [])) if selected is not None else []
        ),
        "selected_final_source_geometry": _safe_geometry(selected) if selected is not None else None,
        "all_engine_histories": histories,
        "interpretation": (
            "Existing transactional emission_event_history serialized after the normal run. "
            "No site count, radius, hazard, RNG, crack path, or constitutive update is altered."
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--audit-out", type=Path, required=True)
    parser.add_argument("driver_args", nargs=argparse.REMAINDER)
    ns = parser.parse_args()
    driver_args = list(ns.driver_args)
    if driver_args and driver_args[0] == "--":
        driver_args = driver_args[1:]
    if not driver_args:
        raise SystemExit("sharp_front driver arguments are required after --")

    Engine, original_reset, registry = _install_registry()
    saved_argv = sys.argv[:]
    try:
        from arrhenius_fracture import sharp_front
        sys.argv = ["sharp_front"] + driver_args
        sharp_front.main()
    except BaseException as exc:
        _dump(ns.audit_out, registry, argv=driver_args, status="failed", error=repr(exc))
        raise
    else:
        _dump(ns.audit_out, registry, argv=driver_args, status="complete")
    finally:
        sys.argv = saved_argv
        Engine.reset = original_reset


if __name__ == "__main__":
    main()
