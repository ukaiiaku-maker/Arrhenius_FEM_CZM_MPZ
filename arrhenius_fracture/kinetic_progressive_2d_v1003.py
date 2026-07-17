"""v10.0.3 delayed progressive transform with live-binding verification.

The transform is constructed only when v9.11 calls ``sharp_front.run_2d``.  At
that point the active v10 engine factory and the mechanics/J/plasticity wrappers
have already been installed.  Their live function objects are synchronized into
the source function globals before the transformed namespace is copied.
"""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from . import fem as _fem
from . import j_integral as _j_integral
from . import mesh as _mesh
from . import plasticity as _plasticity
from . import sharp_front as _sharp_front
from . import kinetic_progressive_2d_v1002 as _v1002
from .kinetic_progressive_2d_v1003_source import (
    build_progressive_run_2d_v1003_source,
)

SCHEMA = "kinetic_campaign_czm_progressive_2d_v10_0_3"

_AUDIT: dict[str, Any] = {}


def reset_progressive_runtime_v1003() -> None:
    _v1002.reset_progressive_runtime_v1002()
    _AUDIT.clear()
    _AUDIT.update({
        "schema": SCHEMA,
        "delayed_transform_entered": False,
        "live_binding_capture_verified": False,
        "engine_factory_called": False,
        "engine_state_model": None,
        "engine_class": None,
        "engine_mro": [],
        "source_budget_total": None,
        "fem_crystal_theta_deg": None,
        "directional_crystal_theta_deg": None,
        "orientation_match": False,
        "binding_ids": {},
    })


def _sync_live_aliases(original_run_2d) -> tuple[dict[str, Any], dict[str, Any]]:
    live = original_run_2d.__globals__
    names = {
        "make_tri_mesh": _mesh.make_tri_mesh,
        "assemble_mechanics": _fem.assemble_mechanics,
        "solve_dirichlet": _fem.solve_dirichlet,
        "compute_J_integral": _j_integral.compute_J_integral,
        "update_plasticity": _plasticity.update_plasticity,
        "build_engine": _sharp_front.build_engine,
    }
    saved = {name: live.get(name) for name in names}
    live.update(names)
    return saved, names


def _restore_live_aliases(original_run_2d, saved: dict[str, Any]) -> None:
    live = original_run_2d.__globals__
    for name, value in saved.items():
        live[name] = value


def build_delayed_progressive_run_2d_v1003(original_run_2d, original_base_build):
    """Return a dispatcher that transforms only after v9.11 patches are live."""

    def delayed(args):
        _AUDIT["delayed_transform_entered"] = True
        saved, live_bindings = _sync_live_aliases(original_run_2d)
        active_build = live_bindings["build_engine"]
        if active_build is original_base_build:
            _restore_live_aliases(original_run_2d, saved)
            raise RuntimeError(
                "v10.0.3 delayed transform reached run_2d before the v10 engine "
                "factory was installed"
            )

        def audited_build(parsed_args, material):
            eng = active_build(parsed_args, material)
            _AUDIT["engine_factory_called"] = True
            _AUDIT["engine_state_model"] = getattr(eng, "state_model", None)
            _AUDIT["engine_class"] = type(eng).__name__
            _AUDIT["engine_mro"] = [c.__name__ for c in type(eng).__mro__]
            capacity = getattr(getattr(eng, "mpz_state", None), "site_capacity", None)
            if capacity is not None:
                _AUDIT["source_budget_total"] = float(np.sum(capacity))
            fem_theta = float(getattr(parsed_args, "crystal_theta_deg", 0.0) or 0.0)
            directional_theta = float(
                getattr(getattr(eng, "_mm", None), "crystal_theta_deg", fem_theta)
            )
            _AUDIT["fem_crystal_theta_deg"] = fem_theta
            _AUDIT["directional_crystal_theta_deg"] = directional_theta
            _AUDIT["orientation_match"] = bool(
                abs(fem_theta - directional_theta) <= 1.0e-12
            )
            if getattr(eng, "state_model", None) != "kinetic_campaign_czm":
                raise RuntimeError(
                    "v10.0.3 live factory returned the wrong front state: "
                    f"{getattr(eng, 'state_model', None)!r}"
                )
            if not bool(getattr(eng, "supports_progressive_kinetic_czm", False)):
                raise RuntimeError(
                    "v10.0.3 engine lacks progressive kinetic-CZM capability"
                )
            if not _AUDIT["orientation_match"]:
                raise RuntimeError(
                    "v10.0.3 FEM/directional crystal orientation mismatch: "
                    f"{fem_theta:g} versus {directional_theta:g} deg"
                )
            return eng

        original_run_2d.__globals__["build_engine"] = audited_build
        try:
            transformed = build_progressive_run_2d_v1003_source(original_run_2d)
            expected = {
                "make_tri_mesh": _mesh.make_tri_mesh,
                "assemble_mechanics": _fem.assemble_mechanics,
                "solve_dirichlet": _fem.solve_dirichlet,
                "compute_J_integral": _j_integral.compute_J_integral,
                "update_plasticity": _plasticity.update_plasticity,
                "build_engine": audited_build,
            }
            mismatches = []
            for name, value in expected.items():
                captured = transformed.__globals__.get(name)
                _AUDIT["binding_ids"][name] = {
                    "captured": id(captured),
                    "expected": id(value),
                    "match": captured is value,
                }
                if captured is not value:
                    mismatches.append(name)
            if mismatches:
                raise RuntimeError(
                    "v10.0.3 transformed namespace captured stale bindings: "
                    + ", ".join(mismatches)
                )
            _AUDIT["live_binding_capture_verified"] = True
        finally:
            _restore_live_aliases(original_run_2d, saved)

        result = transformed(args)
        runtime = _v1002.progressive_runtime_payload_v1002()
        _AUDIT["v1002_runtime_after_run"] = copy.deepcopy(runtime)
        if not runtime.get("full_progressive_trial_loop_active", False):
            raise RuntimeError(
                "v10.0.3 transformed run completed without the dedicated trial "
                "lifecycle; binding audit=" + json.dumps(_AUDIT, default=str)
            )
        return result

    delayed._v1003_delayed_transform = True
    delayed._v1003_original_run_2d = original_run_2d
    return delayed


def progressive_runtime_payload_v1003() -> dict[str, Any]:
    base = _v1002.progressive_runtime_payload_v1002()
    payload = copy.deepcopy(base)
    payload.update(copy.deepcopy(_AUDIT))
    payload["schema"] = SCHEMA
    payload["opening_coupling_env"] = os.environ.get(
        "ARRHENIUS_CZM_OPENING_COUPLING"
    )
    payload["full_progressive_trial_loop_active"] = bool(
        base.get("full_progressive_trial_loop_active", False)
        and _AUDIT.get("delayed_transform_entered", False)
        and _AUDIT.get("live_binding_capture_verified", False)
        and _AUDIT.get("engine_factory_called", False)
        and _AUDIT.get("engine_state_model") == "kinetic_campaign_czm"
        and _AUDIT.get("orientation_match", False)
    )
    return payload


def write_progressive_runtime_audit_v1003(out: str | Path) -> Path:
    root = Path(out)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "kinetic_campaign_czm_progressive_2d_v10_0_3.json"
    path.write_text(json.dumps(
        progressive_runtime_payload_v1003(), indent=2, default=str
    ))
    return path


__all__ = [
    "SCHEMA",
    "reset_progressive_runtime_v1003",
    "build_delayed_progressive_run_2d_v1003",
    "progressive_runtime_payload_v1003",
    "write_progressive_runtime_audit_v1003",
]
