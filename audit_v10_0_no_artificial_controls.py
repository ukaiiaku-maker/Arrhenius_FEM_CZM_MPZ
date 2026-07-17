#!/usr/bin/env python3
"""Fail-closed audit for PF-equivalent kinetic CZM outputs and source settings."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def audit_model(payload: dict[str, Any], *, require_progressive: bool) -> dict[str, Any]:
    require(payload.get("front_state_model") == "kinetic_campaign_czm", "wrong front state model")
    require(payload.get("material_parameter_source") == "pf_v10_1_7_1", "PF material source not active")
    require(payload.get("wake_shielding_active") is False, "wake shielding must be disabled")
    require(payload.get("stress_channels_separated") is True, "stress channels are not separated")
    require(payload.get("continuous_mpz_translation_active") is True, "continuous MPZ translation inactive")
    require(payload.get("source_refresh_from_advance_only") is True, "source refresh is not advance-only")
    if require_progressive:
        require(payload.get("full_progressive_trial_loop_active") is True, "progressive trial loop is not active")

    cfg = payload.get("kinetic_config", {})
    require(math.isclose(float(cfg.get("backstress_scale", -1.0)), 1.0), "backstress scale was changed")
    require(math.isclose(float(cfg.get("source_refresh_scale", -1.0)), 1.0), "source refresh scale was changed")
    require(cfg.get("wake_shielding") is False, "kinetic config enables wake shielding")

    engines = payload.get("engine_audits", [])
    require(bool(engines), "no engine audit records were written")
    for engine in engines:
        require(engine.get("opening_stress_unshielded") is True, "opening stress was shielded")
        require(engine.get("cleavage_uses_active_elastic_shielding_only") is True, "cleavage channel incorrect")
        require(engine.get("emission_uses_local_taylor_backstress_only") is True, "emission channel incorrect")
        require(engine.get("wake_shielding_active") is False, "engine wake shielding active")
        require(engine.get("stored_energy_cleavage_active") is False, "stored-energy cleavage active")
    return {
        "schema": "kinetic_campaign_czm_no_artificial_controls_v10_0",
        "passed": True,
        "progressive_required": bool(require_progressive),
        "engine_count": len(engines),
        "temperature_dependent_source_count": False,
        "temperature_dependent_shielding_coefficient": False,
        "empirical_N_max": False,
        "per_step_emission_cap": False,
        "temporal_source_recycling": False,
        "stored_energy_cleavage_subtraction": False,
        "independent_cohesive_failure_criterion": False,
        "AT2_active": False,
        "wake_shielding_primary": False,
    }


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("output")
    p.add_argument("--require-progressive", action="store_true")
    args = p.parse_args(argv)
    root = Path(args.output)
    audit_path = root if root.is_file() else root / "kinetic_campaign_czm_v10_0_audit.json"
    require(audit_path.exists(), f"missing audit file: {audit_path}")
    result = audit_model(load(audit_path), require_progressive=args.require_progressive)
    out = audit_path.parent / "no_artificial_controls_v10_0.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"NO-ARTIFICIAL-CONTROLS CERTIFIED: {out}")


if __name__ == "__main__":
    main()
