#!/usr/bin/env python3
"""Certify the v10.0.2 anisotropic-elastic, straight-path contract."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def find_one(root: Path, name: str) -> Path:
    direct = root / name
    if direct.exists():
        return direct
    matches = sorted(root.rglob(name))
    require(len(matches) == 1, f"expected exactly one {name}; found {len(matches)}")
    return matches[0]


def audit(root: Path) -> dict[str, Any]:
    root = Path(root)
    model_path = find_one(root, "kinetic_campaign_czm_v10_0_audit.json")
    runtime_path = find_one(
        root, "kinetic_campaign_czm_progressive_2d_v10_0_2.json"
    )
    model = load(model_path)
    runtime = load(runtime_path)

    require(model.get("point_release") == "10.0.2", "wrong point release")
    require(model.get("anisotropic_elasticity_active") is True,
            "model audit does not certify anisotropic elasticity")
    require(model.get("anisotropic_J_active") is True,
            "model audit does not certify anisotropic J")
    require(model.get("anisotropic_path_selection_active") is False,
            "anisotropic path selection was not forced off")
    require(model.get("straight_single_front_mode_I_checkpoint") is True,
            "straight single-front checkpoint contract missing")

    require(runtime.get("anisotropic_elasticity_preserved") is True,
            "runtime did not preserve anisotropic elasticity")
    require(runtime.get("anisotropic_J_preserved") is True,
            "runtime did not preserve anisotropic J")
    require(runtime.get("path_deflection_forced_off") is True,
            "runtime did not force path deflection off")
    require(runtime.get("anisotropic_path_selection_active") is False,
            "runtime reports anisotropic path selection active")
    require(runtime.get("straight_single_front_mode_I_checkpoint") is True,
            "runtime straight checkpoint flag missing")
    require(runtime.get("one_topology_event_per_equilibrium_state") is True,
            "one-topology-event equilibrium contract missing")

    payload = {
        "schema": "v10_0_2_anisotropic_elastic_straight_path_contract",
        "passed": True,
        "anisotropic_elasticity_active": True,
        "anisotropic_J_active": True,
        "path_deflection_active": False,
        "branching_active": False,
        "straight_single_front_mode_I_checkpoint": True,
        "model_audit": str(model_path),
        "runtime_audit": str(runtime_path),
    }
    out = root / "anisotropic_elastic_straight_path_contract_v10_0_2.json"
    out.write_text(json.dumps(payload, indent=2))
    return payload


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("output")
    args = p.parse_args(argv)
    payload = audit(Path(args.output))
    print(
        "V10.0.2 ANISOTROPIC-ELASTIC STRAIGHT-PATH CERTIFIED: "
        f"anisotropic_J={payload['anisotropic_J_active']} "
        f"path_deflection={payload['path_deflection_active']}"
    )


if __name__ == "__main__":
    main()
