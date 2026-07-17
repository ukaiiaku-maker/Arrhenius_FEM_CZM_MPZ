#!/usr/bin/env python3
"""Fail-closed audit for the repaired v10.0.3 one-segment integration."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from arrhenius_fracture.mode_i_first_passage_v10_0_3_progressive import (
    source_population_bound,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def find_one(root: Path, name: str) -> Path:
    direct = root / name
    if direct.exists():
        return direct
    matches = sorted(root.rglob(name))
    require(len(matches) == 1, f"expected exactly one {name}; found {len(matches)}")
    return matches[0]


def read_steps(root: Path) -> tuple[Path, list[dict[str, str]]]:
    paths = sorted(root.glob("steps_*K.csv"))
    require(len(paths) == 1, f"expected one steps CSV; found {len(paths)}")
    with paths[0].open(newline="") as fp:
        rows = list(csv.DictReader(fp))
    require(rows, "steps CSV is empty")
    return paths[0], rows


def f(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    require(math.isfinite(value), f"non-finite {key}")
    return value


def audit(root: Path, target_um: float = 5.0) -> dict[str, Any]:
    root = Path(root)
    runtime = load(find_one(root, "kinetic_campaign_czm_progressive_2d_v10_0_3.json"))
    model = load(find_one(root, "kinetic_campaign_czm_v10_0_3_audit.json"))
    quality = load(find_one(root, "explicit_quality_wrapper_chain_v91856.json"))
    summary = load(find_one(root, "summary.json"))
    results = load(find_one(root, "mode_i_v10_0_3_results.json"))
    _, steps = read_steps(root)

    require(runtime.get("schema") == "kinetic_campaign_czm_progressive_2d_v10_0_3",
            f"wrong runtime schema: {runtime.get('schema')}")
    require(runtime.get("full_progressive_trial_loop_active") is True,
            "dedicated progressive trial lifecycle is inactive")
    require(runtime.get("delayed_transform_entered") is True,
            "delayed transform was not entered")
    require(runtime.get("live_binding_capture_verified") is True,
            "live binding capture was not verified")
    require(runtime.get("opening_coupling_at_transform") == "clock_linear",
            f"wrong coupling at transform: {runtime.get('opening_coupling_at_transform')}")
    require(runtime.get("engine_factory_called") is True,
            "audited v10 engine factory was not called")
    require(runtime.get("engine_state_model") == "kinetic_campaign_czm",
            f"wrong engine state: {runtime.get('engine_state_model')}")
    require(runtime.get("orientation_match") is True,
            "FEM and directional crystal orientations differ")
    require(math.isclose(float(runtime.get("fem_crystal_theta_deg")), 45.0,
                         rel_tol=0.0, abs_tol=1e-12),
            f"unexpected crystal orientation: {runtime.get('fem_crystal_theta_deg')}")

    bindings = runtime.get("binding_ids", {})
    required_bindings = {
        "make_tri_mesh", "assemble_mechanics", "solve_dirichlet",
        "compute_J_integral", "update_plasticity", "build_engine",
    }
    require(required_bindings.issubset(bindings),
            f"missing binding audits: {sorted(required_bindings - set(bindings))}")
    require(all(bool(bindings[name].get("match")) for name in required_bindings),
            "one or more transformed bindings were stale")

    require(int(runtime.get("trial_insertions", 0)) == 1,
            f"expected one trial insertion; got {runtime.get('trial_insertions')}")
    require(int(runtime.get("committed_events", 0)) == 1,
            f"expected one committed event; got {runtime.get('committed_events')}")
    require(int(runtime.get("full_rollbacks", 0)) == 0,
            "unexpected full topology rollback")
    require(int(runtime.get("max_commits_in_outer_interval", 0)) <= 1,
            "more than one commit occurred in one equilibrium state")
    require(float(runtime.get("mpz_advance_on_commit_m", math.nan)) == 0.0,
            "MPZ translated twice at outer commit")

    final = steps[-1]
    extension = f(final, "crack_extension_m")
    require(math.isclose(extension, target_um * 1e-6, rel_tol=0.0, abs_tol=2e-12),
            f"wrong final extension: {extension * 1e6:.12g} um")
    B_final = f(final, "B")
    require(0.0 <= B_final < 1.0 + 1e-12, f"invalid residual B={B_final}")

    budget = float(runtime.get("source_budget_total", math.nan))
    require(math.isfinite(budget) and budget > 0.0, "finite source budget is missing")
    material = model.get("material", {})
    refresh_length = float(material.get("source_refresh_length_m", math.nan))
    require(math.isfinite(refresh_length) and refresh_length > 0.0,
            "source refresh length is missing")
    population_bound = source_population_bound(budget, extension, refresh_length)
    population_tolerance = max(1e-8, 1e-10 * max(population_bound, 1.0))
    max_nem = max(f(row, "N_em") for row in steps)
    require(max_nem <= population_bound + population_tolerance,
            "finite source law violated: "
            f"max N_em={max_nem} > refresh-aware bound={population_bound}")

    require(isinstance(summary, list) and len(summary) == 1, "summary must contain one case")
    require(int(summary[0].get("n_advances", -1)) == 1,
            f"summary advance accounting is wrong: {summary[0].get('n_advances')}")
    require(math.isclose(float(summary[0].get("a_final_mm")), 0.505,
                         rel_tol=0.0, abs_tol=2e-9),
            f"summary final crack position is wrong: {summary[0].get('a_final_mm')}")

    require(isinstance(results, list) and len(results) == 1, "v10.0.3 results missing")
    result = results[0]
    require(result.get("front_state_model") == "kinetic_campaign_czm",
            "result retained legacy front-state label")
    require(result.get("point_release") == "10.0.3",
            "result point release is not 10.0.3")
    require(result.get("B_final") is not None,
            "result B_final is still null")
    require(math.isclose(float(result.get("source_population_bound")), population_bound,
                         rel_tol=1e-12, abs_tol=1e-12),
            "result source-population bound differs from independent audit")

    require(model.get("full_progressive_trial_loop_active") is True,
            "model audit did not certify progressive lifecycle")
    require(model.get("live_binding_capture_verified") is True,
            "model audit did not certify live bindings")
    require(model.get("campaign_dispatch_active") is True,
            "campaign dispatch is not certified")
    require(model.get("stress_channels_separated") is True,
            "stress channels were not certified as separated")
    require(model.get("wake_shielding_active") is False,
            "wake shielding was unexpectedly activated")
    require(model.get("stored_energy_cleavage_active") is False,
            "stored-energy cleavage was unexpectedly activated")
    require(model.get("artificial_sigma_cap_active") is False,
            "artificial stress cap was unexpectedly activated")
    require(model.get("artificial_emission_cap_active") is False,
            "artificial emission cap was unexpectedly activated")
    require(model.get("artificial_N_sat_active") is False,
            "artificial source saturation was unexpectedly activated")
    require(model.get("penalty_convergence_authorized") is False,
            "penalty convergence was prematurely authorized")
    require(model.get("long_progressive_runs_authorized") is False,
            "long progressive growth was prematurely authorized")

    engine_audits = model.get("engine_audits", [])
    require(len(engine_audits) == 1, f"expected one engine audit; found {len(engine_audits)}")
    engine_audit = engine_audits[0]
    require(engine_audit.get("front_state_model") == "kinetic_campaign_czm",
            "engine audit reports the wrong front-state model")
    require(engine_audit.get("stress_channels_separated") is True,
            "engine audit does not certify separated stress channels")
    require(engine_audit.get("continuous_mpz_translation") is True,
            "continuous MPZ translation is not active")
    require(engine_audit.get("source_refresh_from_advance_only") is True,
            "source refresh is not advance-only")
    require(engine_audit.get("wake_shielding_active") is False,
            "engine audit reports wake shielding")
    require(engine_audit.get("stored_energy_cleavage_active") is False,
            "engine audit reports stored-energy cleavage")

    require(not quality.get("quality_vetoes"),
            f"geometry quality vetoes occurred: {quality.get('quality_vetoes')}")
    require(quality.get("consecutive_veto_abort") in (None, False),
            "quality wrapper reported a consecutive-veto abort")

    czm_log_path = find_one(root, "czm_advance_log.json")
    czm = load(czm_log_path)
    czm_rows = czm if isinstance(czm, list) else czm.get("rows", [])
    require(czm_rows, "CZM advance log is empty")
    committed_length = sum(
        math.hypot(float(row["x1"]) - float(row["x0"]),
                   float(row["y1"]) - float(row["y0"]))
        for row in czm_rows
        if float(row.get("damage", 1.0)) >= 1.0 - 1e-10
    )
    require(math.isclose(committed_length, target_um * 1e-6,
                         rel_tol=0.0, abs_tol=2e-12),
            f"committed CZM length is {committed_length * 1e6:.12g} um")

    out = {
        "schema": "v10_0_3_progressive_integration_certification",
        "certified": True,
        "target_um": target_um,
        "B_final": B_final,
        "max_N_em": max_nem,
        "source_budget_total": budget,
        "source_refresh_length_m": refresh_length,
        "source_population_bound": population_bound,
        "runtime": runtime,
    }
    (root / "v10_0_3_progressive_integration_certification.json").write_text(
        json.dumps(out, indent=2, default=str)
    )
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    p.add_argument("--target-um", type=float, default=5.0)
    args = p.parse_args()
    result = audit(args.root, args.target_um)
    print("V10.0.3 PROGRESSIVE INTEGRATION CERTIFIED")
    print(json.dumps({
        "target_um": result["target_um"],
        "B_final": result["B_final"],
        "max_N_em": result["max_N_em"],
        "source_budget_total": result["source_budget_total"],
        "source_population_bound": result["source_population_bound"],
    }, indent=2))


if __name__ == "__main__":
    main()
