#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
from typing import Any


SCHEMA = "state_coupled_material_differentiation_v91858_v1"


def _load_json(path: Path) -> Any:
    if not path.exists():
        raise SystemExit(f"required file missing: {path}")
    return json.loads(path.read_text())


def _float(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if math.isfinite(out) else float(default)


def _read_curve(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise SystemExit(f"R-curve file missing: {path}")
    rows: list[dict[str, float]] = []
    with path.open(newline="") as fp:
        for row in csv.DictReader(fp):
            x = _float(row.get("crack_extension_after_um"), math.nan)
            K = _float(row.get("KJ_MPa_sqrt_m"), math.nan)
            U = _float(row.get("Uapp_m"), math.nan)
            Ksh = _float(row.get("K_shield_MPa_sqrt_m"), 0.0)
            if math.isfinite(x) and math.isfinite(K) and math.isfinite(U) and abs(U) > 0.0:
                rows.append({"x_um": x, "K": K, "U": U, "Kshield": Ksh})
    rows.sort(key=lambda r: r["x_um"])
    if len(rows) < 2:
        raise SystemExit(f"R-curve has fewer than two usable events: {path}")
    return rows


def _source_saturation(case_dir: Path) -> dict[str, float | int | None]:
    path = case_dir / "absolute_hazard_event_audit_v917.json"
    if not path.exists():
        return {"event_count": 0, "inventory_limited_count": 0, "inventory_limited_fraction": None}
    payload = _load_json(path)
    events = list(payload.get("event_summaries", []))
    limited = 0
    usable = 0
    for event in events:
        emitted = _float(event.get("dN_emit_relaxation"), math.nan)
        refreshed = _float(event.get("source_sites_refreshed_during_opening"), math.nan)
        if not math.isfinite(emitted) or not math.isfinite(refreshed) or refreshed <= 0.0:
            continue
        usable += 1
        if math.isclose(emitted, refreshed, rel_tol=1.0e-7, abs_tol=1.0e-10):
            limited += 1
    return {
        "event_count": usable,
        "inventory_limited_count": limited,
        "inventory_limited_fraction": None if usable == 0 else limited / usable,
    }


def _bulk_summary(case_dir: Path) -> dict[str, Any]:
    integration = _load_json(case_dir / "mpz_v9_11_integration_audit.json")
    bulk = dict(integration.get("bulk_PT", {}))
    final = _load_json(case_dir / "bulk_state_v9_11_summary.json")
    mode = str(bulk.get("mode", final.get("bulk_plasticity_mode", "")))
    explicit = bool(
        bulk.get("explicit_mobile_retained_state", final.get("bulk_explicit_mobile_retained_state", False))
    )
    calls = int(_float(final.get("bulk_state_update_calls", final.get("calls", 0)), 0.0))
    accepted_dep = _float(
        final.get("bulk_accepted_dep_mean_acc", final.get("accepted_dep_mean_acc", 0.0)), 0.0
    )
    transfer_count = int(_float(final.get("bulk_remesh_transfer_count", 0), 0.0))
    transfer_status = final.get("bulk_remesh_transfer_status")
    return {
        "mode": mode,
        "explicit_mobile_retained_state": explicit,
        "state_update_calls": calls,
        "accepted_dep_mean_acc": accepted_dep,
        "remesh_transfer_count": transfer_count,
        "remesh_transfer_status": transfer_status,
    }


def _case_payload(row: dict[str, Any]) -> dict[str, Any]:
    case_dir = Path(row["case_dir"])
    bulk = _bulk_summary(case_dir)
    if bulk["mode"] != "bulk_same_pt_km" or not bulk["explicit_mobile_retained_state"]:
        raise SystemExit(
            "material-differentiation campaign used an elastic bulk: "
            f"case={case_dir} mode={bulk['mode']!r} explicit={bulk['explicit_mobile_retained_state']}"
        )
    if bulk["state_update_calls"] <= 0:
        raise SystemExit(f"bulk state was never updated: {case_dir}")
    if bulk["accepted_dep_mean_acc"] <= 0.0:
        raise SystemExit(f"bulk state updated but accumulated no accepted plastic strain: {case_dir}")

    curve = _read_curve(case_dir / "R_curve_topology_events_raw.csv")
    K0 = curve[0]["K"]
    if K0 <= 0.0:
        raise SystemExit(f"nonpositive initial K in {case_dir}")
    for point in curve:
        point["K_normalized"] = point["K"] / K0
        point["geometry_factor"] = point["K"] / point["U"]
        point["shield_fraction"] = abs(point["Kshield"]) / max(abs(point["K"]), 1.0e-30)

    return {
        "class": row.get("class"),
        "T_K": _float(row.get("T_K"), math.nan),
        "case_dir": str(case_dir),
        "bulk": bulk,
        "curve": curve,
        "maximum_direct_shield_fraction": max(p["shield_fraction"] for p in curve),
        "source_inventory": _source_saturation(case_dir),
    }


def _aligned(a: dict[str, Any], b: dict[str, Any]) -> list[tuple[dict[str, float], dict[str, float]]]:
    by_x_a = {round(p["x_um"], 9): p for p in a["curve"]}
    by_x_b = {round(p["x_um"], 9): p for p in b["curve"]}
    common = sorted(set(by_x_a) & set(by_x_b))
    if len(common) < 2:
        raise SystemExit(
            f"fewer than two common crack extensions for {a['class']} and {b['class']}"
        )
    return [(by_x_a[x], by_x_b[x]) for x in common]


def _pair_metrics(
    a: dict[str, Any],
    b: dict[str, Any],
    min_normalized_k_separation: float,
    min_geometry_factor_separation: float,
) -> dict[str, Any]:
    pairs = _aligned(a, b)
    norm_sep = max(abs(x["K_normalized"] - y["K_normalized"]) for x, y in pairs)
    rms = math.sqrt(
        sum((x["K_normalized"] - y["K_normalized"]) ** 2 for x, y in pairs) / len(pairs)
    )
    geom_sep = max(
        abs(x["geometry_factor"] - y["geometry_factor"])
        / max(0.5 * (abs(x["geometry_factor"]) + abs(y["geometry_factor"])), 1.0e-30)
        for x, y in pairs
    )
    collapsed = (
        norm_sep < float(min_normalized_k_separation)
        and geom_sep < float(min_geometry_factor_separation)
    )
    return {
        "class_a": a["class"],
        "class_b": b["class"],
        "T_K": a["T_K"],
        "common_event_count": len(pairs),
        "maximum_normalized_K_separation": norm_sep,
        "rms_normalized_K_separation": rms,
        "maximum_relative_geometry_factor_separation": geom_sep,
        "minimum_required_normalized_K_separation": min_normalized_k_separation,
        "minimum_required_geometry_factor_separation": min_geometry_factor_separation,
        "normalized_shape_collapse": collapsed,
    }


def audit_campaign(
    case_root: Path,
    min_normalized_k_separation: float = 0.02,
    min_geometry_factor_separation: float = 0.01,
    allow_single: bool = False,
) -> dict[str, Any]:
    root = Path(case_root)
    summary = _load_json(root / "v9_13_campaign_summary.json")
    if not isinstance(summary, list) or not summary:
        raise SystemExit(f"campaign summary is empty or invalid: {root}")
    failed = [r for r in summary if int(r.get("subprocess_returncode", r.get("returncode", 1)) or 0) != 0]
    if failed:
        raise SystemExit(f"inner solver failures present: {[(r.get('class'), r.get('returncode')) for r in failed]}")
    incomplete = [r for r in summary if not bool(r.get("target_completed", False))]
    if incomplete:
        raise SystemExit(f"target-incomplete material cases: {[r.get('class') for r in incomplete]}")

    cases = [_case_payload(row) for row in summary]
    pair_rows: list[dict[str, Any]] = []
    for a, b in itertools.combinations(cases, 2):
        if not math.isclose(float(a["T_K"]), float(b["T_K"]), rel_tol=0.0, abs_tol=1.0e-9):
            continue
        pair_rows.append(
            _pair_metrics(
                a,
                b,
                min_normalized_k_separation,
                min_geometry_factor_separation,
            )
        )
    if not pair_rows and not allow_single:
        raise SystemExit("no same-temperature material pair was available for differentiation audit")

    collapsed = [row for row in pair_rows if row["normalized_shape_collapse"]]
    payload = {
        "schema": SCHEMA,
        "case_root": str(root),
        "state_coupled_bulk_required": True,
        "case_count": len(cases),
        "pair_count": len(pair_rows),
        "cases": cases,
        "pairs": pair_rows,
        "collapsed_pair_count": len(collapsed),
        "material_shape_differentiation_passed": len(collapsed) == 0,
        "constitutive_parameter_scaling_added": False,
        "artificial_reload_gate_added": False,
        "artificial_backstress_multiplier_added": False,
    }
    out = root / "state_coupled_material_differentiation_v91858.json"
    out.write_text(json.dumps(payload, indent=2))

    for case in cases:
        sat = case["source_inventory"]["inventory_limited_fraction"]
        sat_text = "n/a" if sat is None else f"{float(sat):.3f}"
        print(
            f"STATE COUPLING {case['class']}: calls={case['bulk']['state_update_calls']} "
            f"accepted_dep={case['bulk']['accepted_dep_mean_acc']:.6e} "
            f"remesh_transfers={case['bulk']['remesh_transfer_count']} "
            f"max_direct_Kshield_over_K={case['maximum_direct_shield_fraction']:.6e} "
            f"inventory_limited_fraction={sat_text}"
        )
    for pair in pair_rows:
        print(
            f"SHAPE AUDIT {pair['class_a']} vs {pair['class_b']}: "
            f"max_normK_sep={pair['maximum_normalized_K_separation']:.6e} "
            f"rms_normK_sep={pair['rms_normalized_K_separation']:.6e} "
            f"max_geometry_factor_sep={pair['maximum_relative_geometry_factor_separation']:.6e} "
            f"collapsed={pair['normalized_shape_collapse']}"
        )
    if collapsed:
        names = [(r["class_a"], r["class_b"]) for r in collapsed]
        raise SystemExit(f"normalized material-response collapse detected: {names}")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_root", type=Path)
    parser.add_argument("--min-normalized-k-separation", type=float, default=0.02)
    parser.add_argument("--min-geometry-factor-separation", type=float, default=0.01)
    parser.add_argument("--allow-single", action="store_true")
    args = parser.parse_args()
    audit_campaign(
        args.case_root,
        min_normalized_k_separation=args.min_normalized_k_separation,
        min_geometry_factor_separation=args.min_geometry_factor_separation,
        allow_single=args.allow_single,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
