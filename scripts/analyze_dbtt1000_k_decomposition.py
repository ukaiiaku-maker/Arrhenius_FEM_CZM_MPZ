#!/usr/bin/env python3
"""DBTT/1000 K PF-vs-FEM/CZM apparent/local K decomposition.

Read-only diagnostic: no fracture kinetics, RNG, mesh, crack path, or
constitutive state is modified.

Three K measures are supported:

KJ_local
    The live J-equivalent K recorded by each solver for its resolved state.

K_app_initial_geometry
    Remote force converted with the shared fresh-crack K/F calibration.  This
    isolates the changing remote load without attributing microscopic tip
    geometry to the measured K.

K_app_SENT
    K_app_initial_geometry multiplied by a conventional projected-crack
    single-edge-tension finite-width correction.  This is an experimental-style
    diagnostic, not an ASTM validity claim at large a/W.

If --reference-calibration is supplied, K_app_reference is also reported from
an independently generated elastic reference-FEM K/F calibration versus
projected crack length.  That is the preferred apparent-K definition for the
final mechanism audit.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

REQUIRED = {
    "step", "Uapp_m", "Ftop_N", "KJ_Pa_sqrtm", "sigma_tip_Pa",
    "sigma_back_Pa", "crack_extension_m", "n_fire",
}


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        missing = REQUIRED.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        rows = []
        for raw in reader:
            row = {}
            for key, value in raw.items():
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    row[key] = value
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no rows")
    return rows


def sent_y(x: np.ndarray) -> np.ndarray:
    """Common finite-width single-edge-tension correction polynomial."""
    return 1.12 - 0.231*x + 10.55*x**2 - 21.72*x**3 + 30.39*x**4


def sent_g(a_m: np.ndarray, width_m: float) -> np.ndarray:
    x = np.asarray(a_m, dtype=float) / float(width_m)
    return np.sqrt(np.pi*np.asarray(a_m, dtype=float)) * sent_y(x)


def load_reference_calibration(path: Path | None):
    if path is None:
        return None
    with path.open(newline="") as fp:
        reader = csv.DictReader(fp)
        needed = {"projected_crack_length_m", "K_over_F_Pa_sqrtm_per_N"}
        missing = needed.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing calibration columns {sorted(missing)}")
        a, kpf = [], []
        for row in reader:
            a.append(float(row["projected_crack_length_m"]))
            kpf.append(float(row["K_over_F_Pa_sqrtm_per_N"]))
    order = np.argsort(a)
    return np.asarray(a)[order], np.asarray(kpf)[order]


def interp_reference(a_m: np.ndarray, calibration):
    if calibration is None:
        return np.full_like(a_m, np.nan, dtype=float)
    x, y = calibration
    if np.min(a_m) < x[0] - 1e-15 or np.max(a_m) > x[-1] + 1e-15:
        raise ValueError("reference calibration does not cover requested crack-length range")
    return np.interp(a_m, x, y)


def first_row_match(a: dict, b: dict) -> dict[str, float]:
    keys = ("Uapp_m", "Ftop_N", "KJ_Pa_sqrtm", "crack_extension_m")
    out = {}
    for key in keys:
        av, bv = float(a[key]), float(b[key])
        scale = max(abs(av), abs(bv), 1e-300)
        out[key] = abs(av-bv)/scale
    return out


def infer_tip_radius(row: dict[str, float]) -> float:
    """Kinematic radius proxy from recorded KJ and sigma_tip.

    It equals the analytic tip radius only if sigma_tip is generated directly
    from KJ/sqrt(2*pi*r) without an active stress cap or unrecorded K offset.
    It is deliberately labelled a proxy in all output.
    """
    k, s = float(row["KJ_Pa_sqrtm"]), float(row["sigma_tip_Pa"])
    if not math.isfinite(k) or not math.isfinite(s) or s <= 0:
        return math.nan
    return (k/s)**2/(2.0*math.pi)


def event_mask(rows: list[dict[str, float]]) -> np.ndarray:
    ext = np.asarray([float(r["crack_extension_m"]) for r in rows])
    fire = np.asarray([float(r.get("n_fire", 0.0)) for r in rows]) > 0
    dext = np.r_[ext[0], np.diff(ext)]
    return fire | (dext > 1e-15)


def process_model(label, rows, *, a0_m, width_m, k_over_f0, reference_calibration):
    ext = np.asarray([float(r["crack_extension_m"]) for r in rows])
    a = a0_m + ext
    force = np.asarray([float(r["Ftop_N"]) for r in rows])
    u = np.asarray([float(r["Uapp_m"]) for r in rows])
    kj = np.asarray([float(r["KJ_Pa_sqrtm"]) for r in rows])
    k_initial = force*k_over_f0
    g0 = float(sent_g(np.asarray([a0_m]), width_m)[0])
    k_sent = k_initial*sent_g(a, width_m)/g0
    kpf_ref = interp_reference(a, reference_calibration)
    k_ref = force*kpf_ref
    is_event = event_mask(rows)

    out = []
    for i, r in enumerate(rows):
        out.append({
            "model": label,
            "source_row": i,
            "step": float(r["step"]),
            "event_row": int(is_event[i]),
            "projected_extension_m": float(ext[i]),
            "projected_crack_length_m": float(a[i]),
            "Uapp_m": float(u[i]),
            "Ftop_N": float(force[i]),
            "stiffness_F_over_U_N_per_m": float(force[i]/u[i]) if u[i] else math.nan,
            "compliance_U_over_F_m_per_N": float(u[i]/force[i]) if force[i] else math.nan,
            "KJ_local_Pa_sqrtm": float(kj[i]),
            "K_app_initial_geometry_Pa_sqrtm": float(k_initial[i]),
            "K_app_SENT_Pa_sqrtm": float(k_sent[i]),
            "K_app_reference_Pa_sqrtm": float(k_ref[i]),
            "deltaK_initial_geometry_minus_local_Pa_sqrtm": float(k_initial[i]-kj[i]),
            "deltaK_SENT_minus_local_Pa_sqrtm": float(k_sent[i]-kj[i]),
            "deltaK_reference_minus_local_Pa_sqrtm": float(k_ref[i]-kj[i]),
            "mapping_ratio_local_over_initial_geometry": float(kj[i]/k_initial[i]) if k_initial[i] else math.nan,
            "mapping_ratio_local_over_SENT": float(kj[i]/k_sent[i]) if k_sent[i] else math.nan,
            "sigma_tip_Pa": float(r["sigma_tip_Pa"]),
            "sigma_back_Pa": float(r["sigma_back_Pa"]),
            "tip_radius_proxy_m": infer_tip_radius(r),
            "N_em": float(r.get("N_em", math.nan)),
            "mpz_K_shield_Pa_sqrtm": float(r.get("mpz_K_shield_Pa_sqrt_m", math.nan)),
            "mpz_mobile_count": float(r.get("mpz_mobile_count", math.nan)),
            "mpz_retained_count": float(r.get("mpz_retained_count", math.nan)),
            "mpz_local_slip_count": float(r.get("mpz_local_slip_count", math.nan)),
        })
    return out


def write_rows(path: Path, rows: Iterable[dict]):
    rows = list(rows)
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)


def closest_event(rows, extension_um):
    ev = [r for r in rows if int(r["event_row"]) == 1] or rows
    return min(ev, key=lambda r: abs(float(r["projected_extension_m"])*1e6-extension_um))


def summary_for(label, rows):
    ev = [r for r in rows if int(r["event_row"]) == 1] or rows
    first, last = ev[0], ev[-1]
    checkpoints = {}
    for ext in (0, 25, 50, 100, 200, 300, 500, 750, 1000):
        r = closest_event(rows, ext)
        checkpoints[str(ext)] = {
            "actual_extension_um": float(r["projected_extension_m"])*1e6,
            "KJ_local_MPa_sqrtm": float(r["KJ_local_Pa_sqrtm"])/1e6,
            "K_app_initial_geometry_MPa_sqrtm": float(r["K_app_initial_geometry_Pa_sqrtm"])/1e6,
            "K_app_SENT_MPa_sqrtm": float(r["K_app_SENT_Pa_sqrtm"])/1e6,
            "K_app_reference_MPa_sqrtm": float(r["K_app_reference_Pa_sqrtm"])/1e6,
            "Ftop_N": float(r["Ftop_N"]),
            "stiffness_F_over_U_N_per_m": float(r["stiffness_F_over_U_N_per_m"]),
            "tip_radius_proxy_um": float(r["tip_radius_proxy_m"])*1e6,
            "sigma_back_GPa": float(r["sigma_back_Pa"])/1e9,
        }
    return {
        "label": label,
        "event_rows": len(ev),
        "first_event_extension_um": float(first["projected_extension_m"])*1e6,
        "last_event_extension_um": float(last["projected_extension_m"])*1e6,
        "KJ_first_MPa_sqrtm": float(first["KJ_local_Pa_sqrtm"])/1e6,
        "KJ_last_MPa_sqrtm": float(last["KJ_local_Pa_sqrtm"])/1e6,
        "remote_force_first_N": float(first["Ftop_N"]),
        "remote_force_last_N": float(last["Ftop_N"]),
        "stiffness_first_N_per_m": float(first["stiffness_F_over_U_N_per_m"]),
        "stiffness_last_N_per_m": float(last["stiffness_F_over_U_N_per_m"]),
        "stiffness_last_over_first": float(last["stiffness_F_over_U_N_per_m"])/float(first["stiffness_F_over_U_N_per_m"]),
        "checkpoints": checkpoints,
    }


def make_plots(outdir: Path, fem, pf):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []
    paths = []
    sets = [("FEM/CZM", [r for r in fem if r["event_row"]]),
            ("PF/sharp-wake", [r for r in pf if r["event_row"]])]

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    for label, rows in sets:
        a = np.asarray([r["projected_extension_m"] for r in rows])*1e6
        ax.plot(a, np.asarray([r["KJ_local_Pa_sqrtm"] for r in rows])/1e6,
                label=f"{label}: live $K_J$")
        ax.plot(a, np.asarray([r["K_app_initial_geometry_Pa_sqrtm"] for r in rows])/1e6,
                linestyle="--", label=f"{label}: remote-load $K_{{app,0}}$")
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel(r"K (MPa$\sqrt{m}$)")
    ax.set_title("DBTT 1000 K: live local K vs remote-load apparent K")
    ax.legend(fontsize=8); ax.grid(alpha=0.25); fig.tight_layout()
    p = outdir/"dbtt1000_K_local_vs_apparent.png"; fig.savefig(p, dpi=180); plt.close(fig); paths.append(str(p))

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    for label, rows in sets:
        a = np.asarray([r["projected_extension_m"] for r in rows])*1e6
        stiff = np.asarray([r["stiffness_F_over_U_N_per_m"] for r in rows])
        ax.plot(a, stiff/stiff[0], label=label)
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel("Global stiffness / first-event stiffness")
    ax.set_title("DBTT 1000 K: global compliance evolution")
    ax.legend(); ax.grid(alpha=0.25); fig.tight_layout()
    p = outdir/"dbtt1000_global_stiffness.png"; fig.savefig(p, dpi=180); plt.close(fig); paths.append(str(p))

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    for label, rows in sets:
        a = np.asarray([r["projected_extension_m"] for r in rows])*1e6
        ax.plot(a, np.asarray([r["tip_radius_proxy_m"] for r in rows])*1e6, label=label)
    ax.set_xlabel("Projected crack extension (µm)")
    ax.set_ylabel("Kinematic tip-radius proxy (µm)")
    ax.set_title("DBTT 1000 K: tip-radius proxy from recorded KJ and tip stress")
    ax.legend(); ax.grid(alpha=0.25); fig.tight_layout()
    p = outdir/"dbtt1000_tip_radius_proxy.png"; fig.savefig(p, dpi=180); plt.close(fig); paths.append(str(p))
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fem-steps", type=Path, required=True)
    ap.add_argument("--pf-steps", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--specimen-width-m", type=float, default=0.002)
    ap.add_argument("--initial-crack-m", type=float, default=0.0005)
    ap.add_argument("--reference-calibration", type=Path)
    ap.add_argument("--first-row-match-rtol", type=float, default=1e-10)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    fem_raw, pf_raw = read_csv(args.fem_steps), read_csv(args.pf_steps)
    match = first_row_match(fem_raw[0], pf_raw[0])
    if max(match.values()) > args.first_row_match_rtol:
        raise RuntimeError(f"fresh-reference first rows do not match: {match}")
    f0, k0 = float(fem_raw[0]["Ftop_N"]), float(fem_raw[0]["KJ_Pa_sqrtm"])
    if f0 == 0:
        raise RuntimeError("fresh-reference force is zero")
    k_over_f0 = k0/f0
    ref = load_reference_calibration(args.reference_calibration)
    fem = process_model("FEM/CZM", fem_raw, a0_m=args.initial_crack_m,
                        width_m=args.specimen_width_m, k_over_f0=k_over_f0,
                        reference_calibration=ref)
    pf = process_model("PF/sharp-wake", pf_raw, a0_m=args.initial_crack_m,
                       width_m=args.specimen_width_m, k_over_f0=k_over_f0,
                       reference_calibration=ref)
    write_rows(args.out/"dbtt1000_k_decomposition.csv", fem+pf)
    summary = {
        "schema": "dbtt1000_pf_fem_k_decomposition_v1",
        "diagnostic_only": True,
        "fresh_reference_first_row_relative_mismatch": match,
        "fresh_reference_K_over_F_Pa_sqrtm_per_N": k_over_f0,
        "initial_crack_m": args.initial_crack_m,
        "specimen_width_m": args.specimen_width_m,
        "apparent_K_definitions": {
            "K_app_initial_geometry": "F_actual * (KJ/F)_fresh_reference",
            "K_app_SENT": "K_app_initial_geometry * g(a_projected)/g(a0), using conventional single-edge-tension Y polynomial",
            "K_app_reference": "F_actual * interpolated independently calibrated elastic reference-FEM K/F versus projected crack length (NaN unless --reference-calibration supplied)",
        },
        "caveat": "K_app_SENT is an experimental-style projected-crack diagnostic, not an ASTM validity claim; prefer K_app_reference when a direct reference-FEM calibration is available.",
        "FEM_CZM": summary_for("FEM/CZM", fem),
        "PF_sharp_wake": summary_for("PF/sharp-wake", pf),
    }
    summary["cross_model"] = {
        "first_event_remote_force_ratio_PF_over_FEM": summary["PF_sharp_wake"]["remote_force_first_N"] / summary["FEM_CZM"]["remote_force_first_N"],
        "last_event_remote_force_ratio_PF_over_FEM": summary["PF_sharp_wake"]["remote_force_last_N"] / summary["FEM_CZM"]["remote_force_last_N"],
        "first_event_KJ_difference_MPa_sqrtm": summary["PF_sharp_wake"]["KJ_first_MPa_sqrtm"] - summary["FEM_CZM"]["KJ_first_MPa_sqrtm"],
        "last_event_KJ_difference_MPa_sqrtm": summary["PF_sharp_wake"]["KJ_last_MPa_sqrtm"] - summary["FEM_CZM"]["KJ_last_MPa_sqrtm"],
    }
    summary["plots"] = make_plots(args.out, fem, pf)
    (args.out/"dbtt1000_k_decomposition_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True)
    )
    print(json.dumps(summary["cross_model"], indent=2))


if __name__ == "__main__":
    main()
