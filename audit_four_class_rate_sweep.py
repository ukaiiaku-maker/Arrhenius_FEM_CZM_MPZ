#!/usr/bin/env python3
"""Audit rate preservation and adaptive time stepping for the 1x/10x/100x CZM sweep."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def rate_label(f: float) -> str:
    return f"rate_{int(f)}x" if float(f).is_integer() else f"rate_{f:g}x"


def find_steps(case_dir: Path) -> Path | None:
    files = sorted(case_dir.glob("steps_*K.csv"))
    return files[0] if files else None


def audit_case(case_dir: Path, factor: float, base_dU: float, base_dt: float,
               adaptive_target: float, adaptive_min_frac: float) -> dict:
    sf = find_steps(case_dir)
    row = {
        "case_dir": str(case_dir),
        "rate_factor": factor,
        "nominal_dU_m": base_dU,
        "nominal_dt_s": base_dt / factor,
        "nominal_opening_rate_m_per_s": base_dU / (base_dt / factor),
    }
    if sf is None:
        row.update(status="missing_steps", n_steps=0)
        return row
    try:
        st = pd.read_csv(sf)
    except Exception as exc:
        row.update(status=f"read_error:{exc}", n_steps=0)
        return row

    row["n_steps"] = len(st)
    if len(st) == 0:
        row["status"] = "empty_steps"
        return row

    required = {"Uapp_m", "dt_cur_s", "adaptive_frac"}
    missing = required.difference(st.columns)
    if missing:
        row.update(status="missing_columns:" + ",".join(sorted(missing)))
        return row

    U = pd.to_numeric(st["Uapp_m"], errors="coerce").to_numpy(float)
    dt = pd.to_numeric(st["dt_cur_s"], errors="coerce").to_numpy(float)
    frac = pd.to_numeric(st["adaptive_frac"], errors="coerce").to_numpy(float)
    dU = np.diff(np.r_[0.0, U])
    valid = np.isfinite(dU) & np.isfinite(dt) & (dt > 0.0) & np.isfinite(frac)
    realized = np.full(len(st), np.nan)
    realized[valid] = dU[valid] / dt[valid]
    expected_rate = base_dU / (base_dt / factor)
    relerr = np.abs(realized[valid] / expected_rate - 1.0) if np.any(valid) else np.array([np.nan])

    dt_expected = (base_dt / factor) * frac
    dt_relerr = np.abs(dt[valid] / dt_expected[valid] - 1.0) if np.any(valid) else np.array([np.nan])
    min_hit = frac <= adaptive_min_frac * (1.0 + 1e-6)

    row.update({
        "status": "ok",
        "realized_rate_median_m_per_s": float(np.nanmedian(realized)),
        "rate_rel_error_max": float(np.nanmax(relerr)),
        "dt_fraction_relation_rel_error_max": float(np.nanmax(dt_relerr)),
        "adaptive_frac_min": float(np.nanmin(frac)),
        "adaptive_frac_p01": float(np.nanpercentile(frac, 1)),
        "adaptive_frac_median": float(np.nanmedian(frac)),
        "adaptive_frac_p99": float(np.nanpercentile(frac, 99)),
        "min_fraction_hit_count": int(np.count_nonzero(min_hit)),
        "min_fraction_hit_fraction": float(np.mean(min_hit)),
        "dt_min_s": float(np.nanmin(dt)),
        "dt_median_s": float(np.nanmedian(dt)),
        "dt_max_s": float(np.nanmax(dt)),
    })
    if "adaptive_dB_total" in st.columns:
        db = pd.to_numeric(st["adaptive_dB_total"], errors="coerce").to_numpy(float)
        row["adaptive_dB_total_max"] = float(np.nanmax(db))
        row["adaptive_dB_target_exceed_fraction"] = float(np.nanmean(db > adaptive_target * 1.05))
    if "crack_extension_m" in st.columns:
        ext = pd.to_numeric(st["crack_extension_m"], errors="coerce").to_numpy(float)
        row["final_crack_extension_um"] = float(np.nanmax(ext) * 1e6)
    row["snapshot_png_exists"] = bool(any(case_dir.glob("field_snapshots_*K.png")))

    # Rate invariance of the adaptive discretization is a hard audit. Hitting
    # the minimum fraction is reported separately because genuinely unstable
    # propagation is allowed to do so by design.
    if not np.isfinite(row["rate_rel_error_max"]) or row["rate_rel_error_max"] > 1e-4:
        row["status"] = "rate_not_preserved"
    elif row["dt_fraction_relation_rel_error_max"] > 1e-7:
        row["status"] = "dt_fraction_mismatch"
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--rate-factors", default="1 10 100")
    ap.add_argument("--base-dU", type=float, default=2e-7)
    ap.add_argument("--base-dt", type=float, default=8.4)
    ap.add_argument("--adaptive-target", type=float, default=0.35)
    ap.add_argument("--adaptive-min-frac", type=float, default=1e-8)
    args = ap.parse_args()

    root = Path(args.root)
    factors = [float(x) for x in args.rate_factors.replace(",", " ").split()]
    rows = []
    for factor in factors:
        rr = root / rate_label(factor)
        if not rr.exists():
            rows.append({"case_dir": str(rr), "rate_factor": factor, "status": "missing_rate_root", "n_steps": 0})
            continue
        for class_dir in sorted(p for p in rr.iterdir() if p.is_dir()):
            for case_dir in sorted(p for p in class_dir.iterdir() if p.is_dir() and p.name.startswith("T")):
                out = audit_case(case_dir, factor, args.base_dU, args.base_dt,
                                 args.adaptive_target, args.adaptive_min_frac)
                out["class"] = class_dir.name
                try:
                    out["T_K"] = int(case_dir.name.split("_", 1)[0][1:])
                except Exception:
                    out["T_K"] = np.nan
                rows.append(out)

    df = pd.DataFrame(rows)
    out_csv = root / "adaptive_timestep_audit.csv"
    df.to_csv(out_csv, index=False)

    summary = {}
    if not df.empty:
        summary = {
            "n_cases_audited": int(len(df)),
            "status_counts": {str(k): int(v) for k, v in df["status"].value_counts(dropna=False).items()},
            "max_rate_relative_error": float(df.get("rate_rel_error_max", pd.Series([np.nan])).max()),
            "total_min_fraction_hits": int(df.get("min_fraction_hit_count", pd.Series([0])).fillna(0).sum()),
            "cases_missing_snapshot_png": int((~df.get("snapshot_png_exists", pd.Series([False] * len(df))).fillna(False)).sum()),
        }
    (root / "adaptive_timestep_audit_summary.json").write_text(json.dumps(summary, indent=2))
    print(df.to_string(index=False) if len(df) <= 30 else df.groupby(["rate_factor", "status"]).size().to_string())
    print(f"WROTE {out_csv}")
    print(f"WROTE {root / 'adaptive_timestep_audit_summary.json'}")

    bad = df[~df["status"].isin(["ok"])] if not df.empty else df
    if len(bad):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
