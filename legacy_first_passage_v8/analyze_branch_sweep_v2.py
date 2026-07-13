#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd


def _safe_read_csv(path: Path):
    try:
        if path.exists() and path.stat().st_size > 0:
            return pd.read_csv(path)
    except Exception:
        pass
    return None


def _path_metrics(path: Path):
    df = _safe_read_csv(path)
    if df is None or len(df) < 2:
        return dict(length_m=np.nan, angle_deg=np.nan, x_final=np.nan, y_final=np.nan)
    x = df['x_m'].to_numpy(float); y = df['y_m'].to_numpy(float)
    dx = np.diff(x); dy = np.diff(y)
    length = float(np.sum(np.sqrt(dx*dx + dy*dy)))
    angle = float(math.degrees(math.atan2(y[-1] - y[0], x[-1] - x[0])))
    return dict(length_m=length, angle_deg=angle, x_final=float(x[-1]), y_final=float(y[-1]))


def _extract(d: dict, key: str, default=np.nan):
    return d.get(key, default) if isinstance(d, dict) else default


def _morphology(row):
    if not np.isfinite(row.get('Kc_first_MPa_sqrt_m', np.nan)) or row.get('n_advances_total', 0) <= 0:
        return 'no_opening'
    if not bool(row.get('branched', False)):
        return 'opened_single_front'
    L1, L2 = row.get('length1_m', np.nan), row.get('length2_m', np.nan)
    a1, a2 = row.get('angle1_path_deg', np.nan), row.get('angle2_path_deg', np.nan)
    nf = row.get('n_advances_total', 0)
    base = 'runaway_' if np.isfinite(nf) and nf > 1e5 else ''
    if not np.isfinite(L1) or not np.isfinite(L2) or min(L1, L2) <= 0:
        return base + 'branched_unresolved'
    ratio = min(L1, L2) / max(L1, L2)
    opposite = (a1 * a2) < 0
    angle_bal = abs(abs(a1) - abs(a2))
    if ratio < 0.2:
        return base + 'daughter_short'
    if opposite and ratio > 0.5 and angle_bal < 12:
        return base + 'symmetric_Y'
    if opposite:
        return base + 'biased_Y'
    return base + 'same_side_or_turning'


def summarize_case(case_dir: Path):
    sfile = case_dir / 'summary.json'
    if not sfile.exists():
        return []
    try:
        summaries = json.loads(sfile.read_text())
    except Exception:
        return []
    try:
        args = json.loads((case_dir / 'run_args.json').read_text())
    except Exception:
        args = {}
    out = []
    for s in summaries:
        T = float(s.get('T', np.nan))
        Ttag = f'{int(round(T)):04d}K'
        diag = _safe_read_csv(case_dir / f'branch_diagnostics_{Ttag}.csv')
        if diag is None:
            cand = list(case_dir.glob(f'branch_diagnostics_*{int(round(T))}K.csv'))
            diag = _safe_read_csv(cand[0]) if cand else None
        dmax = {}
        if diag is not None and len(diag):
            if 'branch_spawned' in diag.columns and diag['branch_spawned'].astype(bool).any():
                r = diag[diag['branch_spawned'].astype(bool)].iloc[0]
            elif 'metric2_over_metric1' in diag.columns:
                arr = diag['metric2_over_metric1'].to_numpy(float)
                r = diag.iloc[int(np.nanargmax(arr))] if np.isfinite(arr).any() else diag.iloc[-1]
            else:
                r = diag.iloc[-1]
            dmax = {c: (float(r[c]) if isinstance(r[c], (np.integer, np.floating)) else r[c]) for c in diag.columns}
        p1 = _path_metrics(case_dir / f'crack_path_{int(round(T))}K.csv')
        p2 = _path_metrics(case_dir / f'crack_path_branch_{int(round(T))}K.csv')
        row = dict(
            case=str(case_dir), T_K=T,
            theta_deg=_extract(args, 'crystal_theta_deg'),
            gamma_aniso=_extract(args, 'cleave_gamma_aniso'),
            branch_overdrive_ratio=_extract(args, 'branch_overdrive_ratio'),
            branch_ratio=_extract(args, 'branch_ratio'),
            material=_extract(args, 'crystal_material_resolved', _extract(args, 'crystal_material', '')),
            Kc_first_MPa_sqrt_m=s.get('Kc_first_MPa_sqrt_m', np.nan),
            branched=bool(s.get('branched', False)),
            deflection_deg=s.get('deflection_deg', np.nan),
            path_span_dy_mm=s.get('path_span_dy_mm', np.nan),
            branch_length_dx_mm=s.get('branch_length_dx_mm', np.nan),
            N_em_init=s.get('N_em_init', np.nan),
            sigma_back_init_GPa=s.get('sigma_back_init_GPa', np.nan),
            r_eff_over_r0_init=s.get('r_eff_over_r0_init', np.nan),
            n_advances_total=s.get('n_advances', np.nan),
            n_advances_primary=s.get('n_advances_primary', np.nan),
            n_advances_branch=s.get('n_advances_branch', np.nan),
            length1_m=p1['length_m'], angle1_path_deg=p1['angle_deg'], x1_final_m=p1['x_final'], y1_final_m=p1['y_final'],
            length2_m=p2['length_m'], angle2_path_deg=p2['angle_deg'], x2_final_m=p2['x_final'], y2_final_m=p2['y_final'],
        )
        for k in ['n_candidates','angle1_deg','angle2_deg','metric1','metric2','metric2_over_metric1',
                  'branch_active','branch_spawned','share_w1','share_w2','lambda_c1','lambda_c2']:
            row[k] = dmax.get(k, np.nan)
        row['morphology'] = _morphology(row)
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    root = Path(args.root)
    rows = []
    for sfile in root.rglob('summary.json'):
        rows.extend(summarize_case(sfile.parent))
    if not rows:
        raise SystemExit(f'No summary.json files found under {root}')
    df = pd.DataFrame(rows).sort_values(['gamma_aniso','branch_overdrive_ratio','theta_deg','T_K'])
    out = Path(args.out) if args.out else root / 'branch_sweep_summary.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} rows)')
    print(df.groupby(['gamma_aniso','branch_overdrive_ratio','morphology']).size().to_string())

if __name__ == '__main__':
    main()
