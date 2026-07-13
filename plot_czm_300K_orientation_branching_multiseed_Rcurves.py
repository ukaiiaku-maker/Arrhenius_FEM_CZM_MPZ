#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def first_existing(df: pd.DataFrame, names: Iterable[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def parse_case(case: Path) -> dict:
    s = str(case)
    seed = np.nan
    theta = np.nan
    branch = False
    klass = case.parts[-4] if len(case.parts) >= 4 else ''
    m = re.search(r'seed(\d+)', s)
    if m:
        seed = int(m.group(1))
    m = re.search(r'theta_([0-9.]+)_(nobranch|weakbranch)', s)
    if m:
        theta = float(m.group(1))
        branch = (m.group(2) == 'weakbranch')
    # project/outroot/class/theta_condition/seed/T...
    parts = case.parts
    for i, part in enumerate(parts):
        if part.startswith('theta_') and i > 0:
            klass = parts[i-1]
            break
    return {'class': klass, 'theta_deg': theta, 'branching': branch, 'seed': seed}


def read_curve(case: Path) -> pd.DataFrame:
    rc = case / 'R_curve_event_sampled.csv'
    if rc.exists():
        df = pd.read_csv(rc)
        kcol = first_existing(df, ['KJ_MPa_sqrt_m', 'K_MPa_sqrt_m', 'Kc_MPa_sqrt_m'])
        xcol = first_existing(df, ['crack_extension_um', 'extension_um', 'delta_a_um'])
        if kcol is None:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        k = pd.to_numeric(df[kcol], errors='coerce')
        x = pd.to_numeric(df[xcol], errors='coerce') if xcol else np.arange(len(df), dtype=float)
        out = pd.DataFrame({'crack_extension_um': x, 'K_MPa_sqrt_m': k})
    else:
        steps = sorted(case.glob('steps_*K.csv'))
        if not steps:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        df = pd.read_csv(steps[0])
        keep = np.zeros(len(df), dtype=bool)
        da_col = first_existing(df, ['da_block_m', 'da_block_um', 'da_block'])
        if da_col:
            vals = pd.to_numeric(df[da_col], errors='coerce').fillna(0).to_numpy(float)
            keep |= vals > 0
        nfire_col = first_existing(df, ['n_fire', 'nfire'])
        if nfire_col:
            vals = pd.to_numeric(df[nfire_col], errors='coerce').fillna(0).to_numpy(float)
            keep |= vals > 0
        ev = df.loc[keep].copy()
        if ev.empty:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        kcol = first_existing(ev, ['KJ_MPa_sqrt_m', 'K_MPa_sqrt_m', 'Kc_MPa_sqrt_m', 'KJ_Pa_sqrtm'])
        if kcol is None:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        k = pd.to_numeric(ev[kcol], errors='coerce') / 1e6 if kcol == 'KJ_Pa_sqrtm' else pd.to_numeric(ev[kcol], errors='coerce')
        xcol_um = first_existing(ev, ['crack_extension_um', 'extension_um', 'delta_a_um'])
        xcol_m = first_existing(ev, ['crack_extension_m'])
        if xcol_um:
            x = pd.to_numeric(ev[xcol_um], errors='coerce')
        elif xcol_m:
            x = pd.to_numeric(ev[xcol_m], errors='coerce') * 1e6
        else:
            x = np.arange(len(ev), dtype=float)
        out = pd.DataFrame({'crack_extension_um': x, 'K_MPa_sqrt_m': k})
    out = out.replace([np.inf, -np.inf], np.nan).dropna()
    out = out[out['K_MPa_sqrt_m'] > 0]
    out = out.sort_values('crack_extension_um').drop_duplicates('crack_extension_um')
    return out.reset_index(drop=True)


def status(case: Path, target: float) -> dict:
    ext = np.nan
    target_hit = False
    log = case / 'run.log'
    if log.exists():
        txt = log.read_text(errors='ignore')
        target_hit = 'reached target crack extension' in txt
        m = re.findall(r'reached target crack extension\s+([0-9.]+)\s+um', txt)
        if m:
            ext = float(m[-1])
    sj = case / 'summary.json'
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            if isinstance(s, list) and s:
                s = s[0]
            if s.get('a_final_mm') is not None:
                val = (float(s['a_final_mm']) - 0.5) * 1000.0
                ext = max(ext, val) if math.isfinite(ext) else val
            for k in ['final_crack_extension_um', 'crack_extension_um', 'extension_um']:
                if s.get(k) is not None:
                    val = float(s[k])
                    ext = max(ext, val) if math.isfinite(ext) else val
        except Exception:
            pass
    return {'final_extension_um': ext, 'target_hit': target_hit, 'complete': bool(math.isfinite(ext) and ext >= 0.98 * target)}


def binned(curve: pd.DataFrame, bin_um: float, max_ext_um: float) -> pd.DataFrame:
    if curve.empty:
        return pd.DataFrame(columns=['bin_center_um', 'K_median_MPa_sqrt_m', 'n'])
    edges = np.arange(0, max_ext_um + bin_um, bin_um)
    x = curve['crack_extension_um'].to_numpy(float)
    y = curve['K_MPa_sqrt_m'].to_numpy(float)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        vals = y[(x >= a) & (x < b)]
        vals = vals[np.isfinite(vals)]
        if len(vals):
            rows.append({'bin_start_um': a, 'bin_end_um': b, 'bin_center_um': 0.5*(a+b), 'K_median_MPa_sqrt_m': float(np.median(vals)), 'n': int(len(vals))})
    return pd.DataFrame(rows)


def mean_from_binned(dfs: list[pd.DataFrame], bin_um: float, max_ext_um: float) -> pd.DataFrame:
    centers = np.arange(0.5*bin_um, max_ext_um + 0.5*bin_um, bin_um)
    rows = []
    for c in centers:
        vals = []
        for df in dfs:
            if df.empty:
                continue
            idx = np.argmin(np.abs(df['bin_center_um'].to_numpy(float) - c))
            if abs(float(df['bin_center_um'].iloc[idx]) - c) < 0.51*bin_um:
                vals.append(float(df['K_median_MPa_sqrt_m'].iloc[idx]))
        vals = np.array(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals):
            rows.append({'bin_center_um': c, 'K_mean_MPa_sqrt_m': float(vals.mean()), 'K_std_MPa_sqrt_m': float(vals.std(ddof=1)) if len(vals)>1 else 0.0, 'K_p10_MPa_sqrt_m': float(np.quantile(vals, 0.1)), 'K_p90_MPa_sqrt_m': float(np.quantile(vals, 0.9)), 'n_seeds': int(len(vals))})
    return pd.DataFrame(rows)


def load_cases(root: Path, target: float, bin_um: float, max_ext_um: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cases = []
    mf = root / 'sweep_manifest_final.json'
    ml = root / 'sweep_manifest_live.json'
    mi = root / 'sweep_manifest_initial.json'
    manifest = mf if mf.exists() else (ml if ml.exists() else (mi if mi.exists() else None))
    if manifest:
        data = json.loads(manifest.read_text())
        for row in data:
            if 'attempt0_incomplete' in str(row.get('case', '')):
                continue
            cases.append(root.parent.parent / row['case'] if not Path(row['case']).is_absolute() else Path(row['case']))
    else:
        cases = sorted(root.glob('*/*/seed*/T*_th*'))

    rows = []
    long = []
    for case in cases:
        if not case.exists() or 'attempt0_incomplete' in str(case):
            continue
        meta = parse_case(case)
        st = status(case, target)
        curve = read_curve(case)
        b = binned(curve, bin_um, max_ext_um)
        rows.append({**meta, **st, 'case_dir': str(case), 'n_points': len(curve), 'n_bins': len(b)})
        if not b.empty:
            bb = b.copy()
            for k, v in meta.items():
                bb[k] = v
            bb['complete'] = st['complete']
            bb['case_dir'] = str(case)
            long.append(bb)
    summary = pd.DataFrame(rows)
    longdf = pd.concat(long, ignore_index=True) if long else pd.DataFrame()

    mean_rows = []
    if not longdf.empty:
        for (theta, branch), group in summary.groupby(['theta_deg', 'branching'], dropna=False):
            dfs = []
            for case in group.loc[group['complete'], 'case_dir']:
                curve = read_curve(Path(case))
                dfs.append(binned(curve, bin_um, max_ext_um))
            mean = mean_from_binned(dfs, bin_um, max_ext_um)
            if not mean.empty:
                mean['theta_deg'] = theta
                mean['branching'] = branch
                mean_rows.append(mean)
    meandf = pd.concat(mean_rows, ignore_index=True) if mean_rows else pd.DataFrame()
    return summary, longdf, meandf


def setup_style():
    plt.rcParams.update({'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 13, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 9, 'lines.linewidth': 2.0, 'savefig.bbox': 'tight'})


def plot_orientation_panel(longdf: pd.DataFrame, meandf: pd.DataFrame, out: Path):
    setup_style()
    no = longdf[longdf['branching'] == False].copy()
    thetas = sorted(no['theta_deg'].dropna().unique())
    if not thetas:
        return
    ncols = 2
    nrows = int(math.ceil(len(thetas) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4.6*nrows), squeeze=False)
    axes = axes.ravel()
    for ax, th in zip(axes, thetas):
        g = no[no['theta_deg'] == th]
        for seed, sg in g.groupby('seed'):
            ax.plot(sg['bin_center_um'], sg['K_median_MPa_sqrt_m'], alpha=0.45, linewidth=1.3)
        mg = meandf[(meandf['branching'] == False) & (meandf['theta_deg'] == th)]
        if not mg.empty:
            ax.plot(mg['bin_center_um'], mg['K_mean_MPa_sqrt_m'], color='black', linewidth=3.0, label='mean')
            ax.fill_between(mg['bin_center_um'], mg['K_mean_MPa_sqrt_m']-mg['K_std_MPa_sqrt_m'], mg['K_mean_MPa_sqrt_m']+mg['K_std_MPa_sqrt_m'], color='black', alpha=0.12, linewidth=0)
        ax.set_title(f'θ = {th:g}°, no branching')
        ax.set_xlabel('Crack extension, Δa (µm)')
        ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
        ax.grid(False)
        ax.legend(frameon=False)
    for ax in axes[len(thetas):]:
        ax.axis('off')
    fig.tight_layout()
    fig.savefig(out / 'orientation_no_branch_seed_Rcurves_panel.png', dpi=320)
    fig.savefig(out / 'orientation_no_branch_seed_Rcurves_panel.svg')
    plt.close(fig)


def plot_mean_overlay(meandf: pd.DataFrame, out: Path):
    setup_style()
    no = meandf[meandf['branching'] == False]
    if no.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5.6))
    for th, g in no.groupby('theta_deg'):
        ax.plot(g['bin_center_um'], g['K_mean_MPa_sqrt_m'], linewidth=2.8, label=f'θ={th:g}°')
    ax.set_xlabel('Crack extension, Δa (µm)')
    ax.set_ylabel(r'Mean $K_J$ (MPa$\sqrt{m}$)')
    ax.set_title('Orientation dependence of mean R-curve, no branching')
    ax.grid(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / 'orientation_no_branch_mean_Rcurves_overlay.png', dpi=320)
    fig.savefig(out / 'orientation_no_branch_mean_Rcurves_overlay.svg')
    plt.close(fig)


def plot_branching(longdf: pd.DataFrame, meandf: pd.DataFrame, out: Path):
    setup_style()
    if longdf.empty:
        return
    ths = sorted(set(longdf.loc[longdf['branching'] == True, 'theta_deg'].dropna()))
    if not ths:
        return
    th = ths[0]
    fig, ax = plt.subplots(figsize=(8, 5.6))
    labels = {False: 'no branch', True: 'weak branch'}
    linestyles = {False: '-', True: '--'}
    for branch in [False, True]:
        g = longdf[(longdf['theta_deg'] == th) & (longdf['branching'] == branch)]
        for seed, sg in g.groupby('seed'):
            ax.plot(sg['bin_center_um'], sg['K_median_MPa_sqrt_m'], alpha=0.25, linewidth=1.0, linestyle=linestyles[branch])
        mg = meandf[(meandf['theta_deg'] == th) & (meandf['branching'] == branch)]
        if not mg.empty:
            ax.plot(mg['bin_center_um'], mg['K_mean_MPa_sqrt_m'], linewidth=3.0, linestyle=linestyles[branch], label=labels[branch])
            ax.fill_between(mg['bin_center_um'], mg['K_mean_MPa_sqrt_m']-mg['K_std_MPa_sqrt_m'], mg['K_mean_MPa_sqrt_m']+mg['K_std_MPa_sqrt_m'], alpha=0.10, linewidth=0)
    ax.set_xlabel('Crack extension, Δa (µm)')
    ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
    ax.set_title(f'Branching comparison at θ = {th:g}°')
    ax.grid(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / 'branching_comparison_seed_Rcurves.png', dpi=320)
    fig.savefig(out / 'branching_comparison_seed_Rcurves.svg')
    plt.close(fig)


def metric_summary(meandf: pd.DataFrame, out: Path):
    if meandf.empty:
        return
    rows = []
    for (theta, branch), g in meandf.groupby(['theta_deg', 'branching']):
        g = g.sort_values('bin_center_um')
        early = g[(g['bin_center_um'] >= 0) & (g['bin_center_um'] <= 50)]['K_mean_MPa_sqrt_m']
        late = g[(g['bin_center_um'] >= 700) & (g['bin_center_um'] <= 1000)]['K_mean_MPa_sqrt_m']
        auc = np.trapezoid(g['K_mean_MPa_sqrt_m'], g['bin_center_um']) / max(float(g['bin_center_um'].max() - g['bin_center_um'].min()), 1.0)
        k0 = float(early.median()) if len(early) else np.nan
        kss = float(late.median()) if len(late) else np.nan
        rows.append({'theta_deg': theta, 'branching': branch, 'K0_mean_0_50': k0, 'Kss_mean_700_1000': kss, 'DeltaK_mean': kss-k0 if np.isfinite(k0) and np.isfinite(kss) else np.nan, 'AUC_mean': auc, 'min_n_seeds_by_bin': int(g['n_seeds'].min()), 'max_n_seeds_by_bin': int(g['n_seeds'].max())})
    pd.DataFrame(rows).to_csv(out / 'condition_Rcurve_metric_summary.csv', index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='runs/czm_Rcurve_300K_orientation_branching_ceramic_multiseed_v2')
    ap.add_argument('--out', default=None)
    ap.add_argument('--bin-um', type=float, default=25.0)
    ap.add_argument('--max-ext-um', type=float, default=1000.0)
    ap.add_argument('--target-ext-um', type=float, default=1000.0)
    args = ap.parse_args()
    root = Path(args.root)
    out = Path(args.out) if args.out else root / 'plots'
    out.mkdir(parents=True, exist_ok=True)
    summary, longdf, meandf = load_cases(root, args.target_ext_um, args.bin_um, args.max_ext_um)
    summary.to_csv(out / 'case_completion_summary.csv', index=False)
    if not longdf.empty:
        longdf.to_csv(out / 'case_binned_Rcurves_long.csv', index=False)
    if not meandf.empty:
        meandf.to_csv(out / 'condition_mean_binned_Rcurves.csv', index=False)
    plot_orientation_panel(longdf, meandf, out)
    plot_mean_overlay(meandf, out)
    plot_branching(longdf, meandf, out)
    metric_summary(meandf, out)
    print(f'WROTE {out}')
    if not summary.empty:
        print(summary.groupby(['theta_deg', 'branching'])['complete'].agg(['sum','count']).to_string())


if __name__ == '__main__':
    main()
