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

plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.labelsize': 15,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'lines.linewidth': 2.0,
    'savefig.bbox': 'tight',
    'figure.dpi': 160,
})

CLASSES = ["ceramic", "peak", "weakT", "DBTT"]


def first_existing(df: pd.DataFrame, names: Iterable[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def read_case_metadata(case_dir: Path) -> dict:
    meta = {
        'case_dir': str(case_dir),
        'seed': np.nan,
        'replicate': np.nan,
        'complete': False,
        'target_hit': False,
        'final_extension_um': np.nan,
    }
    m = re.search(r"seed(\d+)", str(case_dir))
    if m:
        meta['seed'] = int(m.group(1))
    m = re.search(r"replicate_(\d+)", str(case_dir))
    if m:
        meta['replicate'] = int(m.group(1))

    sj = case_dir / 'summary.json'
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            if isinstance(s, list) and s:
                s = s[0]
            if s.get('a_final_mm') is not None:
                meta['final_extension_um'] = (float(s['a_final_mm']) - 0.5) * 1000.0
            for k in ['final_crack_extension_um', 'crack_extension_um', 'extension_um']:
                if s.get(k) is not None:
                    val = float(s[k])
                    if not math.isfinite(meta['final_extension_um']):
                        meta['final_extension_um'] = val
                    else:
                        meta['final_extension_um'] = max(meta['final_extension_um'], val)
        except Exception:
            pass

    log = case_dir / 'run.log'
    if log.exists():
        try:
            txt = log.read_text(errors='ignore')
            meta['target_hit'] = 'reached target crack extension' in txt
            found = re.findall(r"reached target crack extension\s+([0-9.]+)\s+um", txt)
            if found:
                val = float(found[-1])
                if not math.isfinite(meta['final_extension_um']):
                    meta['final_extension_um'] = val
                else:
                    meta['final_extension_um'] = max(meta['final_extension_um'], val)
        except Exception:
            pass

    meta['complete'] = bool(math.isfinite(meta['final_extension_um']) and meta['final_extension_um'] >= 980.0)
    return meta


def read_curve(case_dir: Path) -> pd.DataFrame:
    rc = case_dir / 'R_curve_event_sampled.csv'
    if rc.exists():
        df = pd.read_csv(rc)
        kcol = first_existing(df, ['KJ_MPa_sqrt_m', 'K_MPa_sqrt_m', 'Kc_MPa_sqrt_m'])
        xcol = first_existing(df, ['crack_extension_um', 'extension_um', 'delta_a_um', 'da_cumulative_um'])
        if kcol is None:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        out = pd.DataFrame({
            'crack_extension_um': pd.to_numeric(df[xcol], errors='coerce') if xcol else np.arange(len(df), dtype=float),
            'K_MPa_sqrt_m': pd.to_numeric(df[kcol], errors='coerce'),
        })
        out = out.dropna().sort_values('crack_extension_um')
        out = out[out['K_MPa_sqrt_m'] > 0]
        return out.drop_duplicates('crack_extension_um').reset_index(drop=True)

    steps = case_dir / 'steps_0500K.csv'
    if not steps.exists():
        any_steps = sorted(case_dir.glob('steps_*K.csv'))
        if any_steps:
            steps = any_steps[0]
    if steps.exists():
        df = pd.read_csv(steps)
        keep = np.zeros(len(df), dtype=bool)
        da_col = first_existing(df, ['da_block_m', 'da_block'])
        if da_col:
            keep |= pd.to_numeric(df[da_col], errors='coerce').fillna(0).to_numpy(float) > 0
        nfire_col = first_existing(df, ['n_fire', 'nfire'])
        if nfire_col:
            keep |= pd.to_numeric(df[nfire_col], errors='coerce').fillna(0).to_numpy(float) > 0
        ev = df.loc[keep].copy()
        if ev.empty:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        kcol = first_existing(ev, ['KJ_Pa_sqrtm', 'KJ_MPa_sqrt_m', 'K_MPa_sqrt_m'])
        if kcol is None:
            return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])
        if kcol == 'KJ_Pa_sqrtm':
            kvals = pd.to_numeric(ev[kcol], errors='coerce') / 1e6
        else:
            kvals = pd.to_numeric(ev[kcol], errors='coerce')
        xcol_m = first_existing(ev, ['crack_extension_m'])
        xcol_um = first_existing(ev, ['crack_extension_um', 'extension_um'])
        if xcol_um:
            xvals = pd.to_numeric(ev[xcol_um], errors='coerce')
        elif xcol_m:
            xvals = pd.to_numeric(ev[xcol_m], errors='coerce') * 1e6
        else:
            xvals = np.arange(len(ev), dtype=float)
        out = pd.DataFrame({'crack_extension_um': xvals, 'K_MPa_sqrt_m': kvals})
        out = out.dropna().sort_values('crack_extension_um')
        out = out[out['K_MPa_sqrt_m'] > 0]
        return out.drop_duplicates('crack_extension_um').reset_index(drop=True)

    return pd.DataFrame(columns=['crack_extension_um', 'K_MPa_sqrt_m'])


def build_mean_curve(curves: list[pd.DataFrame], n_grid: int = 300) -> pd.DataFrame:
    valid = [c for c in curves if not c.empty]
    if not valid:
        return pd.DataFrame(columns=['crack_extension_um', 'K_mean_MPa_sqrt_m', 'K_std_MPa_sqrt_m', 'n_contrib'])
    xmax = max(float(c['crack_extension_um'].max()) for c in valid)
    xmin = min(float(c['crack_extension_um'].min()) for c in valid)
    xgrid = np.linspace(xmin, xmax, n_grid)
    mat = np.full((len(valid), len(xgrid)), np.nan)
    for i, c in enumerate(valid):
        x = c['crack_extension_um'].to_numpy(float)
        y = c['K_MPa_sqrt_m'].to_numpy(float)
        if len(x) < 2:
            continue
        order = np.argsort(x)
        x = x[order]
        y = y[order]
        uniq, idx = np.unique(x, return_index=True)
        x = uniq
        y = y[idx]
        mask = (xgrid >= x.min()) & (xgrid <= x.max())
        mat[i, mask] = np.interp(xgrid[mask], x, y)
    n_contrib = np.sum(np.isfinite(mat), axis=0)
    mean = np.nanmean(mat, axis=0)
    std = np.nanstd(mat, axis=0, ddof=1)
    out = pd.DataFrame({
        'crack_extension_um': xgrid,
        'K_mean_MPa_sqrt_m': mean,
        'K_std_MPa_sqrt_m': std,
        'n_contrib': n_contrib,
    })
    out = out[np.isfinite(out['K_mean_MPa_sqrt_m'])].reset_index(drop=True)
    return out


def gather_class(root: Path, klass: str) -> tuple[list[dict], list[pd.DataFrame]]:
    rows: list[dict] = []
    curves: list[pd.DataFrame] = []
    for case in sorted((root / klass).glob('replicate_*_seed*/T500_th45')):
        if 'geometry_veto' in str(case):
            continue
        meta = read_case_metadata(case)
        curve = read_curve(case)
        meta['n_points'] = len(curve)
        rows.append(meta)
        if not curve.empty:
            curve = curve.copy()
            curve['seed'] = meta['seed']
            curve['replicate'] = meta['replicate']
            curve['complete'] = meta['complete']
            curves.append(curve)
    return rows, curves


def plot_class(curves: list[pd.DataFrame], mean_curve: pd.DataFrame, title: str, out_png: Path, out_svg: Path):
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    cmap = plt.get_cmap('tab10')
    for i, curve in enumerate(sorted(curves, key=lambda d: (int(d['replicate'].iloc[0]), int(d['seed'].iloc[0])))):
        seed = int(curve['seed'].iloc[0])
        rep = int(curve['replicate'].iloc[0])
        complete = bool(curve['complete'].iloc[0])
        label = f"rep {rep:02d}, seed {seed}" + ("" if complete else " (incomplete)")
        ls = '-' if complete else '--'
        ax.plot(curve['crack_extension_um'], curve['K_MPa_sqrt_m'], linestyle=ls, alpha=0.9, color=cmap(i % 10), label=label)
    if not mean_curve.empty:
        ax.plot(mean_curve['crack_extension_um'], mean_curve['K_mean_MPa_sqrt_m'], color='black', linewidth=3.0, label='mean across seeds')
    ax.set_xlabel('Crack extension (µm)')
    ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
    ax.set_title(title)
    ax.grid(False)
    ax.legend(frameon=False, ncol=1)
    fig.tight_layout()
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_svg)
    plt.close(fig)


def plot_panel(all_data: dict[str, tuple[list[pd.DataFrame], pd.DataFrame]], out_png: Path, out_svg: Path):
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.4), sharex=False, sharey=False)
    axes = axes.ravel()
    cmap = plt.get_cmap('tab10')
    for ax, klass in zip(axes, CLASSES):
        curves, mean_curve = all_data.get(klass, ([], pd.DataFrame()))
        for i, curve in enumerate(sorted(curves, key=lambda d: (int(d['replicate'].iloc[0]), int(d['seed'].iloc[0])))):
            complete = bool(curve['complete'].iloc[0])
            ls = '-' if complete else '--'
            ax.plot(curve['crack_extension_um'], curve['K_MPa_sqrt_m'], linestyle=ls, alpha=0.85, color=cmap(i % 10))
        if not mean_curve.empty:
            ax.plot(mean_curve['crack_extension_um'], mean_curve['K_mean_MPa_sqrt_m'], color='black', linewidth=3.0)
        ax.set_title(klass)
        ax.set_xlabel('Crack extension (µm)')
        ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
        ax.grid(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_svg)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Plot K versus crack extension for each seed and the average across seeds.')
    ap.add_argument('--root', default='runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45')
    ap.add_argument('--out', default='runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/seed_rcurve_plots')
    ap.add_argument('--classes', default=' '.join(CLASSES))
    ap.add_argument('--n-grid', type=int, default=300)
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    classes = args.classes.replace(',', ' ').split()

    summary_rows = []
    panel_data: dict[str, tuple[list[pd.DataFrame], pd.DataFrame]] = {}

    for klass in classes:
        meta_rows, curves = gather_class(root, klass)
        summary_rows.extend([{**r, 'class': klass} for r in meta_rows])
        mean_curve = build_mean_curve(curves, n_grid=args.n_grid)
        panel_data[klass] = (curves, mean_curve)

        long_rows = []
        for c in curves:
            long_rows.append(c[['replicate', 'seed', 'complete', 'crack_extension_um', 'K_MPa_sqrt_m']].copy())
        long_df = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame(columns=['replicate', 'seed', 'complete', 'crack_extension_um', 'K_MPa_sqrt_m'])
        long_df.insert(0, 'class', klass)
        long_df.to_csv(out / f'{klass}_seed_rcurves_long.csv', index=False)
        mean_curve.to_csv(out / f'{klass}_mean_rcurve.csv', index=False)
        plot_class(curves, mean_curve, f'{klass}: seed R-curves and mean', out / f'{klass}_seed_rcurves_with_mean.png', out / f'{klass}_seed_rcurves_with_mean.svg')

    pd.DataFrame(summary_rows).sort_values(['class', 'replicate', 'seed']).to_csv(out / 'seed_rcurve_case_summary.csv', index=False)
    plot_panel(panel_data, out / 'four_class_seed_rcurves_with_mean_panel.png', out / 'four_class_seed_rcurves_with_mean_panel.svg')

    readme = out / 'README_seed_rcurves_with_mean.txt'
    readme.write_text(
        'For each class, the plot shows K versus crack extension for each seed and a bold black mean curve.\n'
        'Mean curves are built by linear interpolation of each seed curve to a common crack-extension grid;\n'
        'the mean at each x uses all seeds with available data there.\n'
        'Dashed individual curves indicate incomplete seeds.\n'
    )
    print(f'WROTE {out}')


if __name__ == '__main__':
    main()
