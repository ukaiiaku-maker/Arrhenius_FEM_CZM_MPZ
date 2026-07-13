#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def first_existing(df, names):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n in df.columns:
            return n
        if n.lower() in lower:
            return lower[n.lower()]
    return None


def read_curve(case: Path) -> pd.DataFrame:
    rc = case / 'R_curve_event_sampled.csv'
    if rc.exists():
        df = pd.read_csv(rc)
        kcol = first_existing(df, ['KJ_MPa_sqrt_m','K_MPa_sqrt_m','Kc_MPa_sqrt_m'])
        xcol = first_existing(df, ['crack_extension_um','extension_um','delta_a_um'])
        if kcol is None:
            return pd.DataFrame(columns=['crack_extension_um','K_MPa_sqrt_m'])
        x = pd.to_numeric(df[xcol], errors='coerce') if xcol else np.arange(len(df), dtype=float)
        k = pd.to_numeric(df[kcol], errors='coerce')
        out = pd.DataFrame({'crack_extension_um': x, 'K_MPa_sqrt_m': k})
        out = out.dropna().query('K_MPa_sqrt_m > 0').sort_values('crack_extension_um')
        return out.drop_duplicates('crack_extension_um').reset_index(drop=True)
    steps = sorted(case.glob('steps_*K.csv'))
    if not steps:
        return pd.DataFrame(columns=['crack_extension_um','K_MPa_sqrt_m'])
    df = pd.read_csv(steps[0])
    keep = np.zeros(len(df), dtype=bool)
    da = first_existing(df, ['da_block_m','da_block_um','da_block'])
    nf = first_existing(df, ['n_fire','nfire'])
    if da:
        keep |= pd.to_numeric(df[da], errors='coerce').fillna(0).to_numpy(float) > 0
    if nf:
        keep |= pd.to_numeric(df[nf], errors='coerce').fillna(0).to_numpy(float) > 0
    ev = df.loc[keep].copy()
    if ev.empty:
        return pd.DataFrame(columns=['crack_extension_um','K_MPa_sqrt_m'])
    kcol = first_existing(ev, ['KJ_MPa_sqrt_m','K_MPa_sqrt_m','KJ_Pa_sqrtm'])
    if kcol is None:
        return pd.DataFrame(columns=['crack_extension_um','K_MPa_sqrt_m'])
    k = pd.to_numeric(ev[kcol], errors='coerce') / 1e6 if kcol == 'KJ_Pa_sqrtm' else pd.to_numeric(ev[kcol], errors='coerce')
    xum = first_existing(ev, ['crack_extension_um','extension_um'])
    xm = first_existing(ev, ['crack_extension_m'])
    if xum:
        x = pd.to_numeric(ev[xum], errors='coerce')
    elif xm:
        x = pd.to_numeric(ev[xm], errors='coerce') * 1e6
    else:
        x = np.arange(len(ev), dtype=float)
    out = pd.DataFrame({'crack_extension_um': x, 'K_MPa_sqrt_m': k})
    out = out.dropna().query('K_MPa_sqrt_m > 0').sort_values('crack_extension_um')
    return out.drop_duplicates('crack_extension_um').reset_index(drop=True)


def bin_curve(curve: pd.DataFrame, bin_um: float, max_um: float) -> pd.DataFrame:
    if curve.empty:
        return pd.DataFrame(columns=['bin_center_um','K_median_MPa_sqrt_m','n'])
    rows=[]
    x=curve['crack_extension_um'].to_numpy(float); y=curve['K_MPa_sqrt_m'].to_numpy(float)
    edges=np.arange(0, max_um + bin_um, bin_um)
    for a,b in zip(edges[:-1], edges[1:]):
        m=(x>=a)&(x<b)
        vals=y[m]; vals=vals[np.isfinite(vals)]
        if len(vals):
            rows.append({'bin_center_um':0.5*(a+b),'K_median_MPa_sqrt_m':float(np.median(vals)),'n':len(vals)})
    return pd.DataFrame(rows)


def theta_from_manifest_or_path(item):
    if 'theta_deg' in item:
        return float(item['theta_deg'])
    m = re.search(r'theta_([0-9.]+)_', item.get('case',''))
    return float(m.group(1)) if m else np.nan


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root', default='runs/czm_Rcurve_300K_orientation_branching_ceramic_v1')
    ap.add_argument('--out', default=None)
    ap.add_argument('--bin-um', type=float, default=25)
    ap.add_argument('--max-ext-um', type=float, default=1000)
    args=ap.parse_args()
    root=Path(args.root)
    out=Path(args.out) if args.out else root/'plots'
    out.mkdir(parents=True, exist_ok=True)
    manifest_path=root/'sweep_manifest.json'
    if not manifest_path.exists():
        raise SystemExit(f'Missing {manifest_path}')
    manifest=json.loads(manifest_path.read_text())
    rows=[]; curves=[]
    for m in manifest:
        case=Path(m['case'])
        c=read_curve(case)
        b=bin_curve(c, args.bin_um, args.max_ext_um)
        row={**m, 'n_raw':len(c), 'n_bins':len(b)}
        if not c.empty:
            row['max_extension_um']=float(c['crack_extension_um'].max())
            row['K_first_MPa_sqrt_m']=float(c['K_MPa_sqrt_m'].iloc[0])
            row['K_median_MPa_sqrt_m']=float(c['K_MPa_sqrt_m'].median())
        else:
            row['max_extension_um']=np.nan; row['K_first_MPa_sqrt_m']=np.nan; row['K_median_MPa_sqrt_m']=np.nan
        rows.append(row)
        if not b.empty:
            b=b.copy(); b['theta_deg']=float(m['theta_deg']); b['branching']=bool(m['branching']); b['case']=m['case']
            curves.append(b)
    pd.DataFrame(rows).to_csv(out/'orientation_branching_Rcurve_summary.csv', index=False)
    if curves:
        pd.concat(curves, ignore_index=True).to_csv(out/'orientation_branching_binned_Rcurves_long.csv', index=False)

    plt.rcParams.update({'font.size':13, 'axes.labelsize':15, 'axes.titlesize':15, 'legend.fontsize':10})
    # Orientation plot: no-branch cases only
    fig, ax=plt.subplots(figsize=(8,5.6))
    for m in sorted([m for m in manifest if not m['branching']], key=lambda z: z['theta_deg']):
        b=bin_curve(read_curve(Path(m['case'])), args.bin_um, args.max_ext_um)
        if b.empty: continue
        ax.plot(b['bin_center_um'], b['K_median_MPa_sqrt_m'], marker='o', markersize=3, label=f"θ={m['theta_deg']:g}°")
    ax.set_xlabel('Crack extension, Δa (µm)'); ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
    ax.set_title('Orientation sweep, no branching')
    ax.grid(False); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out/'orientation_sweep_no_branch_Rcurves.png', dpi=320); fig.savefig(out/'orientation_sweep_no_branch_Rcurves.svg'); plt.close(fig)

    # Branch comparison at selected theta
    fig, ax=plt.subplots(figsize=(8,5.6))
    branch_cases=sorted(manifest, key=lambda z: (z['theta_deg'], z['branching']))
    # prefer cases at branch theta if present; otherwise all branch/no branch pairs
    btheta = [m['theta_deg'] for m in manifest if m['branching']]
    selected_theta = btheta[0] if btheta else None
    for m in branch_cases:
        if selected_theta is not None and abs(m['theta_deg'] - selected_theta) > 1e-9:
            continue
        b=bin_curve(read_curve(Path(m['case'])), args.bin_um, args.max_ext_um)
        if b.empty: continue
        label=f"θ={m['theta_deg']:g}°, " + ('weak branching' if m['branching'] else 'no branching')
        ax.plot(b['bin_center_um'], b['K_median_MPa_sqrt_m'], marker='o', markersize=3, linestyle='--' if m['branching'] else '-', label=label)
    ax.set_xlabel('Crack extension, Δa (µm)'); ax.set_ylabel(r'$K_J$ (MPa$\sqrt{m}$)')
    ax.set_title('Branching comparison')
    ax.grid(False); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out/'branching_comparison_Rcurves.png', dpi=320); fig.savefig(out/'branching_comparison_Rcurves.svg'); plt.close(fig)
    print(f'WROTE {out}')

if __name__ == '__main__': main()
