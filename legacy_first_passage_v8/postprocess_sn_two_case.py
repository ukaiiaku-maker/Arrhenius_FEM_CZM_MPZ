#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _life_plot(df, out_png: Path, life_col: str, total_col: str, title: str):
    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    for case, g in df.groupby('case', sort=False):
        g = g.sort_values('sigma_a_MPa')
        life = pd.to_numeric(g.get(life_col, np.nan), errors='coerce')
        total = pd.to_numeric(g[total_col], errors='coerce')
        cens = life.isna()
        y = life.where(~cens, total)
        if (~cens).any():
            ax.plot(y[~cens], g.loc[~cens, 'sigma_a_MPa'], 'o-', label=case)
        else:
            # keep legend entry even if all points are censored
            ax.plot([], [], 'o-', label=case)
        if cens.any():
            ax.scatter(y[cens], g.loc[cens, 'sigma_a_MPa'], marker='>', s=55)
    ax.set_xscale('log')
    ax.set_xlabel('cycles')
    ax.set_ylabel(r'nominal stress amplitude $\sigma_a$ (MPa)')
    ax.set_title(title)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()
    fig.tight_layout(); fig.savefig(out_png, dpi=240); plt.close(fig)


def _hydrate_pf_final_state(df: pd.DataFrame, root: Path) -> pd.DataFrame:
    df = df.copy()
    needed = ['B_nuc_final_max','B_emit_final_max','P_final_max','Dloc_final_max','d_final_max']
    for c in needed:
        if c not in df.columns:
            df[c] = np.nan
    for i, r in df.iterrows():
        hist = r.get('history_csv', '')
        candidates = []
        if isinstance(hist, str) and hist:
            p = Path(hist)
            candidates += [p, root / p]
        # canonical local layout
        tag = f"sigmaA_{float(r['sigma_a_MPa']):g}MPa".replace('.', 'p')
        candidates.append(root / str(r['case']) / tag / 'sn_pf2d_history.csv')
        hp = next((p for p in candidates if p.exists()), None)
        if hp is None:
            continue
        try:
            h = pd.read_csv(hp)
        except Exception:
            continue
        if h.empty:
            continue
        last = h.iloc[-1]
        mapping = {
            'B_nuc_final_max':'B_nuc_max',
            'B_emit_final_max':'B_emit_max',
            'P_final_max':'P_max',
            'Dloc_final_max':'Dloc_max',
            'd_final_max':'d_max',
        }
        for dst, src in mapping.items():
            if pd.isna(df.at[i,dst]) and src in h.columns:
                df.at[i,dst] = last[src]
    return df


def _state_plot(df: pd.DataFrame, out_png: Path, col: str, ylabel: str, title: str, logy=True):
    if col not in df.columns or df[col].notna().sum() == 0:
        return
    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    for case, g in df.groupby('case', sort=False):
        g = g.sort_values('sigma_a_MPa')
        ax.plot(g['sigma_a_MPa'], g[col], 'o-', label=case)
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel(r'nominal stress amplitude $\sigma_a$ (MPa)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(); fig.tight_layout(); fig.savefig(out_png, dpi=240); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--v1-root', default='runs/sn_v1_two_case')
    ap.add_argument('--pf2d-root', default='runs/sn_pf2d_two_case')
    args = ap.parse_args()

    vroot = Path(args.v1_root)
    v1 = vroot / 'sn_v1_summary.csv'
    if v1.exists():
        vdf = pd.read_csv(v1)
        _life_plot(vdf, vroot/'SN_V1_two_case.png', 'cycles_to_nucleation', 'cycles_total', 'V1 S-N crack-initiation pilot')
        print(f"wrote {vroot/'SN_V1_two_case.png'}")

    proot = Path(args.pf2d_root)
    p2 = proot / 'sn_pf2d_summary.csv'
    if p2.exists():
        pdf = _hydrate_pf_final_state(pd.read_csv(p2), proot)
        pdf.to_csv(proot/'sn_pf2d_summary_audited.csv', index=False)
        _life_plot(pdf, proot/'SN_PF2D_two_case.png', 'cycles_to_pf_crack', 'cycles_total', '2-D PF S-N blunt-notch pilot — connected crack')
        _life_plot(pdf, proot/'SN_PF2D_nucleation_clock_two_case.png', 'cycles_to_nucleation_clock', 'cycles_total', '2-D PF S-N blunt-notch pilot — nucleation clock')
        _state_plot(pdf, proot/'SN_PF2D_Bnuc_final_vs_stress.png', 'B_nuc_final_max', r'final $\max B_{nuc}$', '2-D crack-opening clock separation', logy=True)
        _state_plot(pdf, proot/'SN_PF2D_dmax_final_vs_stress.png', 'd_final_max', r'final $d_{max}$', '2-D phase-field damage separation', logy=True)
        print(f"wrote {proot/'SN_PF2D_two_case.png'}")
        print(f"wrote {proot/'SN_PF2D_nucleation_clock_two_case.png'}")
        print(f"wrote {proot/'SN_PF2D_Bnuc_final_vs_stress.png'}")
        print(f"wrote {proot/'SN_PF2D_dmax_final_vs_stress.png'}")
        print(f"wrote {proot/'sn_pf2d_summary_audited.csv'}")

if __name__ == '__main__':
    main()
