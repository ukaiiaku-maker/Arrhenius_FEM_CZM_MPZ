#!/usr/bin/env python3
"""Plot and audit the fully Arrhenius V1 two-case S-N initiation sweep."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='runs/sn_v1_arrhenius_two_case')
    args=ap.parse_args()
    root=Path(args.root)
    p=root/'sn_v1_arrhenius_summary.csv'
    if not p.exists():
        raise FileNotFoundError(p)
    df=pd.read_csv(p)

    fig,ax=plt.subplots(figsize=(7.6,5.5))
    for case,g in df.groupby('case',sort=False):
        g=g.sort_values('sigma_a_MPa')
        fail=g[g.status=='failed'].copy()
        cens=g[g.status=='right_censored'].copy()
        other=g[~g.status.isin(['failed','right_censored'])].copy()
        if not fail.empty:
            ax.plot(fail.cycles_to_nucleation,fail.sigma_a_MPa,marker='o',lw=1.8,label=case)
        if not cens.empty:
            ax.scatter(cens.cycles_total,cens.sigma_a_MPa,marker='>',s=55,label=f'{case} censored')
        if not other.empty:
            ax.scatter(other.cycles_total,other.sigma_a_MPa,marker='s',facecolors='none',s=50,label=f'{case} unresolved')
    ax.set_xscale('log')
    ax.set_xlabel('Cycles to crack nucleation / censoring horizon')
    ax.set_ylabel('Stress amplitude (MPa)')
    ax.set_title('V1 fully Arrhenius S–N initiation response')
    ax.grid(True,which='both',alpha=.3)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(root/'SN_V1_Arrhenius_two_case.png',dpi=240); plt.close(fig)

    metrics=[('P_final','Final shielding state P'),('epsp_acc_final','Accumulated plastic strain'),('rho_final_m2','Final dislocation density (m$^{-2}$)'),('B_nuc_final','Final nucleation clock')]
    for col,ylab in metrics:
        if col not in df.columns: continue
        fig,ax=plt.subplots(figsize=(7.2,5.2))
        for case,g in df.groupby('case',sort=False):
            g=g.sort_values('sigma_a_MPa')
            ax.plot(g.sigma_a_MPa,g[col],marker='o',label=case)
        if col in ('rho_final_m2','B_nuc_final'):
            ax.set_yscale('log')
        ax.set_xlabel('Stress amplitude (MPa)'); ax.set_ylabel(ylab)
        ax.set_title(ylab + ' versus stress amplitude')
        ax.grid(True,which='both',alpha=.3); ax.legend(); fig.tight_layout()
        fig.savefig(root/f'SN_V1_{col}_vs_stress.png',dpi=220); plt.close(fig)

    audit=(df.groupby('case').agg(n_points=('case','size'),n_failed=('status',lambda s:(s=='failed').sum()),n_censored=('status',lambda s:(s=='right_censored').sum()),min_failed_cycles=('cycles_to_nucleation','min'),max_failed_cycles=('cycles_to_nucleation','max'),max_P=('P_final','max'),max_rho=('rho_final_m2','max')).reset_index())
    audit.to_csv(root/'SN_V1_two_case_audit.csv',index=False)
    print(f'wrote {root / "SN_V1_Arrhenius_two_case.png"}')
    print(f'wrote {root / "SN_V1_two_case_audit.csv"}')

if __name__=='__main__':
    main()
