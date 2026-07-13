#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _load(root: Path):
    p = root/'sn_pf2d_fullplastic_summary.csv'
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_csv(p)


def _plot_summary(df, root):
    metrics = [
        ('B_nuc_final_max','Final max nucleation clock $B_n$','SN_FULLPLASTIC_Bnuc_final_vs_stress.png','log'),
        ('epsp_acc_final_max','Max accumulated equivalent plastic strain','SN_FULLPLASTIC_epsp_vs_stress.png','log'),
        ('rho_final_max_m2','Max dislocation density (m$^{-2}$)','SN_FULLPLASTIC_rho_vs_stress.png','log'),
        ('residual_sigma1_final_max_Pa','Max residual principal stress (Pa)','SN_FULLPLASTIC_residual_stress_vs_stress.png','linear'),
    ]
    for col,ylabel,name,yscale in metrics:
        if col not in df.columns: continue
        fig,ax=plt.subplots(figsize=(6.8,4.5))
        for case,g in df.groupby('case',sort=False):
            g=g.sort_values('sigma_a_MPa')
            ax.plot(g.sigma_a_MPa,g[col],marker='o',label=case)
        if yscale=='log' and np.nanmax(df[col].to_numpy(float))>0: ax.set_yscale('log')
        ax.set_xlabel('Stress amplitude $\\sigma_a$ (MPa)'); ax.set_ylabel(ylabel); ax.legend(); ax.grid(True,alpha=.25)
        fig.tight_layout(); fig.savefig(root/name,dpi=240); plt.close(fig)

    if {'root_radius_initial_m','root_radius_final_m'}.issubset(df.columns):
        fig,ax=plt.subplots(figsize=(6.8,4.5))
        for case,g in df.groupby('case',sort=False):
            g=g.sort_values('sigma_a_MPa')
            ratio=g.root_radius_final_m/g.root_radius_initial_m
            ax.plot(g.sigma_a_MPa,ratio,marker='o',label=case)
        ax.axhline(1.0,lw=1,ls='--')
        ax.set_xlabel('Stress amplitude $\\sigma_a$ (MPa)'); ax.set_ylabel('Final root radius / initial root radius')
        ax.legend(); ax.grid(True,alpha=.25); fig.tight_layout(); fig.savefig(root/'SN_FULLPLASTIC_root_sharpening_vs_stress.png',dpi=240); plt.close(fig)


def _history_plots(root: Path):
    histories=[]
    for p in root.glob('*/*/sn_pf2d_fullplastic_history.csv'):
        try:
            d=pd.read_csv(p)
            if len(d):
                histories.append((p.parts[-3],p.parts[-2],d))
        except Exception: pass
    stresses=sorted(set(tag for _,tag,_ in histories))
    for tag in stresses:
        subset=[x for x in histories if x[1]==tag]
        if not subset: continue
        fig,axes=plt.subplots(5,1,figsize=(7.2,10.5),sharex=True)
        series=[
            ('B_nuc_max','$B_{nuc,max}$','log'),
            ('epsp_acc_max','max accumulated $\\epsilon_p$','log'),
            ('rho_max_m2','max $\\rho$ (m$^{-2}$)','log'),
            ('residual_sigma1_max_Pa','max residual $\\sigma_1$ (Pa)','linear'),
            ('root_radius_over_initial','$r_{root}/r_{root,0}$','linear'),
        ]
        for case,_,d in subset:
            for ax,(col,ylabel,scale) in zip(axes,series):
                if col in d: ax.plot(d.cycles_total,d[col],label=case)
                ax.set_ylabel(ylabel)
                if scale=='log': ax.set_yscale('log')
                ax.grid(True,alpha=.2)
        axes[-1].set_xlabel('Cycles'); axes[-1].set_xscale('log')
        axes[0].legend(); fig.suptitle(tag.replace('_',' ')); fig.tight_layout(); fig.savefig(root/f'SN_FULLPLASTIC_history_{tag}.png',dpi=220); plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='runs/sn_pf2d_fullplastic_two_case')
    args=ap.parse_args(); root=Path(args.root)
    df=_load(root); _plot_summary(df,root); _history_plots(root)
    print(f'wrote full-plastic S-N diagnostics under {root}')

if __name__=='__main__': main()
