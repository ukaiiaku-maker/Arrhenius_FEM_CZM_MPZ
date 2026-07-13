#!/usr/bin/env python3
"""Reliability- and censor-aware plots for mixed-mode FEM/CZM v2."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);p.add_argument('--out',default=None,type=Path);a=p.parse_args()
    out=a.out if a.out is not None else a.root/'plots_v2';out.mkdir(parents=True,exist_ok=True)
    paths=list(a.root.glob('*/*/mixed_mode_v2_final_summary.csv'))
    if not paths:raise SystemExit(f'no v2 final summaries found under {a.root}')
    df=pd.concat([pd.read_csv(x) for x in paths],ignore_index=True)
    df['event_observed']=df['mode_classification'].eq('brittle')
    df['usable']=df['event_state_control_converged'].astype(bool)&df['projection_reliable'].astype(bool)
    df.to_csv(out/'mixed_mode_v2_plot_data.csv',index=False)
    usable=df[df['usable']];bad=df[~df['usable']]

    fig,ax=plt.subplots(figsize=(8,7))
    for cls,g in usable.groupby('class'):
        ev=g[g['event_observed']];ce=g[~g['event_observed']]
        if len(ev):ax.plot(ev['KI_first_MPa_sqrt_m'],ev['KII_first_MPa_sqrt_m'],'o-',label=f'{cls}: first passage')
        if len(ce):ax.plot(ce['KI_first_MPa_sqrt_m'],ce['KII_first_MPa_sqrt_m'],'^--',label=f'{cls}: censored endpoint')
    if len(bad):ax.scatter(bad['KI_first_MPa_sqrt_m'],bad['KII_first_MPa_sqrt_m'],marker='x',label='unreliable/not converged')
    ax.axhline(0,lw=.8);ax.axvline(0,lw=.8);ax.set_xlabel(r'$K_I$ [MPa$\sqrt{m}$]');ax.set_ylabel(r'$K_{II}$ [MPa$\sqrt{m}$]');ax.set_title('Event-controlled mixed-mode envelope');ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'01_KI_KII_event_controlled.png',dpi=220);plt.close(fig)

    fig,ax=plt.subplots(figsize=(7,6))
    for cls,g in df.groupby('class'):
        ax.plot(g['target_psi_deg'],g['achieved_psi_deg'],'o-',label=cls)
    lim=max(65,float(np.nanmax(np.abs(df[['target_psi_deg','achieved_psi_deg']].to_numpy()))));ax.plot([-lim,lim],[-lim,lim],':',label='target = achieved')
    ax.set_xlabel('Target phase angle [deg]');ax.set_ylabel('Achieved event/endpoint phase angle [deg]');ax.set_title('Closed-loop phase-angle control');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/'02_target_vs_achieved_phase.png',dpi=220);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,5.8))
    for cls,g in usable.groupby('class'):
        ev=g[g['event_observed']];ce=g[~g['event_observed']]
        if len(ev):ax.plot(ev['target_psi_deg'],ev['Kopen_maxhoop_first_MPa_sqrt_m'],'o-',label=f'{cls}: first passage')
        if len(ce):ax.plot(ce['target_psi_deg'],ce['Kopen_maxhoop_first_MPa_sqrt_m'],'^--',label=f'{cls}: censored lower bound')
    ax.set_xlabel('Target phase angle [deg]');ax.set_ylabel(r'$K_{open}$ [MPa$\sqrt{m}$]');ax.set_title('Opening drive versus controlled mixity');ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'03_Kopen_vs_controlled_phase.png',dpi=220);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8.5,5.8))
    for cls,g in df.groupby('class'):
        ax.plot(g['target_psi_deg'],g['projection_psi_spread_deg'],'o-',label=cls)
    ax.axhline(12,ls='--',lw=1,label='reliability limit');ax.set_xlabel('Target phase angle [deg]');ax.set_ylabel('Annulus-to-annulus phase spread [deg]');ax.set_title('Mode-decomposition reliability');ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/'04_projection_reliability.png',dpi=220);plt.close(fig)
    print('wrote',out)
if __name__=='__main__':main()
