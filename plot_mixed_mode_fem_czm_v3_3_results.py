#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args()
    f=a.root/'mixed_mode_v3_3_all_cases.csv'
    if not f.exists():raise SystemExit(f'missing {f}')
    d=pd.read_csv(f);out=a.root/'plots_v3_3';out.mkdir(parents=True,exist_ok=True)
    event=d[d.status=='event'];cens=d[d.status=='right_censored']
    fig,ax=plt.subplots(figsize=(7.3,6.3))
    for cls,g in event.groupby('class'):
        ax.plot(g.KI_first_MPa_sqrt_m,g.KII_first_MPa_sqrt_m,'o-',label=f'{cls}: event')
    for cls,g in cens.groupby('class'):
        ax.plot(g.KI_first_MPa_sqrt_m,g.KII_first_MPa_sqrt_m,'^--',label=f'{cls}: censored endpoint')
    ax.axhline(0,lw=.8);ax.axvline(0,lw=.8);ax.set_xlabel(r'$K_I$ [MPa$\sqrt{m}$]');ax.set_ylabel(r'$K_{II}$ [MPa$\sqrt{m}$]');ax.set_title('J-consistent circular-phase calibrated response');ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'01_KI_KII_J_consistent.png',dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.2,5.6))
    for cls,g in event.groupby('class'):
        ax.plot(g.target_psi_deg,g.Kopen_maxhoop_first_MPa_sqrt_m,'o-',label=f'{cls}: first passage')
    for cls,g in cens.groupby('class'):
        ax.plot(g.target_psi_deg,g.Kopen_maxhoop_first_MPa_sqrt_m,'^--',label=f'{cls}: censored lower bound')
    ax.set_xlabel('Target phase angle [deg]');ax.set_ylabel(r'$K_{open}$ [MPa$\sqrt{m}$]');ax.set_title('Opening drive at event or censoring');ax.grid(alpha=.25)
    if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'02_Kopen_vs_phase.png',dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.2,5.6))
    for cls,g in d.groupby('class'):
        ax.plot(g.target_psi_deg,g.calibrated_loading_angle_deg,'o-',label=cls)
    ax.set_xlabel('Target phase angle [deg]');ax.set_ylabel('Calibrated boundary angle [deg]');ax.set_title('Elastic loading calibration');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'03_loading_angle_calibration.png',dpi=220);plt.close(fig)
    d.to_csv(out/'mixed_mode_v3_3_plot_data.csv',index=False)
    print('wrote',out)
if __name__=='__main__':main()
