#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path
import numpy as np
import pandas as pd

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',default='');a=p.parse_args()
 root=Path(a.root);out=Path(a.out) if a.out else root/'plots';out.mkdir(parents=True,exist_ok=True)
 rows=[]
 for fn in root.rglob('mixed_mode_first_passage_summary.json'):
  d=json.loads(fn.read_text());meta=fn.parent/'campaign_metadata.json'
  if meta.exists():d.update(json.loads(meta.read_text()))
  d['case_dir']=str(fn.parent);rows.append(d)
 if not rows:raise SystemExit(f'no mixed-mode summaries under {root}')
 df=pd.DataFrame(rows);df.to_csv(out/'mixed_mode_first_passage_all_cases.csv',index=False)
 keys=['class','target_psi_deg'];agg=df.groupby(keys).agg(
  n=('solver_seed','count'),KI_mean=('KI_first_MPa_sqrt_m','mean'),KI_sd=('KI_first_MPa_sqrt_m','std'),
  KII_mean=('KII_first_MPa_sqrt_m','mean'),KII_sd=('KII_first_MPa_sqrt_m','std'),
  Kopen_mean=('Kopen_maxhoop_first_MPa_sqrt_m','mean'),Kopen_sd=('Kopen_maxhoop_first_MPa_sqrt_m','std'),
  psi_mean=('mode_phase_first_deg','mean'),psi_sd=('mode_phase_first_deg','std'),kink_mean=('maxhoop_kink_first_deg','mean'),
  N_em_mean=('N_em_final','mean')).reset_index()
 agg.to_csv(out/'mixed_mode_first_passage_grouped.csv',index=False)
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 classes=sorted(df['class'].dropna().unique())
 fig,ax=plt.subplots(figsize=(6.2,5.4))
 for c in classes:
  q=df[df['class']==c];g=agg[agg['class']==c].sort_values('target_psi_deg')
  ax.scatter(q.KI_first_MPa_sqrt_m,q.KII_first_MPa_sqrt_m,alpha=.35,label=f'{c} seeds')
  ax.plot(g.KI_mean,g.KII_mean,'o-',label=f'{c} mean')
 ax.axhline(0,lw=.8);ax.axvline(0,lw=.8);ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.25)
 ax.set(xlabel=r'$K_I$ [MPa$\sqrt{m}$]',ylabel=r'$K_{II}$ [MPa$\sqrt{m}$]',title='Mixed-mode first-passage envelope');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'KI_KII_first_passage_envelope.png',dpi=200);plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4.8))
 for c in classes:
  q=df[df['class']==c];g=agg[agg['class']==c].sort_values('target_psi_deg')
  ax.scatter(q.target_psi_deg,q.Kopen_maxhoop_first_MPa_sqrt_m,alpha=.3)
  ax.errorbar(g.target_psi_deg,g.Kopen_mean,yerr=g.Kopen_sd,marker='o',capsize=3,label=c)
 ax.grid(alpha=.25);ax.set(xlabel='Target phase angle [deg]',ylabel=r'Opening drive at first passage [MPa$\sqrt{m}$]',title='First-passage resistance versus mode mixity');ax.legend();fig.tight_layout();fig.savefig(out/'Kopen_vs_mode_phase.png',dpi=200);plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4.8))
 for c in classes:
  g=agg[agg['class']==c].sort_values('target_psi_deg');ax.plot(g.psi_mean,g.kink_mean,'o-',label=c)
 ax.grid(alpha=.25);ax.set(xlabel='Achieved FEM phase angle [deg]',ylabel='Maximum-hoop kink angle [deg]',title='Predicted first crack direction');ax.legend();fig.tight_layout();fig.savefig(out/'kink_angle_vs_mode_phase.png',dpi=200);plt.close(fig)
 print('wrote',out)
if __name__=='__main__':main()
