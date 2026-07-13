#!/usr/bin/env python3
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();f=a.root/'mixed_mode_v4_anisotropic_all_cases.csv'
if not f.exists():raise SystemExit(f'not found: {f}')
d=pd.read_csv(f);o=a.root/'plots_v4';o.mkdir(parents=True,exist_ok=True)
for col,ylabel,name in [('KJ_reference_first_MPa_sqrt_m',r'$K_{J,ref}$ [MPa$\sqrt{m}$]','01_KJref_vs_phase'),('sigma_cleave_drive_first_GPa','cleavage drive [GPa]','02_cleavage_drive_vs_phase'),('candidate_angle_first_deg','selected crack direction [deg]','03_direction_vs_phase'),('traction_phase_error_first_deg','event phase error [deg]','04_event_phase_error')]:
 fig,ax=plt.subplots(figsize=(8,5.5))
 for k,g in d.groupby('class'):
  g=g.sort_values('target_psi_deg');ax.plot(g.target_psi_deg,g[col],marker='o',label=k)
 ax.set_xlabel('Target process-zone traction phase [deg]');ax.set_ylabel(ylabel);ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(o/(name+'.png'),dpi=220);plt.close(fig)
print('wrote',o)
