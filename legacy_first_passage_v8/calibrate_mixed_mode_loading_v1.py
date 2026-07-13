#!/usr/bin/env python3
"""Elastic calibration of boundary loading angle alpha to FEM SIF phase angle psi."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np


def parse_floats(text):
    return [float(x) for x in str(text).replace(',', ' ').split() if x]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out', default='runs/mixed_mode_fem_czm_v1_calibration')
    p.add_argument('--target-psi-deg', default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--alpha-grid-deg', default='-75 -65 -55 -45 -35 -25 -15 -7.5 0 7.5 15 25 35 45 55 65 75')
    p.add_argument('--U-cal-m', type=float, default=2e-7)
    p.add_argument('--nx', type=int, default=18); p.add_argument('--ny', type=int, default=36)
    p.add_argument('--tip-h-fine', type=float, default=5e-6)
    p.add_argument('--tip-ratio', type=float, default=1.30)
    p.add_argument('--mesh-seed', type=int, default=42)
    p.add_argument('--mode-projection-wedge-deg', type=float, default=18.0)
    p.add_argument('--rJ-outer', type=float, default=25e-6)
    p.add_argument('--crystal-aniso', action='store_true')
    p.add_argument('--crystal-theta-deg', type=float, default=0.0)
    p.add_argument('--crystal-C11-GPa', type=float, default=None)
    p.add_argument('--crystal-C12-GPa', type=float, default=None)
    p.add_argument('--crystal-C44-GPa', type=float, default=None)
    a=p.parse_args()

    from arrhenius_fracture.config import make_emergent_config
    from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
    from arrhenius_fracture.fem import plane_strain_D, assemble_mechanics, solve_dirichlet
    from arrhenius_fracture.j_integral import compute_J_integral, JIntegralConfig
    from arrhenius_fracture.mixed_mode_first_passage_v1 import (
        MixedModeContext, _mixed_solve_factory, project_near_tip_modes,
        maximum_hoop_drive, MODEL_ID)

    cfg=make_emergent_config(); cfg.mesh.nx=a.nx; cfg.mesh.ny=a.ny
    cfg.mesh.tip_h_fine=a.tip_h_fine; cfg.mesh.tip_ratio=a.tip_ratio
    mesh=make_tri_mesh(cfg.geometry,cfg.mesh,seed=a.mesh_seed)
    bnd=make_boundary_data(mesh,cfg.geometry)
    mat=cfg.material
    if a.crystal_aniso:
        from arrhenius_fracture.crystal import cubic_plane_strain_D, W_C11, W_C12, W_C44
        C11=(a.crystal_C11_GPa*1e9 if a.crystal_C11_GPa is not None else W_C11)
        C12=(a.crystal_C12_GPa*1e9 if a.crystal_C12_GPa is not None else W_C12)
        C44=(a.crystal_C44_GPa*1e9 if a.crystal_C44_GPa is not None else W_C44)
        D=cubic_plane_strain_D(C11,C12,C44,a.crystal_theta_deg)
    else:
        D=plane_strain_D(mat)
    u0=np.zeros(mesh.ndof); ep=np.zeros((3,mesh.ne)); rho=np.zeros(mesh.ne)
    d=np.zeros(mesh.nn); d[bnd.notch_nodes]=1.0
    K,R,*_=assemble_mechanics(mesh,u0,ep,rho,d,D,mat)
    jcfg=JIntegralConfig(r_inner_factor=2.0,r_outer_factor=8.0)
    # ell chosen so r_outer equals requested physical radius.
    ell=a.rJ_outer/max(jcfg.r_outer_factor,1e-30)

    def evaluate(alpha):
        ctx=MixedModeContext(alpha,solver_seed=1,stochastic_first_passage=False,
                             projection_wedge_deg=a.mode_projection_wedge_deg)
        solve=_mixed_solve_factory(solve_dirichlet,ctx)
        u,F=solve(K,R,u0.copy(),bnd,0.5*a.U_cal_m,-0.5*a.U_cal_m)
        _K2,_R2,sig,_seq,_s1,psi=assemble_mechanics(mesh,u,ep,rho,d,D,mat)
        J,KJ,ji=compute_J_integral(mesh,u,sig,psi,d,
            np.array([cfg.geometry.a0,0.0]),np.array([1.0,0.0]),mat,ell,cfg=jcfg)
        rmin=max(1.5*(mesh.hbar_tip or mesh.hbar),0.12*a.rJ_outer)
        rmax=max(0.65*a.rJ_outer,1.5*rmin)
        md=project_near_tip_modes(mesh,sig,d,np.array([cfg.geometry.a0,0.0]),
                                  np.array([1.0,0.0]),rmin,rmax,
                                  a.mode_projection_wedge_deg)
        Ko,kink=maximum_hoop_drive(md['KI_Pa_sqrt_m'],md['KII_Pa_sqrt_m'])
        return {'loading_angle_deg':alpha,'achieved_psi_deg':md['mode_phase_deg'],
                'KI_MPa_sqrt_m':md['KI_Pa_sqrt_m']/1e6,
                'KII_MPa_sqrt_m':md['KII_Pa_sqrt_m']/1e6,
                'KJ_MPa_sqrt_m':KJ/1e6,'Kopen_MPa_sqrt_m':Ko/1e6,
                'maxhoop_kink_deg':kink,'generalized_reaction_N':F,
                'projection_n':md['mode_projection_n']}

    samples=[evaluate(x) for x in parse_floats(a.alpha_grid_deg)]
    valid=[r for r in samples if np.isfinite(r['achieved_psi_deg'])]
    valid.sort(key=lambda r:r['achieved_psi_deg'])
    xp=np.array([r['achieved_psi_deg'] for r in valid]); fp=np.array([r['loading_angle_deg'] for r in valid])
    mapping=[]
    for target in parse_floats(a.target_psi_deg):
        alpha=float(np.interp(target,xp,fp))
        r=evaluate(alpha); r['target_psi_deg']=target
        r['psi_error_deg']=r['achieved_psi_deg']-target
        mapping.append(r)
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    for name, rows in [('elastic_alpha_grid.csv',samples),('mixed_mode_loading_calibration.csv',mapping)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'mixed_mode_loading_calibration.json').write_text(json.dumps({
        'model':MODEL_ID,'isotropic_projection':not a.crystal_aniso,
        'settings':vars(a),'mapping':mapping},indent=2))
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.4,4.8))
        ax.plot([r['loading_angle_deg'] for r in valid],[r['achieved_psi_deg'] for r in valid],'o-')
        ax.plot([-75,75],[-75,75],'--',lw=1,label='alpha = psi')
        ax.scatter([r['loading_angle_deg'] for r in mapping],[r['achieved_psi_deg'] for r in mapping],marker='s',label='calibrated targets')
        ax.set(xlabel='Boundary loading angle alpha [deg]',ylabel='FEM phase angle psi [deg]',title='Mixed-mode elastic calibration')
        ax.grid(alpha=.3);ax.legend();fig.tight_layout();fig.savefig(out/'mixed_mode_loading_calibration.png',dpi=180);plt.close(fig)
    except Exception as e: print('plot skipped:',e)
    print('wrote',out/'mixed_mode_loading_calibration.csv')

if __name__=='__main__': main()
