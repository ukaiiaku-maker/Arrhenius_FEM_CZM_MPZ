#!/usr/bin/env python3
"""Robust elastic calibration of boundary loading angle to KI/KII phase angle."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np


def vals(text):
    return [float(x) for x in str(text).replace(',', ' ').split() if x]


def angle_error(a, b):
    return (float(a)-float(b)+180.0) % 360.0 - 180.0


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out', default='runs/mixed_mode_fem_czm_v2_elastic_calibration')
    p.add_argument('--target-psi-deg', default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--alpha-grid-deg', default='-70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70')
    p.add_argument('--U-cal-m', type=float, default=2e-6)
    p.add_argument('--nx', type=int, default=18); p.add_argument('--ny', type=int, default=36)
    p.add_argument('--tip-h-fine', type=float, default=5e-6); p.add_argument('--tip-ratio', type=float, default=1.30)
    p.add_argument('--mesh-seed', type=int, default=42)
    p.add_argument('--mode-projection-angle-deg', type=float, default=105.0)
    p.add_argument('--mode-projection-damage-cutoff', type=float, default=0.85)
    p.add_argument('--crystal-aniso', action='store_true')
    p.add_argument('--crystal-theta-deg', type=float, default=45.0)
    p.add_argument('--crystal-C11-GPa', type=float, default=None)
    p.add_argument('--crystal-C12-GPa', type=float, default=None)
    p.add_argument('--crystal-C44-GPa', type=float, default=None)
    p.add_argument('--psi-tol-deg', type=float, default=0.5)
    p.add_argument('--max-refine-iters', type=int, default=5)
    a=p.parse_args()

    from arrhenius_fracture.config import make_emergent_config
    from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
    from arrhenius_fracture.fem import plane_strain_D, assemble_mechanics, solve_dirichlet
    from arrhenius_fracture.mixed_mode_first_passage_v2 import (
        MixedModeContext, _mixed_solve_factory, project_near_tip_modes, MODEL_ID)

    cfg=make_emergent_config(); cfg.mesh.nx=a.nx; cfg.mesh.ny=a.ny
    cfg.mesh.tip_h_fine=a.tip_h_fine; cfg.mesh.tip_ratio=a.tip_ratio
    mesh=make_tri_mesh(cfg.geometry,cfg.mesh,seed=a.mesh_seed)
    bnd=make_boundary_data(mesh,cfg.geometry); mat=cfg.material
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
    htip=float(mesh.hbar_tip or mesh.hbar)

    cache={}
    def evaluate(alpha):
        key=round(float(alpha),10)
        if key in cache: return dict(cache[key])
        ctx=MixedModeContext(float(alpha),solver_seed=1,stochastic_first_passage=False,
                             projection_angle_deg=a.mode_projection_angle_deg,
                             projection_damage_cutoff=a.mode_projection_damage_cutoff)
        solve=_mixed_solve_factory(solve_dirichlet,ctx)
        u,F=solve(K,R,u0.copy(),bnd,0.5*a.U_cal_m,-0.5*a.U_cal_m)
        _K2,_R2,sig,_seq,_s1,_psi=assemble_mechanics(mesh,u,ep,rho,d,D,mat)
        md=project_near_tip_modes(mesh,sig,d,np.array([cfg.geometry.a0,0.0]),
                                  np.array([1.0,0.0]),1.25*htip,6.0*htip,
                                  a.mode_projection_angle_deg,a.mode_projection_damage_cutoff)
        KI=md['KI_Pa_sqrt_m']; KII=md['KII_Pa_sqrt_m']
        row={'loading_angle_deg':float(alpha),'achieved_psi_deg':md['mode_phase_deg'],
             'KI_MPa_sqrt_m':KI/1e6,'KII_MPa_sqrt_m':KII/1e6,
             'mode_ratio_KII_over_KI':KII/KI if np.isfinite(KI) and abs(KI)>1e-30 else np.nan,
             'generalized_reaction_N':F,'projection_n':md['mode_projection_n'],
             'projection_fit_count':md['mode_projection_fit_count'],
             'projection_rel_rmse':md['mode_projection_rel_rmse'],
             'projection_psi_spread_deg':md['mode_projection_psi_spread_deg'],
             'projection_reliable':md['mode_projection_reliable']}
        cache[key]=dict(row); return row

    grid=[evaluate(x) for x in vals(a.alpha_grid_deg)]
    valid=[r for r in grid if np.isfinite(r['achieved_psi_deg'])]
    if len(valid)<3: raise SystemExit('insufficient valid elastic mode projections')
    valid.sort(key=lambda r:r['loading_angle_deg'])

    mappings=[]
    for target in vals(a.target_psi_deg):
        qt=math.tan(math.radians(target))
        pairs=[]
        for lo,hi in zip(valid[:-1],valid[1:]):
            e0=lo['mode_ratio_KII_over_KI']-qt; e1=hi['mode_ratio_KII_over_KI']-qt
            if np.isfinite(e0) and np.isfinite(e1) and e0*e1<=0 and e1!=e0:
                alpha=lo['loading_angle_deg']-e0*(hi['loading_angle_deg']-lo['loading_angle_deg'])/(e1-e0)
                pairs.append((abs(alpha),alpha))
        if pairs:
            alpha=min(pairs)[1]
        else:
            alpha=min(valid,key=lambda r:abs(r['mode_ratio_KII_over_KI']-qt))['loading_angle_deg']
        hist=[]
        for it in range(a.max_refine_iters):
            r=evaluate(alpha); err=angle_error(r['achieved_psi_deg'],target); hist.append((alpha,err,r))
            if r['projection_reliable'] and abs(err)<=a.psi_tol_deg: break
            da=1.0
            rp=evaluate(min(alpha+da,75.0)); rm=evaluate(max(alpha-da,-75.0))
            slope=angle_error(rp['achieved_psi_deg'],rm['achieved_psi_deg'])/max(rp['loading_angle_deg']-rm['loading_angle_deg'],1e-12)
            if not np.isfinite(slope) or abs(slope)<0.1: slope=np.sign(slope or 1.0)*0.1
            alpha=float(np.clip(alpha-err/slope,-75.0,75.0))
        best=min(hist,key=lambda z:(not bool(z[2]['projection_reliable']),abs(z[1])))
        alpha,err,r=best
        rp=evaluate(min(alpha+1.0,75.0)); rm=evaluate(max(alpha-1.0,-75.0))
        slope=angle_error(rp['achieved_psi_deg'],rm['achieved_psi_deg'])/max(rp['loading_angle_deg']-rm['loading_angle_deg'],1e-12)
        row=dict(r); row.update(target_psi_deg=target,psi_error_deg=err,
                                dpsi_dalpha=slope,elastic_converged=bool(r['projection_reliable'] and abs(err)<=a.psi_tol_deg))
        mappings.append(row)

    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    for name,rows in [('elastic_alpha_grid_v2.csv',grid),('mixed_mode_loading_calibration_v2.csv',mappings)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'mixed_mode_loading_calibration_v2.json').write_text(json.dumps({'model':MODEL_ID,'settings':vars(a),'mapping':mappings},indent=2))
    try:
        import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.8,5.0))
        ax.plot([r['loading_angle_deg'] for r in valid],[r['achieved_psi_deg'] for r in valid],'o-',label='elastic grid')
        ax.scatter([r['loading_angle_deg'] for r in mappings],[r['achieved_psi_deg'] for r in mappings],marker='s',label='refined targets')
        ax.set(xlabel='Boundary loading angle alpha [deg]',ylabel='Projected phase angle psi [deg]',title='v2 robust elastic initialization')
        ax.grid(alpha=.3);ax.legend();fig.tight_layout();fig.savefig(out/'mixed_mode_loading_calibration_v2.png',dpi=200);plt.close(fig)
    except Exception as exc: print('plot skipped:',exc)
    print('wrote',out/'mixed_mode_loading_calibration_v2.csv')

if __name__=='__main__':main()
