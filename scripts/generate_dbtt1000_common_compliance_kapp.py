#!/usr/bin/env python3
"""Common experimental-style K/F calibration from global elastic compliance.

For one fixed specimen mesh, progressively extend a straight, thin crack wake
from the production initial notch.  No MPZ, hazard, plasticity, backstress or
solver-specific evolved tip geometry is present.  At each projected crack
length solve the same small symmetric fixed-grip elastic problem and record
C=U/F.  Then use the compliance method

    G/F^2 = (1/(2 B)) dC/da
    K/F   = sqrt(E'/(2 B) * dC/da)

with unit out-of-plane thickness B=1 m, consistent with the 2-D force ledger.
This defines the common reference used to convert both actual PF and FEM/CZM
remote forces into an apparent experimental-style K_app.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np

from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
from arrhenius_fracture.fem import assemble_mechanics, solve_dirichlet
from arrhenius_fracture.crystal import cubic_plane_strain_D

C11=523.0e9; C12=203.0e9; C44=160.0e9; THETA=0.0
A0=0.5e-3; LX=2.0e-3; LY=4.0e-3; BTHICK=1.0
UAPP=1.0e-6; RHO0=5.0e12
DA=25.0e-6
A_VALUES=np.arange(A0, 1.525e-3+0.5*DA, DA)


def damage_for_a(mesh, geom, a):
    x=mesh.nodes[:,0]; y=mesh.nodes[:,1]
    d=np.zeros(mesh.nn)
    d[(x<=geom.a0)&(np.abs(y)<=geom.notch_half_thickness)]=1.0
    if a<=geom.a0+1e-15:
        return d
    c=mesh.nodes[mesh.elems].mean(axis=1)
    p0=np.array([geom.a0,0.0]); p1=np.array([a,0.0]); seg=p1-p0
    L2=float(seg@seg)+1e-30
    t=np.clip(((c[:,0]-p0[0])*seg[0]+(c[:,1]-p0[1])*seg[1])/L2,0.0,1.0)
    proj=p0[None,:]+t[:,None]*seg[None,:]
    dist2=np.sum((c-proj)**2,axis=1)
    erad=np.sqrt(np.maximum(mesh.area_e,1e-30))
    rad=np.maximum(0.7*erad, 0.5*min(LX/80.0,LY/160.0))
    mask=(dist2<=rad**2)&(c[:,0]>=geom.a0-DA)&(c[:,0]<=a+DA)
    d[mesh.elems[mask].ravel()]=1.0
    return d


def main():
    cfg=make_emergent_config()
    cfg.geometry.Lx=LX; cfg.geometry.Ly=LY; cfg.geometry.a0=A0
    cfg.mesh.nx=80; cfg.mesh.ny=160; cfg.mesh.jitter=0.0; cfg.mesh.tip_h_fine=0.0
    mat=cfg.material
    mesh=make_tri_mesh(cfg.geometry,cfg.mesh,seed=42)
    bnd=make_boundary_data(mesh,cfg.geometry)
    D=cubic_plane_strain_D(C11,C12,C44,THETA)
    ep=np.zeros((3,mesh.ne)); rho=np.full(mesh.ne,RHO0)

    rows=[]
    for a in A_VALUES:
        d=damage_for_a(mesh,cfg.geometry,float(a))
        u=np.zeros(mesh.ndof)
        Kmat,Rint,sigma,seq,s1,psi=assemble_mechanics(mesh,u,ep,rho,d,D,mat,cohesive_network=None)
        u,F=solve_dirichlet(Kmat,Rint,u,bnd,0.5*UAPP,-0.5*UAPP)
        C=UAPP/abs(float(F))
        rows.append({"projected_crack_length_m":float(a),"projected_extension_m":float(a-A0),"Fabs_N":abs(float(F)),"Uapp_m":UAPP,"compliance_m_per_N":C})
        print(f"a={a*1e6:8.3f} um F={abs(float(F)):.9e} C={C:.9e}")

    a=np.asarray([r["projected_crack_length_m"] for r in rows])
    C=np.asarray([r["compliance_m_per_N"] for r in rows])
    if np.any(np.diff(C)<-1e-16*np.max(C)):
        bad=np.where(np.diff(C)<0)[0]
        raise RuntimeError(f"reference compliance is non-monotone at indices {bad.tolist()}")
    dCda=np.gradient(C,a,edge_order=2)
    if np.any(dCda<=0):
        bad=np.where(dCda<=0)[0]
        raise RuntimeError(f"non-positive compliance derivative at indices {bad.tolist()}")
    Eprime=mat.E/(1.0-mat.nu**2)
    kpf=np.sqrt(Eprime*dCda/(2.0*BTHICK))
    for i,r in enumerate(rows):
        r["dC_da_1_per_N"]=float(dCda[i])
        r["K_over_F_Pa_sqrtm_per_N"]=float(kpf[i])
        r["K_reference_at_solved_F_Pa_sqrtm"]=float(kpf[i]*r["Fabs_N"])

    out=Path("dbtt1000_common_compliance_reference"); out.mkdir(exist_ok=True)
    with (out/"dbtt1000_common_compliance_K_over_F.csv").open("w",newline="") as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    manifest={
      "schema":"dbtt1000_common_compliance_reference_v1",
      "diagnostic_only":True,
      "formula":"K_app=F_actual*sqrt(Eprime/(2B)*dC/da)",
      "B_m":BTHICK,
      "Eprime_Pa":Eprime,
      "mesh":{"nx":80,"ny":160,"jitter":0.0,"fixed_for_all_crack_lengths":True},
      "geometry":{"Lx_m":LX,"Ly_m":LY,"a0_m":A0,"da_m":DA},
      "elasticity":{"C11_Pa":C11,"C12_Pa":C12,"C44_Pa":C44,"theta_deg":THETA},
      "state":"linear elastic; no MPZ, no plasticity, no hazard, no evolved microscopic tip state",
      "crack":"production initial notch plus straight thin post-notch wake",
      "compliance_monotone":True,
      "caveat":"Experimental-style common numerical calibration for this model specimen; not an ASTM validity claim."
    }
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
