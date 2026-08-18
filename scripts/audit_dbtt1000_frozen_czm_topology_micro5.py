#!/usr/bin/env python3
"""Five-event frozen DBTT topology microgate using recovered |theta|."""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np

import audit_dbtt1000_frozen_czm_topology as full
import audit_dbtt1000_frozen_czm_topology_short100 as short
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.atomic_path_corridor_adaptive_quality_v10051839 import AdaptiveQualityAtomicPathCorridorCZMBackendV10051839
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import install as install_local_scale_gate, restore as restore_local_scale_gate
from arrhenius_fracture.atomic_path_corridor_czm_v10051839 import reset_audit

OUT=Path('dbtt1000_frozen_czm_topology_micro5')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ledger=json.loads(full.EVENT_LEDGER.read_text())
    lengths=np.asarray(ledger['event_lengths_m'],float)[:5]
    dx=short.DX[:5]
    cfg_a,mesh_a,bnd_a,d_a,u_a=full.initial_state()
    cfg_s,mesh_s,bnd_s,d_s,u_s=full.initial_state()
    startup={'nodes':int(mesh_a.nn),'triangles':int(mesh_a.ne),'hbar_tip_m':float(mesh_a.hbar_tip)}
    if startup['nodes']!=20321 or startup['triangles']!=40502:
        raise RuntimeError(f'production corridor mismatch: {startup}')
    reset_audit(); saved=install_local_scale_gate()
    adaptive=AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(geom=cfg_a.geometry,event_damage=1.0,max_angle_error_deg=35.0)
    sharp=SharpWakeBackend(); rows=[]; path=[]
    p_a=np.array([full.A0_M,0.0]); p_s=p_a.copy()
    def record(iev):
        ext=float((p_a[0]-full.A0_M)*1e6)
        for name,net in [('adaptive_czm_failed',adaptive.cohesive_network),('adaptive_split_no_network',None)]:
            r=full.frozen_solve(mesh_a,bnd_a,d_a,net,cfg_a); rows.append({'event':iev,'representation':name,'extension_um':ext,'failed_interfaces':adaptive.cohesive_network.failed_count(),**r})
        r=full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s); rows.append({'event':iev,'representation':'sharp_wake','extension_um':ext,'failed_interfaces':0,**r})
    record(0); failure=None
    try:
        for iev,(L,ddx) in enumerate(zip(lengths,dx),1):
            rr,p1,sign,trial_failures=short.attempt_signed_event(adaptive,mesh_a,bnd_a,d_a,u_a,p_a,float(ddx),float(L))
            if rr is None:
                failure={'event':iev,'p0_um':(p_a*1e6).tolist(),'dx_um':float(ddx*1e6),'ds_um':float(L*1e6),'trials':trial_failures}; break
            mesh_a,bnd_a,d_a,u_a=rr.mesh,rr.boundary,rr.damage,rr.displacement; p_a=p1
            dy=float(p_a[1]-p_s[1]); p1s=p_s+np.array([float(ddx),dy]); ds=float(np.linalg.norm(p1s-p_s))
            rs=sharp.advance(mesh=mesh_s,boundary=bnd_s,damage=d_s,displacement=u_s,p0=p_s,p1=p1s,direction=(p1s-p_s)/ds,front_id=0,kill_r=max(mesh_s.hbar_tip,0.5e-6))
            if not rs.inserted: raise RuntimeError(f'sharp wake failed event {iev}: {rs.reason}')
            mesh_s,bnd_s,d_s,u_s=rs.mesh,rs.boundary,rs.damage,rs.displacement; p_s=p1s
            path.append({'event':iev,'dx_um':float(ddx*1e6),'ds_um':float(L*1e6),'selected_sign':sign,'x_um':float((p_a[0]-full.A0_M)*1e6),'y_um':float(p_a[1]*1e6),'prior_sign_trial_failures':trial_failures})
            record(iev)
    finally:
        restore_local_scale_gate(saved)
    init={}
    for r in rows:init.setdefault(r['representation'],r['compliance_m2_per_N'])
    for r in rows:r['normalized_compliance']=r['compliance_m2_per_N']/init[r['representation']]
    A=[r for r in rows if r['representation']=='adaptive_czm_failed'];N=[r for r in rows if r['representation']=='adaptive_split_no_network']
    identity=[abs(a['abs_Ftop_N_per_m']-n['abs_Ftop_N_per_m'])/max(a['abs_Ftop_N_per_m'],n['abs_Ftop_N_per_m'],1e-300) for a,n in zip(A,N)]
    with (OUT/'compliance.csv').open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={'schema':'dbtt1000_frozen_topology_micro5_v1','startup':startup,'failure':failure,'path':path,'max_failed_network_identity_rel_diff':max(identity),'rows':rows}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
