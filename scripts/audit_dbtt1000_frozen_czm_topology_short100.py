#!/usr/bin/env python3
"""Fast 100-um gate for the full frozen DBTT topology audit."""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np

import audit_dbtt1000_frozen_czm_topology as full
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.atomic_path_corridor_adaptive_quality_v10051839 import AdaptiveQualityAtomicPathCorridorCZMBackendV10051839
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import install as install_local_scale_gate, restore as restore_local_scale_gate
from arrhenius_fracture.atomic_path_corridor_czm_v10051839 import reset_audit

OUT=Path('dbtt1000_frozen_czm_topology_short100')
TARGETS=(0.0,25.0,50.0,100.0)

def add(rows,name,ext,result,failed=0):
    rows.append({'representation':name,'extension_um':ext*1e6,'failed_interfaces':int(failed),**result})

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ledger=json.loads(full.EVENT_LEDGER.read_text())
    lengths=np.asarray(ledger['event_lengths_m'],float)
    stop=int(np.searchsorted(np.cumsum(lengths)*1e6,100.0)+1)
    lengths=lengths[:stop]
    cfg_a,mesh_a,bnd_a,d_a,u_a=full.initial_state()
    cfg_s,mesh_s,bnd_s,d_s,u_s=full.initial_state()
    startup={'nodes':int(mesh_a.nn),'triangles':int(mesh_a.ne),'hbar_tip_m':float(mesh_a.hbar_tip)}
    if startup['nodes']!=20321 or startup['triangles']!=40502:
        raise RuntimeError(f'production corridor mismatch: {startup}')
    reset_audit(); saved=install_local_scale_gate()
    adaptive=AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(geom=cfg_a.geometry,event_damage=1.0,max_angle_error_deg=35.0)
    sharp=SharpWakeBackend(); rows=[]
    add(rows,'adaptive_czm_failed',0,full.frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a))
    add(rows,'adaptive_split_no_network',0,full.frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a))
    add(rows,'sharp_wake',0,full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s))
    p_a=np.array([full.A0_M,0.0]); p_s=p_a.copy(); ext=0.0; pending=list(TARGETS[1:]); failure=None
    try:
        for iev,L in enumerate(lengths,1):
            prev=ext*1e6; p1=p_a+np.array([float(L),0.0])
            rr=adaptive.advance(mesh=mesh_a,boundary=bnd_a,damage=d_a,displacement=u_a,p0=p_a,p1=p1,direction=np.array([1.0,0.0]),front_id=0)
            if not rr.inserted:
                failure={'event':iev,'extension_before_um':prev,'requested_length_um':float(L*1e6),'reason':rr.reason}; break
            mesh_a,bnd_a,d_a,u_a=rr.mesh,rr.boundary,rr.damage,rr.displacement; p_a=p1; ext=float(p_a[0]-full.A0_M)
            p1s=p_s+np.array([float(L),0.0]); rs=sharp.advance(mesh=mesh_s,boundary=bnd_s,damage=d_s,displacement=u_s,p0=p_s,p1=p1s,direction=np.array([1.0,0.0]),front_id=0,kill_r=max(mesh_s.hbar_tip,0.5e-6))
            mesh_s,bnd_s,d_s,u_s=rs.mesh,rs.boundary,rs.damage,rs.displacement; p_s=p1s
            for cp in [x for x in pending if prev<x<=ext*1e6+1e-9]:
                add(rows,'adaptive_czm_failed',ext,full.frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'adaptive_split_no_network',ext,full.frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'sharp_wake',ext,full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s)); pending.remove(cp)
    finally:
        restore_local_scale_gate(saved)
    initial={}
    for r in rows: initial.setdefault(r['representation'],r['compliance_m2_per_N'])
    for r in rows: r['normalized_compliance']=r['compliance_m2_per_N']/initial[r['representation']]
    A=[r for r in rows if r['representation']=='adaptive_czm_failed']; N=[r for r in rows if r['representation']=='adaptive_split_no_network']
    identity=[abs(a['abs_Ftop_N_per_m']-n['abs_Ftop_N_per_m'])/max(a['abs_Ftop_N_per_m'],n['abs_Ftop_N_per_m'],1e-300) for a,n in zip(A,N)]
    with (OUT/'compliance.csv').open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary={'schema':'dbtt1000_frozen_topology_short100_v1','startup':startup,'events_requested':stop,'failure':failure,'max_failed_network_identity_rel_diff':max(identity),'rows':rows}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
