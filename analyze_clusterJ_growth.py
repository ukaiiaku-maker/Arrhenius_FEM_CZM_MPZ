#!/usr/bin/env python3
import csv, json, os, sys, glob
from pathlib import Path
import numpy as np

def read_summary(case):
    rows = list(Path(case).glob('summary_*.csv'))
    out = {'case': Path(case).name}
    if rows:
        with open(rows[0], newline='') as f:
            r = next(csv.DictReader(f))
            out.update(r)
    fs = list(Path(case).glob('fronts_*.csv'))
    if fs:
        data = np.genfromtxt(fs[0], delimiter=',', names=True)
        if data.size:
            if data.ndim == 0: data = np.array([data])
            out['n_front_ids'] = int(len(set(data['front_id'].astype(int))))
            out['max_x_mm_front_rows'] = float(np.nanmax(data['x_m'])*1e3)
            out['max_branch_len_um'] = float(np.nanmax(data['branch_len_m'])*1e6) if 'branch_len_m' in data.dtype.names else np.nan
            if 'J_source_code' in data.dtype.names:
                out['n_cluster_rows'] = int(np.sum(data['J_source_code']==0))
                out['n_local_rows'] = int(np.sum(data['J_source_code']==1))
                out['n_unresolved_rows'] = int(np.sum(data['J_source_code']==2))
            if 'J_active_elems' in data.dtype.names:
                out['min_active_J_elems'] = int(np.nanmin(data['J_active_elems']))
                out['median_active_J_elems'] = float(np.nanmedian(data['J_active_elems']))
            if 'resolved' in data.dtype.names:
                out['n_resolved_rows'] = int(np.nansum(data['resolved'] == 1))
                out['n_unresolved_rows_by_flag'] = int(np.nansum(data['resolved'] == 0))
            if 'cluster_hold_code' in data.dtype.names:
                # 1 short_branch, 2 parent_overlap, 3 overlap, 4 few_J_elements, 5 probe_failed, 7 resolved
                out['hold_short_rows'] = int(np.nansum(data['cluster_hold_code'] == 1))
                out['hold_parent_overlap_rows'] = int(np.nansum(data['cluster_hold_code'] == 2))
                out['hold_overlap_rows'] = int(np.nansum(data['cluster_hold_code'] == 3))
                out['hold_few_elem_rows'] = int(np.nansum(data['cluster_hold_code'] == 4))
                out['hold_resolved_rows'] = int(np.nansum(data['cluster_hold_code'] == 7))
            if 'stagnant_retired' in data.dtype.names:
                out['n_stagnant_retired_fronts'] = int(len(set(data['front_id'][data['stagnant_retired'] == 1].astype(int))))
            if 'retire_step' in data.dtype.names:
                vals = data['retire_step'][data['retire_step'] >= 0]
                out['first_retire_step'] = int(np.nanmin(vals)) if vals.size else -1
            if 'branch_B_before_spawn' in data.dtype.names:
                out['max_branch_B_before_spawn'] = float(np.nanmax(data['branch_B_before_spawn']))
            if 'branch_lambda_secondary_at_spawn' in data.dtype.names:
                out['max_branch_lambda_secondary_at_spawn'] = float(np.nanmax(data['branch_lambda_secondary_at_spawn']))
            if 'n_geom_adv' in data.dtype.names:
                out['max_geom_adv'] = float(np.nanmax(data['n_geom_adv']))
                # final row per front, summed approximate geometric advances
                ids = sorted(set(data['front_id'].astype(int)))
                total = 0.0; branch_total = 0.0
                for fid in ids:
                    sub = data[data['front_id'].astype(int)==fid]
                    val = float(sub['n_geom_adv'][-1])
                    total += val
                    if fid != 0:
                        branch_total += val
                out['geom_adv_total_final'] = total
                out['geom_adv_branch_final'] = branch_total
    return out

cases = [p for p in glob.glob(os.path.join(sys.argv[1], '*')) if os.path.isdir(p)]
rows = [read_summary(c) for c in sorted(cases)]
fields = sorted({k for r in rows for k in r})
w = csv.DictWriter(sys.stdout, fieldnames=fields)
w.writeheader()
w.writerows(rows)
