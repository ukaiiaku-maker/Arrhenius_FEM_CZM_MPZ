#!/usr/bin/env python3
import argparse, json, os, glob, re
import numpy as np
import pandas as pd

def path_stats(path):
    df = pd.read_csv(path)
    xcol = 'x' if 'x' in df.columns else 'x_m'
    ycol = 'y' if 'y' in df.columns else 'y_m'
    xy = df[[xcol, ycol]].to_numpy(float)
    length = float(np.linalg.norm(np.diff(xy, axis=0), axis=1).sum()) if len(xy) > 1 else 0.0
    return {
        'npts': len(df),
        'length_m': length,
        'dx_m': float(xy[-1,0] - xy[0,0]) if len(xy) else 0.0,
        'dy_m': float(xy[-1,1] - xy[0,1]) if len(xy) else 0.0,
        'xmax_m': float(np.max(xy[:,0])) if len(xy) else np.nan,
        'xmin_m': float(np.min(xy[:,0])) if len(xy) else np.nan,
        'ymax_m': float(np.max(xy[:,1])) if len(xy) else np.nan,
        'ymin_m': float(np.min(xy[:,1])) if len(xy) else np.nan,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    rows = []
    for cdir in sorted(glob.glob(os.path.join(args.root, '**'), recursive=True)):
        if not os.path.isdir(cdir):
            continue
        summ = os.path.join(cdir, 'summary.json')
        ffiles = glob.glob(os.path.join(cdir, 'fronts_*K.csv'))
        if not os.path.exists(summ) and not ffiles:
            continue
        case = os.path.relpath(cdir, args.root)
        row = {'case': case}
        if os.path.exists(summ):
            try:
                s = json.load(open(summ))
                if isinstance(s, list) and s:
                    row.update(s[0])
                elif isinstance(s, dict):
                    row.update(s)
            except Exception:
                pass
        path_rows = []
        for p in glob.glob(os.path.join(cdir, 'crack_path_front*_*.csv')):
            m = re.search(r'front(\d+)_', os.path.basename(p))
            fid = int(m.group(1)) if m else -1
            st = path_stats(p)
            st['front_id'] = fid
            path_rows.append(st)
        if path_rows:
            row['n_path_fronts'] = len(path_rows)
            row['max_path_length_mm'] = max(r['length_m'] for r in path_rows) * 1e3
            row['max_x_reach_mm'] = max(r['xmax_m'] for r in path_rows) * 1e3
            row['sum_path_length_mm'] = sum(r['length_m'] for r in path_rows) * 1e3
            row['fronts_gt_20um'] = sum(r['length_m'] > 20e-6 for r in path_rows)
            row['fronts_gt_50um'] = sum(r['length_m'] > 50e-6 for r in path_rows)
            row['fronts_gt_100um'] = sum(r['length_m'] > 100e-6 for r in path_rows)
        if ffiles:
            f = pd.read_csv(ffiles[0])
            if 'resolved' in f.columns:
                row['resolved_fronts_seen'] = int(f.groupby('front_id')['resolved'].max().sum())
            if 'branch_len_m' in f.columns:
                row['max_branch_len_mm'] = float(f['branch_len_m'].max()) * 1e3
        rows.append(row)
    df = pd.DataFrame(rows)
    if args.out:
        df.to_csv(args.out, index=False)
    if not df.empty:
        cols = [c for c in ['case','T','Kc_first_MPa_sqrt_m','n_advances','n_fronts','max_reach','a_final_mm','max_path_length_mm','max_x_reach_mm','fronts_gt_50um','fronts_gt_100um'] if c in df.columns]
        print(df[cols].to_string(index=False))
    else:
        print('No case outputs found.')
if __name__ == '__main__':
    main()
