#!/usr/bin/env python3
import argparse, csv, json, re
from pathlib import Path

import math


def read_json(path):
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def read_csv_first(path):
    if not path.exists():
        return {}
    try:
        with path.open(newline='') as f:
            rows = list(csv.DictReader(f))
        return rows[0] if rows else {}
    except Exception:
        return {}


def read_csv_last(path):
    if not path.exists():
        return {}
    try:
        with path.open(newline='') as f:
            rows = list(csv.DictReader(f))
        return rows[-1] if rows else {}
    except Exception:
        return {}


def to_float(x):
    try:
        if x is None or x == '' or str(x).lower() == 'nan':
            return float('nan')
        return float(x)
    except Exception:
        return float('nan')


def to_int(x):
    try:
        if x is None or x == '' or str(x).lower() == 'nan':
            return 0
        return int(float(x))
    except Exception:
        return 0


def classify(row):
    kc = to_float(row.get('Kc_first_MPa_sqrt_m'))
    n1 = to_int(row.get('n_advances_primary'))
    n2 = to_int(row.get('n_advances_branch'))
    m21 = to_float(row.get('metric2_over_metric1'))
    a1 = to_float(row.get('angle1_deg'))
    a2 = to_float(row.get('angle2_deg'))
    if not math.isfinite(kc) or n1 <= 0:
        return 'no_opening'
    if n2 <= 0:
        if math.isfinite(m21) and m21 >= 0.75:
            return 'latent_branch_hazard_no_daughter_growth'
        return 'single_front'
    if math.isfinite(a1) and math.isfinite(a2):
        if abs(abs(a1) - abs(a2)) < 7.5 and a1 * a2 < 0:
            return 'symmetric_Y'
        if a1 * a2 < 0:
            return 'biased_Y'
        return 'same_side_split'
    return 'branched'


def infer_params(path, args):
    # Use run_args.json if available; otherwise parse directory name.
    out = {}
    for k in ['temperatures','crystal_theta_deg','cleave_gamma_aniso','branch_overdrive_ratio']:
        if k in args:
            out[k] = args.get(k)
    name = path.name
    pats = {
        'T_K': r'T(\d+)',
        'theta_deg': r'th([0-9.]+)',
        'gamma_aniso': r'g([0-9.]+)',
        'branch_overdrive_ratio': r'r([0-9.]+)',
    }
    for k,p in pats.items():
        m = re.search(p, name)
        if m:
            out[k] = m.group(1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('--out', default=None)
    ns = ap.parse_args()
    root = Path(ns.root)
    rows = []
    for case in sorted([p for p in root.rglob('*') if p.is_dir()]):
        summ = read_json(case/'summary.json')
        if not summ:
            continue
        args = read_json(case/'run_args.json')
        params = infer_params(case, args)
        # Find first temperature-specific diagnostics.
        diag_files = sorted(case.glob('branch_diagnostics_*K.csv'))
        step_files = sorted(case.glob('steps_*K.csv'))
        diag = read_csv_last(diag_files[0]) if diag_files else {}
        steps = read_csv_last(step_files[0]) if step_files else {}
        out = {
            'case': str(case.relative_to(root)),
            'T_K': params.get('T_K', ''),
            'theta_deg': params.get('theta_deg', params.get('crystal_theta_deg','')),
            'gamma_aniso': params.get('gamma_aniso', params.get('cleave_gamma_aniso','')),
            'branch_overdrive_ratio': params.get('branch_overdrive_ratio',''),
            'Kc_first_MPa_sqrt_m': summ.get('Kc_first_MPa_sqrt_m',''),
            'branched': summ.get('branched',''),
            'n_advances_primary': summ.get('n_advances_primary',''),
            'n_advances_branch': summ.get('n_advances_branch',''),
            'N_em_init': summ.get('N_em_init',''),
            'sigma_back_init_GPa': summ.get('sigma_back_init_GPa',''),
            'metric2_over_metric1': diag.get('metric2_over_metric1',''),
            'angle1_deg': diag.get('angle1_deg',''),
            'angle2_deg': diag.get('angle2_deg',''),
            'share_w1': diag.get('share_w1',''),
            'share_w2': diag.get('share_w2',''),
            'step_final': steps.get('step',''),
            'KJ_final': steps.get('KJ_MPa_sqrt_m',''),
        }
        out['morphology'] = classify(out)
        rows.append(out)
    fieldnames = ['case','T_K','theta_deg','gamma_aniso','branch_overdrive_ratio','Kc_first_MPa_sqrt_m','branched','n_advances_primary','n_advances_branch','N_em_init','sigma_back_init_GPa','metric2_over_metric1','angle1_deg','angle2_deg','share_w1','share_w2','step_final','KJ_final','morphology']
    out_path = Path(ns.out) if ns.out else root/'reduced_branch_summary.csv'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'wrote {out_path} ({len(rows)} cases)')

if __name__ == '__main__':
    main()
