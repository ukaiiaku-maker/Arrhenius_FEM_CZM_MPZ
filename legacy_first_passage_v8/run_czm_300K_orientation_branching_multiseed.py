#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable


def split_floats(s: str) -> list[float]:
    return [float(x) for x in s.replace(',', ' ').split() if x.strip()]


def split_ints(s: str) -> list[int]:
    return [int(float(x)) for x in s.replace(',', ' ').split() if x.strip()]


def setopt(cmd: list[str], opt: str, val: str | int | float) -> list[str]:
    cmd = list(cmd)
    if opt in cmd:
        cmd[cmd.index(opt) + 1] = str(val)
    else:
        cmd.extend([opt, str(val)])
    return cmd


def ensure_flag(cmd: list[str], flag: str) -> list[str]:
    cmd = list(cmd)
    if flag not in cmd:
        cmd.append(flag)
    return cmd


def set_seed(cmd: list[str], seed: int) -> list[str]:
    cmd = list(cmd)
    for opt in ["--solver-seed", "--seed", "--random-seed"]:
        if opt in cmd:
            cmd[cmd.index(opt) + 1] = str(seed)
            return cmd
    cmd.extend(["--solver-seed", str(seed)])
    return cmd


def command_from_script(path: Path) -> list[str] | None:
    if not path.exists():
        return None
    txt = path.read_text(errors="ignore")
    for raw in txt.splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if 'run_seeded_sharp_front.py' in line:
            line = re.sub(r"\s*>\s*[^\s]+\s+2>&1\s*$", "", line)
            line = line.split('>', 1)[0].strip()
            try:
                return shlex.split(line)
            except Exception:
                return None
    return None


def class_params(klass: str) -> dict[str, str]:
    # Effective emission values are used here, matching the values passed to run_seeded_sharp_front.py.
    tbl: dict[str, dict[str, str]] = {
        'ceramic': {
            'emit-G00-eV': '1.3922', 'emit-gT-eV-per-K': '0.0091197', 'emit-sigc0-GPa': '2.9802', 'emit-sT-GPa-per-K': '0.0020310', 'emit-exp-a': '0.2383', 'emit-exp-n': '1.1169', 'emit-floor-frac': '0.04100',
            'cleave-G00-eV': '1.7884', 'cleave-gT-eV-per-K': '-0.0002532', 'cleave-sigc0-GPa': '4.7808', 'cleave-sT-GPa-per-K': '0.0006683', 'cleave-exp-a': '1.3908', 'cleave-exp-n': '2.0000', 'cleave-floor-frac': '0.00928', 'cleave-S-hs-kB': '-8.7810', 'cleave-shield-chi': '0.4999', 'n-sat': '655.8',
        },
        'weakT': {
            'emit-G00-eV': '1.1173', 'emit-gT-eV-per-K': '0.0063971', 'emit-sigc0-GPa': '0.9506', 'emit-sT-GPa-per-K': '0.0009817', 'emit-exp-a': '0.5055', 'emit-exp-n': '0.8432', 'emit-floor-frac': '0.03090',
            'cleave-G00-eV': '2.5269', 'cleave-gT-eV-per-K': '0.0016068', 'cleave-sigc0-GPa': '4.0609', 'cleave-sT-GPa-per-K': '0.0026762', 'cleave-exp-a': '0.8360', 'cleave-exp-n': '1.4882', 'cleave-floor-frac': '0.00478', 'cleave-S-hs-kB': '0.1098', 'cleave-shield-chi': '0.1113', 'n-sat': 'inf',
        },
        'peak': {
            'emit-G00-eV': '2.6094', 'emit-gT-eV-per-K': '-0.0009574', 'emit-sigc0-GPa': '5.5176', 'emit-sT-GPa-per-K': '0.0009500', 'emit-exp-a': '0.5254', 'emit-exp-n': '0.8708', 'emit-floor-frac': '0.08442',
            'cleave-G00-eV': '3.8721', 'cleave-gT-eV-per-K': '-0.0026975', 'cleave-sigc0-GPa': '4.2182', 'cleave-sT-GPa-per-K': '-0.0005593', 'cleave-exp-a': '0.8629', 'cleave-exp-n': '1.0271', 'cleave-floor-frac': '0.00528', 'cleave-S-hs-kB': '-0.6268', 'cleave-shield-chi': '0.4163', 'n-sat': '2730.7',
        },
        'DBTT': {
            'emit-G00-eV': '1.695399975', 'emit-gT-eV-per-K': '0.0008853000000000001', 'emit-sigc0-GPa': '3.7521', 'emit-sT-GPa-per-K': '0.0022809', 'emit-exp-a': '0.0944', 'emit-exp-n': '0.835', 'emit-floor-frac': '0.03524',
            'cleave-G00-eV': '2.807', 'cleave-gT-eV-per-K': '0.0045886', 'cleave-sigc0-GPa': '4.1195', 'cleave-sT-GPa-per-K': '0.0000977', 'cleave-exp-a': '0.5925', 'cleave-exp-n': '1.2024', 'cleave-floor-frac': '0.00256', 'cleave-S-hs-kB': '6.9887', 'cleave-shield-chi': '0.79', 'n-sat': '1488.8',
        },
    }
    if klass not in tbl:
        raise SystemExit(f"Unknown class {klass!r}; choose one of {sorted(tbl)}")
    return tbl[klass]


def hardcoded_base(pybin: str, project: Path, klass: str) -> list[str]:
    p = class_params(klass)
    cmd = [
        pybin, '-u', str(project / 'run_seeded_sharp_front.py'),
        '--solver-seed', '1201',
        '--mode', '2d',
        '--nx', '12', '--ny', '24', '--tip-h-fine', '5e-6', '--tip-ratio', '1.3',
        '--dU', '2e-07', '--dt', '8.4', '--steps', '250000', '--n-stagger', '2', '--print-every', '500',
        '--target-crack-extension-um', '1000',
        '--crystal-aniso', '--crystal-compete', '--crystal-material', 'branchy', '--cleave-gamma-aniso', '2',
        '--multihit-m', '3', '--multihit-tau', '1e-6', '--emb-sat-frac', '1',
        '--adaptive-events', '--adaptive-event-target', '0.35', '--adaptive-min-frac', '1e-8', '--adaptive-grow', '4.0',
        '--max-fronts', '1', '--da-phys', '5e-6',
        '--j-decomposition', 'cluster', '--rJ-cluster', '20e-6', '--rJ-outer', '25e-6',
        '--temperatures', '300', '--crystal-theta-deg', '45',
        '--crack-backend', 'adaptive_czm', '--czm-max-angle-error-deg', '60',
        '--emit-barrier-kind', 'exp_floor', '--emit-Tref-K', '300',
        '--cleave-barrier-kind', 'exp_floor', '--cleave-exp-T-mode', 'linear', '--cleave-sigma-S-GPa', '6', '--cleave-S-hs-power', '2', '--cleave-S-hs-Tref-K', '300', '--cleave-Tref-K', '300',
    ]
    for key in ['emit-G00-eV', 'emit-gT-eV-per-K', 'emit-sigc0-GPa', 'emit-sT-GPa-per-K', 'emit-exp-a', 'emit-exp-n', 'emit-floor-frac']:
        cmd += [f'--{key}', p[key]]
    for key in ['cleave-G00-eV', 'cleave-gT-eV-per-K', 'cleave-sigc0-GPa', 'cleave-sT-GPa-per-K', 'cleave-exp-a', 'cleave-exp-n', 'cleave-floor-frac', 'cleave-S-hs-kB', 'cleave-shield-chi', 'n-sat']:
        cmd += [f'--{key}', p[key]]
    return cmd


def find_template(project: Path, pybin: str, klass: str, allow_hardcoded: bool = True) -> list[str]:
    candidates = [
        project / f"runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/{klass}/replicate_01_seed1101/T500_th45/command.txt",
        project / f"runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x/{klass}/T300_th45/command.txt",
        project / f"runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x/{klass}/T500_th45/command.txt",
    ]
    for c in candidates:
        if c.exists():
            try:
                return shlex.split(c.read_text())
            except Exception:
                pass
    script_dirs = [
        project / f"runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/{klass}/replicate_01_seed1101/T500_th45",
        project / f"runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x/{klass}/T300_th45",
    ]
    for d in script_dirs:
        if d.exists():
            for s in sorted(d.glob('run*.sh')):
                cmd = command_from_script(s)
                if cmd:
                    return cmd
    if allow_hardcoded:
        print(f"NOTE: no template found for {klass}; using built-in class parameters.", flush=True)
        return hardcoded_base(pybin, project, klass)
    raise SystemExit(f"Could not find a template command for {klass}")


def normalize_base(cmd: list[str], pybin: str, project: Path) -> list[str]:
    script_idx = None
    for i, item in enumerate(cmd):
        if 'run_seeded_sharp_front.py' in item:
            script_idx = i
            break
    if script_idx is None:
        raise SystemExit('Template command does not contain run_seeded_sharp_front.py')
    tail = list(cmd[script_idx:])
    if not Path(tail[0]).is_absolute():
        tail[0] = str(project / tail[0])
    return [pybin, '-u'] + tail


def apply_common_options(cmd: list[str], args, theta: float, seed: int, branch: bool, case: Path, retry: bool = False) -> list[str]:
    cmd = set_seed(cmd, seed)
    cmd = setopt(cmd, '--temperatures', args.temperature)
    cmd = setopt(cmd, '--crystal-theta-deg', theta)
    cmd = setopt(cmd, '--target-crack-extension-um', args.target_ext_um)
    cmd = setopt(cmd, '--steps', args.retry_steps if retry else args.steps)
    cmd = setopt(cmd, '--print-every', args.print_every)
    cmd = setopt(cmd, '--max-fronts', args.max_fronts_branch if branch else args.max_fronts_base)
    cmd = setopt(cmd, '--out', str(case))
    cmd = setopt(cmd, '--save-snapshots', args.branch_save_snapshots if branch else args.save_snapshots)
    cmd = setopt(cmd, '--nx', args.retry_nx if retry else args.nx)
    cmd = setopt(cmd, '--ny', args.retry_ny if retry else args.ny)
    cmd = setopt(cmd, '--tip-h-fine', args.retry_tip_h_fine if retry else args.tip_h_fine)
    cmd = setopt(cmd, '--tip-ratio', args.retry_tip_ratio if retry else args.tip_ratio)
    cmd = setopt(cmd, '--da-phys', args.retry_da_phys if retry else args.da_phys)
    cmd = setopt(cmd, '--adaptive-event-target', args.retry_event_target if retry else args.event_target)
    cmd = setopt(cmd, '--adaptive-min-frac', args.adaptive_min_frac)
    cmd = setopt(cmd, '--adaptive-grow', args.adaptive_grow)
    cmd = setopt(cmd, '--czm-max-angle-error-deg', args.branch_angle_error_deg if branch else args.angle_error_deg)
    cmd = ensure_flag(cmd, '--no-plots')
    return cmd


def read_text_tail(path: Path, max_bytes: int = 200000) -> str:
    if not path.exists():
        return ''
    try:
        size = path.stat().st_size
        with path.open('rb') as f:
            if size > max_bytes:
                f.seek(size - max_bytes)
            return f.read().decode('utf-8', errors='ignore')
    except Exception:
        return ''


def completion_status(case: Path, target_ext_um: float, initial_a_mm: float = 0.5) -> dict:
    status = {'complete': False, 'extension_um': math.nan, 'target_hit': False, 'n_events': 0, 'geom_veto_count': 0, 'zero_k_tail': False, 'reason': ''}
    log = case / 'run.log'
    txt = read_text_tail(log, max_bytes=2_000_000)
    if txt:
        status['target_hit'] = 'reached target crack extension' in txt
        m = re.findall(r"reached target crack extension\s+([0-9.]+)\s+um", txt)
        if m:
            status['extension_um'] = max(status['extension_um'], float(m[-1])) if math.isfinite(status['extension_um']) else float(m[-1])
        status['geom_veto_count'] = txt.count('GEOMETRY VETO')
        lines = [ln for ln in txt.splitlines() if '[T=' in ln and 'KJ=' in ln]
        if len(lines) >= 8:
            tail = lines[-8:]
            status['zero_k_tail'] = all(re.search(r"KJ=\s*0\.000", ln) for ln in tail)
    for csv in [case / 'R_curve_event_sampled.csv'] + sorted(case.glob('steps_*K.csv')):
        if csv.exists():
            try:
                import pandas as pd
                df = pd.read_csv(csv)
                status['n_events'] = max(status['n_events'], len(df))
                for col, scale in [('crack_extension_um', 1.0), ('extension_um', 1.0), ('crack_extension_m', 1e6)]:
                    if col in df.columns and len(df):
                        ext = float(pd.to_numeric(df[col], errors='coerce').max() * scale)
                        status['extension_um'] = max(status['extension_um'], ext) if math.isfinite(status['extension_um']) else ext
            except Exception:
                pass
    sj = case / 'summary.json'
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            if isinstance(s, list) and s:
                s = s[0]
            if s.get('a_final_mm') is not None:
                val = (float(s['a_final_mm']) - initial_a_mm) * 1000.0
                status['extension_um'] = max(status['extension_um'], val) if math.isfinite(status['extension_um']) else val
            for key in ['final_crack_extension_um', 'crack_extension_um', 'extension_um']:
                if s.get(key) is not None:
                    val = float(s[key])
                    status['extension_um'] = max(status['extension_um'], val) if math.isfinite(status['extension_um']) else val
        except Exception:
            pass
    status['complete'] = bool(math.isfinite(status['extension_um']) and status['extension_um'] >= 0.98 * target_ext_um)
    if status['complete']:
        status['reason'] = 'complete'
    elif status['geom_veto_count']:
        status['reason'] = 'geometry_veto_or_incomplete'
    elif status['zero_k_tail']:
        status['reason'] = 'zero_K_tail_or_incomplete'
    else:
        status['reason'] = 'incomplete_or_missing'
    return status


def monitor_should_kill(case: Path, args, last_ext_seen: float) -> tuple[bool, str, float]:
    txt = read_text_tail(case / 'run.log', max_bytes=500000)
    if not txt:
        return False, '', last_ext_seen
    # update extension from the latest progress lines
    ext = last_ext_seen
    for m in re.findall(r"a=([0-9.]+)mm", txt):
        try:
            # Crack extension from initial 0.5 mm; only a rough live monitor.
            ext = max(ext, (float(m) - 0.5) * 1000.0)
        except Exception:
            pass
    tail_lines = txt.splitlines()[-args.monitor_tail_lines:]
    tail = '\n'.join(tail_lines)
    # Detect repeated topology traps, including both old two-element signatures
    # ([1916, 1964]) and local_hrefine single-element signatures ([759]).
    veto_lines = [ln for ln in tail_lines if 'GEOMETRY VETO' in ln]
    if len(veto_lines) >= args.veto_kill_count:
        sigs = []
        for ln in veto_lines:
            br = re.search(r"\[([^\]]+)\]", ln)
            kind = 'GEOMETRY_VETO'
            km = re.search(r"GEOMETRY VETO front \d+:\s*([^:]+):", ln)
            if km:
                kind = km.group(1).strip()
            if br:
                sigs.append((kind, br.group(1).replace(' ', '')))
            else:
                sigs.append((kind, 'no_element_id'))
        most = max(set(sigs), key=sigs.count)
        if sigs.count(most) >= args.veto_kill_count:
            return True, f"repeated GEOMETRY VETO {most} count={sigs.count(most)}", ext
    prog = [ln for ln in tail_lines if '[T=' in ln and 'KJ=' in ln]
    if len(prog) >= args.zero_k_kill_count:
        z = prog[-args.zero_k_kill_count:]
        if all(re.search(r"KJ=\s*0\.000", ln) for ln in z):
            return True, f"zero-K tail for last {args.zero_k_kill_count} progress lines", ext
    return False, '', ext


def write_run_script(case: Path, project: Path, cmd: list[str]) -> Path:
    case.mkdir(parents=True, exist_ok=True)
    (case / 'command.txt').write_text(shlex.join(cmd) + '\n')
    run = case / 'run_case.sh'
    run.write_text(
        '#!/usr/bin/env bash\n'
        'set -euo pipefail\n'
        f'cd {shlex.quote(str(project))}\n'
        + shlex.join(cmd) + f" > {shlex.quote(str(case / 'run.log'))} 2>&1\n"
    )
    run.chmod(0o755)
    return run


def archive_case(case: Path, suffix: str) -> None:
    if not case.exists():
        return
    stamp = time.strftime('%Y%m%d_%H%M%S')
    dst = case.parent / f"{case.name}_{suffix}_{stamp}"
    case.rename(dst)
    print(f"ARCHIVED {case} -> {dst}", flush=True)


def run_one(case: Path, run_script: Path, args) -> tuple[int, str]:
    print(f"RUN {case}", flush=True)
    proc = subprocess.Popen(['bash', str(run_script)], cwd=str(args.project))
    start = time.time()
    last_ext = 0.0
    kill_reason = ''
    while True:
        rc = proc.poll()
        if rc is not None:
            return rc, kill_reason
        elapsed = time.time() - start
        if args.timeout_hours > 0 and elapsed > args.timeout_hours * 3600:
            kill_reason = f"timeout>{args.timeout_hours}h"
        else:
            should, why, last_ext = monitor_should_kill(case, args, last_ext)
            if should:
                kill_reason = why
        if kill_reason:
            print(f"KILL {case}: {kill_reason}", flush=True)
            try:
                proc.terminate()
                time.sleep(5)
                if proc.poll() is None:
                    proc.kill()
            except Exception:
                pass
            return -9, kill_reason
        time.sleep(args.monitor_interval)


def make_jobs(args) -> list[dict]:
    jobs = []
    for theta in args.thetas:
        for seed in args.seeds:
            rel = Path(args.klass) / f"theta_{theta:05.1f}_nobranch" / f"seed{seed}" / f"T{args.temperature}_th{theta:g}"
            jobs.append({'theta': theta, 'seed': seed, 'branch': False, 'case': args.outroot / rel})
    if args.run_branch:
        for seed in args.seeds:
            theta = args.branch_theta
            rel = Path(args.klass) / f"theta_{theta:05.1f}_weakbranch" / f"seed{seed}" / f"T{args.temperature}_th{theta:g}"
            jobs.append({'theta': theta, 'seed': seed, 'branch': True, 'case': args.outroot / rel})
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser(description='300 K multi-seed orientation + weak-branching FEM/CZM R-curve sweep.')
    ap.add_argument('--project', type=Path, default=Path.cwd())
    ap.add_argument('--python-bin', default=os.environ.get('PYTHON_BIN', 'python'))
    ap.add_argument('--class', dest='klass', default=os.environ.get('CLASS', 'ceramic'))
    ap.add_argument('--temperature', type=int, default=int(os.environ.get('T_K', '300')))
    ap.add_argument('--thetas', type=split_floats, default=split_floats(os.environ.get('THETAS', '0 15 30 45')))
    ap.add_argument('--branch-theta', type=float, default=float(os.environ.get('BRANCH_THETA', '30')))
    ap.add_argument('--seeds', type=split_ints, default=split_ints(os.environ.get('SEEDS', '1201 1202 1203 1204 1205')))
    ap.add_argument('--run-branch', action=argparse.BooleanOptionalAction, default=os.environ.get('RUN_BRANCH', '1') not in {'0','false','False','no','NO'})
    ap.add_argument('--outroot', type=Path, default=Path(os.environ.get('OUTROOT', 'runs/czm_Rcurve_300K_orientation_branching_multiseed_v2')))
    ap.add_argument('--run-now', action=argparse.BooleanOptionalAction, default=os.environ.get('RUN_NOW', '1') not in {'0','false','False','no','NO'})
    ap.add_argument('--allow-hardcoded-template', action=argparse.BooleanOptionalAction, default=True)

    ap.add_argument('--target-ext-um', type=float, default=float(os.environ.get('TARGET_EXT_UM', '1000')))
    ap.add_argument('--steps', type=int, default=int(os.environ.get('STEPS', '80000')))
    ap.add_argument('--retry-steps', type=int, default=int(os.environ.get('RETRY_STEPS', '120000')))
    ap.add_argument('--print-every', type=int, default=int(os.environ.get('PRINT_EVERY', '500')))
    ap.add_argument('--max-fronts-base', type=int, default=int(os.environ.get('MAX_FRONTS_BASE', '1')))
    ap.add_argument('--max-fronts-branch', type=int, default=int(os.environ.get('MAX_FRONTS_BRANCH', '3')))
    ap.add_argument('--save-snapshots', default=os.environ.get('SAVE_SNAPSHOTS', '0'))
    ap.add_argument('--branch-save-snapshots', default=os.environ.get('BRANCH_SAVE_SNAPSHOTS', '10'))

    # Conservative mesh/event defaults to reduce topology traps.
    ap.add_argument('--nx', default=os.environ.get('NX', '12'))
    ap.add_argument('--ny', default=os.environ.get('NY', '24'))
    ap.add_argument('--tip-h-fine', default=os.environ.get('TIP_H_FINE', '5e-6'))
    ap.add_argument('--tip-ratio', default=os.environ.get('TIP_RATIO', '1.3'))
    ap.add_argument('--da-phys', default=os.environ.get('DA_PHYS', '5e-6'))
    ap.add_argument('--event-target', default=os.environ.get('EVENT_TARGET', '0.35'))
    ap.add_argument('--adaptive-min-frac', default=os.environ.get('ADAPTIVE_MIN_FRAC', '1e-8'))
    ap.add_argument('--adaptive-grow', default=os.environ.get('ADAPTIVE_GROW', '4.0'))
    ap.add_argument('--angle-error-deg', default=os.environ.get('CZM_MAX_ANGLE_ERROR_DEG', '60'))
    ap.add_argument('--branch-angle-error-deg', default=os.environ.get('BRANCH_CZM_MAX_ANGLE_ERROR_DEG', '60'))

    # Retry is finer/slower but used only after a detected bad case.
    ap.add_argument('--retry-on-incomplete', action=argparse.BooleanOptionalAction, default=os.environ.get('RETRY_ON_INCOMPLETE', '0') not in {'0','false','False','no','NO'})
    ap.add_argument('--retry-nx', default=os.environ.get('RETRY_NX', '12'))
    ap.add_argument('--retry-ny', default=os.environ.get('RETRY_NY', '24'))
    ap.add_argument('--retry-tip-h-fine', default=os.environ.get('RETRY_TIP_H_FINE', '5e-6'))
    ap.add_argument('--retry-tip-ratio', default=os.environ.get('RETRY_TIP_RATIO', '1.3'))
    ap.add_argument('--retry-da-phys', default=os.environ.get('RETRY_DA_PHYS', '5e-6'))
    ap.add_argument('--retry-event-target', default=os.environ.get('RETRY_EVENT_TARGET', '0.35'))

    # Live monitor / kill guards.
    ap.add_argument('--monitor-interval', type=float, default=float(os.environ.get('MONITOR_INTERVAL', '30')))
    ap.add_argument('--monitor-tail-lines', type=int, default=int(os.environ.get('MONITOR_TAIL_LINES', '250')))
    ap.add_argument('--veto-kill-count', type=int, default=int(os.environ.get('VETO_KILL_COUNT', '8')))
    ap.add_argument('--zero-k-kill-count', type=int, default=int(os.environ.get('ZERO_K_KILL_COUNT', '8')))
    ap.add_argument('--timeout-hours', type=float, default=float(os.environ.get('TIMEOUT_HOURS_PER_CASE', '0')))

    args = ap.parse_args()
    args.project = args.project.resolve()
    if not args.outroot.is_absolute():
        args.outroot = (args.project / args.outroot).resolve()
    args.outroot.mkdir(parents=True, exist_ok=True)

    print(f"python: {args.python_bin}")
    print(f"class:  {args.klass}")
    print(f"T:      {args.temperature} K")
    print(f"thetas: {args.thetas}")
    print(f"seeds:  {args.seeds}")
    print(f"out:    {args.outroot}")
    print(f"mesh:   nx={args.nx} ny={args.ny} tip_h={args.tip_h_fine} da={args.da_phys} event_target={args.event_target}")

    base_template = find_template(args.project, args.python_bin, args.klass, args.allow_hardcoded_template)
    base = normalize_base(base_template, args.python_bin, args.project)
    jobs = make_jobs(args)

    manifest = []
    for job in jobs:
        case = job['case']
        stat = completion_status(case, args.target_ext_um)
        cmd = apply_common_options(base, args, job['theta'], job['seed'], job['branch'], case, retry=False)
        run = write_run_script(case, args.project, cmd)
        row = {**job, 'case': str(case.relative_to(args.project)), 'run_script': str(run.relative_to(args.project)), **stat}
        manifest.append(row)
    (args.outroot / 'sweep_manifest_initial.json').write_text(json.dumps(manifest, indent=2, default=str))

    if not args.run_now:
        print(f"RUN_NOW=0: generated {len(jobs)} cases and scripts only.")
        return 0

    final = []
    for job in jobs:
        case = job['case']
        stat = completion_status(case, args.target_ext_um)
        if stat['complete']:
            print(f"SKIP complete: {case.relative_to(args.project)} ext={stat['extension_um']:.1f} um", flush=True)
            final.append({**job, 'case': str(case.relative_to(args.project)), 'attempt': 0, **stat})
            continue

        cmd = apply_common_options(base, args, job['theta'], job['seed'], job['branch'], case, retry=False)
        run = write_run_script(case, args.project, cmd)
        rc, kill_reason = run_one(case, run, args)
        stat = completion_status(case, args.target_ext_um)
        stat['return_code'] = rc
        stat['kill_reason'] = kill_reason

        if (not stat['complete']) and args.retry_on_incomplete:
            print(f"INCOMPLETE after first attempt: {case.relative_to(args.project)} reason={stat['reason']} ext={stat['extension_um']}", flush=True)
            archive_case(case, 'attempt0_incomplete')
            cmd = apply_common_options(base, args, job['theta'], job['seed'], job['branch'], case, retry=True)
            run = write_run_script(case, args.project, cmd)
            rc2, kill_reason2 = run_one(case, run, args)
            stat2 = completion_status(case, args.target_ext_um)
            stat2['return_code'] = rc2
            stat2['kill_reason'] = kill_reason2
            final.append({**job, 'case': str(case.relative_to(args.project)), 'attempt': 1, **stat2})
            print(f"FINAL retry status: complete={stat2['complete']} ext={stat2['extension_um']} reason={stat2['reason']}", flush=True)
        else:
            final.append({**job, 'case': str(case.relative_to(args.project)), 'attempt': 0, **stat})
            print(f"FINAL status: complete={stat['complete']} ext={stat['extension_um']} reason={stat['reason']}", flush=True)

        (args.outroot / 'sweep_manifest_live.json').write_text(json.dumps(final, indent=2, default=str))

    (args.outroot / 'sweep_manifest_final.json').write_text(json.dumps(final, indent=2, default=str))
    n_complete = sum(1 for r in final if r.get('complete'))
    print(f"DONE campaign: complete {n_complete}/{len(final)}")
    return 0 if n_complete == len(final) else 2


if __name__ == '__main__':
    raise SystemExit(main())
