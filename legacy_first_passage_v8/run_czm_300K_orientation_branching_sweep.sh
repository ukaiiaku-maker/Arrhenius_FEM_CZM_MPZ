#!/usr/bin/env bash
set -euo pipefail

# FEM/CZM R-curve orientation + weak-branching sweep at one temperature and one barrier class.
# Run from the Arrhenius_FEM_CZM project root.

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-/opt/homebrew/Caskroom/miniconda/base/envs/${CONDA_ENV}/bin/python}
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN=$(command -v python)
fi

CLASS=${CLASS:-ceramic}                 # recommended first pass: ceramic; alternative: weakT
T_K=${T_K:-300}
THETAS=${THETAS:-"0 15 30 45"}
BRANCH_THETA=${BRANCH_THETA:-30}
SEED=${SEED:-1201}
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
STEPS=${STEPS:-80000}
PRINT_EVERY=${PRINT_EVERY:-500}

# Keep a version-specific default folder and do not overwrite older campaigns.
OUTROOT=${OUTROOT:-runs/czm_Rcurve_300K_orientation_branching_${CLASS}_v1}

MAX_FRONTS_BASE=${MAX_FRONTS_BASE:-1}
MAX_FRONTS_BRANCH=${MAX_FRONTS_BRANCH:-3}   # weak branching: allow a few fronts but keep the original angle constraints
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-0}
BRANCH_SAVE_SNAPSHOTS=${BRANCH_SAVE_SNAPSHOTS:-12}
RUN_BRANCH=${RUN_BRANCH:-1}
RUN_NOW=${RUN_NOW:-1}

PROJECT_ROOT=$(pwd)
mkdir -p "$OUTROOT"

echo "python: $PYTHON_BIN"
echo "class:  $CLASS"
echo "T_K:    $T_K"
echo "thetas: $THETAS"
echo "out:    $OUTROOT"

"$PYTHON_BIN" - <<'PY'
from __future__ import annotations
import os, shlex, json, math, re
from pathlib import Path

project = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
klass = os.environ["CLASS"]
outroot = (project / os.environ["OUTROOT"]).resolve()
pybin = os.environ["PYTHON_BIN"]
T = os.environ["T_K"]
thetas = os.environ["THETAS"].split()
branch_theta = os.environ["BRANCH_THETA"]
seed = os.environ["SEED"]
target_ext = os.environ["TARGET_EXT_UM"]
steps = os.environ["STEPS"]
print_every = os.environ["PRINT_EVERY"]
max_fronts_base = os.environ["MAX_FRONTS_BASE"]
max_fronts_branch = os.environ["MAX_FRONTS_BRANCH"]
save_snapshots = os.environ["SAVE_SNAPSHOTS"]
branch_save_snapshots = os.environ["BRANCH_SAVE_SNAPSHOTS"]
run_branch = os.environ.get("RUN_BRANCH", "1") not in {"0", "false", "False", "no", "NO"}


def split_command_from_script(path: Path) -> list[str]:
    txt = path.read_text()
    for raw in txt.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "run_seeded_sharp_front.py" in line or "sharp_front" in line:
            line = re.sub(r"\s*>\s*[^\s]+\s+2>&1\s*$", "", line)
            line = line.split(">", 1)[0].strip()
            return shlex.split(line)
    raise RuntimeError(f"Could not find solver command in {path}")


def find_template() -> list[str]:
    candidates = [
        project / f"runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/{klass}/replicate_01_seed1101/T500_th45/command.txt",
        project / f"runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/{klass}/replicate_01_seed1101/T500_th45/run.log",
        project / f"runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x/{klass}/T300_th45/command.txt",
        project / f"runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x/{klass}/T500_th45/command.txt",
    ]
    for c in candidates:
        if c.exists() and c.name == "command.txt":
            return shlex.split(c.read_text())
    script_candidates = list((project / f"runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/{klass}/replicate_01_seed1101/T500_th45").glob("run*.sh"))
    for c in script_candidates:
        try:
            return split_command_from_script(c)
        except Exception:
            pass
    raise SystemExit(
        "Could not find a template command for class " + klass + ".\n"
        "Expected the 500K five-replicate campaign or the rate-sweep campaign to exist."
    )


def setopt(cmd: list[str], opt: str, val: str) -> list[str]:
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


def set_seed(cmd: list[str], val: str) -> list[str]:
    cmd = list(cmd)
    for opt in ["--solver-seed", "--seed", "--random-seed"]:
        if opt in cmd:
            cmd[cmd.index(opt) + 1] = str(val)
            return cmd
    cmd.extend(["--solver-seed", str(val)])
    return cmd


def case_complete(case: Path) -> bool:
    # Simple, robust completion test. Used only for skip/resume.
    log = case / "run.log"
    if log.exists():
        try:
            if "reached target crack extension" in log.read_text(errors="ignore"):
                return True
        except Exception:
            pass
    sj = case / "summary.json"
    if sj.exists():
        try:
            s = json.loads(sj.read_text())
            if isinstance(s, list) and s:
                s = s[0]
            ext = None
            if s.get("a_final_mm") is not None:
                ext = (float(s["a_final_mm"]) - 0.5) * 1000.0
            for k in ["final_crack_extension_um", "crack_extension_um", "extension_um"]:
                if s.get(k) is not None:
                    ext = max(ext or -math.inf, float(s[k]))
            return ext is not None and ext >= 0.98 * float(target_ext)
        except Exception:
            return False
    return False


template = find_template()
# Replace interpreter with requested environment and force unbuffered output.
# Template may already include -u; normalize to pybin -u script ...
if len(template) < 2:
    raise SystemExit("Template command too short")
script_idx = None
for i, item in enumerate(template):
    if item.endswith("run_seeded_sharp_front.py") or "run_seeded_sharp_front.py" in item:
        script_idx = i
        break
if script_idx is None:
    # fallback: keep original after interpreter
    script_idx = 1
base_tail = template[script_idx:]
base = [pybin, "-u"] + base_tail

jobs = []
for th in thetas:
    case_rel = Path(klass) / f"theta_{float(th):05.1f}_nobranch_seed{seed}" / f"T{T}_th{float(th):g}"
    jobs.append({"theta": th, "branch": False, "case": outroot / case_rel, "max_fronts": max_fronts_base, "snapshots": save_snapshots})

if run_branch:
    th = branch_theta
    case_rel = Path(klass) / f"theta_{float(th):05.1f}_weakbranch_seed{seed}" / f"T{T}_th{float(th):g}"
    jobs.append({"theta": th, "branch": True, "case": outroot / case_rel, "max_fronts": max_fronts_branch, "snapshots": branch_save_snapshots})

manifest = []
for job in jobs:
    case = job["case"]
    case.mkdir(parents=True, exist_ok=True)
    cmd = list(base)
    cmd = set_seed(cmd, seed)
    cmd = setopt(cmd, "--temperatures", T)
    cmd = setopt(cmd, "--crystal-theta-deg", job["theta"])
    cmd = setopt(cmd, "--target-crack-extension-um", target_ext)
    cmd = setopt(cmd, "--steps", steps)
    cmd = setopt(cmd, "--print-every", print_every)
    cmd = setopt(cmd, "--max-fronts", job["max_fronts"])
    cmd = setopt(cmd, "--out", str(case))
    cmd = setopt(cmd, "--save-snapshots", str(job["snapshots"]))
    cmd = ensure_flag(cmd, "--no-plots")

    run = case / "run_case.sh"
    run.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"cd {shlex.quote(str(project))}\n"
        + shlex.join(cmd) + f" > {shlex.quote(str(case / 'run.log'))} 2>&1\n"
    )
    run.chmod(0o755)
    manifest.append({
        "case": str(case.relative_to(project)),
        "theta_deg": float(job["theta"]),
        "branching": bool(job["branch"]),
        "max_fronts": int(job["max_fronts"]),
        "seed": int(seed),
        "temperature_K": int(float(T)),
        "run_script": str(run.relative_to(project)),
        "complete_now": case_complete(case),
    })

(outroot / "sweep_manifest.json").write_text(json.dumps(manifest, indent=2))
print(f"WROTE {outroot / 'sweep_manifest.json'}")
for m in manifest:
    status = "SKIP-ready" if m["complete_now"] else "pending"
    print(f"{status:10s} theta={m['theta_deg']:5.1f} branch={m['branching']} max_fronts={m['max_fronts']} -> {m['case']}")
PY

if [[ "$RUN_NOW" == "1" ]]; then
  echo "=== running cases sequentially ==="
  "$PYTHON_BIN" - <<'PY'
import json, subprocess, sys
from pathlib import Path
manifest = json.loads(Path(__import__('os').environ['OUTROOT']).joinpath('sweep_manifest.json').read_text())
for m in manifest:
    case = Path(m['case'])
    run = Path(m['run_script'])
    if m.get('complete_now'):
        print(f"SKIP complete: {case}", flush=True)
        continue
    print(f"START theta={m['theta_deg']} branch={m['branching']} -> {case}", flush=True)
    rc = subprocess.call(['bash', str(run)])
    if rc != 0:
        print(f"FAILED rc={rc}: {case}", flush=True)
        sys.exit(rc)
    print(f"DONE  theta={m['theta_deg']} branch={m['branching']} -> {case}", flush=True)
PY
else
  echo "RUN_NOW=0; generated run_case.sh scripts only."
fi
