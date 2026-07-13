#!/usr/bin/env bash
# =============================================================================
# sweep_T_theta.sh  --  temperature x orientation sweep (continuous anisotropy)
# Run with:   bash sweep_T_theta.sh
# from the directory that CONTAINS the arrhenius_fracture/ package.
#
# Uses the competing-direction model (--crystal-compete): the crack advances
# along argmax of overdrive sigma_nn(phi)/sqrt(gamma(phi)) over a CONTINUUM of
# directions, with a finite cubic cleavage-energy anisotropy delta. Sweeps a few
# temperatures (DBTT axis) x a few orientations (deflection axis) and collates
# Kc(theta,T) and net-deflection(theta,T) into grids you can read at a glance.
# =============================================================================
set -u

# ---------------- CONFIG (tune) ----------------------------------------------
ROOT="Ttheta_$(date +%Y%m%d_%H%M)"
TLIST="300 500 700 900 1100"          # temperatures (DBTT axis)
THETAS="0 15 30 45"                   # orientations  (deflection axis)
DELTA=0.5                             # cleavage-energy anisotropy gamma110/gamma100-1
                                      #   0=isotropic(straight) .. large=lock to {100}; W~0.1-0.3
STEPS=220                             # captures Kc_first at all T + a deflected path
BRANCH=0                             # set to 1 to also allow branching (--crystal-branch)
# -----------------------------------------------------------------------------

mkdir -p "$ROOT/logs"
PKG="python3 -m arrhenius_fracture.sharp_front"
COMMON="--mode 2d --nx 50 --ny 100 --tip-h-fine 0.6e-6 --tip-ratio 1.25 --n-stagger 2 \
--crystal-aniso --crystal-compete --crystal-C44 320e9 \
--emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8 \
--multihit-m 3 --multihit-tau 1e-6 --cleave-H0-eV 3.0 --cleave-shield-chi 0.6 \
--emb-sat-frac 1 --n-sat 2000 --v-rayleigh 1.2e-7 \
--dU 2e-6 --dt 30 --save-snapshots 4 --print-every 9999"
BR=""; [ "$BRANCH" = "1" ] && BR="--crystal-branch"

echo "=== T x theta sweep -> $ROOT   (delta=$DELTA, branch=$BRANCH) ==="
for th in $THETAS; do
  out="$ROOT/theta_${th}"
  echo "[$(date +%H:%M)] >>> theta=$th"
  $PKG $COMMON --crystal-theta-deg $th --cleave-gamma-aniso $DELTA $BR \
    --temperatures $TLIST --steps $STEPS --out "$out" \
    > "$ROOT/logs/theta_${th}.log" 2>&1
  grep -E "Kc_first|BRANCH|no-fracture" "$ROOT/logs/theta_${th}.log" | sed 's/^/      /' || true
done

echo "=== collating grids -> $ROOT/grid_*.csv ==="
python3 - "$ROOT" "$THETAS" "$TLIST" << 'PYEOF'
import os, sys, json, glob
root, thetas, tlist = sys.argv[1], sys.argv[2].split(), sys.argv[3].split()
thetas = [float(t) for t in thetas]; tlist = [float(t) for t in tlist]
# load every summary.json: {(theta,T): rec}
data = {}
for d in sorted(glob.glob(os.path.join(root, "theta_*"))):
    th = float(os.path.basename(d).split("_")[1])
    try:
        recs = json.load(open(os.path.join(d, "summary.json")))
    except Exception:
        continue
    for r in recs:
        data[(th, float(r["T"]))] = r

def fmt(v, w=7, p=1):
    return ("{:>%d}" % w).format("--") if v is None else ("{:>%d.%df}" % (w, p)).format(v)

def grid(field, p=1, transform=lambda x: x):
    lines = ["  theta\\T   " + "".join(f"{int(T):>8}" for T in tlist)]
    for th in thetas:
        row = []
        for T in tlist:
            r = data.get((th, T))
            v = None if (r is None or r.get(field) is None) else transform(r[field])
            row.append(fmt(v, 8, p))
        lines.append(f"  {int(th):>5}    " + "".join(row))
    return "\n".join(lines)

# write a tidy long CSV
import csv
with open(os.path.join(root, "grid_long.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["theta_deg","T_K","Kc_MPa_sqrt_m","deflection_deg",
                                   "path_span_dy_mm","branched","mode"])
    for (th, T), r in sorted(data.items()):
        w.writerow([th, T, r.get("Kc_first_MPa_sqrt_m"), r.get("deflection_deg"),
                    r.get("path_span_dy_mm"), r.get("branched"), r.get("mode")])

print("\n  Kc_first [MPa sqrt(m)]  (DBTT down each column; orientation down each row)")
print(grid("Kc_first_MPa_sqrt_m", 1))
print("\n  net crack deflection [deg]  (0 = straight)")
print(grid("deflection_deg", 1))
print("\n  branched? (1=yes)")
print(grid("branched", 0, transform=lambda b: 1.0 if b else 0.0))
print(f"\n  long-form CSV -> {os.path.join(root,'grid_long.csv')}")
PYEOF

echo "=== DONE.  Also see $ROOT/theta_*/field_snapshots_*.png (geometry)"
echo "             and  $ROOT/theta_*/crack_path_*K.csv  (raw polylines)"
