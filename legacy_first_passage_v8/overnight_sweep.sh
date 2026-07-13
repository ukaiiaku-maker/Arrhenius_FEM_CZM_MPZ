#!/usr/bin/env bash
# =============================================================================
# overnight_sweep.sh  --  arrhenius_fracture exploratory sweep
# Run with:   bash overnight_sweep.sh
# from the directory that CONTAINS the arrhenius_fracture/ package.
# (Use `bash`, not zsh-paste, so inline # comments are safe.)
#
# Five experiments, each isolating ONE axis. All Kc(T)/mode/branch results are
# collated into  <ROOT>/master_summary.csv  at the end, and every run drops
# field_snapshots_<T>K.png (crack geometry) + toughness_vs_temperature.png.
# Tune the CONFIG block, or comment out whole EXP blocks, to fit your night.
# =============================================================================
set -u

# ---------------- CONFIG (tune these) ----------------------------------------
ROOT="sweep_$(date +%Y%m%d_%H%M)"
TFINE="300 400 500 600 700 800 900 1000 1100"   # DBTT temperature grid
TPAIR="300 700"                                  # brittle + near-transition, for orientation
STEPS_DBTT=300                                   # enough for Kc_first + some growth
STEPS_GEOM=320                                   # orientation (deflection) geometry
STEPS_BRANCH=400                                 # branching needs more growth
# -----------------------------------------------------------------------------

mkdir -p "$ROOT/logs"
PKG="python3 -m arrhenius_fracture.sharp_front"

# flags common to EVERY run (resolution, crystal elasticity, cleavage, loading)
COMMON="--mode 2d --nx 50 --ny 100 --tip-h-fine 0.6e-6 --tip-ratio 1.25 --n-stagger 2 \
--crystal-aniso --crystal-C44 320e9 --multihit-m 3 --multihit-tau 1e-6 \
--cleave-H0-eV 3.0 --cleave-shield-chi 0.6 --emb-sat-frac 1 --v-rayleigh 1.2e-7 \
--dU 2e-6 --save-snapshots 6 --print-every 9999"

# metal (W) emission entropy = the DBTT channel (cold-hard / hot-soft)
METAL="--emit-S-T-c0-kB=-20 --emit-S-T-c1=0.02 --emit-S-sigma-max-kB=8"

run () {                       # run LABEL META -- <extra flags...>
  local label="$1"; shift
  local meta="$1"; shift
  local out="$ROOT/$label"
  local log="$ROOT/logs/$label.log"
  echo "#META $meta label=$label" > "$log"
  echo "[$(date +%H:%M)] >>> $label   ($meta)"
  $PKG $COMMON --out "$out" "$@" >> "$log" 2>&1
  grep -E "Kc_first|BRANCH|severed|no-fracture" "$log" | sed 's/^/      /' || true
}

echo "=== overnight sweep -> $ROOT ==="

# === EXP A : DBTT curve x orientation (metal) ================================
# Q: shape of Kc(T); does deflection change the DBTT? (straight th=0 vs deflect th=30)
for th in 0 30; do
  run "A_dbtt_metal_th${th}" "exp=A regime=metal theta=${th}" $METAL \
    --crystal-theta-deg $th --temperatures $TFINE --dt 30 --steps $STEPS_DBTT
done

# === EXP B : orientation sweep -- deflection angle & branching onset =========
# Q: how the crack path changes with theta; where single-deflection -> branch.
# B1 single-front (no branch): climb angle vs theta, at a brittle + a warm T.
for th in 0 10 20 30 35 40 43 45; do
  run "B_orient_th${th}" "exp=B regime=metal theta=${th} branch=0" $METAL \
    --crystal-theta-deg $th --temperatures $TPAIR --dt 30 --steps $STEPS_GEOM
done
# B2 with branching enabled near the degeneracy (theta -> 45).
for th in 40 43 45; do
  run "B_branch_th${th}" "exp=B regime=metal theta=${th} branch=1" $METAL \
    --crystal-branch --crystal-theta-deg $th --temperatures $TPAIR --dt 30 --steps $STEPS_BRANCH
done

# === EXP C : material regime  ceramic <-> metal =============================
# Q: emission strength controls the DBTT. Raise the emission barrier H0 to push
# the ductile transition up/out -> ceramic-like (no upturn). Single axis = emit-H0.
for H0 in 1.8 2.6 3.4 4.0; do
  run "C_emitH0_${H0}" "exp=C theta=30 emitH0=${H0}" \
    --emit-H0-eV $H0 $METAL \
    --crystal-theta-deg 30 --temperatures $TFINE --dt 30 --steps $STEPS_DBTT
done
# explicit ceramic endpoint: emission frozen (no T-entropy ramp) -> pure cleavage.
run "C_ceramic_frozen" "exp=C theta=30 regime=ceramic" \
  --emit-H0-eV 4.0 --emit-S-T-c0-kB=-20 --emit-S-T-c1=0.0 --emit-S-sigma-max-kB=8 \
  --crystal-theta-deg 30 --temperatures $TFINE --dt 30 --steps $STEPS_DBTT

# === EXP D : cap / saturation  (is saturation the DBTT mechanism?) ===========
# Q: does lifting the 30 GPa cap / changing emission saturation move the shelf?
for cap in 30 60 200; do
  run "D_cap_${cap}" "exp=D theta=30 cap=${cap}" $METAL \
    --crystal-theta-deg 30 --sigma-cap-GPa $cap --temperatures $TFINE --dt 30 --steps $STEPS_DBTT
done
for ns in 500 2000 8000; do
  run "D_nsat_${ns}" "exp=D theta=30 nsat=${ns}" $METAL \
    --crystal-theta-deg 30 --n-sat $ns --temperatures $TFINE --dt 30 --steps $STEPS_DBTT
done

# === EXP E : loading-rate effect on the DBTT ================================
# Q: does the transition shift with rate? (Arrhenius: faster load -> higher DBTT)
# NOTE: dt also couples to loading granularity -> read TRENDS, not absolute Kc.
for dt in 10 30 100; do
  run "E_rate_dt${dt}" "exp=E theta=30 dt=${dt}" $METAL \
    --crystal-theta-deg 30 --temperatures $TFINE --dt $dt --steps $STEPS_DBTT
done

echo "=== runs done; collating -> $ROOT/master_summary.csv ==="

# === collation: parse every log into one tidy CSV ===========================
python3 - "$ROOT" << 'PYEOF'
import os, re, sys, glob, csv
root = sys.argv[1]
rows = []
for log in sorted(glob.glob(os.path.join(root, "logs", "*.log"))):
    meta = {}
    text = open(log, errors="ignore").read()
    m = re.search(r'^#META (.*)$', text, re.M)
    if m:
        for kv in m.group(1).split():
            if '=' in kv:
                k, v = kv.split('=', 1); meta[k] = v
    branched = 'BRANCH at' in text
    for line in text.splitlines():
        r = re.search(r'T=(\d+)K:\s*Kc_first=([\d.]+|none).*advances=(\d+).*mode=(\S+)', line)
        if r:
            T, kc, adv, mode = r.groups()
            rows.append({**meta, 'T_K': int(T),
                         'Kc_first': (None if kc == 'none' else float(kc)),
                         'advances': int(adv), 'mode': mode,
                         'branched': int(branched)})
cols = ['exp','label','regime','theta','emitH0','cap','nsat','dt','branch',
        'T_K','Kc_first','advances','mode','branched']
out = os.path.join(root, 'master_summary.csv')
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore'); w.writeheader()
    for r in rows: w.writerow(r)
print(f"  wrote {len(rows)} rows -> {out}")
# compact console digest: Kc(T) per label
from collections import defaultdict
by = defaultdict(dict)
for r in rows: by[r['label']][r['T_K']] = r['Kc_first']
for lab in sorted(by):
    ts = sorted(by[lab])
    s = "  ".join(f"{t}:{('%.1f'%by[lab][t]) if by[lab][t] is not None else '--':>5}" for t in ts)
    print(f"  {lab:22s} {s}")
PYEOF

echo "=== DONE.  Inspect:"
echo "    $ROOT/master_summary.csv            (all Kc(T), modes, branch flags)"
echo "    $ROOT/<label>/field_snapshots_*.png (crack geometry per condition)"
echo "    $ROOT/<label>/toughness_vs_temperature.png"
