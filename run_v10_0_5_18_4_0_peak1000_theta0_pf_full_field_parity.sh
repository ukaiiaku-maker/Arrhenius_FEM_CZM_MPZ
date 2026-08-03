#!/usr/bin/env bash
set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
PF_REFERENCE_ROOT=${PF_REFERENCE_ROOT:-/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1}
PARAMETER_SOURCE_ROOT=${PARAMETER_SOURCE_ROOT:-$ROOT/runtime_inputs/v10_0_5_18_four_class}
OPTION=v913_paper_peak01_0242980_persistent_sites
PF_CASE=${PF_CASE:-$PF_REFERENCE_ROOT/runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1/$OPTION/T1000K_th0_seed8666}
FAMILY_JSON=${FAMILY_JSON:-$PF_REFERENCE_ROOT/runs/v10_2_28_kernel_cache/1447653d199f0b43cb475951092d69444c9b785f6fdf518c723792abb3b1f5e5/family.json}
FAMILY_SHA_REQUIRED=a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a
CAMPAIGN_ROOT=${CAMPAIGN_ROOT:-$ROOT/runs/v10_0_5_18_4_0_peak1000_theta0_seed8666_pf_full_field_parity_v1}
OUT=$CAMPAIGN_ROOT/$OPTION/T1000K
LOG=$OUT/console.log
ANALYSIS_OUT=$CAMPAIGN_ROOT/parity_analysis
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
STEPS=${STEPS:-2000000}
PRINT_EVERY=${PRINT_EVERY:-200}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-20}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
GENERATE_SOLVER_PLOTS=${GENERATE_SOLVER_PLOTS:-1}

for path in \
  "$PARAMETER_SOURCE_ROOT/v10_2_27_v913_four_class_paper_registry.csv" \
  "$PARAMETER_SOURCE_ROOT/four_class_transfer_manifest_v10_0_5_18.json" \
  "$FAMILY_JSON" \
  "$PF_CASE/steps_1000K.csv"
do
  [[ -f "$path" ]] || { echo "ERROR: missing parity input: $path" >&2; exit 1; }
done

FAMILY_JSON=$(cd "$(dirname "$FAMILY_JSON")" && pwd)/$(basename "$FAMILY_JSON")
FAMILY_SHA=$(python - "$FAMILY_JSON" <<'PY'
import hashlib, pathlib, sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
[[ "$FAMILY_SHA" == "$FAMILY_SHA_REQUIRED" ]] || {
  echo "ERROR: PF kernel SHA mismatch" >&2
  echo "Expected: $FAMILY_SHA_REQUIRED" >&2
  echo "Observed: $FAMILY_SHA" >&2
  exit 1
}

GIT_HEAD=$(git -C "$ROOT" rev-parse HEAD)
GIT_BRANCH=$(git -C "$ROOT" branch --show-current)
PACKAGE_VERSION=$(python - <<'PY'
from importlib.metadata import version
print(version('arrhenius-fem-czm'))
PY
)
[[ "$PACKAGE_VERSION" == "10.0.5.18.4.0" ]] || {
  echo "ERROR: expected arrhenius-fem-czm 10.0.5.18.4.0; observed $PACKAGE_VERSION" >&2
  exit 1
}

rm -rf "$OUT" "$ANALYSIS_OUT"
mkdir -p "$OUT"
cat > "$CAMPAIGN_ROOT/campaign_configuration.txt" <<EOF
release=10.0.5.18.4.0
entry=arrhenius_fracture.mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production
git_branch=$GIT_BRANCH
git_head=$GIT_HEAD
package_version=$PACKAGE_VERSION
reference_case=$PF_CASE
parameter_option=$OPTION
temperature_K=1000
theta_deg=0
hazard_seed=8666
target_projected_extension_um=$TARGET_EXT_UM
nx=36
ny=72
dU_m=2e-7
dt_s=8.4
n_stagger=2
tip_h_fine_m=1e-6
tip_ratio=1.20
da_phys_m=5e-6
adaptive_event_target=0.15
signed_kernel_family=$FAMILY_JSON
signed_kernel_family_sha256=$FAMILY_SHA
PF_bulk_plasticity_mode=full_field
FEM_bulk_plasticity_mode=full_field
inherited_solver_bulk_mode_token=bulk_same_pt_km
bulk_model=v10.4.1_bulk_peierls_taylor_detailed_balance_exact_port
PF_crack_backend=sharp_wake
FEM_crack_backend=adaptive_czm_atomic_path_corridor
EOF

# macOS ships Bash 3.2. Under `set -u`, expanding an empty array with
# "${array[@]}" raises an unbound-variable error. Use an initialized scalar
# whose parameter expansion contributes either one flag or no word.
PLOT_FLAG=""
if [[ "$GENERATE_SOLVER_PLOTS" == 0 ]]; then
  PLOT_FLAG="--no-plots"
fi

export ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM="$TARGET_EXT_UM"
export ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY=${ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
export ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO=${ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}
export ARRHENIUS_CORRIDOR_MAX_PATH_SUBSEGMENTS=${ARRHENIUS_CORRIDOR_MAX_PATH_SUBSEGMENTS:-16}
export MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}

echo "[START] option=$OPTION T=1000K theta=0 seed=8666 bulk=full_field"
echo "Git head:     $GIT_HEAD"
echo "PF reference: $PF_CASE"
echo "PF kernel:    $FAMILY_JSON"
echo "Output:       $OUT"

env \
  PYTHONUNBUFFERED=1 \
  CLEAVAGE_HAZARD_MODE=exponential \
  CLEAVAGE_HAZARD_SEED=8666 \
  CLEAVAGE_HAZARD_MIN_THRESHOLD=1e-12 \
  CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled \
  CLEAVAGE_EVENT_MIN_FACTOR=0.5 \
  CLEAVAGE_EVENT_MAX_FACTOR=4.0 \
  EMISSION_HAZARD_SEED=8666 \
  EMISSION_HAZARD_MIN_THRESHOLD=1e-12 \
  EMISSION_MAX_ACTION_SUBSTEP=0.05 \
  EMISSION_MAX_EVENTS_PER_HALF_STEP=10000 \
  EMISSION_RAMP_MAX_LOG_RATE_CHANGE=0.25 \
  EMISSION_RAMP_ACTION_FLOOR=1e-6 \
  EMISSION_INNER_MAX_LOG_HAZARD_CHANGE=0.25 \
  EMISSION_INNER_ACTION_FLOOR=1e-6 \
  EMISSION_EVENT_HORIZON_FACTOR=1.25 \
  EMISSION_INNER_MAX_REFINEMENTS=24 \
  TWO_CHANNEL_PROBE_FALLBACK_MAX_CALLS=8 \
  python -u -m arrhenius_fracture.mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_production \
    --parameter-source-root "$PARAMETER_SOURCE_ROOT" \
    --parameter-option "$OPTION" \
    --signed-kernel-family "$FAMILY_JSON" \
    --tip-refinement-radius-um 330 \
    --selected-cluster-J-outer-um 240 \
    --local-J-outer-um 100 \
    --mode 2d \
    --bulk-plasticity-mode full_field \
    --temperatures 1000 \
    --steps "$STEPS" \
    --nx 36 --ny 72 \
    --tip-h-fine 1e-6 --tip-ratio 1.20 \
    --dU 2e-7 --dt 8.4 \
    --n-stagger 2 \
    --print-every "$PRINT_EVERY" \
    --adaptive-events \
    --adaptive-event-target 0.15 \
    --adaptive-min-frac 1e-8 \
    --adaptive-grow 4 \
    --da-phys 5e-6 \
    --target-crack-extension-um "$TARGET_EXT_UM" \
    --crystal-aniso \
    --crystal-compete \
    --crystal-theta-deg 0 \
    --min-global-forward "$MIN_GLOBAL_FORWARD" \
    --crystal-C11 523e9 --crystal-C12 203e9 --crystal-C44 160e9 \
    --cleave-gamma-aniso 0.3 \
    --crystal-material w \
    --max-fronts 1 \
    --crack-backend adaptive_czm \
    --czm-max-angle-error-deg 35 \
    --j-decomposition cluster \
    --mpz-length-um 50 --mpz-n-bins 80 \
    --save-snapshots "$SAVE_SNAPSHOTS" \
    --snapshot-cols "$SNAPSHOT_COLS" \
    --snapshot-by-crack-extension-um "$SNAPSHOT_BY_EXT_UM" \
    ${PLOT_FLAG:+"$PLOT_FLAG"} \
    --out "$OUT" \
    2>&1 | tee "$LOG"

python "$ROOT/scripts/compare_v10051840_theta0_pf_parity.py" \
  --pf-case "$PF_CASE" \
  --fem-case "$OUT" \
  --out "$ANALYSIS_OUT"

echo "[DONE] option=$OPTION T=1000K theta=0 seed=8666 bulk=full_field"
echo "Parity analysis: $ANALYSIS_OUT"
