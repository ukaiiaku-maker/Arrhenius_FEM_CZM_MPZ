#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
OUTROOT="${OUTROOT:-runs/mpz_v9_verify}"
RUN_PROTOCOL_SMOKES="${RUN_PROTOCOL_SMOKES:-0}"
RUN_2D_SMOKE="${RUN_2D_SMOKE:-0}"

rm -rf "$OUTROOT"
mkdir -p "$OUTROOT"

"$PYTHON_BIN" -m compileall -q \
  arrhenius_fracture \
  audit_legacy_caps_and_ablations.py \
  audit_mpz_three_class_convergence.py \
  build_mpz_analytic_first_passage_atlas.py \
  search_mpz_peierls_taylor_parameters.py \
  fit_mpz_four_classes.py \
  fit_mpz_three_classes.py \
  mpz_run_utils.py \
  run_mpz_dwell.py \
  run_mpz_fatigue_matrix.py \
  run_mpz_fem_czm_validation_matrix.py

PYTHONPATH=. "$PYTHON_BIN" -m pytest -q \
  tests/test_moving_process_zone.py \
  tests/test_mpz_three_class_fit.py \
  tests/test_mpz_analytic_first_passage_atlas.py \
  tests/test_emission_derived_peierls_taylor.py \
  tests/test_bulk_pt_plasticity.py

"$PYTHON_BIN" - <<'PY'
import arrhenius_fracture as af
assert af.__version__ == '0.9.3'
print('package version:', af.__version__)
PY

if [[ "$RUN_PROTOCOL_SMOKES" == "1" ]]; then
  "$PYTHON_BIN" audit_legacy_caps_and_ablations.py \
    --classes ceramic --temperatures 300 --dK-values 0.25 \
    --ablations "baseline no_dN_cap" --Kmax 5 --n-advances 2 \
    --out "$OUTROOT/legacy_audit"

  "$PYTHON_BIN" fit_mpz_four_classes.py --smoke \
    --classes ceramic --temperatures 300 \
    --dK 0.5 --Kmax 5 --n-advances 1 \
    --out "$OUTROOT/fit_smoke"

  "$PYTHON_BIN" build_mpz_analytic_first_passage_atlas.py \
    --initial mpz_three_class_initial_guesses.csv \
    --shape-families "ceramic weakT DBTT" \
    --temperatures "300 700 1200" \
    --Kdot-values "0.005" \
    --samples-per-family 32 \
    --dK 0.25 --refine-dK 0.05 --Kmax 40 \
    --top-per-region 2 \
    --out "$OUTROOT/analytic_atlas_smoke"

  cat > "$OUTROOT/pt_search_smoke_input.csv" <<'CSV'
candidate_id,region,Kdot_MPa_sqrt_m_per_s,ceramic_score,weakT_score,DBTT_precursor_score,refined_Kc_T300,refined_Kc_T1200,emit_G00_eV,emit_gT_eV_per_K,emit_sigc0_GPa,emit_sT_GPa_per_K,emit_exp_a,emit_exp_n,emit_floor_frac
c1,ceramic_intrinsic,0.005,0.1,9,9,18,7,1.3922,0.0091197,2.9802,0.002031,0.2383,1.1169,0.041
w1,weakT_intrinsic,0.005,9,0.1,9,15,15,1.1173,0.0063971,0.9506,0.0009817,0.5055,0.8432,0.0309
d1,DBTT_precursor,0.005,9,9,0.1,15,10,1.6954,0.0008853,3.7521,0.0022809,0.0944,0.835,0.03524
CSV
  "$PYTHON_BIN" search_mpz_peierls_taylor_parameters.py \
    --atlas-shortlist "$OUTROOT/pt_search_smoke_input.csv" \
    --transport-samples 16 --intrinsic-top-per-region 1 \
    --top-per-intrinsic 1 --temperatures "300 700 1200" \
    --strain-rates "1e-5 1e-3" --rho-points 25 \
    --out "$OUTROOT/pt_search_smoke"

  "$PYTHON_BIN" run_mpz_fatigue_matrix.py \
    --classes ceramic --temperatures 300 --Kmax-values 8 \
    --cycles-max 1000 --block-cycles 10 --max-blocks 3 --n-advances 1 \
    --out "$OUTROOT/fatigue_smoke"

  "$PYTHON_BIN" run_mpz_dwell.py \
    --classes ceramic --temperatures 300 --K-MPa-sqrt-m 8 \
    --hold-s 0.001 --dt-initial-s 1e-5 --dt-max-s 1e-4 \
    --max-blocks 10000 --n-advances 1 \
    --out "$OUTROOT/dwell_smoke"
fi

if [[ "$RUN_2D_SMOKE" == "1" ]]; then
  "$PYTHON_BIN" run_mpz_fem_czm_validation_matrix.py \
    --classes ceramic --temperatures 300 --out "$OUTROOT/fem_czm_smoke" \
    --max-jobs 1 --nx 12 --ny 24 --steps 1 --n-stagger 1 \
    --print-every 1 --target-ext-um 1 --save-snapshots 1 \
    --snapshot-by-ext-um 0

  "$PYTHON_BIN" -m arrhenius_fracture.mixed_mode_first_passage_v8 \
    --target-traction-phase-deg 0 --reference-cleavage-shape 1 \
    --mixity-loading-angle-deg 0 --mode 2d --temperatures 300 \
    --out "$OUTROOT/mixed_mode_smoke" --nx 12 --ny 24 \
    --dU 2e-7 --dt 8.4 --steps 1 --n-stagger 1 --print-every 1 \
    --crystal-aniso --crystal-compete --crystal-material w \
    --crystal-theta-deg 45 --cleave-gamma-aniso 0.3 --max-fronts 1 \
    --front-state-model moving_pz --sigma-cap-GPa 0 --save-snapshots 0
fi

echo "MPZ v9.3 verification passed. Outputs: $OUTROOT"
