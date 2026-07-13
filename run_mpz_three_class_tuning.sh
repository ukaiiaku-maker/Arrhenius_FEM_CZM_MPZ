#!/usr/bin/env bash
set -euo pipefail

# Staged three-class moving-process-zone calibration.
# The stages are restartable and write to distinct directories.

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
STAGE=${STAGE:-smoke}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
ROOT=${ROOT:-runs/mpz_v9_1_three_class_tuning}
FIRST_OUT=${FIRST_OUT:-$ROOT/01_first_passage}
RCURVE_OUT=${RCURVE_OUT:-$ROOT/02_rcurve}
JOINT_OUT=${JOINT_OUT:-$ROOT/03_joint}
INITIAL=${INITIAL:-mpz_three_class_initial_guesses.csv}
TARGETS=${TARGETS:-mpz_three_class_design_targets.csv}
PYTHON_BIN=${PYTHON_BIN:-python}

# Ensure Python messages are written promptly when stdout is redirected by nohup.
export PYTHONUNBUFFERED=${PYTHONUNBUFFERED:-1}

stamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

run_py() {
  echo "[$(stamp)] launching: $PYTHON_BIN -u $*"
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$CONDA_ENV" --no-capture-output "$PYTHON_BIN" -u "$@"
  else
    "$PYTHON_BIN" -u "$@"
  fi
}

echo "[$(stamp)] MPZ v9.1 tuning start"
echo "[$(stamp)] stage=$STAGE classes=$CLASSES conda_env=$CONDA_ENV root=$ROOT"
echo "[$(stamp)] initial=$INITIAL targets=$TARGETS"

case "$STAGE" in
  smoke)
    echo "[$(stamp)] smoke output: $ROOT/00_smoke"
    run_py fit_mpz_three_classes.py \
      --smoke \
      --initial "$INITIAL" \
      --targets "$TARGETS" \
      --classes "$CLASSES" \
      --temperatures "300 700 1100" \
      --objective-mode joint \
      --fit-profile reduced \
      --dK 0.25 \
      --da-um 25 \
      --n-advances 41 \
      --early-window-um "25 225" \
      --plateau-window-um "700 1000" \
      --target-dB-substep 0.8 \
      --target-emission-hazard-substep 5 \
      --source-active-fraction-min 1e-3 \
      --out "$ROOT/00_smoke"
    ;;

  first)
    echo "[$(stamp)] first-passage output: $FIRST_OUT popsize=${POPSIZE:-6} maxiter=${MAXITER:-25}"
    run_py fit_mpz_three_classes.py \
      --initial "$INITIAL" \
      --targets "$TARGETS" \
      --classes "$CLASSES" \
      --temperatures "300 500 700 800 900 1100 1200" \
      --objective-mode first \
      --fit-profile first \
      --optimizer de \
      --dK 0.25 \
      --da-um 5 \
      --n-advances 1 \
      --target-dB-substep 0.8 \
      --target-emission-hazard-substep 5 \
      --source-active-fraction-min 1e-3 \
      --popsize "${POPSIZE:-6}" \
      --maxiter "${MAXITER:-25}" \
      --resume --skip-completed \
      --out "$FIRST_OUT"
    ;;

  rcurve)
    test -f "$FIRST_OUT/mpz_three_class_parameters.csv" || {
      echo "Missing Stage 1 parameters: $FIRST_OUT/mpz_three_class_parameters.csv" >&2
      exit 2
    }
    echo "[$(stamp)] R-curve output: $RCURVE_OUT popsize=${POPSIZE:-5} maxiter=${MAXITER:-20}"
    run_py fit_mpz_three_classes.py \
      --initial "$FIRST_OUT/mpz_three_class_parameters.csv" \
      --targets "$TARGETS" \
      --classes "$CLASSES" \
      --temperatures "300 700 900 1100" \
      --objective-mode rcurve \
      --fit-profile rcurve \
      --optimizer de \
      --dK 0.25 \
      --da-um 25 \
      --n-advances 41 \
      --early-window-um "25 225" \
      --plateau-window-um "700 1000" \
      --target-dB-substep 0.8 \
      --target-emission-hazard-substep 5 \
      --source-active-fraction-min 1e-3 \
      --popsize "${POPSIZE:-5}" \
      --maxiter "${MAXITER:-20}" \
      --resume --skip-completed \
      --out "$RCURVE_OUT"
    ;;

  joint)
    test -f "$RCURVE_OUT/mpz_three_class_parameters.csv" || {
      echo "Missing Stage 2 parameters: $RCURVE_OUT/mpz_three_class_parameters.csv" >&2
      exit 2
    }
    echo "[$(stamp)] joint output: $JOINT_OUT maxiter=${MAXITER:-20} maxfev=${MAXFEV:-1200}"
    run_py fit_mpz_three_classes.py \
      --initial "$RCURVE_OUT/mpz_three_class_parameters.csv" \
      --targets "$TARGETS" \
      --classes "$CLASSES" \
      --temperatures "300 500 700 800 900 1100 1200" \
      --objective-mode joint \
      --fit-profile full \
      --optimizer powell \
      --dK 0.20 \
      --da-um 20 \
      --n-advances 51 \
      --early-window-um "20 220" \
      --plateau-window-um "700 1000" \
      --target-dB-substep 0.5 \
      --target-emission-hazard-substep 2 \
      --source-active-fraction-min 3e-4 \
      --maxiter "${MAXITER:-20}" \
      --maxfev "${MAXFEV:-1200}" \
      --resume --skip-completed \
      --out "$JOINT_OUT"
    ;;

  verify)
    test -f "$JOINT_OUT/mpz_three_class_parameters.csv" || {
      echo "Missing Stage 3 parameters: $JOINT_OUT/mpz_three_class_parameters.csv" >&2
      exit 2
    }
    echo "[$(stamp)] convergence output: $ROOT/04_convergence"
    run_py audit_mpz_three_class_convergence.py \
      --parameters "$JOINT_OUT/mpz_three_class_parameters.csv" \
      --classes "$CLASSES" \
      --temperatures "300 400 500 600 700 800 900 1000 1100 1200" \
      --dK-values "0.20 0.10 0.05" \
      --bin-counts "50 100 200" \
      --da-values-um "10 5" \
      --target-extension-um 1000 \
      --out "$ROOT/04_convergence"
    ;;

  *)
    echo "Unknown STAGE=$STAGE; use smoke, first, rcurve, joint, or verify" >&2
    exit 2
    ;;
esac

echo "[$(stamp)] MPZ v9.1 tuning stage completed: $STAGE"
