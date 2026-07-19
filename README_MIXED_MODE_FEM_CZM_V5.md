# Mixed-mode FEM/CZM v5 — anisotropic calibrated-tip first passage

## Purpose

V5 corrects the stress-scale mismatch in v4.

V4 inserted the absolute finite-radius FEM traction at a 10 micrometer probe
radius into Arrhenius barriers that had been calibrated with the sharp-front
process-zone stress.  The resulting cleavage clocks were extremely small, all
cases were right-censored, and different material parameterizations shared the
same elastic endpoint.

V5 restores the validated sharp-front magnitude:

    KJ -> directionally scaled K drive -> existing sigma_tip(K, r_eff)

The anisotropic FEM field is used only to determine dimensionless directional
weights, local mode phase, crystallographic crack direction, and slip-system
resolution.

The implementation is exclusively sharp-interface cohesive mechanics.

## Kinetic coupling

For the current anisotropic stress tensor, v5 computes:

* a cleavage shape measure from maximum crystallographic opening overdrive;
* a slip shape measure from maximum resolved BCC slip shear;
* a Mode-I reference state at the same crystal orientation and probe radius.

The scalar driving amplitudes are

    K_cleave = F_cleave KJ

    K_emit = F_emit KJ

where F_cleave=F_emit=1 in the calibrated zero-phase reference state.  The
existing sharp-front engine then applies its original process-zone radius,
blunting, back stress, shielding, EXP-floor barriers, multi-hit closure, and
adaptive first-passage integration.

Finite-radius traction never replaces the absolute sharp-tip stress scale.

## Event-state mode control

The elastic basis gives an initial boundary loading angle.  Each physical case
then measures the actual process-zone traction phase at first passage or the
right-censored endpoint.  The boundary angle is corrected iteratively until the
phase error is within `EVENT_PSI_TOL_DEG` or `MAX_CONTROL_ITERS` is exhausted.

Each condition writes `mixed_mode_control_history_v5.csv`.

## Material-parameter audit

Every iteration writes `barrier_audit.json`, including a SHA-256 fingerprint of
the class-specific cleavage, emission, shielding, and saturation parameters.
The final campaign CSV carries the same fingerprint.  This makes accidental
class duplication directly visible.

## Installation

Copy the package contents into the project root:

    /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM

The files use new v5 names and do not overwrite v4.

## Verification

    chmod +x verify_mixed_mode_fem_czm_v5.sh
    chmod +x run_mixed_mode_fem_czm_v5_campaign.sh

    CONDA_ENV=arrhenius-fem-czm \
    bash verify_mixed_mode_fem_czm_v5.sh

Expected: 15 tests, all passing, followed by

    MIXED_MODE_V5 verification OK

## Recommended preflight

Run both classes so class separation is tested immediately:

    rm -rf runs/mixed_mode_fem_czm_v5_preflight_cal
    rm -rf runs/mixed_mode_fem_czm_v5_preflight

    CONDA_ENV=arrhenius-fem-czm \
    PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
    CLASSES="ceramic DBTT" \
    TARGET_PSI="-30 0 30" \
    T_K=500 \
    CRYSTAL_THETA_DEG=45 \
    TRACTION_PROBE_RADIUS_M=1e-5 \
    EVENT_PSI_TOL_DEG=2 \
    MAX_CONTROL_ITERS=3 \
    MAX_JOBS=1 \
    RECALIBRATE=1 \
    CALROOT=runs/mixed_mode_fem_czm_v5_preflight_cal \
    OUTROOT=runs/mixed_mode_fem_czm_v5_preflight \
    bash run_mixed_mode_fem_czm_v5_campaign.sh

## Acceptance checks

For the zero-phase calibration row:

    cleavage_factor = 1
    emission_factor = 1

For every selected physical result:

    traction_probe_reliable = True
    event_phase_control_converged = True
    directional_factor_cap_active = False

Ceramic and DBTT must have different `barrier_fingerprint_sha256` values and
should no longer share an endpoint merely because both failed to fracture.
Right-censored cases remain valid lower-bound observations but are never plotted
as first-passage toughnesses.

## Full first-passage campaign

    rm -rf runs/mixed_mode_fem_czm_v5_anisotropic_calibration_500K
    rm -rf runs/mixed_mode_fem_czm_v5_anisotropic_calibrated_tip_500K

    CONDA_ENV=arrhenius-fem-czm \
    PARAMETER_TABLE=four_class_exp_floor_exact_model_inputs.csv \
    CLASSES="ceramic DBTT" \
    TARGET_PSI="-60 -45 -30 -15 0 15 30 45 60" \
    T_K=500 \
    CRYSTAL_THETA_DEG=45 \
    TRACTION_PROBE_RADIUS_M=1e-5 \
    EVENT_PSI_TOL_DEG=2 \
    MAX_CONTROL_ITERS=4 \
    MAX_JOBS=1 \
    RECALIBRATE=1 \
    CALROOT=runs/mixed_mode_fem_czm_v5_anisotropic_calibration_500K \
    OUTROOT=runs/mixed_mode_fem_czm_v5_anisotropic_calibrated_tip_500K \
    bash run_mixed_mode_fem_czm_v5_campaign.sh

## Outputs

Campaign-level:

* `campaign_status_v5.csv`
* `mixed_mode_v5_anisotropic_all_cases.csv`
* `plots_v5/`

Condition-level:

* `mixed_mode_control_history_v5.csv`
* `anisotropic_calibrated_tip_final_summary.json`
* `anisotropic_calibrated_tip_final_summary.csv`

Iteration-level:

* `barrier_audit.json`
* `command.txt`
* `run.log`
* `anisotropic_calibrated_tip_calls.csv`
* `anisotropic_calibrated_tip_first_passage_summary.json`
* normal sharp-front/FEM/CZM CSV outputs

## Scope

This package is for deterministic first-passage screening.  It does not yet run
extended crack growth or threshold quantiles.  The later extension campaign
should retain this anisotropic coupling, run several cumulative-hazard
quantiles, and use the median crack-extension response.
