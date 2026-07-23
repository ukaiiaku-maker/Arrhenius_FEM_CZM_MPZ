# v10.0.5.14 — PF v10.2.22 persistent-site parity for FEM/CZM

## Scope

This version ports the front-local plastic/source closure from
`PF-fracture-fatigue` commit
`198ece3aeb1d193a8c1c4857676fba720c088d27` onto the validated
v10.0.5.13.5 full 2-D FEM/CZM mechanics.

Preserved unchanged:

- plane-strain FEM and cubic anisotropic elasticity;
- domain-integral J and equivalent KJ;
- 330/240/100 micrometre mesh/J policy;
- v10.0.5.13.5 long-corridor Euclidean node deduplication;
- adaptive-CZM topology transaction and quality vetoes;
- stochastic cleavage/first-passage law inherited by the selected entry;
- elastic continuum bulk under `bulk_plasticity_mode=tip_only`.

The changed component is the moving crack-tip MPZ constitutive state.

## Ported physics

- two reduced BCC slip channels;
- separate positive/negative mobile, retained, and accumulated-slip fields;
- persistent areal source density `rho_source0_m2`;
- no source depletion, temporal source recovery, or crack-advance refresh;
- geometry-controlled multiplicity
  `M_s = rho_source0 * c_arc * r_tip * w_eff`;
- physical along-front width with lower bound
  `max(minimum_physical_width, b)`, never MPZ `dx` or cohesive size;
- implicit backward-Euler emission with mechanical-blocking complementarity;
- unsigned mobile+retained density for Taylor backstress;
- signed retained content for shielding (`mobile_shield_fraction=0`);
- Peierls transport, encounter storage, Taylor release, and zero explicit recovery;
- fractional active-to-wake translation of mobile, retained, and slip fields;
- natural resharpening from translated accumulated slip;
- trial/commit state copies and restart serialization.

The exact five v10.2.22 rows are included in
`arrhenius_fracture/data/mpz_v10_0_5_14/`.

## Signed shielding operator

The handoff does not contain the numerical signed shielding influence operator or
activation-to-line conversion. Production therefore fails closed unless
`--signed-shielding-kernel` points to a mechanically derived JSON artifact.

Accepted schemas:

- `v10.2.5_2d_unit_signed_shielding_kernel` from the PF signed-kernel workflow;
- `v10.0.5.14_signed_shielding_kernel` with the same required fields.

Required metadata:

```json
{
  "candidate_independent": true,
  "counts_are_signed_burgers_lines": true,
  "normalization_is_mechanically_derived": true
}
```

The active kernel must have shape `2 x 80`; one positive
`activation_to_line_content_by_system` value is required for each channel. The
FEM/CZM port does not create an analytical or fitted production substitute.

## Install

```bash
cd /Volumes/Data/Data/Nanopillar_calculation

DEST=Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites

git clone \
  --branch v10.0.5.14-persistent-sites-v10222-parity \
  --single-branch \
  https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git \
  "$DEST"

cd "$DEST"
conda activate arrhenius-fem-czm
python -m pip install -e "$PWD" --no-deps
```

## Focused verification

```bash
python -m pytest -q \
  tests/test_persistent_site_v100514.py \
  tests/test_v100513_barrier_only.py \
  tests/test_v1005131_preserved_state.py \
  tests/test_v1005132_startup_resolution_warning.py \
  tests/test_v1005133_tip_only_ramp.py \
  tests/test_v1005134_tip_only_policy_propagation.py \
  tests/test_v1005135_long_corridor.py \
  tests/test_v1005123_phase_c_repairs.py
```

## Candidate 0118 smoke command

```bash
KERNEL=/absolute/path/to/mechanically_derived_signed_kernel.json
OUT=/absolute/path/to/runs/v10_0_5_14_0118_700K_20um_smoke_v1

python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_14_persistent_site \
  --persistent-site-option v912_peak_0118_persistent_sites \
  --signed-shielding-kernel "$KERNEL" \
  --tip-refinement-radius-um 330 \
  --selected-cluster-J-outer-um 240 \
  --local-J-outer-um 100 \
  --mode 2d \
  --bulk-plasticity-mode tip_only \
  --temperatures 700 \
  --steps 20000 \
  --nx 36 --ny 72 \
  --tip-h-fine 2.5e-6 --tip-ratio 1.15 \
  --dU 2e-5 --dt 840 --n-stagger 1 --print-every 100 \
  --adaptive-events --adaptive-event-target 0.05 \
  --adaptive-min-frac 1e-8 --adaptive-grow 4 \
  --da-phys 5e-6 --target-crack-extension-um 20 \
  --crystal-aniso --crystal-compete --crystal-theta-deg 45 \
  --crystal-C11 523e9 --crystal-C12 203e9 --crystal-C44 160e9 \
  --cleave-gamma-aniso 0.3 --crystal-material w \
  --max-fronts 1 --crack-backend adaptive_czm \
  --czm-max-angle-error-deg 35 --j-decomposition cluster \
  --mpz-length-um 50 --mpz-n-bins 80 \
  --save-snapshots 2 --snapshot-cols 2 \
  --snapshot-by-crack-extension-um 10 --no-plots --out "$OUT"
```

## Mandatory runtime invariants

Every accepted increment and production manifest must record:

- `available_site_fraction = 1`;
- `persistent_source_inventory_active = false`;
- `source_depletion_active = false`;
- `source_refresh_active = false`;
- `source_sites_refreshed = 0`;
- `front_width_grid_independent = true`;
- `ahead_of_tip_dx_used_as_front_width_floor = false`;
- `two_channel_drive_reliable = true`.

A failed tensor probe, kernel mismatch, nonzero transferred recovery, or invalid
persistent-source row aborts the run.
