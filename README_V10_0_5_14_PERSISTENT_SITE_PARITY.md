# v10.0.5.14.1 — PF v10.2.22 persistent-site and kernel-family parity

## Scope

This point release ports the front-local plastic/source closure from
`PF-fracture-fatigue` commit
`198ece3aeb1d193a8c1c4857676fba720c088d27` onto the validated
v10.0.5.13.5 full 2-D FEM/CZM mechanics and consumes the actual PF v10.2.14
crack-extension-indexed signed shielding atlas.

Preserved unchanged:

- plane-strain FEM and cubic anisotropic elasticity;
- domain-integral J and equivalent KJ;
- 330/240/100 micrometre mesh/J policy;
- v10.0.5.13.5 long-corridor Euclidean node deduplication;
- adaptive-CZM topology transaction and quality vetoes;
- stochastic cleavage/first-passage law inherited by the selected entry;
- elastic continuum bulk under `bulk_plasticity_mode=tip_only`.

The changed component is the moving crack-tip MPZ constitutive state and its
signed shielding artifact loader.

## Ported persistent-site physics

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

## PF signed-kernel family

The production input is the actual candidate-independent PF atlas:

```text
schema = v10.2.14_active_only_real_signed_2d_shielding_atlas
```

It contains mechanically measured signed-Burgers kernels at cumulative
crack-path extensions of 0, 200, 500, and 800 micrometres. Each source state has
two slip systems on a 40-point physical active grid. The family-level
activation-to-line conversion is used both for line insertion and for the
Taylor-backstress density increment.

The FEM/CZM runtime performs:

1. exact selection at a measured crack-extension state;
2. inverse-distance interpolation between family states, using the atlas
   `neighbors` and `power` values;
3. piecewise-linear projection of the measured physical kernel onto the runtime
   80-bin MPZ coordinates;
4. mode-I signed shielding from the interpolated active kernel;
5. zero wake shielding, matching the active-only PF benchmark.

No crack-extension extrapolation is allowed. No constitutive shielding cap is
applied. A proposed cohesive advance beyond the family envelope is rejected
before the committed moving-frame state is changed.

The static `--signed-shielding-kernel` interface remains available only for
legacy unit tests or an explicitly documented fixed-kernel approximation. PF
parity runs must use `--signed-kernel-family`.

## Update an existing v10.0.5.14 installation

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites

cd "$ROOT"
conda activate arrhenius-fem-czm

git pull --ff-only origin \
  v10.0.5.14-persistent-sites-v10222-parity

python -m pip install -e "$ROOT" --no-deps
```

## Focused verification

```bash
python -m pytest -q \
  tests/test_signed_kernel_family_v1005141.py \
  tests/test_persistent_site_v100514.py \
  tests/test_v100513_barrier_only.py \
  tests/test_v1005131_preserved_state.py \
  tests/test_v1005132_startup_resolution_warning.py \
  tests/test_v1005133_tip_only_ramp.py \
  tests/test_v1005134_tip_only_policy_propagation.py \
  tests/test_v1005135_long_corridor.py \
  tests/test_v1005123_phase_c_repairs.py \
  tests/test_mpz_v9_10_unified_transport.py \
  tests/test_v100510_refinement_support.py \
  tests/test_v100511_same_mesh_energy.py
```

## Validate the actual PF family

```bash
PFROOT=/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1
FAMILY_JSON="$PFROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json"

python - <<'PY'
from pathlib import Path
from arrhenius_fracture.signed_kernel_family_v1005141 import (
    SignedShieldingKernelFamilyV1005141,
)

path = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/"
    "runtime_inputs/v10_2_17/"
    "v10_2_14_active_only_campaign_family.json"
)

family = SignedShieldingKernelFamilyV1005141.from_json(path)
print(family.audit_payload())
PY
```

## Candidate 0118 smoke command

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites
PFROOT=/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1
FAMILY_JSON="$PFROOT/runtime_inputs/v10_2_17/v10_2_14_active_only_campaign_family.json"
OUT="$ROOT/runs/v10_0_5_14_1_0118_700K_20um_family_smoke_v1"

cd "$ROOT"
conda activate arrhenius-fem-czm
rm -rf "$OUT"

python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_14_1_persistent_site_family \
  --persistent-site-option v912_peak_0118_persistent_sites \
  --signed-kernel-family "$FAMILY_JSON" \
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
  --snapshot-by-crack-extension-um 10 --no-plots --out "$OUT" \
  2>&1 | tee "$OUT.console.log"
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
- `two_channel_drive_reliable = true`;
- `kernel_artifact_kind = crack_extension_family`;
- `kernel_interpolation_coordinate = cumulative_crack_path_extension_m`;
- `kernel_extrapolation_allowed = false`;
- `wake_kernel_forced_zero = true`;
- `constitutive_K_shield_cap = false`.

A failed tensor probe, atlas mismatch, nonzero transferred recovery, invalid
persistent-source row, or crack extension outside the atlas envelope aborts the
run.
