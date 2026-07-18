# v10.0.5.3 progressive fatigue loading

## Scope

This point release enables single-front Mode-I fatigue for the v10.0.5.2
parallel opening/tensor-emission kinetic cohesive-zone model. It changes the
loading protocol, not the material model.

The retained v10.0.5.2 physics includes:

- separate opening, cleavage, and tensor-resolved emission drives;
- exact finite source depletion and advance-only source refresh;
- MPZ emission, glide, trapping, release, recovery, and escape state;
- active elastic shielding and local Taylor back stress;
- Strang-split plasticity/cleavage state integration;
- continuous cleavage action and micro-advance;
- progressive `clock_linear` cohesive damage;
- transactional rollback and event localization;
- one committed checkpoint per accepted fatigue outer block, followed by
  re-equilibration before any remaining cycles are reconsidered;
- complete per-channel v10.0.5.2 diagnostics.

No Paris law, crack-growth-per-block law, cyclic degradation law, or new
fatigue-specific constitutive parameter is introduced.

## Numerical loading adapter

The 2-D FEM solve supplies the maximum-load directional `K_open`, `K_cleave`,
and slip-system weights. For each accepted cycle block, midpoint phase
quadrature applies

```
K(phi) = Kmax * [(1 + R)/2 + (1 - R)/2 cos(phi)]
```

with optional zero clipping for compression/closure. Every phase exposure is
advanced by the inherited v10.0.5.2 `integrate_kinetics` method. The existing
adaptive fatigue controller uses a transactional one-cycle prediction to bound
cleavage action and MPZ state increments before accepting a cycle jump.

A topology change stops the accepted outer block after its first committed
checkpoint. Unused block time is not counted as accepted cycles; the next outer
solve re-equilibrates the new crack geometry before selecting another block.

The operator-split cycle jump should be convergence-checked by reducing
`--target-dB`, the `--target-dN-*` limits, and `--max-block-cycles`. The default
production controls are intentionally tighter than the retired fatigue runner.

## Main entry point

```
python -m arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue \
  --v10-material-class DBTT \
  --czm-opening-coupling clock_linear \
  --mode 2d \
  --temperatures 700 \
  --fatigue-cycles --fatigue-hold-load \
  --R 0.1 --frequency-Hz 1000 \
  --cycles-max 1e9 \
  --block-cycles 1e4 --max-block-cycles 1e7 \
  --cycle-block-mode hazard_limited \
  --target-dB 0.01 \
  --target-dN-store 0.01 \
  --target-dN-emit 0.1 \
  --target-dN-mobile 0.1 \
  --n-phase 96 --cyclic-mechanics-phases 16 \
  --steps 4000 --dU 2e-7 --dt 1e-9 \
  --nx 60 --ny 120 \
  --tip-h-fine 2.5e-6 --tip-ratio 1.2 \
  --da-phys 5e-6 --rJ-outer 60e-6 \
  --mpz-length-um 100 --mpz-n-bins 200 \
  --crystal-aniso --crystal-theta-deg 45 \
  --max-fronts 1 --crack-backend adaptive_czm \
  --target-crack-extension-um 50 \
  --out runs/v10_0_5_3_fatigue_pilot
```

`--fatigue-hold-load` is mandatory. The first accepted outer step ramps to the
specified maximum displacement; subsequent steps hold that amplitude while
cycle blocks accumulate.

## Remote stress-range campaign

Use the campaign wrapper when the independent variable is remote stress range:

```
MODE=pilot \
MATERIAL=DBTT \
TEMPERATURES="300 500 700 900" \
DELTA_SIGMA_MPA="200 250 300 350 400" \
R=0.1 \
FREQUENCY_HZ=1000 \
bash run_v10_0_5_3_delta_sigma_fatigue.sh
```

The wrapper performs a model-consistent reaction-force calibration for each
temperature, converts

```
sigma_max = Delta_sigma / (1 - R)
```

to the imposed displacement amplitude, and writes:

- `K_vs_delta_sigma.csv` and `K_vs_delta_sigma.{png,pdf}`;
- `fatigue_growth_points.csv`;
- `da_dN_vs_delta_K.{png,pdf}` when growth is nonzero;
- `remote_stress_calibration.csv`;
- one v10.0.5.3 completion and fatigue audit in every run directory.

`K_vs_delta_sigma.csv` contains requested and reaction-derived stress ranges,
initial/final `KJmax`, inferred `DeltaKJ`, total cycles, crack extension, and
integrated `da/dN`. The block-level file contains local `da/dN`, `KJmax`,
`DeltaKJ`, action increments, and event counts.

## Required convergence gates before production

1. Repeat a representative case with `n_phase = 48, 96, 192`.
2. Halve `target_dB` and all finite `target_dN-*` limits.
3. Reduce `max_block_cycles` by at least a factor of ten.
4. Confirm that `KJmax`, `DeltaKJ`, cycles to first checkpoint, and integrated
   `da/dN` are insensitive within the selected tolerance.
5. Confirm each run's `progressive_fatigue_v10_0_5_3.json` reports the
   progressive trial-CZM lifecycle as active and the legacy fatigue state commit
   as bypassed.
