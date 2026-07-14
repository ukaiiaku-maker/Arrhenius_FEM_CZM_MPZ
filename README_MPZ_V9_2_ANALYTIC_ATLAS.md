# MPZ v9.2 analytical first-passage atlas

## Purpose

The v9.1 numerical `first` stage mixed two different questions:

1. intrinsic virgin-tip cleavage/emission competition;
2. development of a shielding/blunting process-zone state.

A virgin moving-process-zone state contains no mobile defects, retained defects,
or accumulated slip. It therefore cannot be required to possess a developed
DBTT shielding state at its first event.

Version 9.2 replaces that optimization stage with an analytical atlas of the
intrinsic EXP-floor hazards. Transient moving-process-zone simulations are
reserved for candidates whose analytical emission exposure makes later state
development plausible.

## Analytical quantities

For a virgin tip,

```text
sigma_tip(K) = K / sqrt(2*pi*r0)
```

The cleavage first-passage toughness satisfies

```text
B_c(K,T) = integral_0^K lambda_c(K',T) dK'/Kdot = 1.
```

The emission exposure accumulated before cleavage is

```text
H_e(Kc,T) = integral_0^Kc lambda_e(K',T) dK'/Kdot,
p_emit = 1 - exp(-H_e),
N_emit_expected = N_source_total * p_emit.
```

The implementation uses the same EXP-floor free-energy surface, attempt
frequencies, and multihit cleavage rate as the production sharp-front engine.
It does not call `eng.step`, does not translate the process zone, and does not
use the legacy scalar `chi_shield`, `N_sat`, artificial backstress, emission
cap, or stored-energy cleavage offset.

## Sampled space

For each of the ceramic, weak-T, and DBTT shape families, a scrambled Sobol
sequence samples the six intrinsic parameters previously varied by Stage 1:

- `cleave_G00_eV`
- `cleave_gT_eV_per_K`
- `cleave_sigc0_GPa`
- `emit_G00_eV`
- `emit_gT_eV_per_K`
- `emit_sigc0_GPa`

The fixed EXP-floor shape parameters (`a`, `n`, floor fraction) are inherited
from each family row in `mpz_three_class_initial_guesses.csv`. This preserves
the three shape families while mapping their full six-dimensional intrinsic
parameter ranges.

## Region definitions

- `ceramic_intrinsic`: decreasing virgin first-passage toughness and negligible
  emission exposure.
- `weakT_intrinsic`: nearly temperature-independent virgin first passage with
  finite but nonsaturated emission exposure.
- `DBTT_precursor`: cleavage-dominated low temperature, at least a two-decade
  emission crossover, and at least 10% source activation at high temperature.
- `emission_saturated`: the source inventory is largely exhausted before
  cleavage and requires separate rejection/audit.

A `DBTT_precursor` is not yet a DBTT material. It is a candidate whose
intrinsic hazards permit a high-temperature process zone to develop. The next
stage must predict and then verify trapping, retention, shielding, blunting,
and the crack-extension scale for convergence toward the developed state.

## Production run

```bash
nohup env \
SAMPLES_PER_FAMILY=16384 \
KDOT_VALUES="0.005 0.02" \
DK=0.10 \
REFINE_DK=0.01 \
bash run_mpz_analytic_first_passage_atlas.sh \
> runs/mpz_v9_2_analytic_first_passage_atlas.nohup.log 2>&1 &
```

The `0.005` MPa sqrt(m)/s rate matches the completed v9.1 numerical Stage 1.
The `0.02` rate matches the earlier Panel-A workflow. Shortlisted rows are
re-evaluated at `REFINE_DK` after the coarse Sobol atlas.

## Outputs

The default output directory is
`runs/mpz_v9_2_analytic_first_passage_atlas`.

Important files:

- `analytic_first_passage_atlas_candidates.csv.gz`: compressed full, replot-ready Sobol atlas;
- `analytic_first_passage_atlas_regions.csv`: compact region metrics and classifications;
- `analytic_first_passage_atlas_shortlist.csv`: top candidates per region;
- `analytic_first_passage_atlas_shortlist_refined.csv`: fine-increment predictions;
- `mpz_analytic_shortlist_material_rows.csv`: complete material rows for the
  subsequent steady/transient MPZ stage;
- `analytic_first_passage_anchor_predictions.csv`: analytical predictions for
  the initial rows and, when available, the completed v9.1 Stage-1 fitted rows;
- `analytic_first_passage_region_counts.csv`: population counts by family/rate;
- `analytic_atlas_region_map.png`: intrinsic toughness/emission map;
- `analytic_atlas_shortlist_Kc_T.png` and
  `analytic_atlas_shortlist_emission_T.png`.

## Interpretation of the earlier Panel-A code

The attached Panel-A scripts established the useful workflow pattern: preserve
raw tables, compute explicit line metrics, and write replot-ready figures. The
old waterfall itself co-varied `chi_shield` and `N_sat`, so it is retained only
as historical context. Those scalar closure variables are not used in this
atlas.
