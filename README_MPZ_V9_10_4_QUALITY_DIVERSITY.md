# MPZ v9.10.4 isotropic quality-diversity promotion

## Scope

Version 9.10.4 is an isolated parameterization-only layer on top of the existing
v9.10.3 DBTT search. It does not change the FEM/CZM solver, crack geometry,
barriers, transport law, retention law, source law, blunting law, or shielding
kernel.

The reduced v9 moving-process-zone model already uses the intended isotropic
Mode-I-equivalent shielding closure,

```text
K_shield = sum_s,i orientation_s N_s,i G b /
           [(1 - nu) sqrt(2 pi max(x_i, core))]
```

and the v9.10.3 candidates must still pass the unchanged v9.10.2 spatial MPZ and
subsequent 2-D validation.

## Why this layer exists

The v9.10.3 optimizer writes `spatial_promotion_manifest.csv` by sorting accepted
candidates by objective and taking the first `shortlist_count`. This can promote
many near-duplicate candidates from one local basin and spend spatial simulations
without testing qualitatively different parameter/response families.

Version 9.10.4 replaces only that shortlist truncation. It uses:

- objective rank as the quality term;
- robust-scaled distance across candidate parameter columns;
- robust-scaled distance across temperature-resolved initiation, plateau,
  `delta_KR`, early-rise, and plateau-slope responses;
- optional restart-lineage preservation;
- a hard guarantee that every accepted candidate is retained when the number of
  accepted candidates is no larger than the promotion budget.

Failed or near-pass candidates are considered only when spare promotion slots
remain.

## Existing search results

No new global search is required when a completed v9.10.3 result is available.
Run:

```bash
SEARCH_ROOT=runs/mpz_v9_10_3_dbtt_targeted_global_search_v1 \
OUTROOT=runs/mpz_v9_10_4_dbtt_quality_diversity_v1 \
COUNT=10 \
bash run_mpz_v9_10_4_dbtt_quality_diversity.sh
```

This writes:

```text
runs/mpz_v9_10_4_dbtt_quality_diversity_v1/DBTT/
  spatial_promotion_manifest.csv
  quality_diversity_selected.csv
  quality_diversity_temperature_detail.csv
  quality_diversity_selection.json
```

The original v9.10.3 result directory is not modified.

## Spatial promotion

Run the unchanged independent-shape spatial MPZ against the new manifest:

```bash
MANIFEST_ROOT=runs/mpz_v9_10_4_dbtt_quality_diversity_v1 \
OUTROOT=runs/mpz_v9_10_4_dbtt_spatial_screen_50um_v1 \
MAX_PER_CLASS=10 \
TARGET_EXTENSION_UM=50 \
bash run_mpz_v9_10_4_dbtt_spatial_promotion.sh
```

Candidates that retain the finite low-temperature shelf and developed
high-temperature DBTT branch should then be extended to 500 micrometres before
2-D FEM/CZM validation.

## Version boundary

- Repository: `Arrhenius_FEM_CZM_MPZ`
- Branch: `v9_10_4_isotropic_parameterization_local`
- Base: final v9 three-class exporter commit `f44cb61`
- v10 repositories and branches are not inputs and are not modified.
