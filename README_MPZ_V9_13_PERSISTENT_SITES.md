# MPZ v9.13 persistent-site 1-D retest

This additive branch replaces the v9.12 finite source-density reservoir with the
same persistent-site source contract used by the v10.2.21 2-D overlay.

## Constitutive contract

For each reduced slip system,

```text
M_s = rho_source0_m2 * c_arc * r_tip * w_eff
```

where the site multiplicity is persistent and is never consumed by emission.
The effective front width follows the 2-D density-anchored inverse-square-root
law, the tip radius follows accumulated unsigned local slip, and crack advance
resharpens the tip by convecting accumulated slip into the wake.

Emission is solved implicitly against the evolved Taylor backstress. There is
no finite source inventory, no crack-advance source refresh, and no explicit
mobile or retained recovery. Peierls transport, encounter storage, Taylor
release, signed transport, and opposite-sign annihilation remain active.

## Exact 1-D/2-D normalization

The 2-D engine converts one source activation to signed line content using the
mechanically measured kernel field `activation_to_line_content`. Supply the
same signed-kernel family to the 1-D launcher:

```bash
SIGNED_KERNEL_FAMILY_JSON=/path/to/v10_2_14_active_only_campaign_family.json \
MODE=smoke \
bash scripts/run_mpz_v9_13_persistent_top5.sh
```

The launcher extracts the conversion and fails if it is missing, nonpositive,
or state dependent. A source-law-only diagnostic can use the legacy 1-D unit
conversion only through the explicit opt-in:

```bash
ALLOW_UNIT_LINE_CONVERSION=1 \
MODE=smoke \
bash scripts/run_mpz_v9_13_persistent_top5.sh
```

## Runs

Top-1 smoke, 700--1200 K:

```bash
MODE=smoke bash scripts/run_mpz_v9_13_persistent_top5.sh
```

All five candidates, 300--1200 K:

```bash
MODE=full bash scripts/run_mpz_v9_13_persistent_top5.sh
```

Top-two c4 convergence audit, 700--1200 K:

```bash
MODE=convergence bash scripts/run_mpz_v9_13_persistent_top5.sh
```

The main outputs are:

```text
persistent_top5_temperature.csv
persistent_top5_candidate.csv
persistent_top5_summary.json
ranking.csv
```

The runner also stores the pre-correction v9.12 trajectories for direct paired
comparison. The five rows are starting points, not expected to reproduce their
old curves, because the finite source-depletion mechanism has been removed.
