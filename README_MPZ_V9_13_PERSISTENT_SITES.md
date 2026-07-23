# MPZ v9.13 persistent-site 1-D retest

The v10.2.22 top-five 2-D transfer calibration, including the fixed shared
geometry constants, no-escape validated-scalar transport mode, and 49-case
holdout replay, is documented in
[`README_V9_13_V10222_TRANSFER_CALIBRATION.md`](README_V9_13_V10222_TRANSFER_CALIBRATION.md).
The autonomous \(K\)-versus-crack-extension driver and its immutable-candidate
calibration workflow are documented in
[`README_V9_13_AUTONOMOUS_RCURVE_CALIBRATION.md`](README_V9_13_AUTONOMOUS_RCURVE_CALIBRATION.md).

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

## Mesh-independent source width

`w_eff` is an along-front/out-of-plane correlation width.  The MPZ spacing
`dx` resolves distance ahead of the crack tip and is therefore not an
admissible lower bound for `w_eff`.  The persistent-source implementation now
uses

```text
w_min = max(minimum_front_width_m, abs(b))
```

and records the physical minimum, the ratio `w_eff / w_min`, and whether the
physical floor is active.  `minimum_front_width_m` is a required positive
physical parameter.  The supplied baseline is 10 nm; 10, 50, and 100 nm should
be screened before another long 2-D campaign.  The default physical maximum is
the reference width, unless `maximum_front_width_m` is supplied explicitly.

The implicit backstress root also treats the exact mechanical-blocking state as
a valid boundary solution.  This removes false bracket failures caused by
roundoff leaving an infinitesimal positive drive at the blocking density.

## Moving-tip audit mode

The existing production protocol retains discrete accepted-event translation.
An optional coupled audit mode interleaves kinetics and moving-frame
translation so that each translation is no larger than
`moving_tip_cfl * dx`:

```json
{
  "coupled_moving_tip_enabled": true,
  "moving_tip_cfl": 0.25
}
```

This mode is intentionally off by default.  It materially changes the reduced
1-D signed shielding and must be compared with the internal 2-D shielding
kernel before it can replace the accepted-event contract.

The launcher accepts explicit audit overrides:

```bash
PHYSICAL_MIN_FRONT_WIDTH_NM=50 \
COUPLED_MOVING_TIP_AUDIT=0 \
MODE=smoke \
OUT=runs/v9_13_width_50nm \
bash scripts/run_mpz_v9_13_persistent_top5.sh
```

For the initial physical-width screen, use
`PHYSICAL_MIN_FRONT_WIDTH_NM=10`, `50`, and `100` with separate output
directories and the same signed-kernel normalization.

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
