# MPZ v9.7 Peierls–Taylor entropy calibration

## Why this stage is required

The uncapped v9.6 audit removed the artificial density caps, but the historical
Peierls/Taylor energy ratios (`0.005` and `0.02`) produced essentially zero
macroscopic flow stress for most canonical emission surfaces.  In some low-
temperature rows the scaled reference free energy became nearly zero, so the
forward and reverse detailed-balance rates were numerically indistinguishable
at every applied stress.

This is not only excessive temperature softening.  The activation entropy term
vanishes at the reference temperature, so entropy cannot repair an incorrect
reference-strength magnitude.  The calibration is therefore split into two
stages:

1. determine energy-ratio pairs that produce a finite GPa-scale flow stress;
2. hold those energy ratios fixed and vary Peierls and Taylor activation
   entropies independently to control the temperature dependence.

## Entropy-decoupled free energy

The calibration model retains the emission EXP-floor stress shape but uses
independent mechanism entropies:

```text
G0,P(T) = rH,P G00,emit - S*_P kB (T - Tref)
G0,T(T) = rH,T G00,emit - S*_T kB (T - Tref)
```

`S*_P` and `S*_T` are entered in units of `kB`. Positive activation entropy
softens the mechanism with increasing temperature; negative activation entropy
raises the free-energy barrier with temperature. The code also reports the
equivalent slopes `gT = -S*kB` in eV/K to avoid sign ambiguity.

The calibration model retains:

- EXP-floor Peierls and Taylor barriers;
- forward-minus-reverse detailed balance;
- uncapped forest spacing and Taylor stress amplification;
- uncapped correlation hit order;
- separate mobile and forest roles;
- no density, stress, jump-length, mobile-density, or rate cap.

It is not activated in the MPZ or bulk FEM until the calibration output is
reviewed.

## Smoke run

```bash
OUTROOT=runs/mpz_v9_7_pt_entropy_calibration_smoke_v1 \
ENERGY_RATIO_POINTS=7 \
MAGNITUDE_TOP_PER_CLASS=4 \
ENTROPY_SAMPLES=16 \
bash run_mpz_v9_7_pt_entropy_calibration.sh
```

## Production calibration

```bash
nohup env \
OUTROOT=runs/mpz_v9_7_pt_entropy_calibration_v1 \
ENERGY_RATIO_POINTS=17 \
MAGNITUDE_TOP_PER_CLASS=16 \
ENTROPY_SAMPLES=256 \
TARGET_REFERENCE_STRESS_GPA=2.0 \
bash run_mpz_v9_7_pt_entropy_calibration.sh \
> runs/mpz_v9_7_pt_entropy_calibration_v1.nohup.log 2>&1 &
```

## Outputs

- `pt_magnitude_grid.csv`: reference-strength scan at zero activation entropy;
- `pt_magnitude_selected.csv`: energy pairs retained for thermal calibration;
- `pt_entropy_map_all.csv.gz`: complete independent-entropy map;
- `pt_entropy_map_plausible.csv`: finite, rate-sensitive, positive-barrier rows;
- `pt_entropy_shortlist.csv`: diverse thermal-response families;
- `pt_entropy_calibration_summary.json`: counts and search summary;
- `pt_entropy_calibration_config.json`: exact search definition.

The shortlist preserves strong-softening, moderate-softening, near-athermal, and
thermal-hardening families. A common closure is not imposed at this stage.
