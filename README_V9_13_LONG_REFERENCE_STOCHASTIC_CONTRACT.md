# Correct stochastic contract for the v9.13 long reference map

A long v10.2.22 reference case is valid for continuation of the calibrated
52.081839 micrometre loading map only when it preserves the original common
random number and event-length contract:

```text
CLEAVAGE_HAZARD_MODE=exponential
CLEAVAGE_HAZARD_SEED=3621
CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled
CLEAVAGE_EVENT_MIN_FACTOR=0.5
CLEAVAGE_EVENT_MAX_FACTOR=4.0
CLEAVAGE_EVENT_SUBSEGMENT_FRACTION=0.1
ANISOTROPIC_TRANSPORT_MODE=validated_scalar
ANISOTROPIC_USE_AVALANCHE_BACKEND=1
ANISOTROPIC_EMISSION_ENABLED=1
```

The first attempted 110 micrometre reference omitted the explicit hazard and
event-length environment variables. It therefore used deterministic unit
thresholds and fixed event lengths. Its first 16 events did not reproduce the
calibrated map and the resulting 175-case screen is not a valid long-extension
continuation of the prior candidate campaign.

The loading-map extractor now fails closed unless:

1. `stochastic_hazard.mode` is `exponential` with distribution
   `exponential_unit_mean`;
2. `stochastic_avalanche.mode` is `threshold_scaled` with bounds 0.5 and 4.0;
3. the same integrated threshold controls event length;
4. when `--expected-prefix-loading-map` is supplied, all four loading-map arrays
   reproduce the calibrated map event-by-event before any new events are used.

Use a new 2-D output directory and a new 1-D `OUTROOT`; do not resume or mix the
deterministic-reference campaign with the corrected campaign.
