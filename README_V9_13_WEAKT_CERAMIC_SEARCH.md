# v9.13 weak-temperature and ceramic search

## Purpose

This campaign searches only for two material classes:

- weak-temperature/FCC-like response;
- ceramic-like response.

The reused autonomous 1-D runner retains `DBTT` in its historical module and log
names, but no DBTT or peak candidates are selected or scored by this campaign.

## Why a new zero-D population is required

The earlier large zero-D campaign was designed for DBTT and peak behavior. Its
search policy used 85% local samples around peak-family anchors, constrained the
cleavage temperature coefficient to positive values, constrained source density
to at least `1e14 m^-2`, and evaluated only 700--1400 K. Those choices disfavor
near-temperature-flat initiation and negligible process-zone development.

The failed transfers of candidates `0257068` and `0189364` are retained only as
local anchors. The new population is 75% global and permits:

- cleavage and emission temperature coefficients on both sides of zero;
- source density down to `1e11 m^-2`;
- broader emission, Peierls, Taylor, and blunting coordinates;
- physical-surface positivity checks;
- direct scoring of initiation and developed resistance.

## Reduced search hierarchy

1. Vectorized zero-D proxy: 131,072 new Sobol rows to 100 um.
2. Exact zero-D replay: 128 diverse rows per target class to 100 um.
3. Zero-D promotion: eight rows per target class.
4. Autonomous 1-D validation: sixteen total candidates at five temperatures to
   100 um.
5. Final promotion: up to three rows per target class for 2-D validation.

The default 1-D workload is:

```text
16 candidates x 5 temperatures = 80 cases
```

## Temperature grid

The default coarse grid is:

```text
300 600 900 1100 1200 K
```

These points identify the low-temperature baseline, any low-temperature rise, the
midrange response, the response just below the high-temperature endpoint, and the
terminal high-temperature response. A denser grid can be applied later only to
the one or two candidates selected for the paper campaign.

## Class definitions

### Weak-temperature/FCC-like

The strict gate requires:

- initiation resistance at least 8 MPa sqrt(m);
- initiation span no more than 6 MPa sqrt(m);
- developed-resistance span no more than 8 MPa sqrt(m);
- low-temperature initiation rise no more than 3 MPa sqrt(m);
- weak positive median R-curve rise between 1 and 12 MPa sqrt(m);
- no pronounced developed-state peak.

### Ceramic-like

The strict gate requires:

- initiation resistance at least 8 MPa sqrt(m);
- low/mid-temperature initiation span no more than 5 MPa sqrt(m);
- low-temperature initiation rise no more than 2.5 MPa sqrt(m);
- a 0.5--6 MPa sqrt(m) loss at the high-temperature endpoint;
- median absolute R-curve rise no more than 3 MPa sqrt(m);
- maximum absolute R-curve rise no more than 6 MPa sqrt(m);
- no high-temperature rebound or pronounced peak.

## Required inputs

```text
runs/v9_13_weakT_ceramic_paper_handoff_v1/
  v9_13_weakT_ceramic_paper_handoff.csv

runs/v9_13_long_map_exponential_110um_v2/
  v10_2_22_long_rcurve_loading_map_exponential_110um.json
```

No 1000 um loading map is required.

## Production launch

```bash
conda activate arrhenius-fem-czm-v912-gnd

cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_dbtt_temperature_shelf

OUTROOT="$PWD/runs/v9_13_weakT_ceramic_search_100um_v1"

test ! -e "$OUTROOT" || {
  echo "ERROR: output already exists: $OUTROOT"
  exit 1
}

mkdir -p "$OUTROOT"

stty -tostop 2>/dev/null || true

nohup /usr/bin/caffeinate -dimsu \
  /usr/bin/env \
    PYTHON_BIN="$CONDA_PREFIX/bin/python" \
    OUTROOT="$OUTROOT" \
    bash scripts/run_v913_weakT_ceramic_search.sh \
  </dev/null >"$OUTROOT/driver.log" 2>&1 &

PID=$!
echo "$PID" | tee "$OUTROOT/driver.pid"
disown
```

## Final outputs

```text
runs/v9_13_weakT_ceramic_search_100um_v1/
  zero_d/
    zeroD_weakT_ranked.csv
    zeroD_ceramic_ranked.csv
    combined_promoted_registry.csv
  one_d_100um/
    cases/
  analysis_100um/
    case_temperature_metrics.csv
    candidate_metrics.csv
    weakT_ranked.csv
    ceramic_ranked.csv
    promoted_registry.csv
    summary.json
```

A row with `oneD_strict_gate_passed=false` is only the closest available result
and must not be presented as a successful class match.
