# v9.13 weak-temperature and ceramic search

## Why a new zero-D population is required

The earlier large zero-D campaign was designed for DBTT and peak behavior. Its
search policy used 85% local samples around peak-family anchors, constrained the
cleavage temperature coefficient to `cleave_gT_eV_per_K >= 5e-4`, constrained
source density to `rho_source0_m2 >= 1e14`, and evaluated only 700--1400 K.
Those choices are appropriate for the DBTT/peak search but structurally disfavor
near-temperature-flat initiation and negligible process-zone development.

The failed transfers of candidates 0257068 and 0189364 are retained only as local
anchors. The new population is 75% global and permits:

- cleavage and emission temperature coefficients on both sides of zero;
- lower source density down to `1e11 m^-2`;
- broader emission, Peierls, Taylor, and blunting coordinates;
- explicit physical-surface positivity at 300 and 1300 K;
- a 300--1300 K temperature grid.

## Search hierarchy

1. **Vectorized zero-D proxy:** 262,144 new Sobol rows to 100 um.
2. **Exact zero-D replay:** 512 diverse rows per class to 100 um.
3. **Autonomous 1-D gate:** 48 rows per class to 100 um.
4. **Final promotion:** up to five rows per class for 2-D validation.

Both initiation resistance and resistance developed by 100 um enter every class
score. The low-temperature rise from 300 through 700 K is an explicit rejection
term. The 100 um gate is sufficient to detect the unwanted strong early R-curve
while avoiding an unnecessary 1000 um campaign.

## Class definitions

### Weak-temperature/FCC-like

The strict gate requires:

- initiation resistance at least 8 MPa sqrt(m);
- initiation span no more than 6 MPa sqrt(m) over 300--1300 K;
- developed-resistance span no more than 8 MPa sqrt(m);
- 300--700 K initiation rise no more than 3 MPa sqrt(m);
- weak positive median R-curve rise between 1 and 12 MPa sqrt(m);
- no pronounced developed-state peak.

### Ceramic-like

The strict gate requires:

- initiation resistance at least 8 MPa sqrt(m);
- low/mid-temperature initiation span no more than 5 MPa sqrt(m);
- 300--700 K initiation rise no more than 2.5 MPa sqrt(m);
- a 0.5--6 MPa sqrt(m) loss concentrated at high temperature;
- median absolute R-curve rise no more than 3 MPa sqrt(m);
- maximum absolute R-curve rise no more than 6 MPa sqrt(m);
- no high-temperature rebound or pronounced peak.

## Required inputs

The anchor handoff created after the failed transfer assessment:

```text
runs/v9_13_weakT_ceramic_paper_handoff_v1/
  v9_13_weakT_ceramic_paper_handoff.csv
```

The calibrated stochastic loading map with at least 100 um coverage:

```text
runs/v9_13_long_map_exponential_110um_v2/
  v10_2_22_long_rcurve_loading_map_exponential_110um.json
```

No new 1000 um loading map is required.

## Production launch

```bash
conda activate arrhenius-fem-czm-v912-gnd

cd /Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v9_13_dbtt_temperature_shelf

OUTROOT="$PWD/runs/v9_13_weakT_ceramic_search_v1"
mkdir -p "$OUTROOT"

stty -tostop 2>/dev/null || true

nohup /usr/bin/caffeinate -dimsu \
  /usr/bin/env \
    PYTHON_BIN="$CONDA_PREFIX/bin/python" \
    TARGET_EXT_UM=100 \
    MAX_JOBS=4 \
    OUTROOT="$OUTROOT" \
    bash scripts/run_v913_weakT_ceramic_search.sh \
  </dev/null >"$OUTROOT/driver.log" 2>&1 &

PID=$!
echo "$PID" | tee "$OUTROOT/driver.pid"
disown
```

Follow progress with:

```bash
tail -f "$OUTROOT/driver.log"
```

## Final outputs

```text
runs/v9_13_weakT_ceramic_search_v1/
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

The final `analysis_100um/promoted_registry.csv` is the only registry intended for
the next 2-D validation campaign. A row with `oneD_strict_gate_passed=false` is
retained only as a closest candidate and must not be presented as a successful
class match.
