# v10.0.5.14.4 — PF persistent-site long-growth parity

This point release retains the audited PF v10.2.22 persistent-site signed MPZ,
the v10.2.14 crack-extension-indexed shielding atlas, and the validated
v10.0.5.13.5 full 2-D FEM/CZM mechanics.

## Physics retained

- persistent areal source density with no depletion or refresh;
- two signed reduced BCC slip channels;
- unsigned mobile+retained Taylor backstress and signed-retained shielding;
- corrected physical front-width floor;
- implicit emission/backstress complementarity;
- Peierls transport, encounter storage, Taylor release, and zero recovery;
- fractional moving-frame translation and natural resharpening;
- trial/commit and restart preservation;
- kernel-family interpolation by committed crack-path extension.

## v10.0.5.14.4 corrections

### High-temperature transport

The 900–1200 K long-growth campaign exposed false conservation failures in the
v10.0.5.14.3 augmented matrix exponential. The physical mobile/retained state
was mixed with cumulative trapping/release/escape counters, producing poor
scaling for small line populations at high rates.

v10.0.5.14.4 solves only the physical 2N mobile/retained generator with sparse,
L-stable backward Euler. Step doubling controls temporal and state-dependent
rate error. A negligible-active-tail criterion prevents repeated refinement of
a residual population that cannot affect mechanics or the conserved escaped
total. Diagnostic accumulator rows are not part of the solve.

### Adaptive-CZM tip support

The 300–800 K runs reached 80 µm and then repeatedly generated an orphan node
that was also referenced by a cohesive endpoint. The corrected splitter reuses
the backend's authoritative active-tip plus/minus pair, partitions the endpoint
star onto those copies, and repairs a roundoff-induced unsupported side only by
reassigning an orientation-preserving incident triangle. Node numbering and
history arrays are unchanged. The existing quality veto remains fail closed.

### Authoritative accounting

Every step now distinguishes cumulative emitted, active, wake, escaped,
recovered, and discarded line content and reports the signed and relative
balance errors. The inherited `N_em` field remains for compatibility but is
identified as instantaneous active retained content, not cumulative emission.

## Install or update

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites
cd "$ROOT"
conda activate arrhenius-fem-czm
git pull --ff-only origin v10.0.5.14-persistent-sites-v10222-parity
python -m pip install -e "$ROOT" --no-deps
```

Expected package version: `10.0.5.14.4`.

## Verification

```bash
python -m pytest -q \
  tests/test_persistent_site_transport_v1005144.py \
  tests/test_adaptive_czm_tip_support_v1005144.py \
  tests/test_persistent_site_diagnostics_v1005144.py \
  tests/test_persistent_site_transport_v1005143.py \
  tests/test_persistent_site_transport_v1005142.py \
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

## Long-growth campaign

```bash
ROOT=/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_v10_0_5_14_persistent_sites
cd "$ROOT"
conda activate arrhenius-fem-czm
MAX_JOBS=2 bash run_v10_0_5_14_4_0118_300_1200K_200um_family_campaign.sh
```

The default campaign uses candidate 0118, 300–1200 K in 100 K increments,
200 µm extension, 5 µm events, the production PF kernel family, and output root
`v10_0_5_14_4_0118_300_1200K_200um_family_v1`.

## Mandatory runtime invariants

- `available_site_fraction = 1`;
- `source_sites_refreshed = 0`;
- `source_depletion_active = false`;
- `source_refresh_active = false`;
- `front_width_grid_independent = true`;
- `two_channel_drive_reliable = true`;
- `transport_integrator = adaptive_physical_backward_euler_tail_control_v10_0_5_14_4`;
- `transport_cfl_limited = false`;
- all cohesive endpoints have positive bulk incidence;
- final line-content balance is recorded in the production manifest.

The PR remains draft until the corrected 300–1200 K, 200 µm production-atlas
campaign completes and its distance-resolved mechanics and MPZ accounting are
reviewed.
