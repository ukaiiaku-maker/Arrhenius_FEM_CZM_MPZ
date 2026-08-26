# Taylor/Peierls Final Decision

## Decision

Retain both controls: Peak `v913_zeroD_sobol_0242980` and DBTT `v913_zeroD_sobol_0202500`. Within the source-qualified bounds, the ten Taylor/Peierls coordinates alone do not create a credible positive reload-separated R-curve response with the existing Peak/DBTT cleavage and emission barriers.

## Recommended full-precision Taylor/Peierls rows

| Coordinate | Peak control | DBTT control |
|---|---:|---:|
| `peierls_H0_eV` | `6.368386985268444` | `3.525720003992319` |
| `peierls_activation_entropy_kB` | `12.2265643812716` | `-36.33932862430811` |
| `peierls_exp_a` | `1.9097672308795155` | `1.448395521938801` |
| `peierls_exp_n` | `1.1372957159765065` | `1.871259921230376` |
| `taylor_H0_eV` | `1.6986502355337143` | `1.4998139148950578` |
| `taylor_activation_entropy_kB` | `38.45378890633583` | `31.365934647619724` |
| `taylor_exp_a` | `0.4561179580539465` | `0.7048910272121429` |
| `taylor_exp_n` | `1.2506696581840515` | `0.5389550719410181` |
| `taylor_corr_rho_c_m2` | `58780073970578.01` | `22042481205294.297` |
| `taylor_corr_scale` | `0.8638954413677038` | `1.2690231726199324` |

The search changed only mobile/retained partition kinetics. In the reduced model, radius and backstress depend on total state, making fracture response structurally invariant to this partition. Direct PF allows spatial redistribution and shows modest onset changes, but Peak has no reload and DBTT softens on reload. There is no true positive model-form response to promote.

No barrier, memory, lifecycle, provider-specific material row, PF equation, FEM/CZM equation, or production trajectory was modified. No FEM/CZM simulation was launched. The only PF source change is a default-off, no-feedback state observer on branch `codex/oneD-v2-taylor-peierls-pf-transfer`.
