# One-dimensional V2 native baseline validation

## Decision

Only the PF Peak baseline is sufficient for reduced-model use at present: its initial native $K_J$ error is 3.98% and its one-avalanche topology matches. PF DBTT starts within 12.2% but exceeds the qualified 100-µm radius map after three events. FEM/CZM Peak gets the one-avalanche topology but underpredicts qualified onset $G$ by 61.5%; FEM/CZM DBTT underpredicts onset by 60.2% and exceeds the radius map.

| provider | material_class | status | event_count | physical_avalanche_count | first_event_native_KJ_MPa_sqrt_m | qualified_G_onset_min_J_m2 | initial_onset_relative_error | topology_match | terminal_extension_um | max_tip_radius_um |
|---|---|---|---|---|---|---|---|---|---|---|
| PF | Peak | TARGET_RIGHT_CENSORED | 19 | 1 | 56.0203 |  | 0.0397595 | True | 100 | 7.38524 |
| PF | DBTT | RIGHT_CENSORED_DRIVE_MAP_BOUND | 3 | 3 | 63.0328 |  | 0.121816 | False | 7.87308 | 106.699 |
| FEMCZM | Peak | TARGET_RIGHT_CENSORED | 19 | 1 | 35.7555 | 2872.39 | -0.614703 | True | 100 | 2.93631 |
| FEMCZM | DBTT | RIGHT_CENSORED_DRIVE_MAP_BOUND | 16 | 12 | 39.2059 | 3453.5 | -0.60193 | False | 81.6914 | 117.161 |

The 1000-µm extension is credible only for Peak under both mechanics maps. DBTT is right-censored by a qualification boundary, not physically arrested.
