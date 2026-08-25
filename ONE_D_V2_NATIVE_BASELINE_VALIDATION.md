# One-dimensional V2 native baseline validation

## 100 µm gate

All four 1000 K cases reach 100 µm with target-right-censor semantics and exactly reproduce their provider's physical-avalanche count. PF onset error is 4.2% for Peak and 13.1% for DBTT. FEM/CZM qualified structural-G onset error is −7.3% and −9.1%; native KJ remains the kinetic metric.

| provider | material_class | status | first_event_native_KJ_MPa_sqrt_m | onset_error_percent | topology | largest_avalanche_fraction | max_tip_radius_um |
| --- | --- | --- | --- | --- | --- | --- | --- |
| PF | Peak | TARGET_RIGHT_CENSORED | 56.2 | 4.23 | 1/1 | 1 | 7.4 |
| PF | DBTT | TARGET_RIGHT_CENSORED | 63.6 | 13.1 | 2/2 | 0.95 | 49.8 |
| FEMCZM | Peak | TARGET_RIGHT_CENSORED | 55.5 | -7.28 | 1/1 | 1 | 18.6 |
| FEMCZM | DBTT | TARGET_RIGHT_CENSORED | 59.3 | -9.05 | 3/3 | 0.882 | 83.4 |

## 1000 µm continuation

| provider | material_class | status | terminal_extension_um | max_tip_radius_um |
| --- | --- | --- | --- | --- |
| PF | Peak | TARGET_RIGHT_CENSORED | 1e+03 | 7.4 |
| PF | DBTT | RIGHT_CENSORED_DRIVE_MAP_BOUND | 835 | 1.03e+03 |
| FEMCZM | Peak | TARGET_RIGHT_CENSORED | 1e+03 | 292 |
| FEMCZM | DBTT | RIGHT_CENSORED_NUMERICAL_BOUND | 448 | 595 |

Peak completes 1000 µm in both providers. Long DBTT is conditionally useful: PF reaches 835 µm before the radius map bound and FEM/CZM reaches 448 µm before the 500-µm opening bound. These are numerical/domain right-censors, not physical arrests.
